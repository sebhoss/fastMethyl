# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# pipeline.R --- example Illumina methylation analysis pipeline.
#
# Single-file driver: the CONFIG block at the top declares every
# per-study setting; the PIPELINE block below runs one analyze() call
# that orchestrates the whole analysis in the statistically correct
# order:
#
#     build  ->  outlierRemoval  ->  normalization(keep)  ->  analysis
#
# analyze() guarantees that order (outliers are detected on the raw
# betas and excluded from the normalisation reference) and closes the
# GDS afterwards, on normal return AND on error.
#
# Retrieve a writable copy of this template with:
#   file.copy(system.file("scripts", "pipeline.R", package = "fastMethyl"),
#             "pipeline.R")
#
# To run on your own data:
#   1. Edit the CONFIG block below with your real paths + identifiers.
#   2. Rscript pipeline.R
#
# To analyse multiple cohorts, copy this file (`pipeline-cohort1.R`,
# `pipeline-cohort2.R`, ...) and edit the CONFIG block in each copy.

suppressPackageStartupMessages({
    library(fastMethyl)
    library(bigmelon)
})

# ====================================================================
# === CONFIG (edit me) ===============================================
# ====================================================================

# --- Pre-processing inputs (consumed by the build phase) ------------
#
# Directory the IDAT files live under. The samplesheet's `Basename` entries are
# resolved relative to this directory (the slide subdirectory in a Basename is
# preserved).
dataDirectory         <- "/path/to/idat-directory"

# Path to the samplesheet CSV. One row per sample, with a `Basename` column whose
# entries resolve to IDAT basenames under dataDirectory. May live anywhere.
samplesheet           <- "/path/to/samplesheet.csv"

# Path to the cross-reactive (non-specific) probes CSV. Must have a column named
# `TargetID` listing CpG names. Common source: Chen et al. 2013 (PMID 23314698)
# for 450k, McCartney et al. 2016 for EPIC.
crossReactiveProbes   <- "/path/to/cross-reactive-probes.csv"

# Prefix for every output file the pipeline writes (GDS, txt, plots).
# Pick something cohort-specific so multiple analyses do not collide.
datasetName           <- "my_study"

# Directory every output (GDS, txt, plots) is written to. Defaults to the
# current working directory, so a bare `Rscript pipeline.r` keeps writing
# where it always has; set an absolute path to make runs reproducible
# regardless of the directory the pipeline is launched from -- otherwise the
# GDS cache is only found when re-run from the same directory.
outputDir             <- "."

# Annotation R package; must be installed. For 450k use
# IlluminaHumanMethylation450kanno.ilmn12.hg19; for EPIC use
# IlluminaHumanMethylationEPICanno.ilm10b4.hg19; for EPIC v2 use
# IlluminaHumanMethylationEPICv2anno.20a1.hg38.
annotationPackage     <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"

# Sample- and probe-level detection-p QC thresholds. A sample is
# dropped if its mean detection-p >= sampleDetPThreshold; a probe is
# dropped if any surviving sample's detection-p >= probeDetPThreshold.
sampleDetPThreshold   <- 0.01
probeDetPThreshold    <- 0.01

# Staged rebuild control. The pipeline is a linear chain
# (build -> normalize -> cell counts), so naming a stage rebuilds it AND every
# stage after it; anything before it is reused from cache. Values:
#   "none"        reuse every cached stage (the default)
#   "normalize"   reuse the built GDS, redo normalisation + cell counts
#                 (e.g. after changing outlyxPerc / the outlier rule)
#   "build"/"all" redo everything from the IDAT read down
# A logical is also accepted (FALSE = "none", TRUE = "all"). A build-cache MISS
# (changed IDATs / QC thresholds / annotation) rebuilds the GDS regardless and
# cascades downstream automatically, so you only name a stage to force a redo
# whose *inputs* did not change.
rebuild               <- "none"

# IDAT-reader worker count = the full available CPU budget. When a cgroup CPU
# quota is in force (an HPC scheduler, a container, or a systemd CPUQuota=),
# parallel::detectCores() over-subscribes because it reports the host's physical
# cores and ignores the quota, oversizing the parallel fan-out and the streaming
# write bursts. .cgroupCpuBudget() reads the effective core budget from cgroup v2
# (cpu.max: "quota period") or v1 (cpu.cfs_quota_us / cpu.cfs_period_us) and caps
# detectCores() to it; where no quota is set it returns the physical-core count
# unchanged. No core is held back -- when a quota applies the cgroup already
# reserves the rest of the machine for everything else.
.cgroupCpuBudget <- function() {
    physical <- parallel::detectCores(logical = FALSE)
    if (is.na(physical) || physical < 1L) physical <- 1L
    quotaCores <- NA_real_
    v2 <- "/sys/fs/cgroup/cpu.max"
    if (file.exists(v2)) {
        f <- strsplit(readLines(v2, n = 1L, warn = FALSE), "\\s+")[[1L]]
        if (length(f) == 2L && f[[1L]] != "max") {
            q <- suppressWarnings(as.numeric(f[[1L]]))
            p <- suppressWarnings(as.numeric(f[[2L]]))
            if (isTRUE(is.finite(q)) && isTRUE(is.finite(p)) && p > 0) {
                quotaCores <- q / p
            }
        }
    } else {
        qf <- "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
        pf <- "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
        if (file.exists(qf) && file.exists(pf)) {
            q <- suppressWarnings(as.numeric(readLines(qf, n = 1L, warn = FALSE)))
            p <- suppressWarnings(as.numeric(readLines(pf, n = 1L, warn = FALSE)))
            if (isTRUE(is.finite(q)) && q > 0 &&
                  isTRUE(is.finite(p)) && p > 0) {
                quotaCores <- q / p
            }
        }
    }
    if (isTRUE(is.finite(quotaCores))) {
        max(1L, min(physical, as.integer(floor(quotaCores))))
    } else {
        physical
    }
}
readerWorkers         <- .cgroupCpuBudget()

# --- Analysis settings (used by the PIPELINE block below) -----------

# Outlier detection. When TRUE, bigmelon::outlyx() flags multivariate outliers on
# the RAW betas and they are excluded from dasen's quantile reference (the
# outlierRemoval hook below). They are still written to normbetas -- analyze()
# normalises every sample, just against an outlier-free reference -- and reported
# to outliers_<class>.txt / mv_outliers_<class>.txt so you can drop them yourself
# downstream. Flip to FALSE to skip outlier detection entirely.
detectOutliers         <- TRUE

# bigmelon::outlyx() percentile for outlier detection.
outlyxPerc             <- 0.01

# (Normalisation and cell-type caching are driven by `rebuild` above: the
# normalize hook and the cell-count step below skip their work when nothing
# upstream changed -- res$rebuildDownstream -- and their output already exists.
# Set rebuild = "normalize" to force a re-run of both.)

# Passed straight through to bigmelon::estimateCellCounts.gds().
gdsPlatform            <- "EPIC"
compositeCellType      <- "Blood"
cellTypes              <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
referencePlatform      <- "IlluminaHumanMethylationEPIC"
estimateCellCountsPerc <- 1

# ====================================================================
# === PIPELINE (do not edit unless you are changing the pipeline) ====
# ====================================================================

# Phase 0: fail fast on the analysis-phase settings BEFORE the expensive
# pre-processing pass. The build phase validates its own inputs, but a
# typo in referencePlatform or an uninstalled reference-data package would
# otherwise surface only at estimateCellCounts.gds -- after the read + QC +
# GDS write have already burned (potentially) hours. Checking them here
# turns that into a sub-second failure.
stopifnot(
    "`rebuild` must be \"none\"/\"normalize\"/\"build\"/\"all\" or a logical" =
        (is.character(rebuild) && length(rebuild) == 1L &&
             rebuild %in% c("none", "normalize", "build", "all")) ||
        (is.logical(rebuild) && length(rebuild) == 1L && !is.na(rebuild)),
    "`detectOutliers` must be a single TRUE or FALSE" =
        is.logical(detectOutliers) && length(detectOutliers) == 1L &&
        !is.na(detectOutliers),
    "`outlyxPerc` must be a single number in (0, 1]" =
        is.numeric(outlyxPerc) && length(outlyxPerc) == 1L &&
        !is.na(outlyxPerc) && outlyxPerc > 0 && outlyxPerc <= 1,
    "`gdsPlatform` must be one of \"450k\", \"EPIC\"" =
        is.character(gdsPlatform) && length(gdsPlatform) == 1L &&
        gdsPlatform %in% c("450k", "EPIC"),
    "`referencePlatform` must be a recognised Illumina platform" =
        is.character(referencePlatform) && length(referencePlatform) == 1L &&
        referencePlatform %in% c("IlluminaHumanMethylation450k",
                                 "IlluminaHumanMethylationEPIC",
                                 "IlluminaHumanMethylationEPICv2"),
    "`compositeCellType` must be a single non-empty string" =
        is.character(compositeCellType) && length(compositeCellType) == 1L &&
        nzchar(compositeCellType),
    "`cellTypes` must be a non-empty character vector" =
        is.character(cellTypes) && length(cellTypes) >= 1L &&
        all(nzchar(cellTypes)),
    "`estimateCellCountsPerc` must be a single positive number" =
        is.numeric(estimateCellCountsPerc) &&
        length(estimateCellCountsPerc) == 1L &&
        !is.na(estimateCellCountsPerc) && estimateCellCountsPerc > 0)

# The analysis-phase platform knobs (gdsPlatform, referencePlatform) are
# independent variables, but estimateCellCounts.gds picks the wrong reference
# panel -- or fails deep in -- if they disagree with the array the IDATs
# carry. annotationPackage names the array unambiguously, so derive the
# expected platform from it and reject a mismatch now rather than hours later.
arrayFamily <- if (grepl("EPICv2", annotationPackage)) {
    "EPICv2"
} else if (grepl("EPIC", annotationPackage)) {
    "EPIC"
} else if (grepl("450k", annotationPackage)) {
    "450k"
} else {
    stop(sprintf(
        "cannot infer the array platform from annotationPackage = \"%s\"; expected it to contain \"450k\", \"EPIC\", or \"EPICv2\".",
        annotationPackage), call. = FALSE)
}
# estimateCellCounts.gds's gdPlatform/referencePlatform only span 450k and
# EPIC, so an EPIC v2 cohort cannot run the cell-composition step in this
# pipeline; say so up front instead of failing inside estimateCellCounts.gds.
if (arrayFamily == "EPICv2") {
    stop("annotationPackage is an EPIC v2 annotation, but the analysis phase (bigmelon::estimateCellCounts.gds) supports only \"450k\" and \"EPIC\" reference data.\n  Use analyze(..., analysis = NULL, normalization = \"dasen\") to build a normalised EPIC v2 GDS, then supply your own cell-composition step.",
        call. = FALSE)
}
if (gdsPlatform != arrayFamily) {
    stop(sprintf(
        "gdsPlatform = \"%s\" does not match the array implied by annotationPackage = \"%s\" (expected gdsPlatform = \"%s\").",
        gdsPlatform, annotationPackage, arrayFamily), call. = FALSE)
}
.expectedRefPlatform <- paste0("IlluminaHumanMethylation", arrayFamily)
if (referencePlatform != .expectedRefPlatform) {
    stop(sprintf(
        "referencePlatform = \"%s\" does not match the array implied by annotationPackage = \"%s\" (expected referencePlatform = \"%s\").",
        referencePlatform, annotationPackage, .expectedRefPlatform), call. = FALSE)
}

# estimateCellCounts.gds resolves its reference data from a package named
# FlowSorted.<compositeCellType>.<platform>, where <platform> is gdsPlatform
# with EPIC folded to 450k -- there is no EPIC reference set, so the EPIC path
# converts to 450k internally and loads the 450k panel. Mirror that resolution
# exactly so this installed-package check targets the package the analysis
# phase will actually load, and a missing one fails now rather than
# post-preprocessing.
.refPlatform <- if (gdsPlatform == "EPIC") "450k" else gdsPlatform
.referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, .refPlatform)
if (!requireNamespace(.referencePkg, quietly = TRUE)) {
    stop(sprintf(
        "reference-data package \"%s\" (inferred from compositeCellType=\"%s\" + gdsPlatform=\"%s\") is not installed.\n  Install it with BiocManager::install(\"%s\").",
        .referencePkg, compositeCellType, gdsPlatform, .referencePkg),
        call. = FALSE)
}

# All outputs (GDS, txt, plots) are written under outputDir. Resolving every
# path through this one helper keeps a run self-contained in one directory
# regardless of where Rscript was launched from -- including the GDS cache,
# which the build phase keys on the path it is handed.
if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
}
outPath <- function(name) file.path(outputDir, name)

# Outlier-detection hook: run on the RAW betas (before normalisation), write the
# outlier reports, and return the keep mask analyze() feeds to the normalisation
# step so dasen's quantile reference excludes the flagged samples. Returning the
# mask -- rather than rebuilding a smaller GDS -- is all that is needed: analyze()
# still normalises every sample, just against an outlier-free reference, so the
# GDS stays whole. Set detectOutliers = FALSE to skip with "none".
outlierRemoval <- if (detectOutliers) {
    function(gds, res) {
        message("Detecting outliers")
        outliers <- outlyx(gds, plot = FALSE, perc = outlyxPerc)
        write.table(outliers[which(outliers$outliers), ],
                    file = outPath(paste0("outliers_", datasetName, ".txt")),
                    quote = FALSE, row.names = FALSE, col.names = TRUE,
                    sep = "\t")
        write.table(outliers[which(outliers$mv), ],
                    file = outPath(paste0("mv_outliers_", datasetName, ".txt")),
                    quote = FALSE, row.names = FALSE, col.names = TRUE,
                    sep = "\t")
        !outliers$outliers
    }
} else {
    "none"
}

# Normalisation hook: cache-aware streaming dasen. gdsDasen() reproduces
# wateRmelon::dasen exactly while reading the GDS one column block at a time, so
# peak memory stays bounded on large cohorts -- and it honours the `keep` mask
# from the outlier step (reference from survivors, normbetas for every sample).
# Recompute only when something upstream changed (res$rebuildDownstream -- a
# rebuild request or a build-cache miss) or no normbetas node exists yet.
normalization <- function(gds, res, keep) {
    if (res$rebuildDownstream || !("normbetas" %in% ls.gdsn(gds))) {
        message("Normalizing (streaming dasen)")
        gdsDasen(gds, node = "normbetas", keep = keep)
    } else {
        message("Reusing existing normbetas node (set rebuild = \"normalize\" to redo)")
    }
}

# One analyze() call runs build -> outlierRemoval -> normalization -> analysis, then
# closes the GDS (on normal return AND on error). The analysis function receives
# the open, QC'd, outlier-aware-normalised GDS and computes PCA + cell composition.
analyze(
    dataDirectory         = dataDirectory,
    crossReactiveProbes   = crossReactiveProbes,
    samplesheet           = samplesheet,
    gdsOutput             = outPath(paste0(datasetName, ".gds")),
    annotationPackage     = annotationPackage,
    sampleDetPThreshold   = sampleDetPThreshold,
    probeDetPThreshold    = probeDetPThreshold,
    rebuild               = rebuild,
    readerWorkers         = readerWorkers,
    outlierRemoval        = outlierRemoval,
    normalization         = normalization,
    analysis = function(gds, res) {

        # PCA on normalised betas (method = "quick"; "sorted" is much slower).
        # The PCA device is closed via on.exit so a plot() error cannot leave a
        # truncated PDF / dangling device behind (on.exit fires here because
        # FUN is a function).
        message("Performing Principal Components Analysis")
        pca <- prcomp(gds, node.name = "normbetas", method = "quick")
        message("Plotting Principal Components")
        pdf(outPath(paste0("pca_", datasetName, ".pdf")))
        pca_device <- dev.cur()
        on.exit(if (pca_device %in% dev.list()) dev.off(pca_device), add = TRUE)
        plot(pca)
        dev.off(pca_device)

        # Cell-type proportions: the slowest, most memory-hungry step, and the
        # last stage in the chain, so its output is cached -- recomputed only when
        # something upstream changed (res$rebuildDownstream, e.g. a fresh GDS or
        # re-normalisation) or the cellTypes_<class>.txt output does not exist yet.
        cellsPath <- outPath(paste0("cellTypes_", datasetName, ".txt"))
        if (res$rebuildDownstream || !file.exists(cellsPath)) {
            message("Estimating cell type proportions")
            cells <- estimateCellCounts.gds(
                gds,
                gdPlatform        = gdsPlatform,
                bn                = "normbetas",
                perc              = estimateCellCountsPerc,
                compositeCellType = compositeCellType,
                cellTypes         = cellTypes,
                referencePlatform = referencePlatform)
            write.table(cells, file = cellsPath, quote = FALSE,
                        row.names = FALSE, col.names = TRUE, sep = "\t")
        } else {
            message("Reusing existing cell-type estimates (set rebuild = \"normalize\" to recompute)")
        }

        invisible(res)
    })
