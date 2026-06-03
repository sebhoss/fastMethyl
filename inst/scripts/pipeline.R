# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# pipeline.R --- example Illumina methylation analysis pipeline.
#
# Single-file driver: the CONFIG block at the top declares every
# per-study setting; the PIPELINE block below runs pre-processing
# (config in, QC'd GDS out, via fastMethyl) and the analysis on the
# resulting GDS (via bigmelon).
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

# --- Pre-processing inputs (consumed by runPreprocess) -----------
#
# Directory containing IDAT files + the samplesheet CSV. The
# samplesheet must be named `samplesheet_<targetPattern>.csv` and
# must contain a `Basename` column with one row per sample (the slide
# subdirectory in the Basename is preserved).
dataDirectory         <- "/path/to/idat-directory"

# Path to the cross-reactive / non-specific probes CSV. Must have a
# column named `TargetID` listing CpG names. Common source: Chen et
# al. 2013 (PMID 23314698) for 450k, McCartney et al. 2016 for EPIC.
nonSpecificProbesPath <- "/path/to/non-specific-probes.csv"

# Component of the samplesheet filename: targetPattern = "batch1"
# expects samplesheet_batch1.csv under dataDirectory.
targetPattern         <- "batch1"

# Prefix for every output file the pipeline writes (GDS, txt, plots).
# Pick something cohort-specific so multiple analyses do not collide.
datasetClass          <- "my_study"

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

# Cache toggle: when FALSE and a GDS at <datasetClass>.gds already
# exists, runPreprocess skips the read pass and reuses it. Flip to
# TRUE to force a rebuild (e.g. after changing IDAT sources / QC).
forceRebuild          <- FALSE

# IDAT-reader worker count. detectCores(logical = FALSE) returns
# physical cores; leave one free for the OS and other interactive work.
readerWorkers         <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

# --- Analysis settings (used by the PIPELINE block below) -----------
# bigmelon::outlyx() percentile for outlier detection.
outlyxPerc             <- 0.01

# When TRUE, samples outlyx flags as outliers are removed before
# normalisation: a smaller, outlier-free <datasetClass>_clean.gds is built and
# used for dasen / PCA / cell counts. The reasons to enable this are
# statistical (outliers should not skew dasen's quantile reference) and memory
# (fewer columns into the EPIC cell-count step, the most OOM-prone part) -- it
# is NOT a speed optimisation (detecting + subsetting roughly offsets the
# few-percent saved downstream). Default FALSE: outliers are reported but kept,
# leaving the removal decision to you.
dropOutliers          <- FALSE

# Cache toggle for dasen: skip normalisation if `normbetas` already
# exists in the GDS. Flip to TRUE to redo.
forceRenorm            <- FALSE

# Passed straight through to bigmelon::estimateCellCounts.gds().
gdsPlatform            <- "EPIC"
compositeCellType      <- "Blood"
cellTypes              <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
referencePlatform      <- "IlluminaHumanMethylationEPIC"
estimateCellCountsPerc <- 1

# Cache toggle for cell-type estimation: when FALSE and the
# cellTypes_<datasetClass>.txt output already exists, the (slow)
# estimateCellCounts.gds step is skipped and the existing estimates reused.
# Flip to TRUE to recompute (e.g. after changing the reference settings).
forceCellCounts        <- FALSE

# ====================================================================
# === PIPELINE (do not edit unless you are changing the pipeline) ====
# ====================================================================

# Phase 0: fail fast on the analysis-phase settings BEFORE the expensive
# pre-processing pass. runPreprocess validates its own inputs, but a
# typo in referencePlatform or an uninstalled reference-data package would
# otherwise surface only at estimateCellCounts.gds -- after the read + QC +
# GDS write have already burned (potentially) hours. Checking them here
# turns that into a sub-second failure.
stopifnot(
    "`outlyxPerc` must be a single number in (0, 1]" =
        is.numeric(outlyxPerc) && length(outlyxPerc) == 1L &&
        !is.na(outlyxPerc) && outlyxPerc > 0 && outlyxPerc <= 1,
    "`forceRenorm` must be a single TRUE or FALSE" =
        is.logical(forceRenorm) && length(forceRenorm) == 1L &&
        !is.na(forceRenorm),
    "`dropOutliers` must be a single TRUE or FALSE" =
        is.logical(dropOutliers) && length(dropOutliers) == 1L &&
        !is.na(dropOutliers),
    "`forceCellCounts` must be a single TRUE or FALSE" =
        is.logical(forceCellCounts) && length(forceCellCounts) == 1L &&
        !is.na(forceCellCounts),
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
    stop("annotationPackage is an EPIC v2 annotation, but the analysis phase (bigmelon::estimateCellCounts.gds) supports only \"450k\" and \"EPIC\" reference data.\n  Use analyze(..., FUN = NULL) to just build the EPIC v2 GDS, then supply your own cell-composition step.",
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
# which runPreprocess keys on the path it is handed.
if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
}
outPath <- function(name) file.path(outputDir, name)

# Phases 1 + 2 in a single analyze() call. analyze() runs the validated
# pre-processing (config in, QC'd LZ4_RA-compressed GDS out), opens the GDS, and
# hands the open handle to FUN below for the analysis -- guaranteeing the GDS is
# closed afterwards, on normal return AND on error, via on.exit(). A failed run
# therefore never leaves the file locked open. analyze() does no normalisation
# of its own: FUN below runs outlier detection on the raw betas first, then its
# own cache-aware, atomic dasen -- the statistically correct order.
analyze(
    dataDirectory         = dataDirectory,
    nonSpecificProbesPath = nonSpecificProbesPath,
    targetPattern         = targetPattern,
    datasetClass          = outPath(datasetClass),
    annotationPackage     = annotationPackage,
    sampleDetPThreshold   = sampleDetPThreshold,
    probeDetPThreshold    = probeDetPThreshold,
    forceRebuild          = forceRebuild,
    readerWorkers         = readerWorkers,
    FUN = function(gds, res) {

        # Outlier detection on the RAW betas, before normalisation: detecting
        # outliers on un-normalised data flags samples whose technical
        # artefacts are still visible; dasen() may otherwise mask them.
        message("Detecting outliers")
        outliers <- outlyx(gds, plot = FALSE, perc = outlyxPerc)
        write.table(outliers[which(outliers$outliers), ],
                    file = outPath(paste0("outliers_", datasetClass, ".txt")),
                    quote = FALSE, row.names = FALSE, col.names = TRUE,
                    sep = "\t")
        write.table(outliers[which(outliers$mv), ],
                    file = outPath(paste0("mv_outliers_", datasetClass, ".txt")),
                    quote = FALSE, row.names = FALSE, col.names = TRUE,
                    sep = "\t")

        # Optionally drop the flagged outliers before normalisation by building
        # an outlier-free <datasetClass>_clean.gds and working on it. analyze()
        # owns the GDS it opened (`gds`); FUN owns the clean GDS it opens here,
        # so it registers its own on.exit() to close it. Columns are read with
        # readex.gdsn(sel = ...) so dropped samples are never materialised; the
        # clean GDS is cached (rebuilt only on forceRebuild) so the normbetas /
        # cell-count results it carries survive across re-runs.
        workGds <- gds
        if (dropOutliers) {
            flagged <- rownames(outliers)[which(outliers$outliers)]
            cleanPath <- outPath(paste0(datasetClass, "_clean.gds"))
            if (length(flagged) > 0L &&
                  (forceRebuild || !file.exists(cleanPath))) {
                message(sprintf("Dropping %d outlier sample(s) into %s",
                                length(flagged), basename(cleanPath)))
                keep <- !(as.character(
                    read.gdsn(index.gdsn(gds, "pData"))$barcode) %in% flagged)
                if (file.exists(cleanPath)) unlink(cleanPath)
                cleanGds <- createfn.gds(cleanPath)
                for (nd in c("betas", "methylated", "unmethylated", "pvals")) {
                    if (nd %in% ls.gdsn(gds)) {
                        add.gdsn(cleanGds, nd,
                                 val = readex.gdsn(index.gdsn(gds, nd),
                                                   sel = list(NULL, keep)),
                                 compress = "LZ4_RA", closezip = TRUE)
                    }
                }
                add.gdsn(cleanGds, "fData",
                         val = read.gdsn(index.gdsn(gds, "fData")),
                         check = FALSE)
                pd <- read.gdsn(index.gdsn(gds, "pData"))
                add.gdsn(cleanGds, "pData", val = pd[keep, , drop = FALSE])
                for (nd in c("history", "paths")) {
                    if (nd %in% ls.gdsn(gds)) {
                        add.gdsn(cleanGds, nd,
                                 val = read.gdsn(index.gdsn(gds, nd)))
                    }
                }
                closefn.gds(cleanGds)
            }
            if (file.exists(cleanPath)) {
                workGds <- openfn.gds(cleanPath, readonly = FALSE)
                on.exit(closefn.gds(workGds), add = TRUE)
            }
        }

        # Cache-aware, atomic normalisation: dasen writes a scratch node that is
        # promoted to `normbetas` only on success, so a crash cannot leave a
        # half-written node to be reused. A leftover scratch node from a prior
        # crash is discarded unconditionally first.
        scratch <- "normbetas_building"
        if (scratch %in% ls.gdsn(workGds)) {
            delete.gdsn(index.gdsn(workGds, scratch))
        }
        if (forceRenorm || !("normbetas" %in% ls.gdsn(workGds))) {
            message("Normalizing")
            add.gdsn(workGds, scratch)
            dasen(workGds, node = scratch)
            if ("normbetas" %in% ls.gdsn(workGds)) {
                delete.gdsn(index.gdsn(workGds, "normbetas"))
            }
            rename.gdsn(index.gdsn(workGds, scratch), "normbetas")
        } else {
            message("Reusing existing normbetas node (set forceRenorm = TRUE to redo)")
        }

        # PCA on normalised betas (method = "quick"; "sorted" is much slower).
        # The PCA device is closed via on.exit so a plot() error cannot leave a
        # truncated PDF / dangling device behind (on.exit fires here because
        # FUN is a function).
        message("Performing Principal Components Analysis")
        pca <- prcomp(workGds, node.name = "normbetas", method = "quick")
        message("Plotting Principal Components")
        pdf(outPath(paste0("pca_", datasetClass, ".pdf")))
        pca_device <- dev.cur()
        on.exit(if (pca_device %in% dev.list()) dev.off(pca_device), add = TRUE)
        plot(pca)
        dev.off(pca_device)

        # Cell-type proportions: the slowest, most memory-hungry step, so its
        # output is cached -- skipped when cellTypes_<class>.txt exists unless
        # forceCellCounts = TRUE.
        cellsPath <- outPath(paste0("cellTypes_", datasetClass, ".txt"))
        if (forceCellCounts || !file.exists(cellsPath)) {
            message("Estimating cell type proportions")
            cells <- estimateCellCounts.gds(
                workGds,
                gdPlatform        = gdsPlatform,
                bn                = "normbetas",
                perc              = estimateCellCountsPerc,
                compositeCellType = compositeCellType,
                cellTypes         = cellTypes,
                referencePlatform = referencePlatform)
            write.table(cells, file = cellsPath, quote = FALSE,
                        row.names = FALSE, col.names = TRUE, sep = "\t")
        } else {
            message("Reusing existing cell-type estimates (set forceCellCounts = TRUE to recompute)")
        }

        invisible(res)
    })
