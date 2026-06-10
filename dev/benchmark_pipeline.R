# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Per-phase analysis-pipeline benchmark: where does the time + memory in
# pipeline.R's FUN actually go? It times each stage of the shipped pipeline
# separately on a real QC'd GDS, in the order pipeline.R runs them:
#
#   build   analyze(FUN = NULL)        -- the preprocessing/GDS build (baseline)
#   outlyx  outlyx(raw betas)          -- multivariate outlier detection
#   dasen   dasen -> normbetas         -- between-sample quantile normalisation
#   prcomp  prcomp(method = "quick")   -- PCA on normbetas (GDS-native)
#   cells   estimateCellCounts.gds     -- reference-based cell composition
#
# The point is to rank the analysis phases against each other (and against the
# build) so optimisation effort lands where the seconds are -- not on PCA, which
# the PCA benchmark already showed is sub-second.
#
# Each phase opens and closes the GDS *inside* its own measured (forked) call, so
# the node a phase writes (dasen's normbetas) is flushed to disk before the next
# phase's fork opens the file -- the on-disk GDS is the only shared state, never a
# handle inherited across fork.
#
# Run through the wrapper, pointing BENCH_SCRIPT at this file:
#   BENCH_SCRIPT=dev/benchmark_pipeline.R BENCH_N=50,100,200 dev/run-benchmark.sh
# Env knobs:
#   BENCH_N        comma-separated cohort sizes
#   BENCH_WORKERS  reader workers for the build (the analysis phases are serial)
#   BENCH_COMPRESS GDS codec for the build ("" uncompressed default; "LZ4_RA")
#   BENCH_SCRATCH  where cohorts + GDS files are written (default tempdir())

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData)
  library(gdsfmt); library(bigmelon)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
bpCompress <- Sys.getenv("BENCH_COMPRESS", "")
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

# This cohort is built from the 450k minfiData IDATs, so the cell-composition
# reference is the 450k blood panel.
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
hasRef <- requireNamespace("FlowSorted.Blood.450k", quietly = TRUE)

cgCurrent <- "/sys/fs/cgroup/memory.current"
hasCg <- file.exists(cgCurrent)
readCur <- function() as.numeric(readLines(cgCurrent, n = 1L))

# Run fn() in a forked child, returning its value plus elapsed seconds and the
# peak + mean cgroup memory growth (bytes), net of a pre-call baseline. A phase
# that errors (e.g. a reference panel that will not project the cohort) is caught
# so the sweep still reports the phases that ran.
measure <- function(fn) {
  base <- if (hasCg) readCur() else NA_real_
  sf <- tempfile("mem_"); pf <- tempfile("pid_")
  if (hasCg) {
    system(paste0(
      "sh -c 'while true; do cat ", cgCurrent, "; sleep 0.05; done' > ",
      sf, " 2>/dev/null & echo $! > ", pf), wait = FALSE)
    Sys.sleep(0.1)
  }
  out <- NULL
  t <- unname(system.time({
    job <- parallel::mcparallel(
      tryCatch(fn(), error = function(e) structure(NA, err = conditionMessage(e))),
      detached = FALSE)
    out <- parallel::mccollect(job)
  })[3])
  val <- if (length(out)) out[[1L]] else NULL
  err <- if (!is.null(val) && !is.null(attr(val, "err"))) attr(val, "err") else NULL
  if (!is.null(err)) message("    phase error: ", err)
  peak <- NA_real_; avg <- NA_real_
  if (hasCg) {
    Sys.sleep(0.05)
    pid <- tryCatch(readLines(pf, n = 1L), error = function(e) NA_character_)
    if (length(pid) && !is.na(pid)) {
      system(paste("kill", pid), ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    s <- suppressWarnings(as.numeric(readLines(sf, warn = FALSE)))
    s <- s[is.finite(s)]
    if (length(s)) { peak <- max(s) - base; avg <- mean(s) - base }
  }
  unlink(c(sf, pf))
  list(time = t, peak = peak, avg = avg, value = val)
}
mib <- function(x) if (is.na(x)) NA_integer_ else as.integer(round(x / 1048576))

sampleThreshold <- 0.01
probeThreshold  <- 0.01
xreactIds <- c("cg00000029", "cg00000108", "cg00000109",
               "cg00000165", "cg00000236")

mdDir <- system.file("extdata", package = "minfiData")
srcBase <- sub("_Grn\\.idat$", "",
               list.files(mdDir, pattern = "_Grn\\.idat$",
                          recursive = TRUE, full.names = TRUE))
stopifnot(length(srcBase) >= 1L)

makeCohort <- function(n) {
  dir <- tempfile("bench_cohort_", tmpdir = scratch)
  dir.create(dir)
  bns <- sprintf("S%05d", seq_len(n))
  for (i in seq_len(n)) {
    src <- srcBase[(i - 1L) %% length(srcBase) + 1L]
    for (ch in c("_Grn.idat", "_Red.idat")) {
      file.symlink(paste0(src, ch), file.path(dir, paste0(bns[i], ch)))
    }
  }
  write.csv(
    data.frame(Sample_Name = bns, Basename = bns, stringsAsFactors = FALSE),
    file.path(dir, "samplesheet_bench.csv"), row.names = FALSE)
  xreact <- file.path(dir, "xreact.csv")
  write.csv(data.frame(TargetID = xreactIds, stringsAsFactors = FALSE),
            xreact, row.names = FALSE)
  list(dir = dir, bns = bns, xreact = xreact)
}

buildGds <- function(coh, cls) {
  analyze(
    dataDirectory         = coh$dir,
    crossReactiveProbes   = coh$xreact,
    samplesheet           = file.path(coh$dir, "samplesheet_bench.csv"),
    gdsOutput             = paste0(cls, ".gds"),
    annotationPackage     = annoPkg,
    sampleDetPThreshold   = sampleThreshold,
    probeDetPThreshold    = probeThreshold,
    rebuild               = "all",
    readerWorkers         = workers,
    compress              = bpCompress,
    verbose               = 0L,
    analysis              = NULL)
  paste0(cls, ".gds")
}

cat(sprintf("pipeline.R phase breakdown | 450k | build workers=%d | codec=%s\n",
            workers, if (nzchar(bpCompress)) bpCompress else "uncompressed"))
cat(sprintf("memory: cgroup peak/mean over baseline%s | FlowSorted.Blood.450k: %s\n",
            if (hasCg) "" else " UNAVAILABLE", if (hasRef) "yes" else "MISSING"))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))

suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  cls <- tempfile("bench_pl_", tmpdir = scratch)

  message("  build")
  m_build <- measure(function() {
    suppressMessages(buildGds(coh, cls))
    1L
  })
  gp <- paste0(cls, ".gds")

  message("  outlyx")
  m_outlyx <- measure(function() {
    g <- openfn.gds(gp, readonly = TRUE); on.exit(closefn.gds(g))
    o <- outlyx(g, plot = FALSE, perc = 0.01)
    sum(o$outliers)
  })

  message("  dasen")
  m_dasen <- measure(function() {
    g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g))
    if ("normbetas" %in% ls.gdsn(g)) delete.gdsn(index.gdsn(g, "normbetas"))
    add.gdsn(g, "normbetas")
    dasen(g, node = "normbetas")
    TRUE
  })

  message("  prcomp")
  m_prcomp <- measure(function() {
    g <- openfn.gds(gp, readonly = TRUE); on.exit(closefn.gds(g))
    p <- prcomp(g, node.name = "normbetas", method = "quick")
    nrow(p$x)
  })

  message("  cells")
  m_cells <- if (hasRef) measure(function() {
    g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g))
    cc <- estimateCellCounts.gds(
      g, gdPlatform = "450k", bn = "normbetas", perc = 1,
      compositeCellType = "Blood", cellTypes = cellTypes,
      referencePlatform = "IlluminaHumanMethylation450k")
    nrow(cc)
  }) else list(time = NA_real_, peak = NA_real_, avg = NA_real_, value = NULL)

  unlink(coh$dir, recursive = TRUE)
  unlink(paste0(cls, c(".gds", ".gds.buildkey.rds")))
  gc(FALSE)

  rows[[length(rows) + 1L]] <- data.frame(
    N = n,
    build_s = round(m_build$time, 1),
    outlyx_s = round(m_outlyx$time, 1),
    dasen_s = round(m_dasen$time, 1),
    prcomp_s = round(m_prcomp$time, 1),
    cells_s = round(m_cells$time, 1),
    build_peak = mib(m_build$peak),
    outlyx_peak = mib(m_outlyx$peak),
    dasen_peak = mib(m_dasen$peak),
    cells_peak = mib(m_cells$peak))
  # Live per-cohort line (unbuffered stderr) so a long sweep shows progress.
  message(sprintf("  -> N=%d  build=%ss outlyx=%ss dasen=%ss prcomp=%ss cells=%ss",
                  n, format(round(m_build$time, 1)), format(round(m_outlyx$time, 1)),
                  format(round(m_dasen$time, 1)), format(round(m_prcomp$time, 1)),
                  format(round(m_cells$time, 1))))
}

cat("\n=== pipeline.R analysis phases: time (s) + peak memory (MiB), in run order ===\n")
cat("  build = analyze(analysis=NULL); then outlyx -> dasen -> prcomp -> estimateCellCounts.gds\n\n")
print(do.call(rbind, rows), row.names = FALSE)
