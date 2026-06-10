# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Analysis-phase benchmark: how the candidate PCA routines compare on a
# normalised methylation GDS. The preprocessing benchmark (benchmark_analysis.R)
# stops at the QC'd GDS; this one picks up where a real analysis does -- dasen to
# `normbetas`, then PCA -- and times the routes a caller could pick:
#
#   * bigmelon::prcomp(gds, node.name = "normbetas", method = "quick") -- the
#     GDS-native route the shipped pipeline.R uses. Reads the node out of the GDS
#     itself; never materialises a full samples x probes matrix in the caller.
#   * irlba::prcomp_irlba(X, n = ncp) -- truncated SVD over an in-RAM matrix
#     obtained from gdsBetaMatrix(). The hypothesis a fast-PCA swap rests on.
#   * stats::prcomp(X) -- the naive base-R baseline, same in-RAM matrix.
#   * FactoMineR::PCA(X, ncp, graph = FALSE) -- run only if installed.
#
# The matrix routes share one cost the GDS route does not: gdsBetaMatrix() must
# read + transpose the whole node into RAM first. That read is measured on its own
# line so the comparison is honest -- "PCA only" times the decomposition with the
# matrix already resident, "+read" is what a caller actually pays end to end.
#
# Methylation data is p >> n (hundreds of thousands of probes, tens-to-hundreds of
# samples), so the principal subspace lives in at most n dimensions and ncp is
# tiny relative to p. That is exactly the regime irlba is built for -- which is
# why it is worth measuring whether it actually beats the GDS-native prcomp here.
#
# Run it through the wrapper (the only supported entry point), pointing it at this
# script:
#   BENCH_SCRIPT=dev/benchmark_pca.R BENCH_N=50,100,200 dev/run-benchmark.sh
# Env knobs:
#   BENCH_N        comma-separated cohort sizes (shared with the other benchmark)
#   BENCH_WORKERS  reader workers used to build the GDS (PCA itself is serial)
#   BENCH_COMPRESS GDS codec for the build ("" uncompressed, default; "LZ4_RA")
#   BENCH_PCA_NCP  number of components irlba/FactoMineR compute (default 10)
#   BENCH_SCRATCH  where cohorts + GDS files are written (default tempdir())

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData)
  library(gdsfmt); library(bigmelon)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
ncp <- as.integer(Sys.getenv("BENCH_PCA_NCP", "10"))
bpCompress <- Sys.getenv("BENCH_COMPRESS", "")
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

hasIrlba <- requireNamespace("irlba", quietly = TRUE)
hasFactoMineR <- requireNamespace("FactoMineR", quietly = TRUE)

# Peak/average memory from the cgroup (the same approach as benchmark_analysis.R:
# a background poller samples memory.current; figures are net of a pre-call base).
cgCurrent <- "/sys/fs/cgroup/memory.current"
hasCg <- file.exists(cgCurrent)
readCur <- function() as.numeric(readLines(cgCurrent, n = 1L))

# Run fn() in a forked child, returning its value plus elapsed seconds and the
# peak + mean cgroup memory growth (bytes). Forking keeps the parent heap clean
# between routes so each route measures from the same baseline.
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
    job <- parallel::mcparallel(fn(), detached = FALSE)
    out <- parallel::mccollect(job)
  })[3])
  val <- if (length(out)) out[[1L]] else NULL
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

# An N-sample on-disk cohort: symlink the six real 450k IDAT pairs in a cycle,
# plus the samplesheet + cross-reactive CSV analyze() expects.
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

# Build the QC'd GDS (analyze, FUN = NULL) and normalise it (dasen -> normbetas),
# OUTSIDE any measured window -- the build/normalise cost is shared by every PCA
# route, so timing it would only add the same seconds to all of them. Returns the
# GDS path; the file carries a `normbetas` node ready for PCA.
buildNormalisedGds <- function(coh) {
  cls <- tempfile("bench_pca_", tmpdir = scratch)
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
  gp <- paste0(cls, ".gds")
  g <- openfn.gds(gp, readonly = FALSE)
  add.gdsn(g, "normbetas")
  dasen(g, node = "normbetas")
  closefn.gds(g)
  gp
}

cat(sprintf("PCA routes on normalised 450k GDS | ncp=%d | build workers=%d | codec=%s\n",
            ncp, workers, if (nzchar(bpCompress)) bpCompress else "uncompressed"))
cat(sprintf("memory: cgroup peak/mean over baseline%s\n",
            if (hasCg) "" else " UNAVAILABLE"))
cat(sprintf("irlba: %s | FactoMineR: %s\n",
            if (hasIrlba) "yes" else "MISSING", if (hasFactoMineR) "yes" else "absent"))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))

suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  gp <- buildNormalisedGds(coh)

  # GDS-native prcomp: opens its own view of the node, so this is end-to-end.
  m_bm_quick <- measure(function() {
    g <- openfn.gds(gp, readonly = TRUE); on.exit(closefn.gds(g))
    p <- prcomp(g, node.name = "normbetas", method = "quick")
    nrow(p$x)
  })

  # Matrix materialisation, measured on its own so the in-RAM routes can be
  # reported both with and without it.
  m_read <- measure(function() {
    X <- gdsBetaMatrix(gp, node = "normbetas")   # samples x probes
    dim(X)
  })

  m_irlba <- if (hasIrlba) measure(function() {
    X <- gdsBetaMatrix(gp, node = "normbetas")
    p <- irlba::prcomp_irlba(X, n = min(ncp, nrow(X) - 1L),
                             center = TRUE, scale. = FALSE)
    nrow(p$x)
  }) else list(time = NA_real_, peak = NA_real_, avg = NA_real_, value = NULL)

  m_prcomp <- measure(function() {
    X <- gdsBetaMatrix(gp, node = "normbetas")
    p <- stats::prcomp(X, center = TRUE, scale. = FALSE)
    nrow(p$x)
  })

  m_facto <- if (hasFactoMineR) measure(function() {
    X <- gdsBetaMatrix(gp, node = "normbetas")
    p <- FactoMineR::PCA(X, ncp = ncp, scale.unit = FALSE, graph = FALSE)
    nrow(p$ind$coord)
  }) else list(time = NA_real_, peak = NA_real_, avg = NA_real_, value = NULL)

  unlink(coh$dir, recursive = TRUE)
  unlink(paste0(sub("\\.gds$", "", gp), c(".gds", ".gds.buildkey.rds")))
  gc(FALSE)

  rd <- round(m_read$time, 1)
  rows[[length(rows) + 1L]] <- data.frame(
    N = n,
    bm_quick_s = round(m_bm_quick$time, 1),
    read_s = rd,
    irlba_pca_s = round(m_irlba$time - rd, 1),
    irlba_total_s = round(m_irlba$time, 1),
    prcomp_pca_s = round(m_prcomp$time - rd, 1),
    factominer_total_s = round(m_facto$time, 1),
    bm_quick_peak_MiB = mib(m_bm_quick$peak),
    irlba_peak_MiB = mib(m_irlba$peak),
    factominer_peak_MiB = mib(m_facto$peak))
  # Live per-cohort line (unbuffered stderr) so a long sweep shows progress.
  message(sprintf("  -> N=%d  bm_quick=%ss irlba=%ss prcomp=%ss facto=%ss",
                  n, format(round(m_bm_quick$time, 1)), format(round(m_irlba$time, 1)),
                  format(round(m_prcomp$time, 1)), format(round(m_facto$time, 1))))
}

cat("\n=== PCA on normbetas: time (s) + peak memory (MiB) ===\n")
cat("  bm_quick = bigmelon prcomp(method=\"quick\"), GDS-native (reads the node itself)\n")
cat("  read = gdsBetaMatrix() materialise+transpose; *_pca = decomposition with X resident; *_total = read+pca\n\n")
print(do.call(rbind, rows), row.names = FALSE)
