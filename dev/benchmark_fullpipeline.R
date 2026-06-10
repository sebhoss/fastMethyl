# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Apples-to-apples FULL-PIPELINE benchmark: the complete methylation workflow a
# user actually runs, both ways, measuring runtime + peak memory.
#
#   upstream    minfi + bigmelon/wateRmelon, the idiomatic path:
#               read.metharray.exp -> detectionP -> preprocessRaw -> QC -> es2gds
#               (build) ; then outlyx -> dasen -> prcomp -> estimateCellCounts.gds
#               (analysis)
#   fastMethyl  analyze(FUN = NULL) (build) ; then outlyx -> gdsDasen -> prcomp ->
#               estimateCellCounts.gds (analysis)
#
# The comparison is fair by construction:
#   * Both end at the SAME artifacts -- a QC'd bigmelon GDS, a normbetas node, a
#     PCA, and a cell-composition table -- with the SAME QC (sample/probe
#     detection-p, sex-chromosome + cross-reactive drops).
#   * The analysis phase is IDENTICAL code on both sides EXCEPT normalisation:
#     upstream runs bigmelon's in-place dasen.gds, fastMethyl runs the streaming
#     gdsDasen (value-equivalent up to dasen's optional sentrix-position term).
#     outlyx / prcomp / estimateCellCounts.gds are the same call on both, so they
#     add the same cost to both columns -- the difference is the build + dasen.
#   * Both run dasen over ALL samples (gdsDasen keep = NULL) so the normalisation
#     is the same task; fastMethyl can additionally exclude outliers from the
#     reference via a keep mask with no rebuild, but that is left out here to keep
#     the dasen step a like-for-like measurement.
#
# build and analysis are timed as separate segments so the structural win (the
# parallel streaming build) is visible next to the largely-shared analysis.
#
#   BENCH_SCRIPT=dev/benchmark_fullpipeline.R BENCH_N=50,100,200 dev/run-benchmark.sh

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData)
  library(gdsfmt); library(bigmelon)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

cgCurrent <- "/sys/fs/cgroup/memory.current"
hasCg <- file.exists(cgCurrent)
readCur <- function() as.numeric(readLines(cgCurrent, n = 1L))
measure <- function(fn) {
  base <- if (hasCg) readCur() else NA_real_
  sf <- tempfile("mem_"); pf <- tempfile("pid_")
  if (hasCg) {
    system(paste0("sh -c 'while true; do cat ", cgCurrent,
                  "; sleep 0.05; done' > ", sf, " 2>/dev/null & echo $! > ", pf),
           wait = FALSE)
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
  if (!is.null(val) && !is.null(attr(val, "err"))) message("    ERROR: ", attr(val, "err"))
  peak <- NA_real_
  if (hasCg) {
    Sys.sleep(0.05)
    pid <- tryCatch(readLines(pf, n = 1L), error = function(e) NA_character_)
    if (length(pid) && !is.na(pid)) system(paste("kill", pid),
                                            ignore.stdout = TRUE, ignore.stderr = TRUE)
    s <- suppressWarnings(as.numeric(readLines(sf, warn = FALSE)))
    s <- s[is.finite(s)]
    if (length(s)) peak <- max(s) - base
  }
  unlink(c(sf, pf))
  list(time = t, peak = peak)
}
mib <- function(x) if (is.na(x)) NA_integer_ else as.integer(round(x / 1048576))

sampleThreshold <- 0.01
probeThreshold  <- 0.01
xreactIds <- c("cg00000029", "cg00000108", "cg00000109", "cg00000165")

mdDir <- system.file("extdata", package = "minfiData")
srcBase <- sub("_Grn\\.idat$", "",
               list.files(mdDir, pattern = "_Grn\\.idat$",
                          recursive = TRUE, full.names = TRUE))

makeCohort <- function(n) {
  dir <- tempfile("bench_cohort_", tmpdir = scratch); dir.create(dir)
  bns <- sprintf("S%05d", seq_len(n))
  for (i in seq_len(n)) {
    src <- srcBase[(i - 1L) %% length(srcBase) + 1L]
    for (ch in c("_Grn.idat", "_Red.idat")) {
      file.symlink(paste0(src, ch), file.path(dir, paste0(bns[i], ch)))
    }
  }
  write.csv(data.frame(Sample_Name = bns, Basename = bns),
            file.path(dir, "samplesheet_bench.csv"), row.names = FALSE)
  xr <- file.path(dir, "xreact.csv")
  write.csv(data.frame(TargetID = xreactIds), xr, row.names = FALSE)
  list(dir = dir, bns = bns, xr = xr)
}

# Upstream build: minfi read + preprocess + the same QC, then es2gds (bigmelon's
# standard MethylSet -> GDS writer). Returns the GDS path.
upstreamBuild <- function(coh, gp) {
  targets <- data.frame(Basename = file.path(coh$dir, coh$bns))
  rg <- read.metharray.exp(targets = targets)
  detP <- detectionP(rg)
  mset <- preprocessRaw(rg)
  keepS <- colMeans(detP) < sampleThreshold
  detP <- detP[, keepS, drop = FALSE]; mset <- mset[, keepS]
  keepP <- rowSums(detP >= probeThreshold) == 0L
  ann <- getAnnotation(mset)
  sexP <- rownames(mset) %in% ann$Name[ann$chr %in% c("chrX", "chrY")]
  xrP <- rownames(mset) %in% xreactIds
  mset <- mset[keepP & !sexP & !xrP, ]
  if (file.exists(gp)) unlink(gp)
  gds <- es2gds(mset, gp)
  closefn.gds(gds)
  invisible(gp)
}

# fastMethyl build: the streaming analyze(FUN = NULL).
fastMethylBuild <- function(coh, cls) {
  analyze(dataDirectory = coh$dir, crossReactiveProbes = coh$xr,
          samplesheet = file.path(coh$dir, "samplesheet_bench.csv"), gdsOutput = paste0(cls, ".gds"),
          annotationPackage = annoPkg, sampleDetPThreshold = sampleThreshold,
          probeDetPThreshold = probeThreshold, rebuild = "all",
          readerWorkers = workers, compress = "", verbose = 0L, analysis = NULL)
  paste0(cls, ".gds")
}

# Analysis on a built GDS, identical on both sides except the dasen call.
# normFn(gds) writes the normbetas node.
runAnalysis <- function(gp, normFn) {
  g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g))
  invisible(wateRmelon::outlyx(g, plot = FALSE, perc = 0.01))
  normFn(g)
  invisible(prcomp(g, node.name = "normbetas", method = "quick"))
  invisible(estimateCellCounts.gds(
    g, gdPlatform = "450k", bn = "normbetas", perc = 1,
    compositeCellType = "Blood", cellTypes = cellTypes,
    referencePlatform = "IlluminaHumanMethylation450k"))
  TRUE
}
dasenUpstream <- function(g) { add.gdsn(g, "normbetas"); dasen(g, node = "normbetas") }
dasenFastMethyl <- function(g) gdsDasen(g, node = "normbetas")

cat(sprintf("FULL pipeline: upstream (minfi+bigmelon) vs fastMethyl | 450k | workers=%d\n", workers))
cat(sprintf("memory: cgroup peak over baseline%s | cohorts: %s\n\n",
            if (hasCg) "" else " UNAVAILABLE", paste(Ns, collapse = ", ")))
suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  ugp <- tempfile("up_gds_", tmpdir = scratch, fileext = ".gds")
  fcls <- tempfile("fm_gds_", tmpdir = scratch)

  message("  upstream build");    ub <- measure(function() suppressMessages(upstreamBuild(coh, ugp)))
  message("  upstream analysis"); ua <- measure(function() suppressMessages(runAnalysis(ugp, dasenUpstream)))
  message("  fastMethyl build");    fb <- measure(function() suppressMessages(fastMethylBuild(coh, fcls)))
  fgp <- paste0(fcls, ".gds")
  message("  fastMethyl analysis"); fa <- measure(function() suppressMessages(runAnalysis(fgp, dasenFastMethyl)))

  unlink(coh$dir, recursive = TRUE)
  unlink(c(ugp, fgp, paste0(fcls, ".gds.buildkey.rds")))
  gc(FALSE)

  upTot <- ub$time + ua$time
  fmTot <- fb$time + fa$time
  upPeak <- mib(max(ub$peak, ua$peak, na.rm = TRUE))
  fmPeak <- mib(max(fb$peak, fa$peak, na.rm = TRUE))
  rows[[length(rows) + 1L]] <- data.frame(
    N = n,
    up_build_s = round(ub$time, 1), fm_build_s = round(fb$time, 1),
    up_anal_s = round(ua$time, 1),  fm_anal_s = round(fa$time, 1),
    up_total_s = round(upTot, 1),   fm_total_s = round(fmTot, 1),
    total_speedup = round(upTot / fmTot, 2),
    up_build_peak = mib(ub$peak),   fm_build_peak = mib(fb$peak),
    up_peak_MiB = upPeak, fm_peak_MiB = fmPeak)

  # Emit this cohort's result as soon as it is computed (via message -> stderr,
  # which is unbuffered, so a long sweep shows progress instead of only printing
  # the table at the very end).
  message(sprintf(
    paste0("  -> N=%d  total up=%.1fs fm=%.1fs (%.2fx)  ",
           "build up=%.1fs fm=%.1fs  peak up=%dM fm=%dM"),
    n, upTot, fmTot, upTot / fmTot, ub$time, fb$time, upPeak, fmPeak))
}
cat("\n=== Full pipeline: upstream (minfi+bigmelon) vs fastMethyl ===\n")
cat("  build = IDAT->QC'd GDS; anal = outlyx+dasen+prcomp+estimateCellCounts (identical except dasen)\n")
cat("  *_peak_MiB = peak over the whole pipeline (max of build/analysis segments)\n\n")
print(do.call(rbind, rows), row.names = FALSE)
