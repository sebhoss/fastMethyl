# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Prototype + benchmark: a streaming, two-pass between-array quantile
# normalisation ("our own dasenrank") that reproduces wateRmelon::dasen without
# ever holding the full samples x probes matrix and without dasenrank's temp.gds.
#
# Why it can stream: dasen = per-sample dfsfit (dye/background offset) + a
# per-probe-type quantile normalisation whose reference is rowMeans of the sorted
# columns (the mean of order statistics, exactly what limma::normalizeQuantiles
# uses). That reference is *accumulable* -- add each sample's sorted column into a
# running sum -- so:
#   pass 1  accumulate the 4 references (meth/unmeth x Type I/II) over the kept
#           samples (outliers excluded via `keep`), one column block at a time;
#   pass 2  re-read each column, dfsfit it, map its ranks onto the reference
#           (rank + approx, limma's ties branch), compute beta, append.
# dfsfit is done in RAM per column, which is what removes dasenrank's temp.gds.
#
# `keep` is a sample mask so this runs AFTER outlier removal -- the reference is
# built only from survivors, honouring outliers-before-normalisation.
#
# Equivalence oracle: wateRmelon::dasen(M, U, onetwo, fudge, roco = NULL) on the
# full in-RAM matrices (per-sample dfsfit, the well-defined case) must match the
# streaming betas to numerical precision.
#
#   BENCH_SCRIPT=dev/benchmark_streamdasen.R BENCH_N=50,100,200 dev/run-benchmark.sh

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData)
  library(gdsfmt); library(bigmelon); library(wateRmelon); library(limma)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
blockSz <- as.integer(Sys.getenv("BENCH_BLOCK", "64"))
percRef <- as.numeric(Sys.getenv("BENCH_REF_PERC", "1"))  # subsample the reference
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

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
    job <- parallel::mcparallel(fn(), detached = FALSE)
    out <- parallel::mccollect(job)
  })[3])
  val <- if (length(out)) out[[1L]] else NULL
  peak <- NA_real_
  if (hasCg) {
    Sys.sleep(0.05)
    pid <- tryCatch(readLines(pf, n = 1L), error = function(e) NA_character_)
    if (length(pid) && !is.na(pid)) {
      system(paste("kill", pid), ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    s <- suppressWarnings(as.numeric(readLines(sf, warn = FALSE)))
    s <- s[is.finite(s)]
    if (length(s)) peak <- max(s) - base
  }
  unlink(c(sf, pf))
  list(time = t, peak = peak, value = val)
}
mib <- function(x) if (is.na(x)) NA_integer_ else as.integer(round(x / 1048576))

# --- the streaming normalisation under test ---------------------------------

# Per-sample dfsfit offset (roco = NULL case): the mode gap between Type I and
# Type II densities, subtracted from the Type I rows. Verbatim from wateRmelon's
# dfs2 so the adjusted intensities match dasen's.
.dfs2 <- function(x, onetwo) {
  one <- density(x[onetwo == "I"], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  two <- density(x[onetwo == "II"], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  one$x[which.max(one$y)] - two$x[which.max(two$y)]
}
.dfsCol <- function(col, onetwo) {
  col[onetwo == "I"] <- col[onetwo == "I"] - .dfs2(col, onetwo)
  col
}

# Stream dasen over `keep` columns of the GDS into `outNode`. refPerc < 1
# subsamples the probes used to build the reference (the dasenrank speed knob);
# the apply step always maps every probe.
streamDasen <- function(gds, outNode = "normbetas_stream", keep = NULL,
                        fudge = 100, block = 64L, refPerc = 1) {
  meth <- index.gdsn(gds, "methylated")
  unmeth <- index.gdsn(gds, "unmethylated")
  onetwo <- as.character(read.gdsn(index.gdsn(gds, "fData/Type")))
  nprobe <- objdesp.gdsn(meth)$dim[1L]
  nsamp <- objdesp.gdsn(meth)$dim[2L]
  if (is.null(keep)) keep <- rep(TRUE, nsamp)
  kept <- which(keep)
  isI <- onetwo == "I"; isII <- onetwo == "II"
  nI <- sum(isI); nII <- sum(isII)
  # probe subset used for the reference (full or subsampled, deterministic).
  selI <- if (refPerc >= 1) seq_len(nI) else
    seq.int(1L, nI, by = max(1L, as.integer(round(1 / refPerc))))
  selII <- if (refPerc >= 1) seq_len(nII) else
    seq.int(1L, nII, by = max(1L, as.integer(round(1 / refPerc))))
  refMI <- numeric(length(selI)); refUI <- numeric(length(selI))
  refMII <- numeric(length(selII)); refUII <- numeric(length(selII))
  # dfs2 (a density over 2^15 points) is the per-column cost; compute it once in
  # pass 1 and cache the scalar offset so pass 2 only subtracts, never re-densities.
  offM <- numeric(length(kept)); offU <- numeric(length(kept))

  readBlock <- function(node, cols) {
    m <- read.gdsn(node, start = c(1L, cols[1L]),
                   count = c(-1L, length(cols)))
    matrix(as.numeric(m), nrow = nprobe)
  }
  blocks <- split(kept, ceiling(seq_along(kept) / block))

  # pass 1: dfsfit each column, cache its offset, accumulate the references.
  pos <- 0L
  for (cols in blocks) {
    M <- readBlock(meth, cols); U <- readBlock(unmeth, cols)
    for (k in seq_along(cols)) {
      pos <- pos + 1L
      oM <- .dfs2(M[, k], onetwo); oU <- .dfs2(U[, k], onetwo)
      offM[pos] <- oM; offU[pos] <- oU
      mc <- M[, k]; mc[isI] <- mc[isI] - oM
      uc <- U[, k]; uc[isI] <- uc[isI] - oU
      refMI <- refMI + sort(mc[isI])[selI]
      refUI <- refUI + sort(uc[isI])[selI]
      refMII <- refMII + sort(mc[isII])[selII]
      refUII <- refUII + sort(uc[isII])[selII]
    }
  }
  nk <- length(kept)
  refMI <- refMI / nk; refUI <- refUI / nk
  refMII <- refMII / nk; refUII <- refUII / nk

  # extendable output node (probes x 0), appended one block at a time.
  if (outNode %in% ls.gdsn(gds)) delete.gdsn(index.gdsn(gds, outNode))
  out <- add.gdsn(gds, outNode, valdim = c(nprobe, 0L), storage = "double")

  # map fractional ranks onto a reference, limma's ties branch.
  qmap <- function(vals, ref) {
    npos <- length(vals)
    refpos <- (seq_along(ref) - 1) / (length(ref) - 1)
    approx(refpos, ref, (rank(vals) - 1) / (npos - 1),
           ties = list("ordered", mean))$y
  }
  # pass 2: re-read, dfsfit from the cached offsets, quantile-map, beta, append.
  pos <- 0L
  for (cols in blocks) {
    M <- readBlock(meth, cols); U <- readBlock(unmeth, cols)
    betas <- matrix(0, nrow = nprobe, ncol = length(cols))
    for (k in seq_along(cols)) {
      pos <- pos + 1L
      mc <- M[, k]; mc[isI] <- mc[isI] - offM[pos]
      uc <- U[, k]; uc[isI] <- uc[isI] - offU[pos]
      mc[isI] <- qmap(mc[isI], refMI); mc[isII] <- qmap(mc[isII], refMII)
      uc[isI] <- qmap(uc[isI], refUI); uc[isII] <- qmap(uc[isII], refUII)
      betas[, k] <- mc / (mc + uc + fudge)
    }
    append.gdsn(out, betas)
  }
  readmode.gdsn(out)
  invisible(outNode)
}

# --- cohort + build (unmeasured) --------------------------------------------

xreactIds <- c("cg00000029", "cg00000108", "cg00000109")
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
  list(dir = dir, xr = xr)
}
buildGds <- function(coh, cls) {
  analyze(dataDirectory = coh$dir, crossReactiveProbes = coh$xr,
          samplesheet = file.path(coh$dir, "samplesheet_bench.csv"), gdsOutput = paste0(cls, ".gds"),
          annotationPackage = annoPkg, rebuild = "all",
          readerWorkers = workers, compress = "", verbose = 0L, analysis = NULL)
  paste0(cls, ".gds")
}

cat(sprintf("streaming dasen prototype | 450k | block=%d | refPerc=%g | workers=%d\n",
            blockSz, percRef, workers))
cat(sprintf("memory: cgroup peak over baseline%s\n", if (hasCg) "" else " UNAVAILABLE"))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))
suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  cls <- tempfile("bench_sd_", tmpdir = scratch)
  gp <- suppressMessages(buildGds(coh, cls))

  # bigmelon dasen.gds -- the pipeline.R baseline (176 s reference).
  m_bm <- measure(function() {
    g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g))
    if ("normbetas" %in% ls.gdsn(g)) delete.gdsn(index.gdsn(g, "normbetas"))
    add.gdsn(g, "normbetas"); dasen(g, node = "normbetas")
    TRUE
  })
  # streaming dasen.
  m_st <- measure(function() {
    g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g))
    streamDasen(g, "normbetas_stream", block = blockSz, refPerc = percRef)
    TRUE
  })

  # equivalence: streaming vs wateRmelon::dasen(roco = NULL) on the full matrices,
  # and vs bigmelon's normbetas (its dasen.gds output).
  g <- openfn.gds(gp, readonly = TRUE)
  onetwo <- as.character(read.gdsn(index.gdsn(g, "fData/Type")))
  M <- read.gdsn(index.gdsn(g, "methylated"))
  U <- read.gdsn(index.gdsn(g, "unmethylated"))
  ref <- wateRmelon::dasen(M, U, onetwo, fudge = 100, roco = NULL)
  st <- read.gdsn(index.gdsn(g, "normbetas_stream"))
  bm <- read.gdsn(index.gdsn(g, "normbetas"))
  closefn.gds(g)
  okv <- is.finite(ref) & is.finite(st)
  dmax <- max(abs(ref[okv] - st[okv]))
  dq <- quantile(abs(ref[okv] - st[okv]), 0.999)
  okb <- is.finite(bm) & is.finite(st)
  dmax_bm <- max(abs(bm[okb] - st[okb]))

  unlink(coh$dir, recursive = TRUE)
  unlink(paste0(cls, c(".gds", ".gds.buildkey.rds")))
  gc(FALSE)

  rows[[length(rows) + 1L]] <- data.frame(
    N = n,
    bm_dasen_s = round(m_bm$time, 1),
    stream_s = round(m_st$time, 1),
    speedup = round(m_bm$time / m_st$time, 2),
    bm_peak_MiB = mib(m_bm$peak),
    stream_peak_MiB = mib(m_st$peak),
    max_diff_vs_dasenNULL = signif(dmax, 3),
    q999_diff = signif(dq, 3),
    max_diff_vs_bmDasen = signif(dmax_bm, 3))
  # Live per-cohort line (unbuffered stderr) so a long sweep shows progress.
  message(sprintf("  -> N=%d  bm_dasen=%ss stream=%ss (%.2fx)  maxdiff=%.2g",
                  n, format(round(m_bm$time, 1)), format(round(m_st$time, 1)),
                  m_bm$time / m_st$time, dmax))
}
cat("\n=== streaming dasen vs bigmelon dasen.gds: time, memory, equivalence ===\n")
cat("  max_diff_vs_dasenNULL = |streaming - wateRmelon::dasen(roco=NULL)| (the oracle)\n")
cat("  max_diff_vs_bmDasen = |streaming - bigmelon normbetas| (differs by the roco term)\n\n")
print(do.call(rbind, rows), row.names = FALSE)
