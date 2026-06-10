# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Prototype + benchmark: a PARALLEL gdsDasen. The profile (benchmark_gdsdasen_
# profile.R) showed gdsDasen is ~96% compute-bound, dominated by per-column qmap
# (rank+approx, ~57%) and dfs2 (densities, ~31%) -- both embarrassingly parallel.
#
# Design: the MASTER is the only process that touches the GDS (reads each block,
# serial IO ~4%, no gdsfmt concurrent-handle question); the per-column COMPUTE is
# fanned out over MulticoreParam workers. With fork, workers reach the in-master
# block matrices through copy-on-write -- no input is shipped, only results come
# back (pass 1: per-chunk offsets + partial reference sums; pass 2: beta columns).
# The two serial points (summing the reference, the ordered append) are cheap.
#
# Validation: the parallel betas must match the serial gdsDasen AND
# wateRmelon::dasen(roco = NULL) to numerical precision (the reference is summed in
# a different order under fan-out, so expect ~1e-10, not bit-identity).
#
#   BENCH_SCRIPT=dev/benchmark_pargdsdasen.R BENCH_N=50,100,200 dev/run-benchmark.sh

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData)
  library(gdsfmt); library(wateRmelon); library(BiocParallel)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
block <- as.integer(Sys.getenv("BENCH_BLOCK", "64"))
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

cgCurrent <- "/sys/fs/cgroup/memory.current"
hasCg <- file.exists(cgCurrent)
readCur <- function() as.numeric(readLines(cgCurrent, n = 1L))
# Time + peak cgroup growth of fn(), run in the main process (the parallel variant
# forks its own workers, so we avoid the nested mcparallel of the other harnesses).
measureDirect <- function(fn) {
  base <- if (hasCg) readCur() else NA_real_
  sf <- tempfile("mem_"); pf <- tempfile("pid_")
  if (hasCg) {
    system(paste0("sh -c 'while true; do cat ", cgCurrent,
                  "; sleep 0.05; done' > ", sf, " 2>/dev/null & echo $! > ", pf),
           wait = FALSE)
    Sys.sleep(0.1)
  }
  t <- unname(system.time(val <- fn())[3])
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
  list(time = t, peak = peak, value = val)
}
mib <- function(x) if (is.na(x)) NA_integer_ else as.integer(round(x / 1048576))

.dfs2offset <- function(col, isI, isII) {
  one <- density(col[isI], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  two <- density(col[isII], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  one$x[which.max(one$y)] - two$x[which.max(two$y)]
}
.qmapToRef <- function(vals, ref) {
  refpos <- (seq_along(ref) - 1) / (length(ref) - 1)
  approx(refpos, ref, (rank(vals) - 1) / (length(vals) - 1),
         ties = list("ordered", mean))$y
}

# Shared setup: returns the handles/metadata both variants need.
.dasenSetup <- function(g) {
  meth <- index.gdsn(g, "methylated"); unmeth <- index.gdsn(g, "unmethylated")
  onetwo <- as.character(read.gdsn(index.gdsn(g, "fData/Type")))
  dims <- objdesp.gdsn(meth)$dim
  list(meth = meth, unmeth = unmeth, onetwo = onetwo,
       nprobe = dims[1L], nsamp = dims[2L],
       isI = onetwo == "I", isII = onetwo == "II")
}
.readBlock <- function(node_h, cols, nprobe) {
  matrix(as.numeric(read.gdsn(node_h, start = c(1L, cols[1L]),
                              count = c(-1L, length(cols)))), nrow = nprobe)
}

# Serial gdsDasen (faithful copy of R/streaming-normalize.R), for the timing
# baseline and the cross-check.
gdsDasenSerial <- function(gp, node, fudge = 100, blk = 64L) {
  g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g), add = TRUE)
  s <- .dasenSetup(g)
  refMI <- numeric(sum(s$isI)); refUI <- numeric(sum(s$isI))
  refMII <- numeric(sum(s$isII)); refUII <- numeric(sum(s$isII))
  offM <- numeric(s$nsamp); offU <- numeric(s$nsamp)
  blocks <- split(seq_len(s$nsamp), ceiling(seq_len(s$nsamp) / blk))
  for (cols in blocks) {
    mblk <- .readBlock(s$meth, cols, s$nprobe); ublk <- .readBlock(s$unmeth, cols, s$nprobe)
    for (k in seq_along(cols)) {
      j <- cols[[k]]
      oM <- .dfs2offset(mblk[, k], s$isI, s$isII); oU <- .dfs2offset(ublk[, k], s$isI, s$isII)
      offM[j] <- oM; offU[j] <- oU
      mc <- mblk[, k]; mc[s$isI] <- mc[s$isI] - oM
      uc <- ublk[, k]; uc[s$isI] <- uc[s$isI] - oU
      refMI <- refMI + sort(mc[s$isI]); refUI <- refUI + sort(uc[s$isI])
      refMII <- refMII + sort(mc[s$isII]); refUII <- refUII + sort(uc[s$isII])
    }
  }
  refMI <- refMI / s$nsamp; refUI <- refUI / s$nsamp
  refMII <- refMII / s$nsamp; refUII <- refUII / s$nsamp
  if (node %in% ls.gdsn(g)) delete.gdsn(index.gdsn(g, node))
  out <- add.gdsn(g, node, valdim = c(s$nprobe, 0L), storage = "double")
  for (cols in blocks) {
    mblk <- .readBlock(s$meth, cols, s$nprobe); ublk <- .readBlock(s$unmeth, cols, s$nprobe)
    betas <- matrix(0, s$nprobe, length(cols))
    for (k in seq_along(cols)) {
      j <- cols[[k]]
      mc <- mblk[, k]; mc[s$isI] <- mc[s$isI] - offM[j]
      uc <- ublk[, k]; uc[s$isI] <- uc[s$isI] - offU[j]
      mc[s$isI] <- .qmapToRef(mc[s$isI], refMI); mc[s$isII] <- .qmapToRef(mc[s$isII], refMII)
      uc[s$isI] <- .qmapToRef(uc[s$isI], refUI); uc[s$isII] <- .qmapToRef(uc[s$isII], refUII)
      betas[, k] <- mc / (mc + uc + fudge)
    }
    append.gdsn(out, betas)
  }
  readmode.gdsn(out); invisible(node)
}

# Parallel gdsDasen: master reads each block; per-column compute fans out over
# BPPARAM workers (COW access to the block, results returned).
gdsDasenPar <- function(gp, node, fudge = 100, blk = 64L, BPPARAM = SerialParam()) {
  g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g), add = TRUE)
  s <- .dasenSetup(g); isI <- s$isI; isII <- s$isII; nprobe <- s$nprobe
  refMI <- numeric(sum(isI)); refUI <- numeric(sum(isI))
  refMII <- numeric(sum(isII)); refUII <- numeric(sum(isII))
  offM <- numeric(s$nsamp); offU <- numeric(s$nsamp)
  blocks <- split(seq_len(s$nsamp), ceiling(seq_len(s$nsamp) / blk))
  nw <- max(1L, bpnworkers(BPPARAM))

  for (cols in blocks) {                                   # pass 1: reference
    mblk <- .readBlock(s$meth, cols, nprobe); ublk <- .readBlock(s$unmeth, cols, nprobe)
    chunks <- parallel::splitIndices(length(cols), nw)
    parts <- bplapply(chunks, function(idx) {
      pMI <- numeric(sum(isI)); pUI <- numeric(sum(isI))
      pMII <- numeric(sum(isII)); pUII <- numeric(sum(isII))
      oM <- numeric(length(idx)); oU <- numeric(length(idx))
      for (kk in seq_along(idx)) {
        k <- idx[[kk]]
        aM <- .dfs2offset(mblk[, k], isI, isII); aU <- .dfs2offset(ublk[, k], isI, isII)
        oM[kk] <- aM; oU[kk] <- aU
        mc <- mblk[, k]; mc[isI] <- mc[isI] - aM
        uc <- ublk[, k]; uc[isI] <- uc[isI] - aU
        pMI <- pMI + sort(mc[isI]); pUI <- pUI + sort(uc[isI])
        pMII <- pMII + sort(mc[isII]); pUII <- pUII + sort(uc[isII])
      }
      list(idx = idx, oM = oM, oU = oU, pMI = pMI, pUI = pUI, pMII = pMII, pUII = pUII)
    }, BPPARAM = BPPARAM)
    for (p in parts) {
      offM[cols[p$idx]] <- p$oM; offU[cols[p$idx]] <- p$oU
      refMI <- refMI + p$pMI; refUI <- refUI + p$pUI
      refMII <- refMII + p$pMII; refUII <- refUII + p$pUII
    }
  }
  refMI <- refMI / s$nsamp; refUI <- refUI / s$nsamp
  refMII <- refMII / s$nsamp; refUII <- refUII / s$nsamp
  if (node %in% ls.gdsn(g)) delete.gdsn(index.gdsn(g, node))
  out <- add.gdsn(g, node, valdim = c(nprobe, 0L), storage = "double")

  for (cols in blocks) {                                   # pass 2: apply
    mblk <- .readBlock(s$meth, cols, nprobe); ublk <- .readBlock(s$unmeth, cols, nprobe)
    gcols <- cols
    chunks <- parallel::splitIndices(length(cols), nw)
    parts <- bplapply(chunks, function(idx) {
      bc <- matrix(0, nprobe, length(idx))
      for (kk in seq_along(idx)) {
        k <- idx[[kk]]; j <- gcols[[k]]
        mc <- mblk[, k]; mc[isI] <- mc[isI] - offM[j]
        uc <- ublk[, k]; uc[isI] <- uc[isI] - offU[j]
        mc[isI] <- .qmapToRef(mc[isI], refMI); mc[isII] <- .qmapToRef(mc[isII], refMII)
        uc[isI] <- .qmapToRef(uc[isI], refUI); uc[isII] <- .qmapToRef(uc[isII], refUII)
        bc[, kk] <- mc / (mc + uc + fudge)
      }
      list(idx = idx, betas = bc)
    }, BPPARAM = BPPARAM)
    betas <- matrix(0, nprobe, length(cols))
    for (p in parts) betas[, p$idx] <- p$betas
    append.gdsn(out, betas)
  }
  readmode.gdsn(out); invisible(node)
}

# --- cohort + build ----------------------------------------------------------
mdDir <- system.file("extdata", package = "minfiData")
srcBase <- sub("_Grn\\.idat$", "",
               list.files(mdDir, pattern = "_Grn\\.idat$", recursive = TRUE, full.names = TRUE))
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
  xr <- file.path(dir, "xreact.csv"); write.csv(data.frame(TargetID = "cg00000029"), xr, row.names = FALSE)
  list(dir = dir, xr = xr)
}

cat(sprintf("parallel gdsDasen | 450k | block=%d | workers=%d\n", block, workers))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))
suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))
BPP <- MulticoreParam(workers = workers)

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  cls <- tempfile("pd_gds_", tmpdir = scratch)
  suppressMessages(analyze(
    dataDirectory = coh$dir, samplesheet = file.path(coh$dir, "samplesheet_bench.csv"),
    crossReactiveProbes = coh$xr, gdsOutput = paste0(cls, ".gds"),
    annotationPackage = annoPkg, rebuild = "all", readerWorkers = workers,
    compress = "", verbose = 0L, analysis = NULL))
  gp <- paste0(cls, ".gds")

  message("  serial");   mser <- measureDirect(function() gdsDasenSerial(gp, "normbetas_ser", blk = block))
  message("  parallel"); mpar <- measureDirect(function() gdsDasenPar(gp, "normbetas_par", blk = block, BPPARAM = BPP))

  g <- openfn.gds(gp, readonly = TRUE)
  M <- read.gdsn(index.gdsn(g, "methylated")); U <- read.gdsn(index.gdsn(g, "unmethylated"))
  onetwo <- as.character(read.gdsn(index.gdsn(g, "fData/Type")))
  ref <- wateRmelon::dasen(M, U, onetwo, fudge = 100, roco = NULL)
  ser <- read.gdsn(index.gdsn(g, "normbetas_ser")); par <- read.gdsn(index.gdsn(g, "normbetas_par"))
  closefn.gds(g)
  okv <- is.finite(ref) & is.finite(par)
  diff_oracle <- max(abs(ref[okv] - par[okv]))
  diff_serial <- max(abs(ser[is.finite(ser) & is.finite(par)] - par[is.finite(ser) & is.finite(par)]))

  unlink(coh$dir, recursive = TRUE); unlink(c(gp, paste0(gp, ".buildkey.rds")))
  gc(FALSE)
  rows[[length(rows) + 1L]] <- data.frame(
    N = n, serial_s = round(mser$time, 1), parallel_s = round(mpar$time, 1),
    speedup = round(mser$time / mpar$time, 2),
    serial_peak = mib(mser$peak), parallel_peak = mib(mpar$peak),
    max_diff_vs_oracle = signif(diff_oracle, 3), max_diff_vs_serial = signif(diff_serial, 3))
  message(sprintf("  -> N=%d  serial=%.1fs parallel=%.1fs (%.2fx)  maxdiff vs dasen=%.2g vs serial=%.2g",
                  n, mser$time, mpar$time, mser$time / mpar$time, diff_oracle, diff_serial))
}
cat("\n=== parallel gdsDasen vs serial: time, peak, equivalence ===\n\n")
print(do.call(rbind, rows), row.names = FALSE)
