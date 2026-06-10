# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Profile gdsDasen: where does its wall-clock go? It is now the largest analysis
# phase (~150 s at N=200) and the one major phase still serial, so before
# deciding whether to parallelise it we need its compute-vs-IO split.
#
# This instruments a faithful copy of the shipped gdsDasen (R/streaming-normalize.R)
# with per-operation timers, accumulated across both streaming passes:
#   read   read.gdsn of the methylated/unmethylated blocks (IO)        -- IO
#   dfs2   the per-column dye-offset densities (density, n=2^15)        -- compute
#   sort   pass-1 reference accumulation (sort each Type I/II group)    -- compute
#   qmap   pass-2 rank + approx mapping onto the reference              -- compute
#   apply  pass-2 offset subtraction + beta = M/(M+U+fudge)            -- compute
#   write  append.gdsn of the beta blocks                              -- IO
#
# compute = dfs2+sort+qmap+apply (throttle-independent, production-accurate);
# io = read+write (inflated by the wrapper's IOWriteBandwidthMax, so an upper
# bound on the IO share). If compute dominates even under the throttle ->
# parallelise. If io dominates -> the lever is reducing IO (or the perc
# subsample to cut the sort cost), not parallel compute.
#
#   BENCH_SCRIPT=dev/benchmark_gdsdasen_profile.R BENCH_N=50,100,200 dev/run-benchmark.sh

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData); library(gdsfmt)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "50,100,200"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
block <- as.integer(Sys.getenv("BENCH_BLOCK", "64"))
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

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
  write.csv(data.frame(TargetID = "cg00000029"), xr, row.names = FALSE)
  list(dir = dir, xr = xr)
}

# Instrumented copy of gdsDasen's algorithm. Returns the per-op timer totals.
profGdsDasen <- function(gp, fudge = 100, blk = 64L) {
  tm <- new.env()
  for (k in c("read", "dfs2", "sort", "qmap", "apply", "write")) tm[[k]] <- 0
  tic <- function() proc.time()[["elapsed"]]

  g <- openfn.gds(gp, readonly = FALSE); on.exit(closefn.gds(g), add = TRUE)
  meth <- index.gdsn(g, "methylated"); unmeth <- index.gdsn(g, "unmethylated")
  onetwo <- as.character(read.gdsn(index.gdsn(g, "fData/Type")))
  dims <- objdesp.gdsn(meth)$dim; nprobe <- dims[1L]; nsamp <- dims[2L]
  isI <- onetwo == "I"; isII <- onetwo == "II"
  refMI <- numeric(sum(isI)); refUI <- numeric(sum(isI))
  refMII <- numeric(sum(isII)); refUII <- numeric(sum(isII))
  offM <- numeric(nsamp); offU <- numeric(nsamp)

  readBlock <- function(node_h, cols) {
    t0 <- tic()
    m <- matrix(as.numeric(read.gdsn(node_h, start = c(1L, cols[1L]),
                                     count = c(-1L, length(cols)))), nrow = nprobe)
    tm$read <- tm$read + (tic() - t0); m
  }
  dfs2 <- function(col) {
    t0 <- tic()
    one <- density(col[isI], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
    two <- density(col[isII], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
    r <- one$x[which.max(one$y)] - two$x[which.max(two$y)]
    tm$dfs2 <- tm$dfs2 + (tic() - t0); r
  }
  qmap <- function(vals, ref) {
    t0 <- tic()
    refpos <- (seq_along(ref) - 1) / (length(ref) - 1)
    r <- approx(refpos, ref, (rank(vals) - 1) / (length(vals) - 1),
                ties = list("ordered", mean))$y
    tm$qmap <- tm$qmap + (tic() - t0); r
  }
  blocks <- split(seq_len(nsamp), ceiling(seq_len(nsamp) / blk))

  for (cols in blocks) {                                  # pass 1: reference
    mblk <- readBlock(meth, cols); ublk <- readBlock(unmeth, cols)
    for (k in seq_along(cols)) {
      j <- cols[[k]]
      oM <- dfs2(mblk[, k]); oU <- dfs2(ublk[, k]); offM[j] <- oM; offU[j] <- oU
      t0 <- tic()
      mc <- mblk[, k]; mc[isI] <- mc[isI] - oM
      uc <- ublk[, k]; uc[isI] <- uc[isI] - oU
      refMI <- refMI + sort(mc[isI]); refUI <- refUI + sort(uc[isI])
      refMII <- refMII + sort(mc[isII]); refUII <- refUII + sort(uc[isII])
      tm$sort <- tm$sort + (tic() - t0)
    }
  }
  refMI <- refMI / nsamp; refUI <- refUI / nsamp
  refMII <- refMII / nsamp; refUII <- refUII / nsamp

  if ("normbetas_prof" %in% ls.gdsn(g)) delete.gdsn(index.gdsn(g, "normbetas_prof"))
  out <- add.gdsn(g, "normbetas_prof", valdim = c(nprobe, 0L), storage = "double")

  for (cols in blocks) {                                  # pass 2: apply
    mblk <- readBlock(meth, cols); ublk <- readBlock(unmeth, cols)
    betas <- matrix(0, nprobe, length(cols))
    for (k in seq_along(cols)) {
      j <- cols[[k]]
      t0 <- tic()
      mc <- mblk[, k]; mc[isI] <- mc[isI] - offM[j]
      uc <- ublk[, k]; uc[isI] <- uc[isI] - offU[j]
      tm$apply <- tm$apply + (tic() - t0)
      mc[isI] <- qmap(mc[isI], refMI); mc[isII] <- qmap(mc[isII], refMII)
      uc[isI] <- qmap(uc[isI], refUI); uc[isII] <- qmap(uc[isII], refUII)
      t0 <- tic(); betas[, k] <- mc / (mc + uc + fudge); tm$apply <- tm$apply + (tic() - t0)
    }
    t0 <- tic(); append.gdsn(out, betas); tm$write <- tm$write + (tic() - t0)
  }
  readmode.gdsn(out)
  delete.gdsn(index.gdsn(g, "normbetas_prof"))
  as.list(tm)
}

cat(sprintf("gdsDasen profile | 450k | block=%d | build workers=%d\n", block, workers))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))
suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  cls <- tempfile("prof_gds_", tmpdir = scratch)
  message("  building GDS")
  suppressMessages(analyze(
    dataDirectory = coh$dir, samplesheet = file.path(coh$dir, "samplesheet_bench.csv"),
    crossReactiveProbes = coh$xr, gdsOutput = paste0(cls, ".gds"),
    annotationPackage = annoPkg, rebuild = "all", readerWorkers = workers,
    compress = "", verbose = 0L, analysis = NULL))
  gp <- paste0(cls, ".gds")

  message("  profiling gdsDasen")
  t0 <- proc.time()[["elapsed"]]
  tm <- profGdsDasen(gp, blk = block)
  total <- proc.time()[["elapsed"]] - t0

  unlink(coh$dir, recursive = TRUE); unlink(c(gp, paste0(gp, ".buildkey.rds")))
  gc(FALSE)
  compute <- tm$dfs2 + tm$sort + tm$qmap + tm$apply
  io <- tm$read + tm$write
  rows[[length(rows) + 1L]] <- data.frame(
    N = n, total_s = round(total, 1),
    read_s = round(tm$read, 1), dfs2_s = round(tm$dfs2, 1),
    sort_s = round(tm$sort, 1), qmap_s = round(tm$qmap, 1),
    apply_s = round(tm$apply, 1), write_s = round(tm$write, 1),
    compute_s = round(compute, 1), io_s = round(io, 1),
    compute_pct = round(100 * compute / (compute + io)))
  message(sprintf("  -> N=%d  total=%.1fs  compute=%.1fs io=%.1fs (compute %d%%)  [read=%.1f dfs2=%.1f sort=%.1f qmap=%.1f]",
                  n, total, compute, io, round(100 * compute / (compute + io)),
                  tm$read, tm$dfs2, tm$sort, tm$qmap))
}
cat("\n=== gdsDasen profile: per-operation seconds (compute = dfs2+sort+qmap+apply, io = read+write) ===\n\n")
print(do.call(rbind, rows), row.names = FALSE)
