# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# gdsDasen(): a streaming, two-pass between-array quantile normalisation that
# reproduces wateRmelon::dasen (per-sample dfsfit, roco = NULL) bit-for-bit while
# reading the GDS one column block at a time -- it never assembles the full
# samples x probes matrix and writes no scratch GDS (bigmelon's dasenrank spills a
# temp.gds for the dfsfit'd intensities; we do dfsfit in RAM per column instead).
#
# Why a between-array method can stream at all: dasen's quantile-normalisation
# reference is rowMeans of the sorted columns -- the mean of order statistics --
# which is accumulable one sample at a time. So:
#   pass 1  read each column, dfsfit it (caching the per-sample dfs offset so
#           pass 2 never recomputes the density), and add its sorted Type I / II
#           values into the four running references (methylated/unmethylated x
#           Type I/II), accumulating ONLY over the `keep` samples;
#   pass 2  re-read every column, dfsfit from the cached offset, map each
#           column's ranks onto the reference (rank + approx, limma's ties
#           branch), form beta = M / (M + U + fudge), and append.
#
# `keep` is the outlier mask: the reference is built from survivors only (the QN
# reference must exclude outliers), but normbetas is written for EVERY sample, so
# the GDS stays rectangular (normbetas matches the methylated/pData sample count)
# and a flagged-but-kept outlier still gets a normalised column -- just one
# normalised against an outlier-free reference. This is what lets the analyze()
# orchestration normalise correctly without a clean-GDS rebuild.
#
# It depends only on stats (density/approx) + base (sort/rank) + gdsfmt, so it
# adds no package dependency; the value-identity to wateRmelon::dasen is pinned by
# tests, not by calling wateRmelon.

# Per-sample dfs2 offset: the gap between the modes of the Type I and Type II
# intensity densities (verbatim from wateRmelon's dfs2, so the adjusted
# intensities match dasen's). roco-based sentrix-position smoothing is not
# applied -- this is the per-sample (roco = NULL) dasen, which is what the
# equivalence tests pin against.
.dfs2offset <- function(col, isI, isII) {
  one <- stats::density(col[isI], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  two <- stats::density(col[isII], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
  one$x[which.max(one$y)] - two$x[which.max(two$y)]
}

# Map a column's values onto a quantile reference exactly as
# limma::normalizeQuantiles does in its ties = TRUE branch: rank (averaging
# ties), rescale to [0, 1], and interpolate the reference at those positions.
.qmapToRef <- function(vals, ref) {
  refpos <- (seq_along(ref) - 1) / (length(ref) - 1)
  stats::approx(refpos, ref, (rank(vals) - 1) / (length(vals) - 1),
                ties = list("ordered", mean))$y
}

# Split seq_len(n) into up to k contiguous chunks of near-equal size (>= 1 each),
# the per-worker work units. Avoids a parallel/base dependency for the split.
.colChunks <- function(n, k) {
  if (n <= 0L) {
    return(list())
  }
  k <- max(1L, min(as.integer(k), n))
  starts <- floor(seq(0, n, length.out = k + 1L))
  lapply(seq_len(k), function(i) (starts[[i]] + 1L):starts[[i + 1L]])
}

# Exported -----------------------------------------------------------------

# Stream dasen normalisation into `node`, reproducing wateRmelon::dasen bit-for-
# bit. `gds` is the open read-write handle analyze()'s FUN receives or a path
# (opened read-write, closed on return). `keep` is a logical mask over samples
# selecting which build the quantile reference (default: all); normbetas is
# always written for every sample. `block` is the column-block size that bounds
# peak memory -- larger reads fewer times, smaller holds less at once.
#
# `BPPARAM` is the OPT-IN parallel backend. The default (NULL -> SerialParam)
# runs serial and lean -- the master holds one block at a time and nothing is
# forked. Pass a BiocParallel backend (e.g. MulticoreParam(4)) to fan the
# per-column compute (the cost is ~90% per-column: dfs2 densities + the rank/map)
# across workers: the master still does all GDS IO and forks workers that reach
# the in-master block via copy-on-write, so it is value-identical either way
# (the parallel reference is summed per worker then combined; differences are at
# floating-point precision). It trades ~2x speed for higher -- but block-bound,
# not cohort-bound -- peak memory (the COW of the forked workers), so it is opt-in
# rather than the default: the lean serial path keeps fastMethyl's memory
# guarantee, and parallelism is there when RAM is not the constraint.
gdsDasen <- function(gds, node = "normbetas", keep = NULL, fudge = 100,
                     block = 64L, BPPARAM = NULL) {
  .check(is.character(node) && length(node) == 1L && !is.na(node) && nzchar(node),
         "`node` must be a single non-empty string.")
  .check(is.numeric(fudge) && length(fudge) == 1L && !is.na(fudge) && fudge >= 0,
         "`fudge` must be a single non-negative number.")
  .check(is.numeric(block) && length(block) == 1L && !is.na(block) && block >= 1L,
         "`block` must be a single positive integer.")
  block <- as.integer(block)
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  } else if (!methods::is(BPPARAM, "BiocParallelParam")) {
    .userStop("`BPPARAM` must be a BiocParallelParam object (e.g. ",
              "BiocParallel::MulticoreParam()), or NULL for the serial default.")
  }
  handle <- .gdsHandle(gds, readonly = FALSE)
  on.exit(handle$close(), add = TRUE)
  g <- handle$gds

  present <- gdsfmt::ls.gdsn(g)
  for (req in c("methylated", "unmethylated", "fData/Type")) {
    base <- strsplit(req, "/", fixed = TRUE)[[1L]][[1L]]
    if (!(base %in% present)) {
      .userStopf("gdsDasen needs the `%s` node; the GDS has: %s.",
                 base, paste(present, collapse = ", "))
    }
  }
  meth <- gdsfmt::index.gdsn(g, "methylated")
  unmeth <- gdsfmt::index.gdsn(g, "unmethylated")
  onetwo <- as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "fData/Type")))
  dims <- gdsfmt::objdesp.gdsn(meth)$dim
  nprobe <- dims[1L]
  nsamp <- dims[2L]
  isI <- onetwo == "I"
  isII <- onetwo == "II"
  if (!any(isI) || !any(isII)) {
    .userStop("fData/Type must contain both \"I\" and \"II\" probes for dasen.")
  }

  if (is.null(keep)) {
    keep <- rep(TRUE, nsamp)
  }
  .check(is.logical(keep) && length(keep) == nsamp,
         "`keep` must be a logical vector with one entry per sample.")
  .check(any(keep & !is.na(keep)),
         "`keep` must select at least one sample for the reference.")
  keep[is.na(keep)] <- FALSE

  nI <- sum(isI)
  nII <- sum(isII)
  nworkers <- max(1L, bpnworkers(BPPARAM))
  refMI <- numeric(nI)
  refUI <- numeric(nI)
  refMII <- numeric(nII)
  refUII <- numeric(nII)
  offM <- numeric(nsamp)
  offU <- numeric(nsamp)

  readBlock <- function(node_h, cols) {
    matrix(as.numeric(gdsfmt::read.gdsn(
      node_h, start = c(1L, cols[1L]), count = c(-1L, length(cols))
    )), nrow = nprobe)
  }
  allCols <- seq_len(nsamp)
  blocks <- split(allCols, ceiling(allCols / block))

  # Pass 1: dfsfit every column (cache its offset), accumulate the references
  # over the `keep` samples only. The master reads each block; the per-column
  # work is fanned over BPPARAM workers (each handles a contiguous column chunk,
  # reaching the block via copy-on-write under a forking backend) and returns its
  # offsets plus partial reference sums, which the master reduces. With the
  # default SerialParam this is one in-process chunk == the plain serial loop.
  for (cols in blocks) {
    mblk <- readBlock(meth, cols)
    ublk <- readBlock(unmeth, cols)
    parts <- bplapply(.colChunks(length(cols), nworkers), function(idx) {
      oM <- numeric(length(idx))
      oU <- numeric(length(idx))
      pMI <- numeric(nI)
      pUI <- numeric(nI)
      pMII <- numeric(nII)
      pUII <- numeric(nII)
      for (kk in seq_along(idx)) {
        k <- idx[[kk]]
        aM <- .dfs2offset(mblk[, k], isI, isII)
        aU <- .dfs2offset(ublk[, k], isI, isII)
        oM[kk] <- aM
        oU[kk] <- aU
        if (keep[cols[[k]]]) {
          mc <- mblk[, k]
          mc[isI] <- mc[isI] - aM
          uc <- ublk[, k]
          uc[isI] <- uc[isI] - aU
          pMI <- pMI + sort(mc[isI])
          pUI <- pUI + sort(uc[isI])
          pMII <- pMII + sort(mc[isII])
          pUII <- pUII + sort(uc[isII])
        }
      }
      list(idx = idx, oM = oM, oU = oU,
           pMI = pMI, pUI = pUI, pMII = pMII, pUII = pUII)
    }, BPPARAM = BPPARAM)
    for (p in parts) {
      offM[cols[p$idx]] <- p$oM
      offU[cols[p$idx]] <- p$oU
      refMI <- refMI + p$pMI
      refUI <- refUI + p$pUI
      refMII <- refMII + p$pMII
      refUII <- refUII + p$pUII
    }
  }
  nk <- sum(keep)
  refMI <- refMI / nk
  refUI <- refUI / nk
  refMII <- refMII / nk
  refUII <- refUII / nk

  if (node %in% gdsfmt::ls.gdsn(g)) {
    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(g, node))
  }
  out <- gdsfmt::add.gdsn(g, node, valdim = c(nprobe, 0L), storage = "double")

  # Pass 2: dfsfit from the cached offset, map ranks onto the reference, form
  # beta. Workers compute their column chunk's betas; the master appends each
  # block in order. Written for EVERY sample so normbetas stays rectangular.
  for (cols in blocks) {
    mblk <- readBlock(meth, cols)
    ublk <- readBlock(unmeth, cols)
    parts <- bplapply(.colChunks(length(cols), nworkers), function(idx) {
      bc <- matrix(0, nrow = nprobe, ncol = length(idx))
      for (kk in seq_along(idx)) {
        k <- idx[[kk]]
        j <- cols[[k]]
        mc <- mblk[, k]
        mc[isI] <- mc[isI] - offM[j]
        uc <- ublk[, k]
        uc[isI] <- uc[isI] - offU[j]
        mc[isI] <- .qmapToRef(mc[isI], refMI)
        mc[isII] <- .qmapToRef(mc[isII], refMII)
        uc[isI] <- .qmapToRef(uc[isI], refUI)
        uc[isII] <- .qmapToRef(uc[isII], refUII)
        bc[, kk] <- mc / (mc + uc + fudge)
      }
      list(idx = idx, betas = bc)
    }, BPPARAM = BPPARAM)
    if (length(parts) == 1L) {
      betas <- parts[[1L]]$betas
    } else {
      betas <- matrix(0, nrow = nprobe, ncol = length(cols))
      for (p in parts) betas[, p$idx] <- p$betas
    }
    gdsfmt::append.gdsn(out, betas)
  }
  gdsfmt::readmode.gdsn(out)
  invisible(node)
}
