# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Read-side adapters: project a finished QC'd GDS into the container shape a
# downstream package expects. The GDS that analyze() builds is the canonical,
# low-memory store -- it is written by streaming so a cohort-sized matrix never
# lives in RAM. These adapters are the opposite direction: they materialise a
# representation on demand, sized to the caller's needs, for tools that cannot
# read a GDS directly.
#
# Two consumers are covered:
#   * base-R / matrix tools (FactoMineR, stats::prcomp, irlba, umap, Rtsne) want a
#     dense `individuals x variables` matrix -- gdsBetaMatrix() with the default
#     transpose = TRUE delivers samples-as-rows;
#   * the SummarizedExperiment ecosystem (limma, sva) wants an SE whose assays are
#     `features x samples` -- gdsSummarizedExperiment() builds one with rowData
#     from fData and colData from pData.
#
# GDS-native analysis (bigmelon's dasen / outlyx / prcomp / estimateCellCounts.gds)
# needs none of this -- it operates on the open handle directly. An adapter is
# only worth its memory cost when the destination package is not GDS-aware; that
# trade-off (materialise the matrix vs. stream over the GDS) is the caller's to
# make, which is why these are explicit calls rather than an analyze() output mode.

# Resolve `gds` to an open handle plus a matching close action. A gds.class is
# used as-is and left open (its owner -- typically analyze()'s loan pattern --
# closes it); a path is opened here (read-only by default; read-write for
# callers that add a node) and closed when the caller returns, so a standalone
# call (e.g. on analyze(FUN = NULL)$gds_path) is self-contained and never leaves
# a dangling handle.
.gdsHandle <- function(gds, readonly = TRUE) {
  if (methods::is(gds, "gds.class")) {
    list(gds = gds, close = function() invisible(NULL))
  } else if (is.character(gds) && length(gds) == 1L && !is.na(gds) &&
               nzchar(gds)) {
    if (!file.exists(gds)) {
      stop(sprintf("GDS file \"%s\" does not exist.", gds), call. = FALSE)
    }
    handle <- gdsfmt::openfn.gds(gds, readonly = readonly)
    list(gds = handle, close = function() gdsfmt::closefn.gds(handle))
  } else {
    stop("`gds` must be an open gds.class handle or a single GDS file path.",
         call. = FALSE)
  }
}

# Probe (row) and sample (column) names for the matrix nodes. bigmelon-style GDS
# files store a two-element `paths` node naming the fData/pData columns that hold
# the canonical row/column identifiers; honour it when present so the names match
# what bigmelon itself would report, and fall back to the conventional
# fData/Name + pData/barcode otherwise. Either lookup degrades to NULL (unnamed
# dimension) rather than erroring, because a name vector of the wrong length is
# worse than none.
.gdsDimNames <- function(gds) {
  nodes <- gdsfmt::ls.gdsn(gds)
  paths <- if ("paths" %in% nodes) {
    tryCatch(as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "paths"))),
             error = function(e) NULL)
  }
  rowPath <- if (length(paths) >= 1L && nzchar(paths[[1L]])) {
    paths[[1L]]
  } else {
    "fData/Name"
  }
  colPath <- if (length(paths) >= 2L && nzchar(paths[[2L]])) {
    paths[[2L]]
  } else {
    "pData/barcode"
  }
  readVec <- function(path) {
    tryCatch(
      as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, path))),
      error = function(e) NULL
    )
  }
  list(probes = readVec(rowPath), samples = readVec(colPath))
}

# Attach probe/sample names to a probes x samples matrix when the name vectors
# line up with its dimensions. A length mismatch (a hand-edited GDS, a node that
# was compacted without its names) silently leaves that dimension unnamed.
.applyGdsDimNames <- function(m, dn) {
  if (!is.null(dn$probes) && length(dn$probes) == nrow(m)) {
    rownames(m) <- dn$probes
  }
  if (!is.null(dn$samples) && length(dn$samples) == ncol(m)) {
    colnames(m) <- dn$samples
  }
  m
}

# Exported -----------------------------------------------------------------

# Extract a 2-D GDS node (a beta / methylated / unmethylated / p-value matrix, or
# a normalised node such as `normbetas`) as a base R matrix.
#
# The GDS stores these nodes as probes x samples. transpose = TRUE (the default)
# returns samples x probes -- the individuals x variables orientation FactoMineR,
# stats::prcomp, irlba, umap and Rtsne all expect -- so the typical call feeds a
# matrix tool directly. transpose = FALSE keeps the on-disk probes x samples
# layout for code that wants features in rows.
#
# `gds` is either the open handle analyze()'s FUN receives or a path to a GDS on
# disk; a path is opened read-only and closed before returning.
gdsBetaMatrix <- function(gds, node = "betas", transpose = TRUE) {
  stopifnot(
    "`node` must be a single non-empty string" =
      is.character(node) && length(node) == 1L && !is.na(node) && nzchar(node),
    "`transpose` must be a single TRUE or FALSE" =
      is.logical(transpose) && length(transpose) == 1L && !is.na(transpose)
  )
  handle <- .gdsHandle(gds)
  on.exit(handle$close(), add = TRUE)
  g <- handle$gds

  if (!(node %in% gdsfmt::ls.gdsn(g))) {
    stop(sprintf(
      "node \"%s\" is not present in the GDS (available: %s).",
      node, paste(gdsfmt::ls.gdsn(g), collapse = ", ")
    ), call. = FALSE)
  }
  m <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, node))
  if (!is.matrix(m) || length(dim(m)) != 2L) {
    stop(sprintf("node \"%s\" is not a 2-D matrix.", node), call. = FALSE)
  }
  m <- .applyGdsDimNames(m, .gdsDimNames(g))
  if (transpose) t(m) else m
}

# Build a SummarizedExperiment from the GDS matrix nodes, with rowData taken from
# fData (probe annotation) and colData from pData (sample annotation). The assays
# keep the GDS probes x samples orientation SE expects, so no transpose happens.
#
# `assays` names the nodes to load; any not present in the file are dropped (so a
# GDS that omits, say, `pvals`, or one carrying a `normbetas` node you ask for,
# both work). rowData / colData are attached only when their row count matches the
# assay dimensions, so a GDS without fData/pData still yields a valid SE.
gdsSummarizedExperiment <- function(
    gds,
    assays = c("betas", "methylated", "unmethylated", "pvals")) {
  stopifnot(
    "`assays` must be a non-empty character vector" =
      is.character(assays) && length(assays) >= 1L && all(nzchar(assays))
  )
  handle <- .gdsHandle(gds)
  on.exit(handle$close(), add = TRUE)
  g <- handle$gds

  present <- gdsfmt::ls.gdsn(g)
  wanted <- assays[assays %in% present]
  if (length(wanted) == 0L) {
    stop(sprintf(
      "none of the requested assay nodes (%s) are present in the GDS (available: %s).",
      paste(assays, collapse = ", "), paste(present, collapse = ", ")
    ), call. = FALSE)
  }

  dn <- .gdsDimNames(g)
  mats <- lapply(wanted, function(nd) {
    .applyGdsDimNames(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, nd)), dn)
  })
  names(mats) <- wanted

  se <- SummarizedExperiment::SummarizedExperiment(assays = mats)

  # rowData/colData carry the GDS dim identifiers as their own row names, so
  # assigning them does not blank the SE's dimnames (a DataFrame built without
  # row.names would reset rownames(se) / colnames(se) to NULL).
  fData <- if ("fData" %in% present) {
    tryCatch(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "fData")),
             error = function(e) NULL)
  }
  if (is.data.frame(fData) && nrow(fData) == nrow(se)) {
    SummarizedExperiment::rowData(se) <-
      S4Vectors::DataFrame(fData, row.names = rownames(se))
  }
  pData <- if ("pData" %in% present) {
    tryCatch(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "pData")),
             error = function(e) NULL)
  }
  if (is.data.frame(pData) && nrow(pData) == ncol(se)) {
    SummarizedExperiment::colData(se) <-
      S4Vectors::DataFrame(pData, row.names = colnames(se))
  }
  se
}
