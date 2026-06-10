# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Adapter tests: gdsBetaMatrix() / gdsSummarizedExperiment() project a finished
# GDS into a base matrix and a SummarizedExperiment. They touch only gdsfmt (+
# SummarizedExperiment), never minfi, so they live in the pure bucket. A small
# bigmelon-shaped GDS is built by hand -- probes x samples matrix nodes plus
# fData / pData / paths -- so the orientation, naming, and metadata round-trip
# can be asserted without reading any IDAT.

# Write a minimal bigmelon-style GDS: each matrix node is probes (rows) x samples
# (cols), fData is the probe annotation, pData carries the sample barcodes, and
# `paths` names the fData/pData columns holding the canonical dim identifiers.
.makeTinyGds <- function(path, nProbe = 5L, nSample = 3L) {
  probes <- sprintf("cg%08d", seq_len(nProbe))
  samples <- sprintf("S%03d", seq_len(nSample))
  mat <- function(seed) {
    matrix(seed + seq_len(nProbe * nSample), nrow = nProbe, ncol = nSample)
  }
  gf <- gdsfmt::createfn.gds(path)
  on.exit(gdsfmt::closefn.gds(gf), add = TRUE)
  gdsfmt::add.gdsn(gf, "betas", val = mat(0) / 100, closezip = TRUE)
  gdsfmt::add.gdsn(gf, "methylated", val = mat(1000L), closezip = TRUE)
  gdsfmt::add.gdsn(gf, "unmethylated", val = mat(2000L), closezip = TRUE)
  gdsfmt::add.gdsn(gf, "pvals", val = mat(0) / 1000, closezip = TRUE)

  fnode <- gdsfmt::addfolder.gdsn(gf, "fData")
  gdsfmt::put.attr.gdsn(fnode, "R.class", "data.frame")
  gdsfmt::add.gdsn(fnode, "Name", val = probes, check = FALSE)
  gdsfmt::add.gdsn(fnode, "chr",
                   val = rep(c("chr1", "chr2"), length.out = nProbe),
                   check = FALSE)

  gdsfmt::add.gdsn(gf, "pData",
                   val = data.frame(barcode = samples,
                                    Sample_Group = rep(c("A", "B"),
                                                       length.out = nSample),
                                    stringsAsFactors = FALSE))
  gdsfmt::add.gdsn(gf, "paths", val = c("fData/Name", "pData/barcode"))
  list(probes = probes, samples = samples)
}

test_that("gdsBetaMatrix returns samples x probes with names by default", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  meta <- .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  m <- gdsBetaMatrix(path)
  expect_equal(dim(m), c(3L, 5L))                  # samples x probes
  expect_identical(rownames(m), meta$samples)
  expect_identical(colnames(m), meta$probes)
})

test_that("gdsBetaMatrix transpose = FALSE keeps the on-disk probes x samples", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  meta <- .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  m <- gdsBetaMatrix(path, transpose = FALSE)
  expect_equal(dim(m), c(5L, 3L))                  # probes x samples
  expect_identical(rownames(m), meta$probes)
  expect_identical(colnames(m), meta$samples)
  # transpose = TRUE is the exact transpose of transpose = FALSE.
  expect_identical(gdsBetaMatrix(path), t(m))
})

test_that("gdsBetaMatrix can read a non-default node", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  expect_equal(dim(gdsBetaMatrix(path, node = "methylated")), c(3L, 5L))
})

test_that("gdsBetaMatrix errors on a missing node, naming what is available", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  expect_error(gdsBetaMatrix(path, node = "normbetas"), "not present")
})

test_that("gdsBetaMatrix accepts an open handle and does not close it", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  g <- gdsfmt::openfn.gds(path, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  m <- gdsBetaMatrix(g)
  expect_equal(dim(m), c(3L, 5L))
  # The handle the caller owns is still usable -- the adapter left it open.
  expect_true("betas" %in% gdsfmt::ls.gdsn(g))
})

test_that("gdsBetaMatrix rejects a non-existent path and a bad argument", {
  expect_error(gdsBetaMatrix("/no/such/file.gds"), "does not exist")
  expect_error(gdsBetaMatrix(42), "open gds.class handle or a single GDS file path")
})

test_that("gdsSummarizedExperiment carries assays, rowData and colData", {
  skip_if_not_installed("gdsfmt")
  skip_if_not_installed("SummarizedExperiment")
  path <- tempfile(fileext = ".gds")
  meta <- .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  se <- gdsSummarizedExperiment(path)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(5L, 3L))                 # features x samples
  expect_identical(
    sort(SummarizedExperiment::assayNames(se)),
    sort(c("betas", "methylated", "unmethylated", "pvals"))
  )
  expect_identical(rownames(se), meta$probes)
  expect_identical(colnames(se), meta$samples)
  expect_true("chr" %in% names(SummarizedExperiment::rowData(se)))
  expect_true("Sample_Group" %in% names(SummarizedExperiment::colData(se)))
})

test_that("gdsSummarizedExperiment drops absent nodes and errors when none remain", {
  skip_if_not_installed("gdsfmt")
  skip_if_not_installed("SummarizedExperiment")
  path <- tempfile(fileext = ".gds")
  .makeTinyGds(path)
  on.exit(unlink(path), add = TRUE)

  # `normbetas` is absent; only the two present nodes survive.
  se <- gdsSummarizedExperiment(path, assays = c("betas", "normbetas", "pvals"))
  expect_identical(sort(SummarizedExperiment::assayNames(se)),
                   sort(c("betas", "pvals")))
  expect_error(gdsSummarizedExperiment(path, assays = "normbetas"),
               "none of the requested assay nodes")
})
