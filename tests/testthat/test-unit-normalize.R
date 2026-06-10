# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# gdsDasen() equivalence + behaviour. The contract is that the streaming, two-
# pass normaliser reproduces wateRmelon::dasen (per-sample dfsfit, roco = NULL)
# bit-for-bit -- so the oracle here is wateRmelon::dasen on the same intensities.
# A small synthetic GDS (methylated/unmethylated/fData$Type) is built directly
# with gdsfmt, so no IDAT read or minfi reader is involved; only the equivalence
# of the maths is under test. wateRmelon is a Suggests, so the file skips when it
# (or gdsfmt) is absent.

# Build a tiny bigmelon-shaped GDS with realistic-ish intensities: Type I and
# Type II probes drawn from offset distributions so dfsfit has a real mode gap to
# correct, positive and within dfs2's [0, 5000] density window.
.makeIntensityGds <- function(path, nProbeEach = 800L, nSample = 8L, seed = 1L) {
  set.seed(seed)
  nprobe <- 2L * nProbeEach
  onetwo <- rep(c("I", "II"), each = nProbeEach)
  draw <- function(shift) {
    m <- matrix(
      stats::rgamma(nprobe * nSample, shape = 3, scale = 300) + 100,
      nrow = nprobe, ncol = nSample
    )
    m[onetwo == "I", ] <- m[onetwo == "I", ] + shift
    m
  }
  meth <- draw(150)
  unmeth <- draw(120)
  gf <- gdsfmt::createfn.gds(path)
  on.exit(gdsfmt::closefn.gds(gf), add = TRUE)
  gdsfmt::add.gdsn(gf, "methylated", val = meth, closezip = TRUE)
  gdsfmt::add.gdsn(gf, "unmethylated", val = unmeth, closezip = TRUE)
  fnode <- gdsfmt::addfolder.gdsn(gf, "fData")
  gdsfmt::put.attr.gdsn(fnode, "R.class", "data.frame")
  gdsfmt::add.gdsn(fnode, "Name", val = sprintf("cg%08d", seq_len(nprobe)),
                   check = FALSE)
  gdsfmt::add.gdsn(fnode, "Type", val = onetwo, check = FALSE)
  list(meth = meth, unmeth = unmeth, onetwo = onetwo)
}

test_that("gdsDasen reproduces wateRmelon::dasen bit-for-bit (keep = all)", {
  skip_if_not_installed("gdsfmt")
  skip_if_not_installed("wateRmelon")
  path <- tempfile(fileext = ".gds")
  dat <- .makeIntensityGds(path)
  on.exit(unlink(path), add = TRUE)

  gdsDasen(path, node = "normbetas")

  g <- gdsfmt::openfn.gds(path, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  got <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "normbetas"))
  # Oracle on the exact intensities the GDS holds, per-sample dfsfit (roco=NULL).
  ref <- wateRmelon::dasen(dat$meth, dat$unmeth, dat$onetwo,
                           fudge = 100, roco = NULL)

  expect_equal(dim(got), dim(ref))
  expect_equal(got, unname(ref), tolerance = 1e-10)
})

test_that("gdsDasen writes normbetas for every sample even with a keep mask", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  dat <- .makeIntensityGds(path, nSample = 6L)
  on.exit(unlink(path), add = TRUE)

  keep <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)  # drop sample 6 from the reference
  gdsDasen(path, node = "normbetas", keep = keep)

  g <- gdsfmt::openfn.gds(path, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  got <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "normbetas"))
  # normbetas stays rectangular: a column for every sample, not just survivors.
  expect_equal(ncol(got), 6L)
  expect_equal(nrow(got), 1600L)
  expect_true(all(is.finite(got)))
})

test_that("excluding a sample from the reference changes the result", {
  skip_if_not_installed("gdsfmt")
  skip_if_not_installed("wateRmelon")
  path <- tempfile(fileext = ".gds")
  .makeIntensityGds(path, nSample = 6L)
  on.exit(unlink(path), add = TRUE)

  gdsDasen(path, node = "nb_all")                                   # keep all
  gdsDasen(path, node = "nb_sub",
           keep = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))           # drop one

  g <- gdsfmt::openfn.gds(path, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  a <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "nb_all"))
  b <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "nb_sub"))
  # A different reference set yields a different (still valid) normalisation.
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("the built-in dasen hook caches: skips unless rebuildDownstream", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  .makeIntensityGds(path, nSample = 6L)
  on.exit(unlink(path), add = TRUE)

  hook <- .resolveNormalize("dasen")
  g <- gdsfmt::openfn.gds(path, readonly = FALSE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  dims <- gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(g, "methylated"))$dim
  # Seed a sentinel normbetas distinct from any real dasen output (all -1).
  gdsfmt::add.gdsn(g, "normbetas", val = matrix(-1, dims[1L], dims[2L]))

  # Node exists and nothing rebuilt upstream -> skip, sentinel preserved.
  hook(g, list(rebuildDownstream = FALSE), NULL)
  kept <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "normbetas"))
  expect_true(all(kept == -1))

  # rebuildDownstream = TRUE -> recompute, sentinel overwritten by real betas.
  hook(g, list(rebuildDownstream = TRUE), NULL)
  done <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "normbetas"))
  expect_false(any(done == -1))
  expect_true(all(done >= 0 & done <= 1))
})

test_that("gdsDasen validates arguments and required nodes", {
  skip_if_not_installed("gdsfmt")
  path <- tempfile(fileext = ".gds")
  .makeIntensityGds(path)
  on.exit(unlink(path), add = TRUE)

  expect_error(gdsDasen(path, node = ""), "non-empty string")
  expect_error(gdsDasen(path, fudge = -1), "non-negative")
  expect_error(gdsDasen(path, keep = c(TRUE, FALSE)), "one entry per sample")

  empty <- tempfile(fileext = ".gds")
  gf <- gdsfmt::createfn.gds(empty); gdsfmt::closefn.gds(gf)
  on.exit(unlink(empty), add = TRUE)
  expect_error(gdsDasen(empty), "methylated")
})
