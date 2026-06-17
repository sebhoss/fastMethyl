# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Tests for the centralised error helpers (R/errors.R) and analyze()'s phase
# wrapping. The contract: an unexpected failure anywhere in the pipeline
# surfaces with what analyze() was doing, the underlying cause, a remedy to try,
# and the bug-report URL; an actionable fastMethyl_user_error passes through
# verbatim, never relabelled as a bug. All of this is exercised without IDATs,
# minfi, gdsfmt, or bigmelon by injecting analyze()'s heavy operations.

# ---- the helpers -----------------------------------------------------

test_that(".userStop raises a classed, actionable condition", {
  err <- tryCatch(.userStop("fix your ", "input"), error = identity)
  expect_s3_class(err, "fastMethyl_user_error")
  expect_s3_class(err, "error")
  expect_identical(conditionMessage(err), "fix your input")
})

test_that(".check passes a true condition and fails a false one as a user error", {
  expect_true(.check(TRUE, "unused message"))
  err <- tryCatch(.check(1 > 2, "must be greater"), error = identity)
  expect_s3_class(err, "fastMethyl_user_error")
  expect_match(conditionMessage(err), "must be greater")
})

test_that(".withPhase returns the value on success", {
  expect_identical(.withPhase("doing x", 1 + 1), 2)
})

test_that(".withPhase wraps an unexpected error with phase, cause, hint, bug URL", {
  err <- tryCatch(
    .withPhase("opening the file", stop("disk gone"), hint = "try rebuilding"),
    error = identity)
  expect_s3_class(err, "fastMethyl_error")
  msg <- conditionMessage(err)
  expect_match(msg, "while opening the file")
  expect_match(msg, "disk gone")
  expect_match(msg, "try rebuilding")
  expect_match(msg, "github.com/sebhoss/fastMethyl/issues", fixed = TRUE)
  expect_match(msg, "sessionInfo", fixed = TRUE)
})

test_that(".withPhase relays a fastMethyl_user_error unchanged (no bug footer)", {
  err <- tryCatch(
    .withPhase("opening the file", .userStop("you passed a bad path")),
    error = identity)
  expect_s3_class(err, "fastMethyl_user_error")
  expect_identical(conditionMessage(err), "you passed a bad path")
  expect_false(grepl("issues", conditionMessage(err)))
})

test_that(".withPhase adds context exactly once (nested phases do not double-wrap)", {
  inner <- function() .withPhase("the inner step", stop("root cause"))
  err <- tryCatch(.withPhase("the outer step", inner()), error = identity)
  msg <- conditionMessage(err)
  expect_match(msg, "the inner step")
  # The outer wrapper relays the already-wrapped fastMethyl_error rather than
  # wrapping it again, so the outer phase name and a second footer never appear.
  expect_false(grepl("the outer step", msg))
  expect_equal(length(gregexpr("issues", msg)[[1L]]), 1L)
})

test_that(".phaseMessage omits the hint block when no hint is given", {
  msg <- .phaseMessage("doing y", "the cause")
  expect_match(msg, "doing y")
  expect_match(msg, "the cause")
  expect_false(grepl("What you can try", msg))
})

# ---- analyze() phase wrapping (injected deps) ------------------------

test_that("analyze wraps a preprocessing failure with context and the bug URL", {
  err <- tryCatch(
    analyze(dataDirectory = "D", annotationPackage = "A",
            analysis = function(gds, res) "never reached",
            .preprocess = function(...) stop("minfi blew up"),
            .open = function(p) "GDS", .close = function(g) invisible()),
    error = identity)
  expect_s3_class(err, "fastMethyl_error")
  msg <- conditionMessage(err)
  expect_match(msg, "building the QC'd GDS")
  expect_match(msg, "minfi blew up")
  expect_match(msg, "issues", fixed = TRUE)
})

test_that("analyze relays an actionable preprocessing user error unchanged", {
  err <- tryCatch(
    analyze(dataDirectory = "D", annotationPackage = "A",
            analysis = function(gds, res) "never reached",
            .preprocess = function(...) .userStop("samplesheet missing Basename"),
            .open = function(p) "GDS", .close = function(g) invisible()),
    error = identity)
  expect_s3_class(err, "fastMethyl_user_error")
  expect_identical(conditionMessage(err), "samplesheet missing Basename")
})

test_that("analyze wraps an analysis-function failure and still closes the GDS", {
  closed <- 0L
  err <- tryCatch(
    analyze(dataDirectory = "D", annotationPackage = "A",
            analysis = function(gds, res) stop("kaboom in FUN"),
            .preprocess = function(...) list(gds_path = "/tmp/x.gds",
                                             targets = "T", rebuilt = FALSE),
            .open = function(p) "GDS",
            .close = function(g) closed <<- closed + 1L),
    error = identity)
  expect_s3_class(err, "fastMethyl_error")
  msg <- conditionMessage(err)
  expect_match(msg, "running your analysis function")
  expect_match(msg, "kaboom in FUN")
  # The loan-pattern close runs on error too, so a failed analysis never leaks
  # the open handle.
  expect_equal(closed, 1L)
})

test_that("analyze wraps a GDS-open failure with the rebuild hint", {
  err <- tryCatch(
    analyze(dataDirectory = "D", annotationPackage = "A",
            analysis = function(gds, res) "never reached",
            .preprocess = function(...) list(gds_path = "/tmp/x.gds",
                                             targets = "T", rebuilt = FALSE),
            .open = function(p) stop("cannot open"),
            .close = function(g) invisible()),
    error = identity)
  expect_s3_class(err, "fastMethyl_error")
  msg <- conditionMessage(err)
  expect_match(msg, "opening the QC'd GDS file")
  expect_match(msg, "cannot open")
  expect_match(msg, "rebuild")
})
