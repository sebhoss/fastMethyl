# Pure tests for analyze() and its resource-management core .withGds(). The
# heavy operations are injectable, so the control flow -- argument validation,
# preprocess forwarding, optional normalisation, the FUN callback, the return
# value, and the close-even-on-error guarantee -- is exercised here without any
# IDATs, minfi, gdsfmt, or bigmelon.

# ---- .withGds (the loan-pattern core) --------------------------------

test_that(".withGds opens, runs fun, returns its value, and closes once", {
    log <- new.env()
    log$closed <- 0L
    out <- .withGds(
        "path/to.gds",
        function(gds) paste0("used:", gds),
        .open  = function(p) { log$opened <- p; "GDS" },
        .close = function(g) { log$closed <- log$closed + 1L; log$closedg <- g })
    expect_identical(out, "used:GDS")
    expect_identical(log$opened, "path/to.gds")
    expect_identical(log$closed, 1L)
    expect_identical(log$closedg, "GDS")
})

test_that(".withGds closes the GDS even when fun errors, and propagates", {
    closed <- 0L
    expect_error(
        .withGds(
            "p",
            function(gds) stop("boom in user code"),
            .open  = function(p) "GDS",
            .close = function(g) closed <<- closed + 1L),
        regexp = "boom in user code")
    expect_equal(closed, 1L)
})

# ---- .dasenNormalize guard -------------------------------------------

test_that(".dasenNormalize errors actionably when bigmelon is unavailable", {
    expect_error(
        .dasenNormalize(NULL, "normbetas", .have = function(p) FALSE),
        regexp = "bigmelon")
})

# ---- analyze() argument validation -----------------------------------

test_that("analyze rejects a non-function FUN", {
    expect_error(analyze(FUN = "not a function"), regexp = "FUN")
})

test_that("analyze rejects a non-scalar/NA normalize", {
    expect_error(analyze(FUN = identity, normalize = "yes"), regexp = "normalize")
    expect_error(analyze(FUN = identity, normalize = c(TRUE, FALSE)),
                 regexp = "normalize")
    expect_error(analyze(FUN = identity, normalize = NA), regexp = "normalize")
})

test_that("analyze rejects an invalid normNode", {
    expect_error(analyze(FUN = identity, normNode = ""), regexp = "normNode")
    expect_error(analyze(FUN = identity, normNode = c("a", "b")),
                 regexp = "normNode")
})

# ---- analyze() control flow (all hooks injected) ---------------------

# A recorder set of hooks; returns the hooks plus the shared log env.
.analyze_hooks <- function(gds_path = "/tmp/x.gds") {
    log <- new.env()
    log$closed <- 0L
    log$normalized <- NULL
    res <- list(gds_path = gds_path, targets = "TARGETS", keepSamples = TRUE)
    list(
        log = log, res = res,
        pre  = function(...) { log$ppargs <- list(...); res },
        open = function(p) { log$opened <- p; "GDS" },
        clos = function(g) { log$closed <- log$closed + 1L },
        norm = function(gds, node) { log$normalized <- node; invisible(node) })
}

test_that("analyze forwards ... to preprocess, normalises, calls FUN(gds,res), returns FUN value, closes", {
    h <- .analyze_hooks()
    out <- analyze(
        dataDirectory = "D", annotationPackage = "A",
        FUN = function(gds, res) list(g = gds, r = res),
        .preprocess = h$pre, .open = h$open, .close = h$clos,
        .normalize = h$norm)
    expect_identical(h$log$ppargs,
                     list(dataDirectory = "D", annotationPackage = "A"))
    expect_identical(h$log$opened, "/tmp/x.gds")
    expect_identical(h$log$normalized, "normbetas")     # default normNode
    expect_identical(out$g, "GDS")
    expect_identical(out$r, h$res)
    expect_identical(h$log$closed, 1L)
})

test_that("analyze skips normalisation when normalize = FALSE", {
    h <- .analyze_hooks()
    analyze(FUN = function(gds, res) NULL, normalize = FALSE,
            .preprocess = h$pre, .open = h$open, .close = h$clos,
            .normalize = h$norm)
    expect_null(h$log$normalized)
    expect_identical(h$log$closed, 1L)
})

test_that("analyze passes a custom normNode to the normaliser", {
    h <- .analyze_hooks()
    analyze(FUN = function(gds, res) NULL, normNode = "myNode",
            .preprocess = h$pre, .open = h$open, .close = h$clos,
            .normalize = h$norm)
    expect_identical(h$log$normalized, "myNode")
})

test_that("analyze closes the GDS even when FUN errors", {
    h <- .analyze_hooks()
    expect_error(
        analyze(FUN = function(gds, res) stop("user code failed"),
                .preprocess = h$pre, .open = h$open, .close = h$clos,
                .normalize = h$norm),
        regexp = "user code failed")
    expect_identical(h$log$closed, 1L)
})

test_that("analyze closes the GDS even when normalisation errors", {
    h <- .analyze_hooks()
    expect_error(
        analyze(FUN = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos,
                .normalize = function(gds, node) stop("dasen failed")),
        regexp = "dasen failed")
    expect_identical(h$log$closed, 1L)
})

# ---- analyze() verbose diagnostics -----------------------------------

test_that("analyze(verbose = 2L) logs each external-data phase", {
    h <- .analyze_hooks()
    msgs <- character(0)
    withCallingHandlers(
        analyze(dataDirectory = "D", verbose = 2L,
                FUN = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos,
                .normalize = h$norm),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl("\\[analyze\\] starting", msgs)))
    expect_true(any(grepl("preprocessing complete", msgs)))
    expect_true(any(grepl("opening GDS", msgs)))
    expect_true(any(grepl("normalising node 'normbetas'", msgs)))
    expect_true(any(grepl("analysis function returned", msgs)))
    # verbose still reached preprocess (it is read from ..., not consumed).
    expect_identical(h$log$ppargs$verbose, 2L)
})

test_that("analyze is silent at the default verbose level", {
    h <- .analyze_hooks()
    expect_silent(
        analyze(dataDirectory = "D",
                FUN = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos,
                .normalize = h$norm))
})
