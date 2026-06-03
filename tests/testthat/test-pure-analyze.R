# Pure tests for analyze() and its resource-management core .withGds(). The
# heavy operations are injectable, so the control flow -- argument validation,
# preprocess forwarding, the optional FUN callback, the return value, and the
# close-even-on-error guarantee -- is exercised here without any IDATs, minfi,
# gdsfmt, or bigmelon.

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

# ---- analyze() argument validation -----------------------------------

test_that("analyze rejects a non-NULL FUN that is not a function", {
    expect_error(analyze(FUN = "not a function"), regexp = "FUN")
})

# ---- analyze() control flow (all hooks injected) ---------------------

# A recorder set of hooks; returns the hooks plus the shared log env.
.analyze_hooks <- function(gds_path = "/tmp/x.gds") {
    log <- new.env()
    log$closed <- 0L
    res <- list(gds_path = gds_path, targets = "TARGETS", keepSamples = TRUE)
    list(
        log = log, res = res,
        pre  = function(...) { log$ppargs <- list(...); res },
        open = function(p) { log$opened <- p; "GDS" },
        clos = function(g) { log$closed <- log$closed + 1L })
}

test_that("analyze forwards ... to preprocess, calls FUN(gds, res), returns its value, and closes", {
    h <- .analyze_hooks()
    out <- analyze(
        dataDirectory = "D", annotationPackage = "A",
        FUN = function(gds, res) list(g = gds, r = res),
        .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(h$log$ppargs,
                     list(dataDirectory = "D", annotationPackage = "A"))
    expect_identical(h$log$opened, "/tmp/x.gds")
    expect_identical(out$g, "GDS")
    expect_identical(out$r, h$res)
    expect_identical(h$log$closed, 1L)
})

test_that("analyze closes the GDS even when FUN errors", {
    h <- .analyze_hooks()
    expect_error(
        analyze(FUN = function(gds, res) stop("user code failed"),
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "user code failed")
    expect_identical(h$log$closed, 1L)
})

# ---- analyze() without a closure (just build the GDS) ----------------

test_that("analyze without FUN returns the preprocess result without opening the GDS", {
    h <- .analyze_hooks()
    out <- analyze(dataDirectory = "D",
                   .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(out, h$res)       # preprocessing result, not a FUN value
    expect_null(h$log$opened)          # preprocessing already built+closed it
    expect_identical(h$log$closed, 0L)
})

# ---- analyze() verbose diagnostics -----------------------------------

test_that("analyze(verbose = 2L) logs each external-data phase", {
    h <- .analyze_hooks()
    msgs <- character(0)
    withCallingHandlers(
        analyze(dataDirectory = "D", verbose = 2L,
                FUN = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl("\\[analyze\\] starting", msgs)))
    expect_true(any(grepl("preprocessing complete", msgs)))
    expect_true(any(grepl("opening GDS", msgs)))
    expect_true(any(grepl("analysis function returned", msgs)))
    # verbose still reached preprocess (it is read from ..., not consumed).
    expect_identical(h$log$ppargs$verbose, 2L)
})

test_that("analyze is silent at the default verbose level", {
    h <- .analyze_hooks()
    expect_silent(
        analyze(dataDirectory = "D",
                FUN = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos))
})
