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
    expect_error(analyze(analysis = "not a function"), regexp = "analysis")
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
        analysis = function(gds, res) list(g = gds, r = res),
        .preprocess = h$pre, .open = h$open, .close = h$clos)
    # The user's named arguments reach .preprocess unchanged; analyze also
    # appends the resolved verbose (0L here -- unspecified means quiet) and the
    # build-stage forceRebuild translated from `rebuild` (FALSE by default).
    expect_identical(h$log$ppargs,
                     list(dataDirectory = "D", annotationPackage = "A",
                          verbose = 0L, forceRebuild = FALSE))
    expect_identical(h$log$opened, "/tmp/x.gds")
    expect_identical(out$g, "GDS")
    # FUN sees the preprocess result, augmented with the rebuild cascade flag
    # (FALSE here -- no rebuild requested and the mock reports no cache miss).
    expect_identical(out$r[names(h$res)], h$res)
    expect_false(out$r$rebuildDownstream)
    expect_identical(h$log$closed, 1L)
})

# ---- analyze() is the single authority on verbose --------------------
#
# The contract: an analyze() call that does not name verbose forwards 0L (quiet)
# to .preprocess rather than letting .preprocess apply its own standalone
# default, and a user-supplied level is normalised once (legacy logical ->
# integer, the same mapping analyze applies to its own diagnostics) so
# .preprocess sees exactly that normalised value.

test_that("analyze forwards verbose = 0L to preprocess when verbose is unspecified", {
    h <- .analyze_hooks()
    analyze(dataDirectory = "D",
            .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(h$log$ppargs$verbose, 0L)
})

test_that("analyze forwards an explicit integer verbose unchanged to preprocess", {
    for (lvl in c(0L, 1L, 2L)) {
        h <- .analyze_hooks()
        analyze(dataDirectory = "D", verbose = lvl,
                .preprocess = h$pre, .open = h$open, .close = h$clos)
        expect_identical(h$log$ppargs$verbose, lvl)
    }
})

test_that("analyze normalises a legacy logical verbose before forwarding it", {
    h <- .analyze_hooks()
    analyze(dataDirectory = "D", verbose = TRUE,
            .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(h$log$ppargs$verbose, 1L)

    h <- .analyze_hooks()
    analyze(dataDirectory = "D", verbose = FALSE,
            .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(h$log$ppargs$verbose, 0L)
})

test_that("analyze forwards the resolved verbose on the no-FUN path too", {
    h <- .analyze_hooks()
    analyze(dataDirectory = "D", verbose = 2L,
            .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(h$log$ppargs$verbose, 2L)
    expect_null(h$log$opened)          # no FUN: the GDS is never opened
})

test_that("analyze never passes verbose to preprocess twice", {
    # do.call would error with a duplicated formal if analyze both forwarded the
    # raw ... verbose and appended a resolved one; a clean run proves it writes
    # the normalised value back over the raw one instead of adding a second.
    h <- .analyze_hooks()
    # verbose = 0L so analyze emits no diagnostics of its own; the duplicate
    # check in .preprocess is independent of the level.
    expect_silent(
        analyze(dataDirectory = "D", verbose = 0L,
                analysis = function(gds, res) NULL,
                .preprocess = function(...) {
                    args <- list(...)
                    stopifnot(sum(names(args) == "verbose") == 1L)
                    h$res
                },
                .open = h$open, .close = h$clos))
})

test_that("analyze rejects an out-of-range or malformed verbose before preprocessing", {
    h <- .analyze_hooks()
    # .normalize_verbose runs in analyze, so a bad level fails fast and
    # .preprocess is never reached (no recorded args).
    expect_error(
        analyze(dataDirectory = "D", verbose = 3L,
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "verbose")
    expect_null(h$log$ppargs)

    expect_error(
        analyze(dataDirectory = "D", verbose = NA,
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "verbose")
    expect_error(
        analyze(dataDirectory = "D", verbose = c(1L, 2L),
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "verbose")
    expect_error(
        analyze(dataDirectory = "D", verbose = "loud",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "verbose")
})

test_that("analyze closes the GDS even when FUN errors", {
    h <- .analyze_hooks()
    expect_error(
        analyze(analysis = function(gds, res) stop("user code failed"),
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

# ---- analyze() hook orchestration ------------------------------------
#
# The contract: analyze runs build -> outlierRemoval -> normalize(keep) -> FUN.
# outlierRemoval(gds, res) returns a keep mask; that mask reaches normalize as its
# third argument AND FUN via res$outlierKeep; both hooks accept a method name, a
# function, or "none". The fake .open returns the string "GDS", which the hooks
# and FUN receive as `gds`, so the wiring is observable without a real GDS.

test_that("analyze runs outlierRemoval -> normalize(keep) -> FUN in order, threading keep", {
    h <- .analyze_hooks()
    order <- character(0)
    keepMask <- c(TRUE, FALSE, TRUE)
    out <- analyze(
        dataDirectory = "D",
        outlierRemoval = function(gds, res) { order <<- c(order, "outlier"); keepMask },
        normalization = function(gds, res, keep) {
            order <<- c(order, "normalize")
            h$log$normKeep <- keep; h$log$normGds <- gds
        },
        analysis = function(gds, res) { order <<- c(order, "fun"); res$outlierKeep },
        .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_identical(order, c("outlier", "normalize", "fun"))
    expect_identical(h$log$normKeep, keepMask)       # keep reaches normalize
    expect_identical(h$log$normGds, "GDS")           # on the open handle
    expect_identical(out, keepMask)                  # FUN saw res$outlierKeep
    expect_identical(h$log$closed, 1L)               # GDS closed once
})

test_that("analyze opens the GDS to run hooks even without FUN, and returns res", {
    h <- .analyze_hooks()
    out <- analyze(
        dataDirectory = "D",
        normalization = function(gds, res, keep) h$log$normRan <- TRUE,
        .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_true(isTRUE(h$log$normRan))               # hook ran...
    expect_identical(h$log$opened, "/tmp/x.gds")     # ...so the GDS was opened
    expect_identical(out$gds_path, "/tmp/x.gds")     # and res is returned
    expect_identical(h$log$closed, 1L)
})

test_that("analyze with both hooks \"none\" and no FUN never opens the GDS", {
    h <- .analyze_hooks()
    analyze(dataDirectory = "D", outlierRemoval = "none", normalization = "none",
            .preprocess = h$pre, .open = h$open, .close = h$clos)
    expect_null(h$log$opened)
    expect_identical(h$log$closed, 0L)
})

test_that("analyze rejects unknown hook method names before preprocessing", {
    h <- .analyze_hooks()
    expect_error(
        analyze(dataDirectory = "D", outlierRemoval = "bogus",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "unknown outlierRemoval method")
    expect_null(h$log$ppargs)                        # failed before preprocessing

    expect_error(
        analyze(dataDirectory = "D", normalization = "bogus",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "unknown normalization method")
})

test_that("analyze closes the GDS even when a hook errors", {
    h <- .analyze_hooks()
    expect_error(
        analyze(dataDirectory = "D",
                outlierRemoval = function(gds, res) stop("outlier hook failed"),
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        regexp = "outlier hook failed")
    expect_identical(h$log$closed, 1L)
})

# ---- analyze() deprecated argument aliases ---------------------------
#
# Renamed arguments keep working via a deprecation warning + remap. Two of them
# also changed semantics, so the alias reconstructs the old behaviour:
# targetPattern (a filename fragment) -> a samplesheet path under dataDirectory,
# and datasetClass (a prefix) -> a gdsOutput path with the .gds extension.

test_that("analyze remaps deprecated preprocessing aliases with a warning", {
    h <- .analyze_hooks()
    expect_warning(
        analyze(dataDirectory = "D", targetPattern = "batch1",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        "targetPattern.*deprecated")
    expect_identical(h$log$ppargs$samplesheet,
                     file.path("D", "samplesheet_batch1.csv"))
    expect_null(h$log$ppargs$targetPattern)

    h <- .analyze_hooks()
    expect_warning(
        analyze(dataDirectory = "D", datasetClass = "my_study",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        "datasetClass.*deprecated")
    expect_identical(h$log$ppargs$gdsOutput, "my_study.gds")

    h <- .analyze_hooks()
    expect_warning(
        analyze(dataDirectory = "D", nonSpecificProbesPath = "/x.csv",
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        "nonSpecificProbesPath.*deprecated")
    expect_identical(h$log$ppargs$crossReactiveProbes, "/x.csv")
})

test_that("analyze remaps the deprecated FUN and normalize aliases", {
    h <- .analyze_hooks()
    expect_warning(
        out <- analyze(dataDirectory = "D",
                       FUN = function(gds, res) "ran",
                       .preprocess = h$pre, .open = h$open, .close = h$clos),
        "FUN.*deprecated")
    expect_identical(out, "ran")           # the deprecated FUN ran as `analysis`

    h <- .analyze_hooks(); ran <- FALSE
    expect_warning(
        analyze(dataDirectory = "D",
                normalize = function(gds, res, keep) ran <<- TRUE,
                .preprocess = h$pre, .open = h$open, .close = h$clos),
        "normalize.*deprecated")
    expect_true(ran)                       # the deprecated normalize ran
})

test_that("analyze strips a deprecated alias passed as NULL (legacy FUN = NULL)", {
    h <- .analyze_hooks()
    expect_warning(
        out <- analyze(dataDirectory = "D", FUN = NULL,
                       .preprocess = h$pre, .open = h$open, .close = h$clos),
        "FUN.*deprecated")
    # FUN must not be forwarded to .preprocess: a fixed-formal runPreprocess
    # would reject it ("unused argument"). Name-based stripping removes it even
    # though the value is NULL (a value-based `!is.null` check would not).
    expect_false("FUN" %in% names(h$log$ppargs))
    # FUN = NULL means "just build the GDS", so the build result is returned and
    # the GDS is never opened.
    expect_identical(out, h$res)
    expect_null(h$log$opened)
})

# ---- analyze() unknown-argument guard --------------------------------
#
# The real .preprocess (runPreprocess) has fixed formals and no `...`, so a
# stray/misspelled name would surface as R's opaque "unused argument". analyze
# catches it first with an actionable message. A fixed-formal mock (no `...`)
# reproduces that surface here; the `function(...)` hooks elsewhere bypass it.

test_that("analyze rejects an unknown forwarded argument, listing valid names", {
    pre <- function(dataDirectory, annotationPackage, verbose, forceRebuild) {
        list(gds_path = "/tmp/x.gds")
    }
    err <- tryCatch(
        analyze(dataDirectory = "D", annotationPackage = "A", bogusArg = 1,
                .preprocess = pre, .open = function(p) "G",
                .close = function(g) NULL),
        error = identity)
    expect_s3_class(err, "fastMethyl_user_error")
    expect_match(conditionMessage(err), "bogusArg")
    expect_match(conditionMessage(err), "Valid preprocessing arguments")
})

test_that("analyze forwards analyze-injected args past the unknown-argument guard", {
    # verbose and forceRebuild are injected by analyze; the guard must not flag
    # them, so a fixed-formal .preprocess that declares them runs cleanly.
    pre <- function(dataDirectory, annotationPackage, verbose, forceRebuild) {
        list(gds_path = "/tmp/x.gds", targets = "T", rebuilt = FALSE)
    }
    out <- analyze(dataDirectory = "D", annotationPackage = "A",
                   .preprocess = pre, .open = function(p) "G",
                   .close = function(g) NULL)
    expect_identical(out$gds_path, "/tmp/x.gds")
})

# ---- analyze() staged rebuild ----------------------------------------
#
# `rebuild` names a stage; that stage and every stage after it rebuild (the
# pipeline is a linear chain build -> normalize -> the FUN's outputs). analyze
# realises the build stage through .preprocess's forceRebuild and cascades the
# rest via res$rebuildDownstream, which the cache-aware "dasen" hook and the
# user's FUN read.

test_that(".resolveRebuild maps stage names, none/all, and logicals", {
    expect_identical(.resolveRebuild(NULL), list(build = FALSE, normalize = FALSE))
    expect_identical(.resolveRebuild("none"), list(build = FALSE, normalize = FALSE))
    expect_identical(.resolveRebuild(FALSE), list(build = FALSE, normalize = FALSE))
    expect_identical(.resolveRebuild("normalize"), list(build = FALSE, normalize = TRUE))
    expect_identical(.resolveRebuild("build"), list(build = TRUE, normalize = TRUE))
    expect_identical(.resolveRebuild("all"), list(build = TRUE, normalize = TRUE))
    expect_identical(.resolveRebuild(TRUE), list(build = TRUE, normalize = TRUE))
})

test_that(".resolveRebuild treats forceRebuild as a deprecated alias", {
    expect_warning(r <- .resolveRebuild(NULL, TRUE), "deprecated")
    expect_identical(r, list(build = TRUE, normalize = TRUE))
    expect_warning(r <- .resolveRebuild(NULL, FALSE), "deprecated")
    expect_identical(r, list(build = FALSE, normalize = FALSE))
    # An explicit `rebuild` wins over the alias and skips the deprecation path.
    expect_silent(r <- .resolveRebuild("normalize", TRUE))
    expect_identical(r, list(build = FALSE, normalize = TRUE))
})

test_that(".resolveRebuild rejects an unknown stage or malformed value", {
    expect_error(.resolveRebuild("cellcounts"), "unknown rebuild stage")
    expect_error(.resolveRebuild(c("build", "normalize")), "must be a stage name")
})

test_that("analyze translates rebuild to forceRebuild and cascades downstream", {
    capture <- function(rebuildArg, rebuilt) {
        seen <- new.env()
        analyze(
            dataDirectory = "D", rebuild = rebuildArg,
            analysis = function(gds, res) { seen$down <- res$rebuildDownstream; NULL },
            .preprocess = function(...) {
                seen$fr <- list(...)$forceRebuild
                list(gds_path = "/tmp/x.gds", rebuilt = rebuilt)
            },
            .open = function(p) "GDS", .close = function(g) NULL)
        seen
    }
    # rebuild = "build": forceRebuild TRUE to preprocess, downstream rebuilds.
    s <- capture("build", FALSE); expect_true(s$fr); expect_true(s$down)
    # rebuild = "normalize": build reused (forceRebuild FALSE) but downstream rebuilds.
    s <- capture("normalize", FALSE); expect_false(s$fr); expect_true(s$down)
    # rebuild = "none" with a cache hit: nothing downstream rebuilds.
    s <- capture("none", FALSE); expect_false(s$fr); expect_false(s$down)
    # rebuild = "none" but the build was a cache MISS (rebuilt = TRUE): cascades
    # so normbetas / FUN outputs are not silently reused against a fresh GDS.
    s <- capture("none", TRUE); expect_false(s$fr); expect_true(s$down)
})

# ---- analyze() verbose diagnostics -----------------------------------

test_that("analyze(verbose = 2L) logs each external-data phase", {
    h <- .analyze_hooks()
    msgs <- character(0)
    withCallingHandlers(
        analyze(dataDirectory = "D", verbose = 2L,
                analysis = function(gds, res) NULL,
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
                analysis = function(gds, res) NULL,
                .preprocess = h$pre, .open = h$open, .close = h$clos))
})
