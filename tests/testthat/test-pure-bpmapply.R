# Pure tests for the fault-tolerant bpmapply wrappers in R/zzz.R.
#
# These exercise the failure/retry machinery that the integration suite
# never reaches (it only ever reads valid IDATs, so a worker never
# fails). Everything here runs on BiocParallel::SerialParam(), so no
# minfi load and no real I/O is required -- BiocParallel is attached for
# bare-name resolution of bpmapply/bptry/bpok in the pure runner.

suppressPackageStartupMessages(library(BiocParallel))

test_that(".bpmapply_with_progress rejects malformed iter_args", {
    expect_error(
        .bpmapply_with_progress(FUN = identity, iter_args = 1,
                                BPPARAM = SerialParam(), verbose = 0L,
                                label = "t"))
    expect_error(
        .bpmapply_with_progress(FUN = identity, iter_args = list(),
                                BPPARAM = SerialParam(), verbose = 0L,
                                label = "t"))
    expect_error(
        .bpmapply_with_progress(FUN = `+`, iter_args = list(x = 1:3, y = 1:2),
                                BPPARAM = SerialParam(), verbose = 0L,
                                label = "t"))
})

test_that(".bpmapply_with_progress (verbose < 2) delegates and keeps names", {
    res <- .bpmapply_with_progress(
        FUN = function(x, y) x + y,
        iter_args = list(x = c(a = 1, b = 2), y = c(10, 20)),
        BPPARAM = SerialParam(), verbose = 0L, label = "t")
    expect_equal(res, list(a = 11, b = 22))
})

test_that(".bpmapply_with_progress (verbose >= 2) chunks and reports progress", {
    # Unnamed iterating arg exercises the names-less result_names branch;
    # n = 5 with one SerialParam worker yields >1 chunk -> >1 progress line.
    expect_message(
        res <- .bpmapply_with_progress(
            FUN = function(x) x * 2L,
            iter_args = list(x = 1:5),
            BPPARAM = SerialParam(), verbose = 2L, label = "prog"),
        "progress")
    expect_equal(unlist(res, use.names = FALSE), c(2L, 4L, 6L, 8L, 10L))
    # No-names branch: result_names is character(n), i.e. all empty.
    expect_true(all(names(res) == ""))
})

test_that(".bpmapply_try returns the plain result on full success", {
    res <- .bpmapply_try(
        FUN = function(x) x + 1L,
        iter_args = list(x = c(a = 1L, b = 2L)),
        MoreArgs = list(), BPPARAM = SerialParam(), label = "ok")
    expect_equal(res, list(a = 2L, b = 3L))
})

test_that(".bpmapply_try retries failed elements once and recovers", {
    # A FUN that fails the first time element 2 is seen and succeeds on the
    # retry. The per-key counter lives in an env so SerialParam's in-process
    # re-run (via BPREDO) observes the incremented count.
    calls <- new.env(parent = emptyenv())
    flaky <- function(x) {
        key <- as.character(x)
        seen <- if (is.null(calls[[key]])) 0L else calls[[key]]
        calls[[key]] <- seen + 1L
        if (x == 2L && seen == 0L) stop("transient failure")
        x * 10L
    }
    expect_message(
        res <- .bpmapply_try(
            FUN = flaky, iter_args = list(x = c(1L, 2L, 3L)),
            MoreArgs = list(), BPPARAM = SerialParam(), label = "retry"),
        "retrying once")
    expect_equal(unlist(res, use.names = FALSE), c(10L, 20L, 30L))
})

test_that(".bpmapply_try stops, naming inputs, when a failure persists", {
    # Named inputs so the error message uses the names branch (nm[failed]).
    boom <- function(x) if (x == 2L) stop("always boom") else x
    expect_error(
        .bpmapply_try(
            FUN = boom, iter_args = list(x = c(a = 1L, b = 2L, c = 3L)),
            MoreArgs = list(), BPPARAM = SerialParam(), label = "boom"),
        "failed after one retry")
})

test_that(".bpmapply_try truncates the failure list past ten elements", {
    # 12 unnamed failing inputs: hits the as.character(failed) label branch
    # and the '... and N more' truncation arm (n_show = 10).
    allbad <- function(x) stop("nope")
    expect_error(
        .bpmapply_try(
            FUN = allbad, iter_args = list(x = 1:12),
            MoreArgs = list(), BPPARAM = SerialParam(), label = "many"),
        "and 2 more")
})
