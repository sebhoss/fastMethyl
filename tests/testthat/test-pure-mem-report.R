# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Pure tests for the verbose-diagnostics helpers .mem_report() and .path_size().
# No IDATs, minfi, or GDS needed -- just the filesystem and /proc (or the gc
# fallback).

test_that(".mem_report returns a single non-empty string", {
    r <- .mem_report()
    expect_type(r, "character")
    expect_length(r, 1L)
    expect_true(nzchar(r))
    # On Linux it reports RSS; everywhere it is at least a heap figure.
    expect_match(r, "RSS|heap")
})

test_that(".path_size reports a human size for a real file", {
    tmp <- tempfile()
    on.exit(unlink(tmp), add = TRUE)
    writeBin(raw(2048L), tmp)
    s <- .path_size(tmp)
    expect_type(s, "character")
    expect_match(s, "B")            # e.g. "2 kB"
})

test_that(".path_size sums a directory's contents", {
    d <- tempfile()
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    writeBin(raw(1024L), file.path(d, "a.bin"))
    writeBin(raw(1024L), file.path(d, "b.bin"))
    expect_match(.path_size(d), "B")
})

test_that(".path_size returns sentinels for missing or empty paths", {
    expect_identical(.path_size(file.path(tempdir(), "does-not-exist-xyz")),
                     "missing")
    expect_identical(.path_size(NA_character_), "unknown")
    expect_identical(.path_size(""), "unknown")
    expect_identical(.path_size(c("a", "b")), "unknown")
})
