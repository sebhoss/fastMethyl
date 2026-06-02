# Unit tests for the small helpers in R/zzz.R that the fork relies on
# everywhere. These have no fixtures, so they run unconditionally and
# act as the suite's cheap canary if R/zzz.R is broken (every reader
# test depends on .normalize_verbose).

test_that(".normalize_verbose accepts integers in [0, 2]", {
    expect_identical(.normalize_verbose(0L), 0L)
    expect_identical(.normalize_verbose(1L), 1L)
    expect_identical(.normalize_verbose(2L), 2L)
})

test_that(".normalize_verbose accepts numeric coercible to [0, 2]", {
    expect_identical(.normalize_verbose(0), 0L)
    expect_identical(.normalize_verbose(2), 2L)
})

test_that(".normalize_verbose coerces legacy TRUE/FALSE", {
    expect_identical(.normalize_verbose(TRUE),  1L)
    expect_identical(.normalize_verbose(FALSE), 0L)
})

test_that(".normalize_verbose rejects out-of-range and malformed inputs", {
    expect_error(.normalize_verbose(3L),         regexp = "0, 1, or 2")
    expect_error(.normalize_verbose(-1L),        regexp = "0, 1, or 2")
    expect_error(.normalize_verbose(NA),         regexp = "verbose")
    expect_error(.normalize_verbose(c(1L, 2L)),  regexp = "verbose")
    expect_error(.normalize_verbose("1"),        regexp = "verbose")
})
