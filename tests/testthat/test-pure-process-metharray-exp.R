# Input-validation tests for processMethArrayExp(). These cover the
# early validation block that fires before any IDAT read or annotation
# lookup -- all pure base R, so they run in the pure runner without
# library(minfi).

test_that("processMethArrayExp errors when targets is missing", {
    expect_error(processMethArrayExp(gds_path = "x.gds",
                                       annotation_package = "pkg"),
                 regexp = "data\\.frame.+Basename")
})

test_that("processMethArrayExp errors when targets is not a data.frame", {
    expect_error(processMethArrayExp(targets = list(Basename = "x"),
                                       gds_path = "x.gds",
                                       annotation_package = "pkg"),
                 regexp = "data\\.frame.+Basename")
})

test_that("processMethArrayExp errors when targets lacks a Basename column", {
    expect_error(processMethArrayExp(
        targets = data.frame(Sample_Name = "S1"),
        gds_path = "x.gds",
        annotation_package = "pkg"),
        regexp = "missing.+Basename")
})

test_that("processMethArrayExp errors when gds_path is missing or malformed", {
    df <- data.frame(Basename = "x", stringsAsFactors = FALSE)
    expect_error(processMethArrayExp(targets = df,
                                       annotation_package = "pkg"),
                 regexp = "gds_path")
    expect_error(processMethArrayExp(targets = df, gds_path = "",
                                       annotation_package = "pkg"),
                 regexp = "gds_path")
    expect_error(processMethArrayExp(targets = df, gds_path = c("a", "b"),
                                       annotation_package = "pkg"),
                 regexp = "gds_path")
})

test_that("processMethArrayExp errors when annotation_package is missing or malformed", {
    df <- data.frame(Basename = "x", stringsAsFactors = FALSE)
    expect_error(processMethArrayExp(targets = df, gds_path = "x.gds"),
                 regexp = "annotation_package")
    expect_error(processMethArrayExp(targets = df, gds_path = "x.gds",
                                       annotation_package = ""),
                 regexp = "annotation_package")
})

test_that("processMethArrayExp errors when drop_sex_chromosomes is not a single logical", {
    df <- data.frame(Basename = "x", stringsAsFactors = FALSE)
    expect_error(processMethArrayExp(targets = df, gds_path = "x.gds",
                                       annotation_package = "pkg",
                                       drop_sex_chromosomes = "yes"),
                 regexp = "drop_sex_chromosomes")
    expect_error(processMethArrayExp(targets = df, gds_path = "x.gds",
                                       annotation_package = "pkg",
                                       drop_sex_chromosomes = NA),
                 regexp = "drop_sex_chromosomes")
})
