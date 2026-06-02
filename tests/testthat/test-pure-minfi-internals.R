# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Pure tests for the minfi internal helper reimplemented in
# R/minfi-internals.R (.getAnnotationString). These need neither minfi S4
# objects nor I/O.

test_that(".getAnnotationString resolves scalar and named annotations", {
    expect_equal(.getAnnotationString("ilmn12.hg19"), "ilmn12.hg19anno")
    expect_equal(
        unname(.getAnnotationString(c(array = "IlluminaHumanMethylation450k",
                                      annotation = "ilmn12.hg19"))),
        "IlluminaHumanMethylation450kanno.ilmn12.hg19")
})

test_that(".getAnnotationString rejects Unknown and unrecognised shapes", {
    expect_error(.getAnnotationString("Unknown"), "Unknown")
    expect_error(.getAnnotationString(c(array = "Unknown", annotation = "x")),
                 "Unknown")
    expect_error(.getAnnotationString(c("a", "b")),
                 "unable to get the annotation string")
})
