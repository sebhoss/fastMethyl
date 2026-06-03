# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0

library(testthat)
# Integration tests assert on minfi / SummarizedExperiment generics
# (annotation, getGreen, assay, colData, ...), so those packages are attached.
# fastMethyl is attached last: its readMethArray* deliberately mask minfi's.
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(minfi)
  library(fastMethyl)
})

test_check("fastMethyl")
