# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Installation-contract tests. fastMethyl is a thin performance layer that
# imports a fixed surface from minfi and never vendors it -- so it installs and
# loads only if the minfi it is resolved against exports exactly that surface,
# and only if every package its NAMESPACE pulls from is also declared in
# DESCRIPTION. These tests fail loudly, naming the offending symbol or package,
# instead of letting a mismatch surface as a cryptic load error or mid-pipeline
# "object not found".
#
# The surface is small precisely so it stays stable across minfi versions: both
# the Bioconductor release and the sebhoss/minfi fork export all of it. (The
# fork carries a higher version -- 1.59.0.9000 vs the release's 1.59.0 -- so a
# user who installed the fork first and later installs fastMethyl can have
# BiocManager try to "correct" the too-new minfi back to the release; whichever
# minfi ends up installed, these tests assert it still satisfies the contract.)

# Values fastMethyl imports from minfi (mirrors importFrom(minfi, ...) in
# NAMESPACE). The integer-indexed reader and the streaming pipeline reach minfi
# only through these.
minfi_value_contract <- c(
    "getProbeInfo", "getManifest", "getManifestInfo", "getControlAddress",
    "getAnnotation", "RGChannelSet", "RGChannelSetExtended"
)

# S4 classes fastMethyl imports from minfi (importClassesFrom(minfi, ...)). The
# readers construct and dispatch on these.
minfi_class_contract <- c("RGChannelSet", "RGChannelSetExtended")

test_that("the installed minfi exports every value fastMethyl imports", {
    skip_if_not_installed("minfi")
    missing <- setdiff(minfi_value_contract, getNamespaceExports("minfi"))
    expect_identical(
        missing, character(0),
        info = paste("minfi no longer exports:", paste(missing, collapse = ", "))
    )
})

test_that("the installed minfi defines every class fastMethyl imports", {
    skip_if_not_installed("minfi")
    for (cl in minfi_class_contract) {
        expect_false(
            is.null(methods::getClassDef(cl, package = "minfi")),
            info = paste("minfi no longer defines class:", cl)
        )
    }
})

test_that("every package fastMethyl's NAMESPACE imports is declared in DESCRIPTION", {
    # An importFrom() / importClassesFrom() whose package is absent from
    # DESCRIPTION's Depends/Imports installs locally (the dep happens to be
    # present) but fails for a clean-room user with "there is no package called
    # '<pkg>'". Pin the two lists together.
    imported_pkgs <- setdiff(names(getNamespaceImports("fastMethyl")), "base")

    desc <- packageDescription("fastMethyl")
    declared <- unlist(lapply(c("Depends", "Imports", "LinkingTo"), function(field) {
        if (is.null(desc[[field]])) {
            return(character(0))
        }
        entries <- strsplit(desc[[field]], ",")[[1]]
        trimws(sub("\\s*\\(.*\\)\\s*$", "", entries)) # drop version constraints
    }))
    declared <- setdiff(declared, "R")

    undeclared <- setdiff(imported_pkgs, declared)
    expect_identical(
        undeclared, character(0),
        info = paste("imported in NAMESPACE but not declared in DESCRIPTION:",
                     paste(undeclared, collapse = ", "))
    )
})

test_that("analyze is the exported public entry point with the documented signature", {
    expect_true("analyze" %in% getNamespaceExports("fastMethyl"))
    expect_true(is.function(analyze))

    fmls <- names(formals(analyze))
    expect_true("..." %in% fmls)
    expect_true(all(c("analysis", "outlierRemoval", "normalization",
                      ".preprocess", ".open", ".close") %in% fmls))
})
