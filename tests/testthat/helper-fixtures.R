# Shared fixtures for the fork's testthat suite. Every reader test in this
# directory runs against the minfiData 450k IDATs and (optionally) the
# minfiDataEPIC IDATs; these helpers locate them and skip the calling
# test gracefully if the dataset package is not installed in the
# environment the tests are being run in.

minfidata_basenames <- function() {
    skip_if_not_installed("minfiData")
    baseDir <- system.file("extdata", package = "minfiData")
    files <- list.files(baseDir, pattern = "_Grn\\.idat$",
                        recursive = TRUE, full.names = TRUE)
    sub("_Grn\\.idat$", "", files)
}

minfidataepic_basenames <- function() {
    skip_if_not_installed("minfiDataEPIC")
    baseDir <- system.file("extdata", package = "minfiDataEPIC")
    files <- list.files(baseDir, pattern = "_Grn\\.idat$",
                        recursive = TRUE, full.names = TRUE)
    sub("_Grn\\.idat$", "", files)
}

# Gate for tests that read real IDATs, write GDS files, or otherwise
# perform non-trivial I/O. Unit tests should not call this; integration
# files call it at the top of every test_that block. Set
# FASTMETHYL_RUN_INTEGRATION=1 to run them.
skip_if_slow <- function() {
    if (!nzchar(Sys.getenv("FASTMETHYL_RUN_INTEGRATION"))) {
        skip("integration test (real IDAT I/O); set FASTMETHYL_RUN_INTEGRATION=1 to run")
    }
}

# Gate for tests that drive processMethArray(): on top of the minfiData IDATs
# they need bigmelon + gdsfmt to write and read the GDS file the worked
# pipeline.r produces.
skip_if_no_gds <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
    skip_if_not_installed("bigmelon")
    skip_if_not_installed("gdsfmt")
}

# processMethArray() resolves its gdsfmt writers through the package namespace,
# but the worked pipeline.r attaches bigmelon + gdsfmt, and the downstream GDS
# assertions here read the file back through them. Attach them for any test
# that drives the GDS path.
attach_gds_pkgs <- function() {
    suppressPackageStartupMessages({
        library(bigmelon)
        library(gdsfmt)
    })
}

# Materialize a GDS written by processMethArray() into the
# plain in-memory matrices + identifiers the tests assert on. gdsfmt stores
# the bulk matrices without dimnames; the locus and sample identifiers live
# in the fData/pData nodes (the bigmelon convention), so we read those
# separately.
.read_gds <- function(path) {
    gf <- gdsfmt::openfn.gds(path)
    on.exit(gdsfmt::closefn.gds(gf), add = TRUE)
    fData <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "fData"))
    pData <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "pData"))
    list(
        meth     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "methylated")),
        unmeth   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "unmethylated")),
        detP     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "pvals")),
        features = fData$Name,
        samples  = pData$barcode)
}

# Run the fused reader into a throwaway GDS and return its materialized
# contents. The on-disk file is deleted once it has been read back into
# memory, so callers get a self-contained snapshot.
.run_gds <- function(basenames, ...) {
    path <- tempfile(fileext = ".gds")
    on.exit(unlink(path), add = TRUE)
    processMethArray(basenames, gds_path = path, verbose = 0L, ...)
    .read_gds(path)
}
