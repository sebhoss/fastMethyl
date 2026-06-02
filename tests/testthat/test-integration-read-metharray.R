# Tests for the fork's readMethArray() and readMethArrayExp().
# The fork adds BPPARAM (parallel IDAT reading) and an integer-indexed
# assembly hot path; these tests pin the dims, class, and content
# invariants that the example pipeline depends on, plus the BPPARAM
# parallel/serial equivalence that the integer-indexed worker path
# could plausibly break.

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

test_that("readMethArray returns an RGChannelSet with the expected dims for minfiData IDATs", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()
    rg <- readMethArray(basenames, verbose = 0L)

    expect_s4_class(rg, "RGChannelSet")
    expect_equal(ncol(rg), length(basenames))
    expect_gt(nrow(rg), 600000L)
    expect_setequal(colnames(rg), basename(basenames))
    expect_equal(annotation(rg)[["array"]], "IlluminaHumanMethylation450k")
})

test_that("readMethArray(extended = TRUE) returns SDs and NBeads matrices", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()
    rge <- readMethArray(basenames, extended = TRUE, verbose = 0L)

    expect_s4_class(rge, "RGChannelSetExtended")
    expect_true(all(c("Green", "Red", "GreenSD", "RedSD", "NBeads")
                    %in% assayNames(rge)))
    expect_equal(dim(assay(rge, "GreenSD")), dim(assay(rge, "Green")))
    expect_equal(dim(assay(rge, "NBeads")),  dim(assay(rge, "Green")))
})

test_that("readMethArray output is invariant under BPPARAM choice", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()

    rg_serial <- readMethArray(basenames,
                                BPPARAM = BiocParallel::SerialParam(),
                                verbose = 0L)
    rg_par <- readMethArray(basenames,
                             BPPARAM = BiocParallel::MulticoreParam(2L),
                             verbose = 0L)

    # The integer-indexed worker hot path runs on both backends; if any
    # character-rowname lookup re-enters silently on one backend, these
    # matrices diverge.
    expect_identical(getGreen(rg_par),  getGreen(rg_serial))
    expect_identical(getRed(rg_par),    getRed(rg_serial))
    expect_identical(rownames(rg_par),  rownames(rg_serial))
    expect_identical(colnames(rg_par),  colnames(rg_serial))
    expect_identical(annotation(rg_par), annotation(rg_serial))
})

test_that("readMethArray errors loudly on missing IDAT files", {
    skip_if_slow()
    skip_if_not_installed("BiocParallel")
    expect_error(
        readMethArray(file.path(tempdir(), "definitely_not_a_real_basename"),
                       verbose = 0L),
        regexp = "do not exist")
})

test_that("readMethArrayExp(targets = ...) wires Basename into pData and matches readMethArray", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()
    targets <- data.frame(
        Sample_Name = paste0("S", seq_along(basenames)),
        Basename = basenames,
        stringsAsFactors = FALSE)

    rg_exp <- readMethArrayExp(targets = targets, verbose = 0L)
    rg_ref <- readMethArray(basenames, verbose = 0L)

    expect_s4_class(rg_exp, "RGChannelSet")
    expect_identical(getGreen(rg_exp), getGreen(rg_ref))
    expect_identical(getRed(rg_exp),   getRed(rg_ref))
    expect_true("Sample_Name" %in% colnames(colData(rg_exp)))
    expect_equal(colData(rg_exp)$Sample_Name, targets$Sample_Name)
})

test_that("readMethArrayExp(targets) requires a Basename column", {
    skip_if_slow()
    skip_if_not_installed("BiocParallel")
    bad <- data.frame(Sample_Name = "S1", path = "x", stringsAsFactors = FALSE)
    expect_error(readMethArrayExp(targets = bad), regexp = "Basename")
})
