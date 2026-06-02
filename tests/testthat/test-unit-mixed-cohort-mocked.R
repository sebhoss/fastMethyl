# Mixed-size cohort (same array type, different probe counts -- e.g.
# early-access EPIC vs final-release EPIC) is a real but rare config we
# cannot exercise from stock minfiData / minfiDataEPIC (they have either
# different array types or uniform sizes within type). Rather than
# fabricate IDAT binaries, mock illuminaio::readIDAT at the minfi
# namespace binding so the master sees two synthetic samples whose
# probe counts differ but both map to the same array via
# .guessArrayTypes. This exercises the verbose "mixed probe counts"
# message in readMethArray and the sampleIdx branch of pull().

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

# Synthetic IDAT-shape data: a Quants matrix with addresses "1".."N" as
# rownames and Mean / SD / NBeads columns. Two samples in the 450k band
# (622000..623000) with different probe counts; intersection is the
# smaller address set. The mock dispatches on filename so master and
# worker calls produce consistent shapes.
.mock_idat_factory <- function(path_to_n) {
    function(file, what = NULL, ...) {
        n <- NULL
        for (key in names(path_to_n)) {
            if (grepl(key, file, fixed = TRUE)) { n <- path_to_n[[key]]; break }
        }
        if (is.null(n)) stop("mock: unknown file ", file)
        if (identical(what, "nSNPsRead")) return(n)
        addrs <- as.character(seq_len(n))
        Quants <- matrix(as.integer(seq_len(n) %% 1000L + 100L),
                         nrow = n, ncol = 3L,
                         dimnames = list(addrs, c("Mean", "SD", "NBeads")))
        list(Quants = Quants)
    }
}

test_that("readMethArray prints the mixed-probe-count message and indexes per-sample addresses via sampleIdx", {
    # Two synthetic samples in the 450k band -- sample A at the canonical
    # 450k count, sample B at a slightly larger count (still inside the
    # 622000-623000 band). Intersection over their integer addresses is
    # the smaller set.
    tmpdir <- tempfile("mixed_mocked_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    bases <- c(file.path(tmpdir, "sampleA"), file.path(tmpdir, "sampleB"))
    for (b in bases) {
        file.create(paste0(b, "_Grn.idat"))
        file.create(paste0(b, "_Red.idat"))
    }

    mock <- .mock_idat_factory(list(sampleA = 622399L, sampleB = 622500L))
    # readIDAT is imported into the fastMethyl namespace; mock the consumer
    # binding so calls inside readMethArray's worker resolve to the mock.
    testthat::local_mocked_bindings(readIDAT = mock, .package = "fastMethyl")

    msgs <- character(0)
    rg <- withCallingHandlers(
        readMethArray(bases, verbose = 1L),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_s4_class(rg, "RGChannelSet")
    # Intersection of "1".."622399" and "1".."622500" is "1".."622399":
    # 622399 common addresses.
    expect_equal(nrow(rg), 622399L)
    expect_equal(ncol(rg), 2L)
    expect_true(any(grepl("mixed probe counts auto-detected", msgs)))
})
