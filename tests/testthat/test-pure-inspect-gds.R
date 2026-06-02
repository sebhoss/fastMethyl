# Pure tests for .inspectExistingGDS() -- the cache-validation helper that
# decides whether an existing GDS file can be trusted and reused instead
# of rebuilt. It only ever touches gdsfmt (namespace-qualified) and base
# R, so it runs in the pure runner: we fabricate real GDS files with
# gdsfmt directly, no minfi and no IDATs.

skip_if_not_installed("gdsfmt")

# Build a GDS that mimics the streaming writer's node layout. Toggles let
# each test omit a node (or mis-size betas) to drive a specific branch.
make_gds <- function(path, barcodes = c("bc1", "bc2"), betas_ncol = NULL,
                     with_betas = TRUE, with_pdata = TRUE,
                     with_barcode = TRUE) {
    g <- gdsfmt::createfn.gds(path)
    on.exit(gdsfmt::closefn.gds(g), add = TRUE)
    if (with_betas) {
        ncol <- if (is.null(betas_ncol)) length(barcodes) else betas_ncol
        gdsfmt::add.gdsn(g, "betas", val = matrix(0.5, nrow = 3L, ncol = ncol))
    }
    if (with_pdata) {
        pd <- if (with_barcode) {
            data.frame(barcode = barcodes, stringsAsFactors = FALSE)
        } else {
            data.frame(other = barcodes, stringsAsFactors = FALSE)
        }
        gdsfmt::add.gdsn(g, "pData", val = pd)
    }
    invisible(path)
}

sheet <- function(barcodes) {
    data.frame(Basename = file.path("slideX", barcodes),
               Sample_Name = paste0("S", seq_along(barcodes)),
               stringsAsFactors = FALSE)
}

test_that(".inspectExistingGDS returns NULL for a non-GDS file", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    writeLines("not a gds file", f)
    expect_null(.inspectExistingGDS(f, sheet(c("bc1", "bc2"))))
})

test_that(".inspectExistingGDS returns NULL when required nodes are absent", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    make_gds(f, with_pdata = FALSE)         # betas present, pData missing
    expect_null(.inspectExistingGDS(f, sheet(c("bc1", "bc2"))))
})

test_that(".inspectExistingGDS returns NULL when the barcode node is unreadable", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    make_gds(f, with_barcode = FALSE)       # pData present but no barcode child
    expect_null(.inspectExistingGDS(f, sheet(c("bc1", "bc2"))))
})

test_that(".inspectExistingGDS returns NULL when betas columns != barcodes", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    make_gds(f, barcodes = c("bc1", "bc2"), betas_ncol = 3L)
    expect_null(.inspectExistingGDS(f, sheet(c("bc1", "bc2"))))
})

test_that(".inspectExistingGDS reconstructs targets for a trustworthy GDS", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    make_gds(f, barcodes = c("bc1", "bc2"))
    res <- .inspectExistingGDS(f, sheet(c("bc1", "bc2")))
    expect_type(res, "list")
    expect_identical(rownames(res$targets), c("bc1", "bc2"))
    expect_identical(res$keepSamples, c(TRUE, TRUE))
    expect_null(res$keepProbes)
    expect_identical(res$unmatched, 0L)
})

test_that(".inspectExistingGDS flags GDS columns absent from the samplesheet", {
    f <- tempfile(fileext = ".gds")
    on.exit(unlink(f), add = TRUE)
    # GDS has bc1 + bc9; samplesheet only knows bc1 -> bc9 is unmatched (NA row).
    make_gds(f, barcodes = c("bc1", "bc9"))
    res <- .inspectExistingGDS(f, sheet(c("bc1", "bc2")))
    expect_identical(res$unmatched, 1L)
    expect_identical(res$keepSamples, c(TRUE, FALSE))
    expect_true(is.na(res$targets["bc9", "Sample_Name"]))
})
