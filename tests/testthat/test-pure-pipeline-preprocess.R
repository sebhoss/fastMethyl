# Pure tests for runPreprocess()'s own guards -- the logic that
# runs after .validateArgs but before any IDAT byte is read, so it needs
# neither minfi nor real IDAT files. Behavioural coverage of the full
# read+QC+GDS pass lives in the integration suite.
#
# runPreprocess lives in the minfi namespace; the pure runner
# sources it (see helper-pure.R).

test_that("runPreprocess rejects a BPPARAM that is not a BiocParallelParam", {
    # Minimal valid on-disk fixture so .validateArgs passes and execution
    # reaches the BPPARAM resolution block.
    root <- tempfile("pp_")
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    data_dir <- file.path(root, "data")
    slide <- file.path(data_dir, "5723646052")
    dir.create(slide, recursive = TRUE)
    base <- file.path(slide, "5723646052_R02C02")
    file.create(paste0(base, "_Grn.idat"), paste0(base, "_Red.idat"))
    write.csv(
        data.frame(Sample_Name = "S1",
                   Basename = "5723646052/5723646052_R02C02",
                   stringsAsFactors = FALSE),
        file.path(data_dir, "samplesheet_real.csv"), row.names = FALSE)
    xreact_path <- file.path(root, "xreact.csv")
    write.csv(data.frame(TargetID = "cg00000029", stringsAsFactors = FALSE),
              xreact_path, row.names = FALSE)

    # A plain string fails the methods::is(..., "BiocParallelParam") guard
    # regardless of whether BiocParallel's classes are loaded, so this
    # errors before any IDAT/minfi work begins.
    expect_error(
        runPreprocess(
            dataDirectory       = data_dir,
            crossReactiveProbes = xreact_path,
            samplesheet         = file.path(data_dir, "samplesheet_real.csv"),
            gdsOutput           = "pp_test.gds",
            annotationPackage   = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            forceRebuild        = TRUE,
            BPPARAM             = "not-a-param",
            .check_annotation   = function(pkg) TRUE),
        regexp = "BPPARAM.*BiocParallelParam")
})

# ---- cache-reuse path (no IDAT read) ---------------------------------
#
# These build an on-disk fixture plus a real GDS at the gdsOutput path,
# so runPreprocess takes the .inspectExistingGDS cache-hit branch and
# never reaches processMethArrayExp. gdsfmt is namespace-qualified inside
# the helper, so only an install (not an attach, not minfi) is needed.

# Two slide samples; barcode is basename(Basename).
.pp_fixture <- function(barcodes = c("5723646052_R02C02", "5723646052_R04C01")) {
    root <- tempfile("ppcache_")
    data_dir <- file.path(root, "data")
    slide <- file.path(data_dir, "5723646052")
    dir.create(slide, recursive = TRUE)
    for (bc in barcodes) {
        base <- file.path(slide, bc)
        file.create(paste0(base, "_Grn.idat"), paste0(base, "_Red.idat"))
    }
    write.csv(
        data.frame(Sample_Name = paste0("S", seq_along(barcodes)),
                   Basename = file.path("5723646052", barcodes),
                   stringsAsFactors = FALSE),
        file.path(data_dir, "samplesheet_real.csv"), row.names = FALSE)
    xreact <- file.path(root, "xreact.csv")
    write.csv(data.frame(TargetID = "cg00000029", stringsAsFactors = FALSE),
              xreact, row.names = FALSE)
    list(root = root, data_dir = data_dir, xreact = xreact,
         sheet = file.path(data_dir, "samplesheet_real.csv"),
         gdsOutput = file.path(root, "cohort.gds"),
         gdsPath = file.path(root, "cohort.gds"))
}

# Minimal bigmelon-shaped GDS: betas matrix + pData/barcode.
.write_cache_gds <- function(path, barcodes) {
    g <- gdsfmt::createfn.gds(path)
    on.exit(gdsfmt::closefn.gds(g), add = TRUE)
    gdsfmt::add.gdsn(g, "betas",
                     val = matrix(0.5, nrow = 3L, ncol = length(barcodes)))
    gdsfmt::add.gdsn(g, "pData",
                     val = data.frame(barcode = barcodes,
                                      stringsAsFactors = FALSE))
    invisible(path)
}

# Write a build-key sidecar matching the fixture + the given config, so the
# cache-reuse path treats the GDS as up to date (the "match" branch). Defaults
# mirror runPreprocess's own defaults for the threshold arguments.
.write_matching_buildkey <- function(fx, annotationPackage = "AnnoPkg",
                                     sampleDetPThreshold = 0.01,
                                     probeDetPThreshold = 0.01) {
    key <- .currentBuildKey(
        annotationPackage, sampleDetPThreshold, probeDetPThreshold,
        fx$sheet, fx$xreact)
    saveRDS(key, .buildKeyPath(fx$gdsPath))
    invisible(key)
}

test_that("runPreprocess reuses a valid existing GDS instead of rebuilding", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    .write_cache_gds(fx$gdsPath, c("5723646052_R02C02", "5723646052_R04C01"))
    .write_matching_buildkey(fx)

    # BPPARAM = NULL (default) exercises the MulticoreParam derivation too. A
    # matching build-key sidecar means the reuse is silent -- no staleness
    # warning -- so expect_message alone captures the happy path.
    expect_message(
        res <- runPreprocess(
            dataDirectory         = fx$data_dir,
            crossReactiveProbes   = fx$xreact,
            samplesheet           = fx$sheet,
            gdsOutput             = fx$gdsOutput,
            annotationPackage     = "AnnoPkg",
            .check_annotation     = function(pkg) TRUE),
        "Reusing existing GDS")
    expect_equal(res$gds_path, fx$gdsPath)
    expect_equal(nrow(res$targets), 2L)
    expect_identical(res$keepSamples, c(TRUE, TRUE))
    expect_null(res$keepProbes)
})

test_that("runPreprocess warns but reuses a GDS that has no build-key sidecar", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    # No sidecar (e.g. a GDS built before build-key tracking) -> cannot verify
    # the config matches, so the file is reused with a staleness warning.
    .write_cache_gds(fx$gdsPath, c("5723646052_R02C02", "5723646052_R04C01"))

    expect_warning(
        res <- suppressMessages(runPreprocess(
            dataDirectory         = fx$data_dir,
            crossReactiveProbes   = fx$xreact,
            samplesheet           = fx$sheet,
            gdsOutput             = fx$gdsOutput,
            annotationPackage     = "AnnoPkg",
            .check_annotation     = function(pkg) TRUE)),
        regexp = "no build-key sidecar")
    expect_equal(res$gds_path, fx$gdsPath)
})

test_that("runPreprocess rebuilds when the build key no longer matches", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    .write_cache_gds(fx$gdsPath, c("5723646052_R02C02", "5723646052_R04C01"))
    # Sidecar records a different detection-p threshold than the call uses, so
    # the cached GDS is stale and the rebuild branch runs -- which reaches
    # processMethArrayExp and stops on the (fake) uninstalled annotation
    # package, confirming the cache was bypassed rather than reused.
    .write_matching_buildkey(fx, sampleDetPThreshold = 0.05)

    expect_message(
        expect_error(
            runPreprocess(
                dataDirectory         = fx$data_dir,
                crossReactiveProbes   = fx$xreact,
                samplesheet           = fx$sheet,
                gdsOutput             = fx$gdsOutput,
                annotationPackage     = "AnnoPkg",
                sampleDetPThreshold   = 0.01,
                .check_annotation     = function(pkg) TRUE),
            regexp = "annotation package not installed"),
        regexp = "different configuration")
})

test_that("runPreprocess warns that a supplied BPPARAM ignores readerWorkers", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    .write_cache_gds(fx$gdsPath, c("5723646052_R02C02", "5723646052_R04C01"))
    .write_matching_buildkey(fx)

    expect_warning(
        runPreprocess(
            dataDirectory         = fx$data_dir,
            crossReactiveProbes   = fx$xreact,
            samplesheet           = fx$sheet,
            gdsOutput             = fx$gdsOutput,
            annotationPackage     = "AnnoPkg",
            readerWorkers         = 2L,
            BPPARAM               = BiocParallel::SerialParam(),
            .check_annotation     = function(pkg) TRUE),
        regexp = "readerWorkers.*ignored")
})

test_that("runPreprocess warns when the cached GDS has samples absent from the sheet", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    # GDS holds one barcode the samplesheet does not -> unmatched > 0.
    .write_cache_gds(fx$gdsPath, c("5723646052_R02C02", "9999999999_R01C01"))
    .write_matching_buildkey(fx)

    expect_warning(
        suppressMessages(runPreprocess(
            dataDirectory         = fx$data_dir,
            crossReactiveProbes   = fx$xreact,
            samplesheet           = fx$sheet,
            gdsOutput             = fx$gdsOutput,
            annotationPackage     = "AnnoPkg",
            .check_annotation     = function(pkg) TRUE)),
        regexp = "not present in the current")
})

test_that("runPreprocess rebuilds when an existing GDS is unreadable", {
    skip_if_not_installed("gdsfmt")
    fx <- .pp_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    # A corrupt (non-GDS) file at the cache path -> .inspectExistingGDS
    # returns NULL -> the rebuild branch runs. processMethArrayExp then
    # stops because the (fake) annotation package is not installed, which is
    # exactly the rebuild path we want to exercise without real IDAT reads.
    writeLines("corrupt", fx$gdsPath)
    expect_error(
        suppressMessages(runPreprocess(
            dataDirectory         = fx$data_dir,
            crossReactiveProbes   = fx$xreact,
            samplesheet           = fx$sheet,
            gdsOutput             = fx$gdsOutput,
            annotationPackage     = "AnnoPkg",
            .check_annotation     = function(pkg) TRUE)),
        regexp = "annotation package not installed")
})

test_that(".compareBuildKey treats a corrupt or non-list sidecar as absent", {
    tmp <- tempfile(fileext = ".gds")
    on.exit(unlink(c(tmp, .buildKeyPath(tmp))), add = TRUE)
    key <- list(a = 1)
    # readRDS succeeds but yields a non-list -> cannot trust it -> absent.
    saveRDS(42L, .buildKeyPath(tmp))
    expect_identical(.compareBuildKey(tmp, key)$status, "absent")
    # garbage bytes -> readRDS errors -> tryCatch NULL -> absent.
    writeLines("not an rds file", .buildKeyPath(tmp))
    expect_identical(.compareBuildKey(tmp, key)$status, "absent")
})
