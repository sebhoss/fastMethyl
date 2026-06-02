# Coverage-driven edge case tests: branches in readMethArray() and
# processMethArray() that the happy-path tests do not exercise.
# Each test pins one specific branch and is paired with a comment
# identifying what it covers.

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

# Copies a real IDAT to a gzipped path under `dir_base`. Used by the
# .gz-fallback tests below: the readers' path-resolution code appends
# ".gz" when "<base>_<channel>.idat" is missing on disk.
.copy_as_gz <- function(src_path, dst_gz_path) {
    con_in <- file(src_path, "rb")
    on.exit(close(con_in), add = TRUE)
    con_out <- gzfile(dst_gz_path, "wb")
    on.exit(close(con_out), add = TRUE)
    repeat {
        buf <- readBin(con_in, what = "raw", n = 1e6)
        if (length(buf) == 0L) break
        writeBin(buf, con_out)
    }
}

# Fabricates a workspace with each sample's IDATs gzipped (no plain
# .idat present), returns the base names that callers should pass to
# readMethArray() / processMethArray().
.gz_only_workspace <- function() {
    src <- minfidata_basenames()
    dir <- tempfile("gz_idat_")
    dir.create(dir)
    out_bases <- character(length(src))
    for (i in seq_along(src)) {
        base_out <- file.path(dir, sprintf("s%02d", i))
        .copy_as_gz(paste0(src[i], "_Grn.idat"), paste0(base_out, "_Grn.idat.gz"))
        .copy_as_gz(paste0(src[i], "_Red.idat"), paste0(base_out, "_Red.idat.gz"))
        out_bases[i] <- base_out
    }
    attr(out_bases, "dir") <- dir
    out_bases
}

test_that("readMethArray follows the .gz fallback when only gzipped IDATs are on disk", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- .gz_only_workspace()
    on.exit(unlink(attr(bases, "dir"), recursive = TRUE), add = TRUE)
    rg <- readMethArray(bases, verbose = 0L)
    expect_s4_class(rg, "RGChannelSet")
    expect_equal(ncol(rg), length(bases))
})

test_that("processMethArray follows the .gz fallback when only gzipped IDATs are on disk", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    bases <- .gz_only_workspace()
    on.exit(unlink(attr(bases, "dir"), recursive = TRUE), add = TRUE)
    out <- .run_gds(bases)
    expect_equal(length(out$samples), length(bases))
})

test_that("readMethArray refuses to merge IDATs from different array types", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("minfiDataEPIC")
    # 450k + EPIC have different probe counts AND different array names,
    # which hits the !sameLength && !sameArray branch.
    mixed <- c(minfidata_basenames()[1L], minfidataepic_basenames()[1L])
    expect_error(readMethArray(mixed, verbose = 0L),
                 regexp = "different size and type")
})

test_that("processMethArray refuses when address intersection across versions is empty", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("minfiDataEPIC")
    mixed <- c(minfidata_basenames()[1L], minfidataepic_basenames()[1L])
    # The 2nd-version address read fires (it is a different probe count),
    # then .streamingBuildIndices intersects the two address spaces and
    # finds zero shared loci -- triggering the "no loci left" stop, which
    # fires before any GDS write.
    expect_error(
        processMethArray(mixed, gds_path = tempfile(fileext = ".gds"),
                          verbose = 0L),
        regexp = "no loci left")
})

test_that("processMethArray refuses when drop_probes empties the manifest", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
    bases <- minfidata_basenames()
    # Dropping every locus in the manifest empties the output, tripping the
    # "no loci left" stop before any GDS write. Pull the locus universe from
    # the resolved 450k manifest rather than a full read.
    manifest <- fastMethyl:::.streamingResolveManifest(622399L)$manifest
    drop_all <- getManifestInfo(manifest, "locusNames")
    expect_error(
        processMethArray(bases, gds_path = tempfile(fileext = ".gds"),
                          drop_probes = drop_all, verbose = 0L),
        regexp = "no loci left")
})

test_that("processMethArray(pData = ...) writes the supplied phenotype columns into the GDS pData node", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("bigmelon")
    skip_if_not_installed("gdsfmt")
    suppressPackageStartupMessages({
        library(bigmelon)
        library(gdsfmt)
    })
    bases <- minfidata_basenames()
    pd <- data.frame(
        Sample_Name = paste0("S", seq_along(bases)),
        Group = rep(c("A", "B"), length.out = length(bases)),
        stringsAsFactors = FALSE)
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    processMethArray(bases, gds_path = gds_path,
                         pData = pd, verbose = 0L)
    gf <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(gf), add = TRUE)
    pd_node <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "pData"))
    expect_true("Group" %in% names(pd_node))
    expect_equal(pd_node$Group, pd$Group)
})

test_that("readMethArrayExp(base = dir, recursive = TRUE) scans IDATs by directory listing", {
    skip_if_slow()
    skip_if_no_minfidata()
    baseDir <- system.file("extdata", package = "minfiData")
    rg <- readMethArrayExp(base = baseDir, recursive = TRUE, verbose = 0L)
    expect_s4_class(rg, "RGChannelSet")
    expect_gt(ncol(rg), 0L)
})

test_that("readMethArraySheet parses a CSV next to IDATs and emits a Basename column", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    # Build a minimal samplesheet next to the IDATs in a tempdir; reuse
    # the IDATs via symlink to keep the test self-contained.
    workdir <- tempfile("sheet_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    for (b in bases) {
        for (ch in c("_Grn.idat", "_Red.idat")) {
            file.symlink(paste0(b, ch),
                         file.path(workdir, paste0(basename(b), ch)))
        }
    }
    # Sentrix_ID -> Slide, Sentrix_Position -> Array; readMethArraySheet
    # combines them to look up <Slide>_<Array>_Grn.idat in the directory.
    sheet <- data.frame(
        Sample_Name = paste0("S", seq_along(bases)),
        Sentrix_ID = vapply(strsplit(basename(bases), "_"), `[[`, character(1L), 1L),
        Sentrix_Position = vapply(strsplit(basename(bases), "_"), `[[`, character(1L), 2L),
        stringsAsFactors = FALSE)
    write.csv(sheet, file.path(workdir, "sheet.csv"), row.names = FALSE)

    out <- readMethArraySheet(base = workdir, pattern = "sheet\\.csv$",
                                recursive = FALSE, verbose = FALSE)
    expect_true("Basename" %in% names(out))
    expect_equal(nrow(out), length(bases))
})

test_that("readMethArray reports Red-channel-only file misses (not silently emits the Green-channel mask)", {
    skip_if_slow()
    skip_if_no_minfidata()
    src <- minfidata_basenames()
    dir <- tempfile("redmissing_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bases <- character(length(src))
    for (i in seq_along(src)) {
        b <- file.path(dir, sprintf("s%02d", i))
        file.symlink(paste0(src[i], "_Grn.idat"), paste0(b, "_Grn.idat"))
        # Only symlink Red for samples 1..(N-1); leave the last Red missing.
        if (i < length(src)) {
            file.symlink(paste0(src[i], "_Red.idat"), paste0(b, "_Red.idat"))
        }
        bases[i] <- b
    }
    err <- tryCatch(readMethArray(bases, verbose = 0L),
                    error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "do not exist")
    # The fix here is that the missing file actually appears in the
    # message -- before the fix the mask referenced the wrong vector and
    # the list was empty.
    expect_match(conditionMessage(err), sprintf("s%02d_Red\\.idat", length(src)))
})

test_that("processMethArray reports Red-channel-only file misses", {
    skip_if_slow()
    skip_if_no_minfidata()
    src <- minfidata_basenames()
    dir <- tempfile("redmissing_stream_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    bases <- character(length(src))
    for (i in seq_along(src)) {
        b <- file.path(dir, sprintf("s%02d", i))
        file.symlink(paste0(src[i], "_Grn.idat"), paste0(b, "_Grn.idat"))
        if (i < length(src)) {
            file.symlink(paste0(src[i], "_Red.idat"), paste0(b, "_Red.idat"))
        }
        bases[i] <- b
    }
    err <- tryCatch(
        processMethArray(bases, gds_path = tempfile(fileext = ".gds"),
                          verbose = 0L),
        error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), sprintf("s%02d_Red\\.idat", length(src)))
})

test_that("processMethArray errors when Green-channel IDATs cannot be resolved", {
    skip_if_slow()
    expect_error(
        processMethArray(file.path(tempdir(), "not_a_real_basename"),
                          gds_path = tempfile(fileext = ".gds"),
                          verbose = 0L),
        regexp = "do not exist")
})

test_that("processMethArray(verbose = 2L) prints the mixed-cohort detection message when probe counts differ", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("minfiDataEPIC")
    mixed <- c(minfidata_basenames()[1L], minfidataepic_basenames()[1L])
    # The verbose >= 2 mixed-cohort line fires before .streamingBuildIndices
    # hits the empty-intersection stop, so capture the messages emitted
    # ahead of the error.
    msgs <- character(0)
    err <- tryCatch(
        withCallingHandlers(
            processMethArray(mixed, gds_path = tempfile(fileext = ".gds"),
                              verbose = 2L),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            }),
        error = identity)
    expect_s3_class(err, "error")
    expect_true(any(grepl("mixed-version cohort detected", msgs)))
})

test_that("processMethArray(verbose = 2L) prints the dropped-samples list when QC removes any", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    bases <- minfidata_basenames()
    # An impossibly strict sample threshold drops every sample, the only
    # deterministic way to trigger the QC-block path on minfiData's uniform
    # colMeans. The empty-cohort GDS write may then fail; the messages we
    # assert on are emitted before that, so tolerate a downstream error.
    msgs <- character(0)
    tryCatch(
        withCallingHandlers(
            processMethArray(bases, gds_path = tempfile(fileext = ".gds"),
                              sample_detP_threshold = 1e-300,
                              verbose = 2L),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            }),
        error = identity)
    expect_true(any(grepl("dropped samples:", msgs)))
    expect_true(any(grepl("QC drops .+ samples", msgs)))
})

test_that("readMethArray realizes Green/Red matrices through a non-NULL DelayedArray backend when one is set", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("HDF5Array")
    skip_if_not_installed("DelayedArray")
    prior <- DelayedArray::getAutoRealizationBackend()
    DelayedArray::setAutoRealizationBackend("HDF5Array")
    on.exit(DelayedArray::setAutoRealizationBackend(prior), add = TRUE)
    bases <- minfidata_basenames()
    rg <- readMethArray(bases, verbose = 0L)
    # When a backend is set, getGreen()/getRed() return DelayedMatrix-backed
    # storage instead of plain integer matrices. The realize() call inside
    # pull() is the only path that produces that shape.
    expect_s4_class(getGreen(rg), "DelayedMatrix")
})
