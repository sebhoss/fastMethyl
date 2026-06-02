# verbose = 0L is the silent default; 1L surfaces phase summaries +
# per-sample worker reads; 2L additionally surfaces master-side progress
# (every chunk completes), array-resolved messages, memory footprints,
# and QC drop details. These tests pin that each level emits the
# messages the example pipeline relies on for long-running runs (the
# user passes verbose = 2L for big cohorts so the terminal does not look
# frozen). They also exercise the per-chunk dispatch in
# .bpmapply_with_progress, which only runs at verbose >= 2.

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

test_that("readMethArray(verbose = 1L) prints worker reads + phase summary", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()
    msgs <- capture_messages(readMethArray(basenames, verbose = 1L))
    expect_true(any(grepl("Reading .+_Grn\\.idat", msgs)))
    expect_true(any(grepl("Read idat files in", msgs)))
    expect_true(any(grepl("Creating data matrices", msgs)))
    expect_true(any(grepl("Instantiating final object", msgs)))
})

test_that("readMethArray(verbose = 2L) additionally prints worker count + Array resolution + chunk progress", {
    skip_if_slow()
    skip_if_no_minfidata()
    basenames <- minfidata_basenames()
    msgs <- capture_messages(readMethArray(basenames, verbose = 2L))
    expect_true(any(grepl("IDAT pair\\(s\\) found", msgs)))
    expect_true(any(grepl("Reference probe count", msgs)))
    expect_true(any(grepl("Array: IlluminaHumanMethylation450k", msgs)))
    # The chunked-progress dispatch in .bpmapply_with_progress
    expect_true(any(grepl("\\[readMethArray\\] progress", msgs)))
})

test_that("processMethArray(verbose = 1L) prints phase summaries", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    msgs <- capture_messages(
        processMethArray(basenames, gds_path = gds_path, verbose = 1L))
    expect_true(any(grepl("Reading .+_Grn\\.idat", msgs)))
    expect_true(any(grepl("Streaming .+ to GDS", msgs)))
    expect_true(any(grepl("Streamed in", msgs)))
})

test_that("processMethArray(verbose = 2L) prints header scan, array, progress, streaming, and GDS write", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    msgs <- capture_messages(
        processMethArray(basenames, gds_path = gds_path, verbose = 2L))
    expect_true(any(grepl("scanning IDAT headers", msgs)))
    expect_true(any(grepl("Array resolved", msgs)))
    expect_true(any(grepl("IDAT pair\\(s\\) found", msgs)))
    expect_true(any(grepl("\\[processMethArray\\] progress", msgs)))
    expect_true(any(grepl("Streaming .+ to GDS", msgs)))
    expect_true(any(grepl("building annotation", msgs)))
    expect_true(any(grepl("annotation built", msgs)))
    expect_true(any(grepl("GDS file size", msgs)))
})

test_that("processMethArray(verbose = 1L) prints loci-kept message when drop_probes shrinks the manifest", {
    skip_if_slow()
    skip_if_no_gds()
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    # Drop the first 1000 manifest loci -- triggers the "loci kept after
    # drop_probes" message. The locus universe comes from the resolved 450k
    # manifest, so no reference read is needed.
    manifest <- fastMethyl:::.streamingResolveManifest(622399L)$manifest
    drop <- getManifestInfo(manifest, "locusNames")[seq_len(1000L)]
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    msgs <- capture_messages(
        processMethArray(basenames, gds_path = gds_path,
                          drop_probes = drop, verbose = 1L))
    expect_true(any(grepl("loci kept after drop_probes", msgs)))
})

test_that("processMethArray(verbose = 1L) prints QC-drops message when probe_detP_threshold drops probes", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    msgs <- capture_messages(
        processMethArray(basenames, gds_path = gds_path,
                          probe_detP_threshold = 0.001,
                          verbose = 1L))
    expect_true(any(grepl("QC drops .+ probes", msgs)))
})
