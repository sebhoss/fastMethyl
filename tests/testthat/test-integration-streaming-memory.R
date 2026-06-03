# Regression tests for the streaming write's chunking + per-block memory
# release. processMethArray() processes samples in column blocks, frees each
# block's matrices (and the manifest / annotation scaffolding) before the next
# block, and appends to the GDS so a full cohort-sized matrix is never held.
#
# The dramatic failure mode -- forked MulticoreParam workers copy-on-write-
# multiplying retained per-block matrices into an OOM on a large real cohort --
# is fork- and machine-specific and not robustly assertable in CI. These tests
# instead pin the two deterministic properties that the fix must preserve:
#   1. output is identical no matter how the cohort is split into blocks
#      (also the only coverage of the parallel/fork code path), and
#   2. peak heap does not grow with cohort size -- a regression that reverted to
#      assembling whole-cohort matrices would scale peak heap with sample count.

# Symlink the minfiData IDAT pairs `reps` times into `dir`, returning the
# replicated basenames. Replication is by symlink, so it costs no disk.
.replicate_minfidata <- function(dir, reps) {
    src <- minfidata_basenames()
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    stems <- character(0)
    for (r in seq_len(reps)) {
        for (b in src) {
            stem <- file.path(dir, sprintf("rep%03d_%s", r, basename(b)))
            for (ch in c("_Grn.idat", "_Red.idat")) {
                file.symlink(normalizePath(paste0(b, ch)), paste0(stem, ch))
            }
            stems <- c(stems, stem)
        }
    }
    stems
}

# Peak R vector heap (bytes) reached while evaluating `expr`, measured from a
# reset high-water mark so only `expr`'s allocations count.
.peak_heap_bytes <- function(expr) {
    gc(reset = TRUE)
    force(expr)
    unname(gc()["Vcells", "max used"]) * 8
}

test_that("streamed GDS is identical regardless of how the cohort is blocked", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    basenames <- minfidata_basenames()        # 6 samples
    g_serial <- tempfile(fileext = ".gds")
    g_par <- tempfile(fileext = ".gds")
    on.exit(unlink(c(g_serial, g_par)), add = TRUE)

    # SerialParam: block_size = min(6, 1*4, cap) = 4 -> 2 blocks.
    # MulticoreParam(3): block_size = min(6, 3*4, cap) = 6 -> 1 block, and the
    # workers actually fork (on Linux), exercising the parallel path the rest of
    # the suite leaves on SerialParam.
    processMethArray(basenames, gds_path = g_serial, verbose = 0L,
                     BPPARAM = BiocParallel::SerialParam())
    processMethArray(basenames, gds_path = g_par, verbose = 0L,
                     BPPARAM = BiocParallel::MulticoreParam(3))

    s <- .read_gds(g_serial)
    p <- .read_gds(g_par)
    expect_identical(s$meth, p$meth)
    expect_identical(s$unmeth, p$unmeth)
    expect_identical(s$detP, p$detP)
    expect_identical(s$features, p$features)
    expect_identical(s$samples, p$samples)
})

test_that("streaming peak heap does not grow with cohort size", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    small_dir <- tempfile("mem_small_")
    big_dir <- tempfile("mem_big_")
    g_small <- tempfile(fileext = ".gds")
    g_big <- tempfile(fileext = ".gds")
    on.exit(unlink(c(small_dir, big_dir, g_small, g_big), recursive = TRUE),
            add = TRUE)
    small <- .replicate_minfidata(small_dir, 5L)     # 30 samples
    big <- .replicate_minfidata(big_dir, 20L)        # 120 samples

    peak_small <- .peak_heap_bytes(
        processMethArray(small, gds_path = g_small, verbose = 0L,
                         BPPARAM = BiocParallel::SerialParam()))
    peak_big <- .peak_heap_bytes(
        processMethArray(big, gds_path = g_big, verbose = 0L,
                         BPPARAM = BiocParallel::SerialParam()))

    nLoci <- length(.read_gds(g_small)$features)
    extra_samples <- length(big) - length(small)     # 90
    # If the stream were holding whole-cohort matrices, the 90 extra samples
    # would add at least one cohort-sized double matrix (90 x nLoci x 8 B) to the
    # peak. Block-bounded streaming keeps the growth far below that -- the fixed
    # annotation-build cost is identical in both runs and cancels in the
    # difference, leaving only per-sample scaling. (Generous bound: a true
    # assemble-everything regression grows by several such matrices.)
    leak_floor <- as.numeric(extra_samples) * nLoci * 8
    expect_lt(peak_big - peak_small, leak_floor)
})
