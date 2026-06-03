# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Pure tests for the streaming QC primitives: the per-block sample/probe masks
# (.streamingChunkKeepSamples, .streamingChunkProbeFails) and the final report +
# all-dropped guards (.streamingQCfinalize). Synthetic detection-p blocks and
# mask vectors exercise every branch without reading any IDATs.

test_that(".streamingChunkKeepSamples keeps every column when no threshold", {
    detP <- matrix(0.5, nrow = 4L, ncol = 3L)
    expect_identical(.streamingChunkKeepSamples(detP, NULL), rep(TRUE, 3L))
})

test_that(".streamingChunkKeepSamples drops columns whose mean detP is too high", {
    # column 3's mean (0.9) is above 0.5; columns 1-2 (mean 0.01) survive.
    detP <- cbind(c(0.01, 0.01), c(0.01, 0.01), c(0.9, 0.9))
    expect_identical(.streamingChunkKeepSamples(detP, 0.5),
                     c(TRUE, TRUE, FALSE))
})

test_that(".streamingChunkProbeFails flags probes failing in any column", {
    # probe 2 has detP 0.9 in column 1 -> fails at threshold 0.5; probe 1 clean.
    detP <- rbind(c(0.01, 0.02), c(0.9, 0.01))
    expect_identical(.streamingChunkProbeFails(detP, 0.5), c(FALSE, TRUE))
})

test_that(".streamingAnnoKeep drops probes absent from the annotation", {
    out_loci <- c("cg1", "cg2", "cg3", "cg4")
    # cg3 is not in the annotation -> must be folded out even though QC kept it.
    anno <- c("cg1", "cg2", "cg4", "cg9")
    keepProbes <- c(TRUE, TRUE, TRUE, FALSE)   # cg4 already dropped by QC
    expect_identical(.streamingAnnoKeep(keepProbes, out_loci, anno),
                     c(TRUE, TRUE, FALSE, FALSE))
})

test_that(".streamingAnnoKeep is a no-op when the annotation covers every locus", {
    out_loci <- c("cg1", "cg2", "cg3")
    keepProbes <- c(TRUE, FALSE, TRUE)
    expect_identical(
        .streamingAnnoKeep(keepProbes, out_loci, c("cg1", "cg2", "cg3", "cgX")),
        keepProbes)
})

test_that(".streamingQCfinalize passes silently when nothing is dropped", {
    expect_silent(
        .streamingQCfinalize(c(TRUE, TRUE), c(TRUE, TRUE, TRUE),
                             c("a", "b"), NULL, NULL, verbose = 2L))
})

test_that(".streamingQCfinalize errors when every sample is dropped", {
    expect_error(
        .streamingQCfinalize(c(FALSE, FALSE), c(TRUE, TRUE),
                             c("a", "b"), sample_detP_threshold = 0.1,
                             probe_detP_threshold = NULL, verbose = 0L),
        regexp = "all 2 samples were dropped")
})

test_that(".streamingQCfinalize errors when every probe is dropped", {
    expect_error(
        .streamingQCfinalize(c(TRUE, TRUE), c(FALSE, FALSE, FALSE, FALSE),
                             c("a", "b"), sample_detP_threshold = NULL,
                             probe_detP_threshold = 0.1, verbose = 0L),
        regexp = "all 4 probes were dropped")
})

test_that(".streamingQCfinalize prints the QC summary and dropped-sample list when verbose", {
    msgs <- character(0)
    withCallingHandlers(
        .streamingQCfinalize(c(TRUE, FALSE), c(TRUE, TRUE),
                             c("keepme", "dropme"),
                             sample_detP_threshold = 0.5,
                             probe_detP_threshold = NULL, verbose = 2L),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl("QC drops 1/2 samples", msgs)))
    expect_true(any(grepl("dropped samples: dropme", msgs)))
})
