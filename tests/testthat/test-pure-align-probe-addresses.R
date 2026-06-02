# .alignProbeAddresses() resolves the mixed-cohort branch of readMethArray
# without needing fabricated IDATs. It takes the worker output list (one
# entry per sample with $n and $addr) plus the master-side reference
# address vector, and decides whether the cohort is uniform (sampleIdx
# NULL, all rows fixed to refAddr) or mixed (sampleIdx populated with
# per-sample integer remaps into the intersection).

test_that("uniform cohort returns NULL sampleIdx and refAddr-aligned commonAddresses", {
    refAddr <- c("a", "b", "c", "d")
    sampleData <- list(
        list(n = 4L, addr = NULL, gMean = 1:4, rMean = 1:4),
        list(n = 4L, addr = NULL, gMean = 5:8, rMean = 5:8))
    aligned <- .alignProbeAddresses(sampleData, refAddr)

    expect_null(aligned$sampleIdx)
    expect_identical(aligned$commonAddresses, refAddr)
})

test_that("mixed-size cohort intersects address vectors and emits per-sample integer remaps", {
    # Two "versions": refAddr (4 probes) and a "smaller version" missing
    # "b". The intersection is {a, c, d}, so each sample's vector must be
    # re-indexed to those positions in its own address order.
    refAddr     <- c("a", "b", "c", "d")
    smallerAddr <- c("a", "c", "d")
    sampleData <- list(
        list(n = 4L, addr = NULL,        gMean = 11:14, rMean = 21:24),
        list(n = 3L, addr = smallerAddr, gMean = 31:33, rMean = 41:43))
    aligned <- .alignProbeAddresses(sampleData, refAddr)

    expect_identical(aligned$commonAddresses, c("a", "c", "d"))
    expect_type(aligned$sampleIdx, "list")
    expect_equal(length(aligned$sampleIdx), 2L)
    # sample 1's row order is refAddr (a,b,c,d): rows 1,3,4 are the intersection
    expect_identical(aligned$sampleIdx[[1L]], c(1L, 3L, 4L))
    # sample 2's row order is smallerAddr (a,c,d): rows 1,2,3 are the intersection
    expect_identical(aligned$sampleIdx[[2L]], c(1L, 2L, 3L))
})

test_that("mixed-size cohort: pulling values via sampleIdx reproduces the address-aligned matrix", {
    refAddr     <- c("a", "b", "c", "d")
    smallerAddr <- c("a", "c", "d")
    sampleData <- list(
        list(n = 4L, addr = NULL,        gMean = c(a = 11, b = 12, c = 13, d = 14)),
        list(n = 3L, addr = smallerAddr, gMean = c(a = 31,         c = 33, d = 34)))
    aligned <- .alignProbeAddresses(sampleData, refAddr)
    # Reproduce the matrix-pull step that readMethArray() performs after
    # alignment: each sample contributes its $gMean[sampleIdx[[j]]] column.
    out <- vapply(seq_along(sampleData),
                  function(j) sampleData[[j]]$gMean[aligned$sampleIdx[[j]]],
                  numeric(length(aligned$commonAddresses)))
    expect_equal(out[, 1L], c(a = 11, c = 13, d = 14))
    expect_equal(out[, 2L], c(a = 31, c = 33, d = 34))
})
