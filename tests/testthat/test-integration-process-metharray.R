# Tests for the fork's processMethArray() — the fused-pipeline reader
# that does readMethArray + detectionP + preprocessRaw + es2gds in one
# streaming pass straight to a GDS.
#
# The unfiltered GDS run doubles as the equivalence target: its
# methylated/unmethylated/pvals nodes must equal preprocessRaw(readMethArray(...))
# + detectionP(...) on the same IDATs. The QC arg paths (drop_probes,
# sample_detP_threshold, probe_detP_threshold) are checked against post-hoc
# subsets of that unfiltered GDS so any divergence means the in-function
# shortcut diverged from the manual equivalent.

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

# Cache the heavy fixtures (the upstream MethylSet/detectionP reference and
# a no-QC GDS run read back into plain matrices) at file scope so the
# dozen-odd tests below do not each pay for the ~1 s read + assembly.
# testthat::local_* would cap the cache at one test_that block; this caches
# across the whole file.
.fork_fixtures <- new.env(parent = emptyenv())
.fixtures <- function() {
    if (!is.null(.fork_fixtures$ready)) return(.fork_fixtures)
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    rgSet <- readMethArray(basenames, verbose = 0L)
    .fork_fixtures$basenames <- basenames
    .fork_fixtures$mSet_ref  <- preprocessRaw(rgSet)
    .fork_fixtures$detP_ref  <- detectionP(rgSet)
    .fork_fixtures$full      <- .run_gds(basenames)
    .fork_fixtures$ready <- TRUE
    .fork_fixtures
}

test_that("GDS output is byte-identical to preprocessRaw(readMethArray()) + detectionP()", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()

    expect_identical(f$full$meth,     unname(getMeth(f$mSet_ref)))
    expect_identical(f$full$unmeth,   unname(getUnmeth(f$mSet_ref)))
    expect_identical(f$full$detP,     unname(f$detP_ref))
    expect_identical(f$full$features, featureNames(f$mSet_ref))
    expect_identical(f$full$samples,  colnames(f$mSet_ref))
})

test_that("the streamed GDS has the expected pipeline nodes", {
    skip_if_slow()
    skip_if_no_minfidata()
    skip_if_not_installed("bigmelon")
    skip_if_not_installed("gdsfmt")
    # The worked pipeline.r attaches bigmelon + gdsfmt; the node assertions
    # below read the GDS back through gdsfmt.
    suppressPackageStartupMessages({
        library(bigmelon)
        library(gdsfmt)
    })
    basenames <- minfidata_basenames()

    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)
    res <- processMethArray(basenames, gds_path = gds_path, verbose = 0L)

    expect_true(file.exists(gds_path))
    gf <- gdsfmt::openfn.gds(gds_path)
    on.exit(gdsfmt::closefn.gds(gf), add = TRUE)
    node_names <- gdsfmt::ls.gdsn(gf)
    # bigmelon-compatible pipeline nodes; the raw signal layer
    # (methylated / unmethylated / pvals) plus the derived `betas` and
    # the manifest / phenotype / history scaffolding.
    for (expected in c("betas", "methylated", "unmethylated", "pvals",
                       "fData", "pData", "history")) {
        expect_true(expected %in% node_names,
                    info = sprintf("GDS missing node '%s'", expected))
    }
})

test_that("drop_probes shrinks output to (loci minus dropped) and preserves values", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()

    set.seed(1L)
    drop <- sample(f$full$features,
                   size = floor(length(f$full$features) * 0.05))
    dropped <- .run_gds(f$basenames, drop_probes = drop)

    kept <- setdiff(f$full$features, drop)
    expect_equal(length(dropped$features), length(kept))
    expect_identical(dropped$features, f$full$features[f$full$features %in% kept])
    sub_idx <- match(dropped$features, f$full$features)
    expect_identical(dropped$meth,   f$full$meth[sub_idx, , drop = FALSE])
    expect_identical(dropped$unmeth, f$full$unmeth[sub_idx, , drop = FALSE])
    expect_identical(dropped$detP,   f$full$detP[sub_idx, , drop = FALSE])
})

test_that("probe_detP_threshold drops only rows whose detP ever crosses the threshold", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()
    probe_thr <- 0.01

    pf <- .run_gds(f$basenames, probe_detP_threshold = probe_thr)
    keep_p <- which(rowSums(f$full$detP >= probe_thr) == 0L)
    expect_equal(length(pf$features), length(keep_p))
    expect_identical(pf$features, f$full$features[keep_p])
    expect_identical(pf$meth, f$full$meth[keep_p, , drop = FALSE])
})

test_that("probe-QC row compaction preserves every matrix node (betas included)", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()
    probe_thr <- 0.01
    # Write an unfiltered GDS and a probe-QC'd GDS (the latter exercises the
    # post-stream row-compaction rewrite), then read all four matrix nodes back
    # directly so the comparison includes betas, not just meth/unmeth/detP.
    full_path <- tempfile(fileext = ".gds")
    qc_path <- tempfile(fileext = ".gds")
    on.exit(unlink(c(full_path, qc_path)), add = TRUE)
    processMethArray(f$basenames, gds_path = full_path, verbose = 0L)
    processMethArray(f$basenames, gds_path = qc_path,
                     probe_detP_threshold = probe_thr, verbose = 0L)
    rd <- function(path, node) {
        g <- gdsfmt::openfn.gds(path)
        on.exit(gdsfmt::closefn.gds(g), add = TRUE)
        gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, node))
    }
    keep_p <- which(rowSums(f$full$detP >= probe_thr) == 0L)
    expect_lt(length(keep_p), length(f$full$features))   # compaction did fire
    for (node in c("betas", "methylated", "unmethylated", "pvals")) {
        full_m <- rd(full_path, node)
        qc_m <- rd(qc_path, node)
        expect_identical(qc_m, full_m[keep_p, , drop = FALSE])
    }
})

test_that("sample_detP_threshold drops columns whose colMeans(detP) clears the threshold", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()
    cm <- colMeans(f$full$detP)
    if (length(unique(cm)) == 1L) {
        skip("All samples share identical colMeans(detP); no discriminating threshold")
    }
    sample_thr <- median(cm)
    sf <- .run_gds(f$basenames, sample_detP_threshold = sample_thr)
    kept_samples <- which(cm < sample_thr)
    expect_equal(length(sf$samples), length(kept_samples))
    expect_identical(sf$meth, f$full$meth[, kept_samples, drop = FALSE])
    expect_identical(sf$detP, f$full$detP[, kept_samples, drop = FALSE])
})

test_that("drop_probes + probe_detP_threshold combine to a single manual mask", {
    skip_if_slow()
    skip_if_no_gds()
    f <- .fixtures()
    set.seed(2L)
    probe_thr <- 0.01
    drop <- sample(f$full$features,
                   size = floor(length(f$full$features) * 0.05))
    both <- .run_gds(f$basenames, drop_probes = drop,
                     probe_detP_threshold = probe_thr)
    keep_mask <- (!(f$full$features %in% drop)) &
        (rowSums(f$full$detP >= probe_thr) == 0L)
    expect_equal(length(both$features), sum(keep_mask))
    expect_identical(both$meth,   f$full$meth[keep_mask, , drop = FALSE])
    expect_identical(both$unmeth, f$full$unmeth[keep_mask, , drop = FALSE])
})

test_that("processMethArray refuses to run without a gds_path", {
    skip_if_slow()
    expect_error(processMethArray(character(0)),
                 regexp = "gds_path")
})

test_that("verbose argument rejects out-of-range integers", {
    skip_if_slow()
    expect_error(processMethArray(character(0), verbose = 3L),
                 regexp = "verbose")
})

test_that("multi-version helper resolves to single-version on uniform input and intersects on mixed input", {
    skip_if_slow()
    # .streamingBuildIndices is the internal that drives mixed-array
    # cohort dispatch. The minfiData fixtures are uniform 450k so the
    # multi-version path is exercised here via a synthetic "smaller
    # version" that drops 1% of addresses from the same manifest --
    # this is the same construction dev/audit_streaming.R uses.
    skip_if_not_installed("minfiData")
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")

    # Use 450k's refN so .streamingResolveManifest picks the installed
    # 450k manifest deterministically.
    refN_450k <- 622399L
    mInfo <- fastMethyl:::.streamingResolveManifest(refN_450k)
    manifest <- mInfo$manifest

    full_addrs <- as.character(sort(unique(c(
        getProbeInfo(manifest, type = "I-Red")$AddressA,
        getProbeInfo(manifest, type = "I-Red")$AddressB,
        getProbeInfo(manifest, type = "I-Green")$AddressA,
        getProbeInfo(manifest, type = "I-Green")$AddressB,
        getProbeInfo(manifest, type = "II")$AddressA,
        getControlAddress(manifest, controlType = "NEGATIVE")))))
    set.seed(7L)
    dropped <- sample(full_addrs, size = floor(length(full_addrs) * 0.01))
    smaller <- setdiff(full_addrs, dropped)

    idx_single <- fastMethyl:::.streamingBuildIndices(manifest,
        list(full_addrs), drop_probes = NULL)
    locusNames <- getManifestInfo(manifest, "locusNames")
    expect_equal(length(idx_single$out_loci), length(locusNames))
    expect_equal(idx_single$n_versions, 1L)

    idx_multi <- fastMethyl:::.streamingBuildIndices(manifest,
        list(full_addrs, smaller), drop_probes = NULL)
    expect_equal(idx_multi$n_versions, 2L)
    expect_true(all(idx_multi$out_loci %in% idx_single$out_loci))
    expect_lt(length(idx_multi$out_loci), length(idx_single$out_loci))
    for (per in idx_multi$per_version) {
        expect_false(anyNA(unlist(per[setdiff(names(per), "refN")])))
    }
})

test_that(".streamingBuildIndices drops probes whose address is absent from a single version", {
    # Regression for real-world IDATs (e.g. GSE261142) that carry fewer
    # addresses than the installed manifest references: the single-version path
    # must filter probes to the addresses actually present, not assume the IDAT
    # holds every manifest address. The old fast path skipped the membership
    # test, so an absent address became an NA index and tripped the internal
    # "address present in intersection but match failed" invariant.
    skip_if_slow()
    skip_if_not_installed("minfiData")
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
    manifest <- fastMethyl:::.streamingResolveManifest(622399L)$manifest
    full_addrs <- as.character(sort(unique(c(
        getProbeInfo(manifest, type = "I-Red")$AddressA,
        getProbeInfo(manifest, type = "I-Red")$AddressB,
        getProbeInfo(manifest, type = "I-Green")$AddressA,
        getProbeInfo(manifest, type = "I-Green")$AddressB,
        getProbeInfo(manifest, type = "II")$AddressA,
        getControlAddress(manifest, controlType = "NEGATIVE")))))

    full <- fastMethyl:::.streamingBuildIndices(manifest, list(full_addrs),
                                                drop_probes = NULL)
    # A single version that is *missing* some manifest addresses, exactly as a
    # real IDAT lacking a few hundred of them would be.
    set.seed(11L)
    partial <- setdiff(full_addrs, sample(full_addrs, 300L))
    expect_error(
        idx <- fastMethyl:::.streamingBuildIndices(manifest, list(partial),
                                                   drop_probes = NULL),
        NA)   # must NOT error (the old fast path did)
    expect_equal(idx$n_versions, 1L)
    # The probes referencing the absent addresses are dropped, not kept as NA.
    expect_lt(length(idx$out_loci), length(full$out_loci))
    expect_true(all(idx$out_loci %in% full$out_loci))
    for (per in idx$per_version) {
        expect_false(anyNA(unlist(per[setdiff(names(per), "refN")])))
    }
})
