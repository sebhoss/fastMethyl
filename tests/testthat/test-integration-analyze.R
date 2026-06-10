# End-to-end analyze() on minfiData via its production defaults
# (runPreprocess + openfn.gds + closefn.gds), which the pure tests stub out.
# Normalisation runs inside FUN (dasen), and the test asserts the GDS is
# genuinely closed when analyze returns.

test_that("analyze runs end-to-end with real defaults and closes the GDS", {
    skip_if_slow()
    skip_if_no_gds()
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    attach_gds_pkgs()

    bns <- minfidata_basenames()
    base <- system.file("extdata", package = "minfiData")
    rel <- sub(paste0(base, "/"), "", bns, fixed = TRUE)   # "<slide>/<array>"

    fixture <- tempfile("analyze_e2e_")
    dir.create(fixture)
    on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
    # Symlink the minfiData IDAT pairs into the fixture's <slide>/<array> layout
    # so runPreprocess resolves them from the samplesheet Basename.
    for (i in seq_along(bns)) {
        dir.create(file.path(fixture, dirname(rel[i])), recursive = TRUE,
                   showWarnings = FALSE)
        for (ch in c("_Grn.idat", "_Red.idat")) {
            file.symlink(paste0(bns[i], ch),
                         file.path(fixture, paste0(rel[i], ch)))
        }
    }
    write.csv(
        data.frame(Sample_Name = paste0("S", seq_along(rel)),
                   Basename = rel, stringsAsFactors = FALSE),
        file.path(fixture, "samplesheet_test.csv"), row.names = FALSE)
    xreact <- file.path(fixture, "xreact.csv")
    write.csv(data.frame(TargetID = "cg00000029", stringsAsFactors = FALSE),
              xreact, row.names = FALSE)

    out <- analyze(
        dataDirectory         = fixture,
        crossReactiveProbes   = xreact,
        samplesheet           = file.path(fixture, "samplesheet_test.csv"),
        gdsOutput             = file.path(fixture, "ana.gds"),
        annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        verbose               = 0L,
        analysis = function(gds, res) {
            dasen(gds, node = "normbetas")   # normalisation lives in FUN now
            list(nodes = ls.gdsn(gds),
                 nsamp = objdesp.gdsn(index.gdsn(gds, "betas"))$dim[[2]],
                 targets = res$targets)
        })

    expect_true("normbetas" %in% out$nodes)        # FUN's dasen created it
    expect_equal(out$nsamp, length(bns))
    expect_equal(nrow(out$targets), length(bns))

    # The GDS must be closed: re-opening it must succeed (an open handle would
    # make openfn.gds error).
    g2 <- openfn.gds(file.path(fixture, "ana.gds"), readonly = TRUE)
    expect_s3_class(g2, "gds.class")
    closefn.gds(g2)
})

# The hook orchestration on real data, via the BUILT-IN method names: analyze
# resolves outlierRemoval = "outlyx" and normalization = "dasen" (gdsDasen) to real
# bigmelon/fastMethyl calls, runs them in order, threads the keep mask, and the
# normalised node is usable downstream (prcomp). The pure tests cover the wiring
# with custom-function spies; this pins the built-in resolution end to end.
test_that("analyze drives the outlyx + dasen hooks end-to-end on real data", {
    skip_if_slow()
    skip_if_no_gds()
    skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    attach_gds_pkgs()

    bns <- minfidata_basenames()
    base <- system.file("extdata", package = "minfiData")
    rel <- sub(paste0(base, "/"), "", bns, fixed = TRUE)

    fixture <- tempfile("analyze_hooks_")
    dir.create(fixture)
    on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
    for (i in seq_along(bns)) {
        dir.create(file.path(fixture, dirname(rel[i])), recursive = TRUE,
                   showWarnings = FALSE)
        for (ch in c("_Grn.idat", "_Red.idat")) {
            file.symlink(paste0(bns[i], ch),
                         file.path(fixture, paste0(rel[i], ch)))
        }
    }
    write.csv(
        data.frame(Sample_Name = paste0("S", seq_along(rel)),
                   Basename = rel, stringsAsFactors = FALSE),
        file.path(fixture, "samplesheet_test.csv"), row.names = FALSE)
    xreact <- file.path(fixture, "xreact.csv")
    write.csv(data.frame(TargetID = "cg00000029", stringsAsFactors = FALSE),
              xreact, row.names = FALSE)

    out <- analyze(
        dataDirectory         = fixture,
        crossReactiveProbes   = xreact,
        samplesheet           = file.path(fixture, "samplesheet_test.csv"),
        gdsOutput             = file.path(fixture, "hooks.gds"),
        annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        verbose               = 0L,
        outlierRemoval        = "outlyx",   # built-in -> bigmelon::outlyx
        normalization          = "dasen",    # built-in -> fastMethyl::gdsDasen
        analysis = function(gds, res) {
            # normalize ran before FUN: normbetas exists and is rectangular.
            stopifnot("normbetas" %in% ls.gdsn(gds))
            list(
                hasNorm = "normbetas" %in% ls.gdsn(gds),
                nbDim   = objdesp.gdsn(index.gdsn(gds, "normbetas"))$dim,
                keepLen = length(res$outlierKeep),
                npcs    = length(prcomp(gds, node.name = "normbetas",
                                        method = "quick")$sdev))
        })

    expect_true(out$hasNorm)                       # the dasen hook created it
    expect_equal(out$nbDim[[2]], length(bns))      # one normbetas column per sample
    expect_equal(out$keepLen, length(bns))         # outlyx keep mask, one per sample
    expect_gt(out$npcs, 0L)                        # prcomp consumes the node

    g2 <- openfn.gds(file.path(fixture, "hooks.gds"), readonly = TRUE)
    expect_s3_class(g2, "gds.class")
    closefn.gds(g2)
})
