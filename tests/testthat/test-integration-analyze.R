# End-to-end analyze() on minfiData via its production defaults
# (runPreprocess + openfn.gds + dasen + closefn.gds), which the pure tests
# stub out. Also asserts the GDS is genuinely closed when analyze returns.

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
        nonSpecificProbesPath = xreact,
        targetPattern         = "test",
        datasetClass          = file.path(fixture, "ana"),
        annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        verbose               = 0L,
        FUN = function(gds, res) {
            list(nodes = ls.gdsn(gds),
                 nsamp = objdesp.gdsn(index.gdsn(gds, "betas"))$dim[[2]],
                 targets = res$targets)
        })

    expect_true("normbetas" %in% out$nodes)        # normalize = TRUE ran dasen
    expect_equal(out$nsamp, length(bns))
    expect_equal(nrow(out$targets), length(bns))

    # The GDS must be closed: re-opening it must succeed (an open handle would
    # make openfn.gds error).
    g2 <- openfn.gds(file.path(fixture, "ana.gds"), readonly = TRUE)
    expect_s3_class(g2, "gds.class")
    closefn.gds(g2)
})
