# Integration coverage for the streaming reader paths the rest of the
# suite leaves untouched:
#   * processMethArrayExp()        -- the experiment wrapper
#   * runPreprocess() end-to-end   -- the rebuild branch
#
# All need real minfiData IDATs plus the bigmelon/gdsfmt stack + the 450k
# annotation package, so they are gated behind MINFI_RUN_INTEGRATION like the
# rest of the integration tier.

# Relative "<slide>/<array>" basenames under the minfiData extdata dir.
.minfidata_rel <- function(basenames) {
    base <- system.file("extdata", package = "minfiData")
    sub("^/", "", substring(basenames, nchar(base) + 1L))
}

# ---- processMethArrayExp() -----------------------------------------

test_that("processMethArrayExp drops sex-chromosome probes and filters targets", {
    skip_if_slow()
    skip_if_no_gds()
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    attach_gds_pkgs()
    base <- system.file("extdata", package = "minfiData")
    basenames <- minfidata_basenames()
    targets <- data.frame(Sample_Name = paste0("S", seq_along(basenames)),
                          Basename = .minfidata_rel(basenames),
                          stringsAsFactors = FALSE)
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)

    res <- processMethArrayExp(
        base = base, targets = targets, gds_path = gds_path,
        annotation_package = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        verbose = 0L)

    expect_true(file.exists(res$gds_path))
    expect_equal(nrow(res$targets), sum(res$keepSamples))
    # Sex-chromosome probes were removed, so the GDS holds fewer loci than
    # an unfiltered run on the same IDATs.
    full <- .run_gds(basenames)
    got <- .read_gds(gds_path)
    expect_lt(length(got$features), length(full$features))
})

test_that("processMethArrayExp accepts base = NULL with pre-resolved Basenames", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    # base = NULL -> targets$Basename is used verbatim as full paths.
    targets <- data.frame(Basename = basenames, stringsAsFactors = FALSE)
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)

    res <- processMethArrayExp(
        base = NULL, targets = targets, gds_path = gds_path,
        annotation_package = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        drop_sex_chromosomes = FALSE, verbose = 0L)

    expect_true(file.exists(res$gds_path))
    expect_equal(nrow(res$targets), length(basenames))
})

test_that("processMethArrayExp with drop_sex_chromosomes = FALSE keeps every locus", {
    skip_if_slow()
    skip_if_no_gds()
    attach_gds_pkgs()
    base <- system.file("extdata", package = "minfiData")
    basenames <- minfidata_basenames()
    targets <- data.frame(Basename = .minfidata_rel(basenames),
                          stringsAsFactors = FALSE)
    gds_path <- tempfile(fileext = ".gds")
    on.exit(unlink(gds_path), add = TRUE)

    res <- processMethArrayExp(
        base = base, targets = targets, gds_path = gds_path,
        annotation_package = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        drop_sex_chromosomes = FALSE, verbose = 0L)

    got <- .read_gds(gds_path)
    full <- .run_gds(basenames)
    expect_identical(length(got$features), length(full$features))
})

# ---- runPreprocess() end-to-end --------------------------------

test_that("runPreprocess rebuilds a GDS from a samplesheet end-to-end", {
    skip_if_slow()
    skip_if_no_gds()
    skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    attach_gds_pkgs()
    basenames <- minfidata_basenames()
    rel <- .minfidata_rel(basenames)

    root <- tempfile("pp_e2e_")
    data_dir <- file.path(root, "data")
    dir.create(data_dir, recursive = TRUE)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    # Symlink the real IDAT pairs into the fixture's <slide>/<array> layout.
    for (i in seq_along(basenames)) {
        dst <- file.path(data_dir, rel[i])
        dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
        file.symlink(paste0(basenames[i], "_Grn.idat"), paste0(dst, "_Grn.idat"))
        file.symlink(paste0(basenames[i], "_Red.idat"), paste0(dst, "_Red.idat"))
    }
    write.csv(data.frame(Sample_Name = paste0("S", seq_along(rel)),
                         Basename = rel, stringsAsFactors = FALSE),
              file.path(data_dir, "samplesheet_e2e.csv"), row.names = FALSE)
    xreact <- file.path(root, "xreact.csv")
    write.csv(data.frame(TargetID = "cg00000029", stringsAsFactors = FALSE),
              xreact, row.names = FALSE)

    res <- runPreprocess(
        dataDirectory         = data_dir,
        crossReactiveProbes   = xreact,
        samplesheet           = file.path(data_dir, "samplesheet_e2e.csv"),
        gdsOutput             = file.path(root, "cohort.gds"),
        annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        forceRebuild          = TRUE,
        BPPARAM               = BiocParallel::SerialParam(),
        verbose               = 0L)

    expect_equal(res$gds_path, file.path(root, "cohort.gds"))
    expect_true(file.exists(res$gds_path))
    expect_equal(nrow(res$targets), sum(res$keepSamples))
})
