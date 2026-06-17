# Tests for .validateArgs() -- the inline-args replacement
# for the old file-config parser. Each test builds a tiny synthetic
# fixture in tempdir (samplesheet + empty-but-existent IDAT files +
# cross-reactive CSV) and then calls the function with tweaked args
# to assert an actionable error pointing at the broken value.
#
# Empty IDAT files are deliberate: file.exists() returns TRUE, which
# is all .validateArgs checks. readIDAT() is never
# called here -- that happens later inside processMethArrayExp().

# The pure runner avoids loading minfi (and its annotation packages,
# which transitively pull minfi in). Tests therefore stub the
# annotation check with a function that just returns TRUE -- the real
# annotation behaviour is exercised in the integration suite.
.ann_ok <- function(pkg) TRUE
.ann_no <- function(pkg) FALSE

# Build a complete valid fixture and return paths the caller can twist.
.make_fixture <- function() {
    root <- tempfile("validate_args_")
    dir.create(root)
    data_dir <- file.path(root, "data")
    dir.create(data_dir)

    # Two slide subdirs with one sample each, mirroring the minfiData
    # layout. Empty .idat files are enough: file.exists is the only
    # check this function performs against them.
    slide_a <- file.path(data_dir, "5723646052")
    slide_b <- file.path(data_dir, "5723646053")
    dir.create(slide_a); dir.create(slide_b)
    for (b in c(file.path(slide_a, "5723646052_R02C02"),
                file.path(slide_b, "5723646053_R04C02"))) {
        file.create(paste0(b, "_Grn.idat"), paste0(b, "_Red.idat"))
    }

    # Samplesheet matching the cohort.
    sheet <- data.frame(
        Sample_Name = c("S1", "S2"),
        Basename = c("5723646052/5723646052_R02C02",
                     "5723646053/5723646053_R04C02"),
        stringsAsFactors = FALSE)
    sheet_path <- file.path(data_dir, "samplesheet_real.csv")
    write.csv(sheet, sheet_path, row.names = FALSE)

    # Cross-reactive probes CSV with the required column.
    xreact <- data.frame(TargetID = c("cg00000029", "cg00000108"),
                         stringsAsFactors = FALSE)
    xreact_path <- file.path(root, "xreact.csv")
    write.csv(xreact, xreact_path, row.names = FALSE)

    list(root = root,
         sheet_path = sheet_path,
         data_dir = data_dir,
         xreact_path = xreact_path)
}

# Default args that match the fixture; tests pass `...` to override
# any of them per scenario.
.default_args <- function(fx, ...) {
    base <- list(
        dataDirectory         = fx$data_dir,
        crossReactiveProbes   = fx$xreact_path,
        samplesheet           = fx$sheet_path,
        gdsOutput             = file.path(fx$root, "test.gds"),
        annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        sampleDetPThreshold   = 0.01,
        probeDetPThreshold    = 0.01,
        forceRebuild          = TRUE,
        readerWorkers         = 1L,
        .check_annotation     = .ann_ok)
    overrides <- list(...)
    base[names(overrides)] <- overrides
    base
}

.call_validate <- function(fx, ...) {
    do.call(.validateArgs, .default_args(fx, ...))
}

# ---- happy path ------------------------------------------------------

test_that(".validateArgs accepts a complete valid fixture", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    parsed <- .call_validate(fx)

    expect_named(parsed, c("sheetPath", "sampleData",
                           "basenames", "xReactiveProbes"))
    expect_equal(parsed$sheetPath, fx$sheet_path)
    expect_equal(nrow(parsed$sampleData), 2L)
    expect_length(parsed$basenames, 2L)
    expect_true("TargetID" %in% names(parsed$xReactiveProbes))
})

# ---- type / range checks (one per .require_* used internally) -------

test_that(".validateArgs errors on non-string dataDirectory", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, dataDirectory = 42),
                 regexp = "dataDirectory.*single non-empty string")
})

test_that(".validateArgs errors on out-of-range sampleDetPThreshold", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, sampleDetPThreshold = 1.5),
                 regexp = "sampleDetPThreshold.*single number")
})

test_that(".validateArgs errors on out-of-range probeDetPThreshold", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, probeDetPThreshold = -0.1),
                 regexp = "probeDetPThreshold.*single number")
})

test_that(".validateArgs errors on non-logical forceRebuild", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, forceRebuild = "yes"),
                 regexp = "forceRebuild.*single TRUE or FALSE")
})

test_that(".validateArgs errors on non-positive-integer readerWorkers", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, readerWorkers = 0),
                 regexp = "readerWorkers.*positive integer")
})

# ---- filesystem existence checks ------------------------------------

test_that(".validateArgs errors when dataDirectory does not exist", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(
        .call_validate(fx, dataDirectory = file.path(fx$root, "no_such_data")),
        regexp = "dataDirectory.*existing directory")
})

test_that(".validateArgs errors when crossReactiveProbes does not exist", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(
        .call_validate(fx, crossReactiveProbes = file.path(fx$root, "missing.csv")),
        regexp = "crossReactiveProbes.*existing file")
})

test_that(".validateArgs errors when samplesheet does not point to an existing file", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, samplesheet = file.path(fx$root, "no_such_sheet.csv")),
                 regexp = "samplesheet.*existing file")
})

# ---- gdsOutput destination ------------------------------------------

test_that(".validateArgs errors when gdsOutput's directory does not exist", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    bad <- file.path(fx$root, "no_such_dir", "out.gds")
    err <- tryCatch(.call_validate(fx, gdsOutput = bad), error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "gdsOutput")
    expect_match(conditionMessage(err), "does not exist")
})

test_that(".validateArgs accepts a relative gdsOutput whose parent exists", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    # A bare filename resolves against the working directory (dirname "."),
    # which exists and is writable during the test run, so this must pass.
    expect_silent(.call_validate(fx, gdsOutput = "relative_out.gds"))
})

test_that(".validateArgs errors when gdsOutput's directory is not writable", {
    skip_on_os("windows")  # POSIX permission bits; root ignores them
    skip_if(unname(Sys.info()["user"]) == "root")
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    ro <- file.path(fx$root, "readonly")
    dir.create(ro)
    Sys.chmod(ro, mode = "0500")
    on.exit(Sys.chmod(ro, mode = "0700"), add = TRUE)  # restore so unlink works
    err <- tryCatch(.call_validate(fx, gdsOutput = file.path(ro, "out.gds")),
                    error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "not writable")
})

test_that(".validateArgs errors when an input file is unreadable", {
    skip_on_os("windows")
    skip_if(unname(Sys.info()["user"]) == "root")
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    Sys.chmod(fx$xreact_path, mode = "0000")
    on.exit(Sys.chmod(fx$xreact_path, mode = "0600"), add = TRUE)
    err <- tryCatch(.call_validate(fx), error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "crossReactiveProbes")
    expect_match(conditionMessage(err), "not readable")
})

# ---- annotation package check ---------------------------------------

test_that(".validateArgs errors when the annotation package is not installed", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    expect_error(.call_validate(fx, .check_annotation = .ann_no),
                 regexp = "annotationPackage.*not installed")
})

# ---- samplesheet structure ------------------------------------------

test_that(".validateArgs errors when the samplesheet lacks a Basename column", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    write.csv(data.frame(Sample_Name = "S1"), fx$sheet_path, row.names = FALSE)
    expect_error(.call_validate(fx),
                 regexp = "missing the required `Basename` column")
})

test_that(".validateArgs errors when a Basename entry is empty / NA", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    write.csv(data.frame(Sample_Name = c("S1", "S2"),
                         Basename = c("foo/bar", NA),
                         stringsAsFactors = FALSE),
              fx$sheet_path, row.names = FALSE)
    expect_error(.call_validate(fx),
                 regexp = "empty or NA Basename entries")
})

# ---- IDAT existence -------------------------------------------------

test_that(".validateArgs errors when an IDAT pair is missing", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    unlink(file.path(fx$data_dir, "5723646052",
                     "5723646052_R02C02_Grn.idat"))
    err <- tryCatch(.call_validate(fx), error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "missing IDAT files")
    expect_match(conditionMessage(err), "5723646052_R02C02")
})

test_that(".validateArgs truncates the missing-IDAT list past five samples", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    # Seven samplesheet rows, none with IDATs on disk: n_bad (7) exceeds the
    # five-sample display cap, exercising the '... and N more' arm.
    write.csv(
        data.frame(Sample_Name = paste0("S", 1:7),
                   Basename = sprintf("slide/sample_%d", 1:7),
                   stringsAsFactors = FALSE),
        fx$sheet_path, row.names = FALSE)
    err <- tryCatch(.call_validate(fx), error = identity)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "7 sample\\(s\\)")
    expect_match(conditionMessage(err), "and 2 more")
})

# ---- cross-reactive probes CSV structure ----------------------------

test_that(".validateArgs errors when the cross-reactive CSV lacks TargetID", {
    fx <- .make_fixture()
    on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
    write.csv(data.frame(Probe_ID = c("a", "b")),
              fx$xreact_path, row.names = FALSE)
    expect_error(.call_validate(fx),
                 regexp = "missing the required `TargetID` column")
})
