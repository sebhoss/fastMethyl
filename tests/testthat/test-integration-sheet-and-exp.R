# Coverage for readMethArraySheet() column-detection variants and
# readMethArrayExp() directory-scan branches. readMethArraySheet has
# four column families (Sentrix_*, Array_ID, Slide_ID, Plate_ID, plus
# coercion-only Pool_ID / Sample_Plate / Sample_Well) and two base-mode
# branches (directory of CSVs / list of CSV files). readMethArrayExp
# has a targets path and a directory-scan path; the directory-scan path
# additionally has warnings for orphan Grn / Red files.

skip_if_no_minfidata <- function() {
    skip_if_not_installed("minfiData")
    skip_if_not_installed("BiocParallel")
}

.write_sheet <- function(workdir, sheet, name = "sheet.csv") {
    write.csv(sheet, file.path(workdir, name), row.names = FALSE)
    file.path(workdir, name)
}

.symlink_idats <- function(workdir, basenames) {
    for (b in basenames) {
        for (ch in c("_Grn.idat", "_Red.idat")) {
            file.symlink(paste0(b, ch),
                         file.path(workdir, paste0(basename(b), ch)))
        }
    }
}

test_that("readMethArraySheet detects Array_ID + Slide_ID columns (alternative naming)", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    workdir <- tempfile("sheet_alt_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    .symlink_idats(workdir, bases)

    parts <- do.call(rbind, strsplit(basename(bases), "_"))
    sheet <- data.frame(
        Sample_Name = paste0("S", seq_along(bases)),
        Slide_ID  = parts[, 1L],
        Array_ID  = parts[, 2L],
        Plate_ID  = "P1",
        Pool_ID   = "Pool1",
        Sample_Plate = "Sp1",
        Sample_Well  = "A01",
        stringsAsFactors = FALSE)
    .write_sheet(workdir, sheet)
    out <- readMethArraySheet(base = workdir, pattern = "sheet\\.csv$",
                                recursive = FALSE, verbose = FALSE)
    expect_true(all(c("Array", "Slide", "Plate", "Basename") %in% names(out)))
    expect_true("Pool_ID" %in% names(out))
})

test_that("readMethArraySheet warns when Array or Slide cannot be inferred", {
    skip_if_slow()
    workdir <- tempfile("sheet_nowarn_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    # No Array column source, no Slide column source -- sheet readers
    # should warn for each missing identifier rather than fail outright.
    sheet <- data.frame(
        Sample_Name = c("S1", "S2"),
        Notes = c("a", "b"),
        stringsAsFactors = FALSE)
    .write_sheet(workdir, sheet)
    out <- suppressWarnings(
        readMethArraySheet(base = workdir, pattern = "sheet\\.csv$",
                             recursive = FALSE, verbose = FALSE))
    expect_warning(
        readMethArraySheet(base = workdir, pattern = "sheet\\.csv$",
                             recursive = FALSE, verbose = FALSE),
        regexp = "Could not infer (array|slide) name")
    expect_false("Basename" %in% names(out))
})

test_that("readMethArraySheet(verbose = TRUE) prints the CSV files it finds", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    workdir <- tempfile("sheet_verbose_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    .symlink_idats(workdir, bases)
    parts <- do.call(rbind, strsplit(basename(bases), "_"))
    sheet <- data.frame(
        Sample_Name = paste0("S", seq_along(bases)),
        Sentrix_ID  = parts[, 1L],
        Sentrix_Position = parts[, 2L],
        stringsAsFactors = FALSE)
    .write_sheet(workdir, sheet)
    msgs <- capture_messages(
        readMethArraySheet(base = workdir, pattern = "sheet\\.csv$",
                             recursive = FALSE, verbose = TRUE))
    expect_true(any(grepl("Found the following CSV files", msgs)))
})

test_that("readMethArraySheet errors when base does not exist", {
    skip_if_slow()
    expect_error(
        readMethArraySheet(base = file.path(tempdir(), "definitely_no_dir"),
                             verbose = FALSE),
        regexp = "does not exists")
})

test_that("readMethArraySheet errors when base mixes files and directories", {
    skip_if_slow()
    skip_if_no_minfidata()
    workdir <- tempfile("sheet_mixed_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    sheet <- data.frame(Sample_Name = "S1", stringsAsFactors = FALSE)
    csv_path <- .write_sheet(workdir, sheet, name = "single.csv")
    # base = c(<dir>, <file>) -- info$isdir is mixed.
    expect_error(
        readMethArraySheet(base = c(workdir, csv_path), verbose = FALSE),
        regexp = "needs to be either directories or files")
})

test_that("readMethArraySheet exercises the 'base is files (not directories)' branch", {
    skip_if_slow()
    skip_if_no_minfidata()
    workdir <- tempfile("sheet_filelist_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    sheet <- data.frame(Sample_Name = "S1", stringsAsFactors = FALSE)
    csv_path <- .write_sheet(workdir, sheet, name = "single.csv")
    # When `base` is a file rather than a directory, the upstream code
    # takes the else branch (csvfiles <- list.files(base, ...)). The
    # branch is exercised here for coverage; the output shape is not
    # meaningful in this mode and is not asserted.
    expect_no_error(
        suppressWarnings(readMethArraySheet(base = csv_path, verbose = FALSE)))
})

test_that("readMethArrayExp(base, targets) joins base path with each Basename", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    # Same base dir resolved two different ways: with `base + targets`,
    # the function builds `file.path(base, basename(targets$Basename))`.
    # Use bases pointing at minfiData's extdata directly so the join
    # resolves to real IDAT paths.
    parent_dirs <- unique(dirname(bases))
    expect_length(parent_dirs, 2L)  # minfiData splits IDATs across two slide dirs

    # Use just one slide so the file.path(base, ...) call is unambiguous.
    slide <- parent_dirs[[1L]]
    same_slide_bases <- bases[dirname(bases) == slide]
    targets <- data.frame(
        Sample_Name = paste0("S", seq_along(same_slide_bases)),
        Basename = basename(same_slide_bases),
        stringsAsFactors = FALSE)
    rg <- readMethArrayExp(base = slide, targets = targets, verbose = 0L)
    expect_s4_class(rg, "RGChannelSet")
    expect_equal(ncol(rg), nrow(targets))
})

test_that("readMethArrayExp errors when no IDAT files are found under base", {
    skip_if_slow()
    workdir <- tempfile("exp_empty_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    expect_error(readMethArrayExp(base = workdir, verbose = 0L),
                 regexp = "No IDAT files were found")
})

test_that("readMethArrayExp errors when only single-channel IDATs are present", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    workdir <- tempfile("exp_grn_only_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    # Symlink Grn for one set of basenames and Red for an entirely
    # different set so the intersect comes up empty.
    file.symlink(paste0(bases[[1L]], "_Grn.idat"),
                 file.path(workdir, "a_Grn.idat"))
    file.symlink(paste0(bases[[2L]], "_Red.idat"),
                 file.path(workdir, "b_Red.idat"))
    expect_error(readMethArrayExp(base = workdir, verbose = 0L),
                 regexp = "No IDAT files with both Red and Green")
})

test_that("readMethArrayExp warns when extra Grn or Red files lack a sibling channel", {
    skip_if_slow()
    skip_if_no_minfidata()
    bases <- minfidata_basenames()
    workdir <- tempfile("exp_orphans_")
    dir.create(workdir)
    on.exit(unlink(workdir, recursive = TRUE), add = TRUE)
    .symlink_idats(workdir, bases)
    # Add one orphan Grn (no matching Red) and one orphan Red (no matching Grn).
    file.symlink(paste0(bases[[1L]], "_Grn.idat"),
                 file.path(workdir, "orphan_grn_Grn.idat"))
    file.symlink(paste0(bases[[1L]], "_Red.idat"),
                 file.path(workdir, "orphan_red_Red.idat"))
    expect_warning(
        expect_warning(readMethArrayExp(base = workdir, verbose = 0L),
                       regexp = "green channel"),
        regexp = "red channel")
})
