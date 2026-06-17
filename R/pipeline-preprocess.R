# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Pipeline-scaffolding function: take the inputs every methylation
# pre-processing pipeline needs as plain named arguments, validate
# every input we can check without touching IDAT bytes, run the fused
# reader, and return a result object the caller's analysis code plugs
# into.
#
# .validateArgs (internal) does all validation and the
# preliminary file I/O (samplesheet load, x-reactive CSV load, IDAT
# file.exists checks) without ever calling processMethArrayExp.
# Tests target it directly with synthetic fixtures so the full
# validation surface can be exercised in <2 s without loading minfi.
# runPreprocess (exported) composes .validateArgs +
# processMethArrayExp.
#
# Type-and-range checks use base R's `stopifnot()` with named
# expressions (idiomatic since R 4.0): the name becomes the error
# message verbatim. Value-bearing errors (path missing, samplesheet
# malformed, etc.) use inline `stop(sprintf(...))` so the offending
# value appears in the message.

# Internal: an existing input path must also be readable, not merely present.
# file.access(mode = 4) tests read permission; a non-zero result means the
# current user cannot read it (e.g. a samplesheet owned by another user).
.assertReadable <- function(path, label) {
  if (file.access(path, mode = 4L) != 0L) {
    stop(
      sprintf("`%s` (\"%s\") exists but is not readable by the current user.",
              label, path),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Internal: validate that `path` is a usable destination for a GDS we will
# create (and overwrite on a rebuild). The only real failure modes are a parent
# directory that does not exist or is not writable, so check both up front and
# turn them into an actionable error rather than a cryptic gdsfmt write failure
# raised deep inside the build, after IDATs have already been read.
#
# Relative paths are fully supported: gdsfmt resolves them against the working
# directory like any other file, so a relative `path` is fine as long as its
# parent exists. The error names the *resolved absolute* directory so a relative
# path makes clear exactly where the write would have landed.
.assertWritableTarget <- function(path, label) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop(
      sprintf(
        paste0(
          "`%s` (\"%s\") points into a directory that does not exist:\n  %s\n",
          "  Create the directory first, or pass a path under one that exists ",
          "(an absolute path removes any ambiguity about the working directory)."
        ),
        label, path, normalizePath(parent, mustWork = FALSE)
      ),
      call. = FALSE
    )
  }
  if (file.access(parent, mode = 2L) != 0L) {
    stop(
      sprintf(
        paste0("`%s` (\"%s\"): the directory %s is not writable, so the output ",
               "file cannot be created there."),
        label, path, normalizePath(parent)
      ),
      call. = FALSE
    )
  }
  if (file.exists(path) && file.access(path, mode = 2L) != 0L) {
    stop(
      sprintf(
        paste0("`%s` (\"%s\") already exists but is not writable, so it cannot ",
               "be overwritten."),
        label, path
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Internal: validate every pre-processing input + load the supporting
# CSV files. Returns a list with the loaded samplesheet, resolved
# basenames, sheet path, and the cross-reactive probes CSV. Does NOT
# read any IDAT bytes -- it only checks that the files exist on disk.
.validateArgs <- function(
  dataDirectory,
  crossReactiveProbes,
  samplesheet,
  gdsOutput,
  annotationPackage,
  sampleDetPThreshold,
  probeDetPThreshold,
  forceRebuild,
  readerWorkers,
  .check_annotation = function(pkg) requireNamespace(pkg, quietly = TRUE)
) {
  # --- types + ranges -----------------------------------------------
  stopifnot(
    "`dataDirectory` must be a single non-empty string" =
      is.character(dataDirectory) && length(dataDirectory) == 1L &&
      !is.na(dataDirectory) && nzchar(dataDirectory),
    "`crossReactiveProbes` must be a single non-empty string" =
      is.character(crossReactiveProbes) &&
      length(crossReactiveProbes) == 1L &&
      !is.na(crossReactiveProbes) && nzchar(crossReactiveProbes),
    "`samplesheet` must be a single non-empty string" =
      is.character(samplesheet) && length(samplesheet) == 1L &&
      !is.na(samplesheet) && nzchar(samplesheet),
    "`gdsOutput` must be a single non-empty string" =
      is.character(gdsOutput) && length(gdsOutput) == 1L &&
      !is.na(gdsOutput) && nzchar(gdsOutput),
    "`annotationPackage` must be a single non-empty string" =
      is.character(annotationPackage) &&
      length(annotationPackage) == 1L &&
      !is.na(annotationPackage) && nzchar(annotationPackage),
    "`sampleDetPThreshold` must be a single number in (0, 1)" =
      is.numeric(sampleDetPThreshold) &&
      length(sampleDetPThreshold) == 1L &&
      !is.na(sampleDetPThreshold) &&
      sampleDetPThreshold > 0 && sampleDetPThreshold < 1,
    "`probeDetPThreshold` must be a single number in (0, 1)" =
      is.numeric(probeDetPThreshold) &&
      length(probeDetPThreshold) == 1L &&
      !is.na(probeDetPThreshold) &&
      probeDetPThreshold > 0 && probeDetPThreshold < 1,
    "`forceRebuild` must be a single TRUE or FALSE" =
      is.logical(forceRebuild) && length(forceRebuild) == 1L &&
      !is.na(forceRebuild),
    "`readerWorkers` must be a positive integer" =
      is.numeric(readerWorkers) && length(readerWorkers) == 1L &&
      !is.na(readerWorkers) && readerWorkers >= 1L &&
      readerWorkers == as.integer(readerWorkers)
  )

  # --- filesystem existence (value-bearing errors) ------------------
  if (!dir.exists(dataDirectory)) {
    stop(
      sprintf(
        "`dataDirectory` (\"%s\") does not point to an existing directory.",
        dataDirectory
      ),
      call. = FALSE
    )
  }
  .assertReadable(dataDirectory, "dataDirectory")
  if (!file.exists(crossReactiveProbes)) {
    stop(
      sprintf(
        "`crossReactiveProbes` (\"%s\") does not point to an existing file.",
        crossReactiveProbes
      ),
      call. = FALSE
    )
  }
  .assertReadable(crossReactiveProbes, "crossReactiveProbes")

  # `gdsOutput` is created (and overwritten on a rebuild), so validate its
  # destination up front -- before any IDAT byte is read -- so an unwritable or
  # non-existent output directory fails immediately with an actionable message
  # instead of a cryptic gdsfmt error deep inside the build.
  .assertWritableTarget(gdsOutput, "gdsOutput")

  # `samplesheet` is the path to the CSV itself (one row per sample, with a
  # `Basename` column). The IDAT files those basenames name are resolved under
  # `dataDirectory`, so the samplesheet may live anywhere.
  sheetPath <- samplesheet
  if (!file.exists(sheetPath)) {
    stop(
      sprintf(
        "`samplesheet` (\"%s\") does not point to an existing file.",
        sheetPath
      ),
      call. = FALSE
    )
  }
  .assertReadable(sheetPath, "samplesheet")

  # --- annotation package installed ---------------------------------
  # .check_annotation is injectable so the pure tests can stub it
  # (calling requireNamespace on a real Illumina annotation package
  # transitively loads minfi, which defeats the point of the pure
  # runner). Default in production is requireNamespace.
  if (!.check_annotation(annotationPackage)) {
    stop(
      sprintf(
        "`annotationPackage` (\"%s\") is not installed. Run BiocManager::install(\"%s\") and re-run the pipeline.",
        annotationPackage, annotationPackage
      ),
      call. = FALSE
    )
  }

  # --- samplesheet structure ----------------------------------------
  sampleData <- utils::read.csv(sheetPath, stringsAsFactors = FALSE)
  if (!"Basename" %in% names(sampleData)) {
    stop(
      sprintf(
        "samplesheet at %s is missing the required `Basename` column.\n  Columns found: %s.",
        sheetPath, paste(names(sampleData), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(is.na(sampleData$Basename)) || any(!nzchar(sampleData$Basename))) {
    stop(
      sprintf(
        "samplesheet at %s has empty or NA Basename entries.",
        sheetPath
      ),
      call. = FALSE
    )
  }

  # --- IDAT files exist for every sample ----------------------------
  basenames <- paste0(dataDirectory, "/", sampleData$Basename)
  grn_paths <- paste0(basenames, "_Grn.idat")
  red_paths <- paste0(basenames, "_Red.idat")
  grn_missing <- !file.exists(grn_paths) & !file.exists(paste0(grn_paths, ".gz"))
  red_missing <- !file.exists(red_paths) & !file.exists(paste0(red_paths, ".gz"))
  if (any(grn_missing) || any(red_missing)) {
    bad <- basenames[grn_missing | red_missing]
    n_bad <- length(bad)
    n_show <- min(n_bad, 5L)
    stop(n_bad,
      " sample(s) from the samplesheet have missing IDAT files:\n  ",
      paste(bad[seq_len(n_show)], collapse = "\n  "),
      if (n_bad > n_show) {
        sprintf("\n  ... and %d more", n_bad - n_show)
      } else {
        ""
      },
      "\n  (Looked for both .idat and .idat.gz under ",
      dataDirectory, ".)",
      call. = FALSE
    )
  }

  # --- cross-reactive probes CSV structure --------------------------
  xReactiveProbes <- utils::read.csv(
    file = crossReactiveProbes,
    stringsAsFactors = FALSE
  )
  if (!"TargetID" %in% names(xReactiveProbes)) {
    stop(
      sprintf(
        "CSV at %s is missing the required `TargetID` column.\n  Columns found: %s.",
        crossReactiveProbes,
        paste(names(xReactiveProbes), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  list(
    sheetPath = sheetPath,
    sampleData = sampleData,
    basenames = basenames,
    xReactiveProbes = xReactiveProbes
  )
}

# Internal: inspect an already-built GDS so the cache-reuse path can report
# the *true* surviving samples instead of pretending every input row made it
# through QC. Opens the file read-only, asserts the nodes a streaming build
# always writes are present, and reconstructs one targets row per surviving
# GDS column (matched to the input samplesheet by sample barcode, which the
# writer stores as basename(Basename) in pData/barcode).
#
# Returns NULL when the file cannot be trusted -- unreadable, truncated by a
# killed run, or missing the expected nodes. The caller treats NULL as a
# cache miss and rebuilds, so a half-written GDS from an OOM-killed run can
# never be silently reused.
#
# keepProbes is deliberately NULL: which probes were dropped at build time is
# not recoverable from the GDS alone (the original pre-drop manifest is not
# stored), so reporting anything else would be a lie. keepSamples is a mask
# over the *input* samplesheet rows, matching the fresh-build contract.
.inspectExistingGDS <- function(gdsPath, sampleData) {
  # nocov start
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    return(NULL)
  }
  # nocov end

  gds <- tryCatch(gdsfmt::openfn.gds(gdsPath, readonly = TRUE),
    error = function(e) NULL
  )
  if (is.null(gds)) {
    return(NULL)
  }
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  nodes <- tryCatch(gdsfmt::ls.gdsn(gds), error = function(e) character(0L))
  if (!all(c("betas", "pData") %in% nodes)) {
    return(NULL)
  }

  gdsBarcodes <- tryCatch(
    as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "pData/barcode"))),
    error = function(e) NULL
  )
  if (is.null(gdsBarcodes) || length(gdsBarcodes) == 0L) {
    return(NULL)
  }

  # betas must be a 2-D node whose column count matches the barcode vector;
  # a mismatch means the write was interrupted between nodes.
  bdim <- tryCatch(gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(gds, "betas"))$dim,
    error = function(e) NULL
  )
  if (is.null(bdim) || length(bdim) != 2L ||
        bdim[2L] != length(gdsBarcodes)) {
    return(NULL)
  }

  inBarcodes <- basename(sampleData$Basename)
  # One targets row per GDS column, in GDS column order. Rows for GDS
  # samples absent from the current samplesheet come back as NA but keep
  # their barcode as the row name, so the result still describes the GDS.
  idx <- match(gdsBarcodes, inBarcodes)
  targets <- sampleData[idx, , drop = FALSE]
  rownames(targets) <- gdsBarcodes
  keepSamples <- inBarcodes %in% gdsBarcodes

  list(
    targets = targets,
    keepSamples = keepSamples,
    keepProbes = NULL,
    unmatched = sum(is.na(idx))
  )
}

# Build-key sidecar (written next to the GDS as `<gds>.buildkey.rds`):
# records the inputs that determine the GDS's *content*. On a later run the
# stored key is compared against the current configuration; any difference
# means the cached GDS no longer reflects the requested pre-processing, so it
# is rebuilt rather than reused. Without this, editing a detection-p threshold
# (or swapping the annotation package / cross-reactive list) and re-running
# with forceRebuild = FALSE would silently reuse a stale GDS and the new
# settings would have no effect.
#
# The samplesheet and cross-reactive CSVs are captured by md5 rather than by
# value: a content edit (different samples, different probes) changes the
# digest and so invalidates the cache, while their absolute paths do not enter
# the key.
.buildKeyPath <- function(gdsPath) paste0(gdsPath, ".buildkey.rds")

.currentBuildKey <- function(annotationPackage, sampleDetPThreshold,
                             probeDetPThreshold, sheetPath,
                             crossReactiveProbes) {
  list(
    annotationPackage   = annotationPackage,
    sampleDetPThreshold = sampleDetPThreshold,
    probeDetPThreshold  = probeDetPThreshold,
    samplesheet_md5     = unname(tools::md5sum(sheetPath)),
    xreactive_md5       = unname(tools::md5sum(crossReactiveProbes))
  )
}

# Compare the sidecar against the current key. Returns a list with `status`
# ("match", "mismatch", or "absent") and, for a mismatch, a `detail` string
# naming the differing fields. A sidecar that is missing or unreadable is
# "absent" -- the caller still trusts the GDS (an older or hand-built file may
# simply predate the sidecar) but warns that the match cannot be verified.
.compareBuildKey <- function(gdsPath, key) {
  keyPath <- .buildKeyPath(gdsPath)
  if (!file.exists(keyPath)) {
    return(list(status = "absent", detail = ""))
  }
  prev <- tryCatch(readRDS(keyPath), error = function(e) NULL)
  if (!is.list(prev)) {
    return(list(status = "absent", detail = ""))
  }
  fields <- union(names(prev), names(key))
  diffs <- Filter(function(f) !identical(prev[[f]], key[[f]]), fields)
  if (length(diffs) == 0L) {
    list(status = "match", detail = "")
  } else {
    list(status = "mismatch", detail = paste(diffs, collapse = ", "))
  }
}

runPreprocess <- function(
  dataDirectory,
  crossReactiveProbes,
  samplesheet,
  gdsOutput,
  annotationPackage,
  sampleDetPThreshold = 0.01,
  probeDetPThreshold = 0.01,
  forceRebuild = FALSE,
  readerWorkers = 1L,
  verbose = 2L,
  BPPARAM = NULL,
  compress = "",
  .check_annotation = function(pkg) requireNamespace(pkg, quietly = TRUE)
) {
  # `compress` is the gdsfmt codec for the matrix nodes. The default "" (no
  # compression) writes ~2x faster -- compression is the single dominant cost in
  # the whole pipeline -- for a file only ~1.3x larger, so it is the right default
  # for a working GDS. Pass "LZ4_RA" to shrink the file for archiving/sharing at
  # the cost of runtime. It does not change the values, so it is deliberately NOT
  # part of the build key -- a cached GDS reads back identically whatever codec it
  # was written with.
  stopifnot(
    "`compress` must be a single string (e.g. \"LZ4_RA\" or \"\")" =
      is.character(compress) && length(compress) == 1L && !is.na(compress)
  )
  parsed <- .validateArgs(
    dataDirectory         = dataDirectory,
    crossReactiveProbes   = crossReactiveProbes,
    samplesheet           = samplesheet,
    gdsOutput           = gdsOutput,
    annotationPackage     = annotationPackage,
    sampleDetPThreshold   = sampleDetPThreshold,
    probeDetPThreshold    = probeDetPThreshold,
    forceRebuild          = forceRebuild,
    readerWorkers         = readerWorkers,
    .check_annotation     = .check_annotation
  )

  message(sprintf(
    "Config validated: %d samples, %d cross-reactive probes.",
    nrow(parsed$sampleData), nrow(parsed$xReactiveProbes)
  ))

  # Resolve the parallel backend. `readerWorkers` is the simple knob -- an
  # integer worker count realised as a MulticoreParam, which is all most
  # users need. Advanced users who require a different backend (SnowParam
  # on Windows or a socket cluster, BatchtoolsParam under an HPC
  # scheduler) pass a ready-made `BPPARAM`, in which case readerWorkers is
  # ignored. Validating it here keeps a bad backend a sub-second failure
  # rather than one that surfaces after IDAT reading has started.
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = readerWorkers)
  } else {
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
      stop("`BPPARAM` must be a BiocParallelParam object (e.g. ",
        "BiocParallel::SnowParam()), or NULL to derive one from ",
        "readerWorkers.",
        call. = FALSE
      )
    }
    if (readerWorkers != 1L) {
      warning("`BPPARAM` was supplied, so `readerWorkers` is ignored; ",
        "set the worker count on the BPPARAM object instead.",
        call. = FALSE
      )
    }
  }

  # `gdsOutput` is the literal path to the output GDS file (e.g. ".../x.gds");
  # the build-key sidecar sits next to it as "<gdsOutput>.buildkey.rds".
  gdsPath <- gdsOutput
  buildKey <- .currentBuildKey(
    annotationPackage, sampleDetPThreshold, probeDetPThreshold,
    parsed$sheetPath, crossReactiveProbes
  )

  cached <- NULL
  rebuildReason <- NULL
  if (!forceRebuild && file.exists(gdsPath)) {
    keyCmp <- .compareBuildKey(gdsPath, buildKey)
    if (keyCmp$status == "mismatch") {
      rebuildReason <- sprintf(
        "was built with a different configuration (%s changed)", keyCmp$detail
      )
    } else {
      cached <- .inspectExistingGDS(gdsPath, parsed$sampleData)
      if (!is.null(cached) && keyCmp$status == "absent") {
        warning(sprintf(
          paste0(
            "reusing %s but cannot verify it was built with the current ",
            "configuration (no build-key sidecar). Set forceRebuild = TRUE ",
            "if you have changed the detection-p thresholds, the annotation ",
            "package, the samplesheet, or the cross-reactive probe list."
          ),
          gdsPath
        ), call. = FALSE)
      }
    }
  }
  if (!is.null(cached)) {
    message(sprintf(
      "Reusing existing GDS file %s (%d samples; set forceRebuild = TRUE to overwrite)",
      gdsPath, nrow(cached$targets)
    ))
    if (cached$unmatched > 0L) {
      warning(sprintf(
        paste0(
          "%d sample(s) in %s are not present in the current ",
          "samplesheet; their targets rows are NA. Set ",
          "forceRebuild = TRUE if the GDS is stale."
        ),
        cached$unmatched, gdsPath
      ), call. = FALSE)
    }
    targets <- cached$targets
    keepSamples <- cached$keepSamples
    keepProbes <- cached$keepProbes
  } else {
    if (!forceRebuild && file.exists(gdsPath)) {
      message(sprintf(
        "Existing GDS file %s %s; rebuilding.",
        gdsPath,
        if (!is.null(rebuildReason)) {
          rebuildReason
        } else {
          "is unreadable or incomplete"
        }
      ))
    }
    res <- processMethArrayExp(
      base                  = dataDirectory,
      targets               = parsed$sampleData,
      gds_path              = gdsPath,
      annotation_package    = annotationPackage,
      drop_probes           = parsed$xReactiveProbes$TargetID,
      sample_detP_threshold = sampleDetPThreshold,
      probe_detP_threshold  = probeDetPThreshold,
      compress              = compress,
      verbose               = verbose,
      BPPARAM               = BPPARAM
    )
    # Record the build key only after the GDS write succeeds, so an
    # interrupted build never leaves a sidecar claiming a stale file is valid.
    saveRDS(buildKey, .buildKeyPath(gdsPath))
    targets <- res$targets
    keepSamples <- res$keepSamples
    keepProbes <- res$keepProbes
  }

  message(sprintf(
    "Pre-processing finished: GDS at %s, %d samples surviving.",
    gdsPath, nrow(targets)
  ))

  # `rebuilt` records whether this call actually (re)wrote the GDS -- TRUE on a
  # forced rebuild OR a cache miss/mismatch, FALSE when a valid cache was reused.
  # analyze() cascades it downstream: a fresh GDS makes any cached normbetas /
  # FUN outputs stale, so they must be recomputed even when `rebuild` did not name
  # them.
  invisible(list(
    gds_path = gdsPath,
    targets = targets,
    keepSamples = keepSamples,
    keepProbes = keepProbes,
    rebuilt = is.null(cached)
  ))
}
