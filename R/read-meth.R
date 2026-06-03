# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Parallel IDAT reader. Public API mirrors minfi::readMethArray(.exp) with an
# added BPPARAM= argument and an integer-indexed worker assembly hot path.

# Internal functions -----------------------------------------------------------

# Resolves the row-alignment between a list of per-sample IDAT outputs and
# the reference address vector pre-read in master. For a uniform cohort
# (every sample's probe count == refN) the alignment degenerates to the
# identity and sampleIdx is NULL. For a mixed cohort (same array type,
# different probe counts -- e.g. early-access vs final-release EPIC) it
# computes the intersection over the distinct address vectors and a
# per-sample integer index that maps the intersection back into each
# version's row order.
#
# Returns list(commonAddresses, sampleIdx). sampleIdx is NULL for the
# uniform case; otherwise a list of integer vectors, one per sample.
.alignProbeAddresses <- function(sampleData, refAddr) {
  ns <- vapply(sampleData, function(s) s$n, integer(1L))
  if (length(unique(ns)) == 1L) {
    return(list(
      commonAddresses = as.character(refAddr),
      sampleIdx = NULL
    ))
  }
  distinct_ns <- sort(unique(ns))
  version_addrs <- lapply(distinct_ns, function(target_n) {
    idx <- match(target_n, ns)
    a <- sampleData[[idx]]$addr
    if (is.null(a)) refAddr else a
  })
  commonAddresses <- as.character(Reduce("intersect", version_addrs))
  version_match <- lapply(
    version_addrs,
    function(a) match(commonAddresses, as.character(a))
  )
  sample_version_id <- match(ns, distinct_ns)
  sampleIdx <- lapply(
    sample_version_id,
    function(v) version_match[[v]]
  )
  list(commonAddresses = commonAddresses, sampleIdx = sampleIdx)
}

.guessArrayTypes <- function(nProbes) {
  if (nProbes >= 622000 && nProbes <= 623000) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylation450k",
      annotation = .default.450k.annotation
    )
  } else if (nProbes >= 1050000 && nProbes <= 1053000) {
    # NOTE: "Current EPIC scan type"
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationEPIC",
      annotation = .default.epic.annotation
    )
  } else if (nProbes >= 1032000 && nProbes <= 1033000) {
    # NOTE: "Old EPIC scan type"
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationEPIC",
      annotation = .default.epic.annotation
    )
  } else if (nProbes >= 1105000 && nProbes <= 1105300) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationEPICv2",
      annotation = .default.epicv2.annotation
    )
  } else if (nProbes >= 55200 && nProbes <= 55400) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylation27k",
      annotation = .default.27k.annotation
    )
  } else if (nProbes >= 54700 && nProbes <= 54800) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationAllergy",
      annotation = .default.allergy.annotation
    )
  } else if (nProbes >= 41000 && nProbes <= 41100) {
    arrayAnnotation <- c(
      array = "HorvathMammalMethylChip40",
      annotation = "test.unknown"
    )
  } else if (nProbes >= 43650 && nProbes <= 43680) {
    arrayAnnotation <- c(
      array = "IlluminaHumanMethylationAllergy",
      annotation = .default.allergy.annotation
    )
  } else {
    arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
  }
  arrayAnnotation
}

# Exported functions -----------------------------------------------------------

readMethArray <- function(basenames, extended = FALSE, verbose = 0L,
                          BPPARAM = SerialParam()) {
  # `verbose` accepts either a logical (back-compat with the historical
  # TRUE/FALSE API) or an integer in {0, 1, 2}; see ?readMethArray. v >= 2
  # turns on master-side progress lines so long-running parallel reads
  # always show activity instead of going silent between worker batches.
  verbose <- .normalize_verbose(verbose)

  # TODO: Need to think about API. Currently, if `BACKEND` is NULL then the
  #       RGSet is matrix-backed. However, within the DelayedArray package,
  #       BACKEND = NULL means to create a DelayedArray with an ordinary
  #       array as the seed. This is a source of slight tension. However, in
  #       general there is no advantage to using a DelayedMatrix with an
  #       ordinary matrix as the seed over directly using the ordinary matrix
  #       (at least within minfi).
  BACKEND <- getAutoRealizationBackend()

  # Check files exist
  basenames <- sub("_Grn\\.idat.*", "", basenames)
  basenames <- sub("_Red\\.idat.*", "", basenames)
  stopifnot(!anyDuplicated(basenames))
  G.files <- paste(basenames, "_Grn.idat", sep = "")
  names(G.files) <- basename(basenames)
  these.dont.exists <- !file.exists(G.files)
  if (any(these.dont.exists)) {
    G.files[these.dont.exists] <- paste0(G.files[these.dont.exists], ".gz")
  }
  if (!all(file.exists(G.files))) {
    noexits <- sub("\\.gz", "", G.files[!file.exists(G.files)])
    stop(
      "The following specified files do not exist:",
      paste(noexits, collapse = ", ")
    )
  }
  R.files <- paste(basenames, "_Red.idat", sep = "")
  names(R.files) <- basename(basenames)
  these.dont.exists <- !file.exists(R.files)
  if (any(these.dont.exists)) {
    R.files[these.dont.exists] <- paste0(R.files[these.dont.exists], ".gz")
  }
  if (!all(file.exists(R.files))) {
    noexits <- sub("\\.gz", "", R.files[!file.exists(R.files)])
    stop(
      "The following specified files do not exist:",
      paste(noexits, collapse = ", ")
    )
  }

  if (verbose >= 2L) {
    message(sprintf(
      "[readMethArray] %d IDAT pair(s) found; using %d worker(s) via %s",
      length(G.files), bpnworkers(BPPARAM), class(BPPARAM)[1L]
    ))
  }

  # Pre-read the first IDAT to obtain a reference set of probe addresses;
  # every worker is then handed `refN` and skips returning its own `addr`
  # whenever its array has the same probe count. For a single-array study
  # this collapses per-sample IPC by ~5-15MB (the serialised rownames),
  # which at N>=1000 is the dominant overhead in the parallel path.
  refQ <- readIDAT(G.files[[1L]])[["Quants"]]
  refAddr <- rownames(refQ)
  refN <- length(refAddr)
  rm(refQ)

  if (verbose >= 2L) {
    message(sprintf(
      "[readMethArray] Reference probe count from first IDAT: %d", refN
    ))
  }

  # Each worker reads both channels of one sample and returns only the
  # vectors needed for the final matrices. addr is sent back only when
  # the probe count differs from the reference -- the master then needs it
  # to compute the address intersection for the mixed-cohort path.
  read_one <- function(g_file, r_file, refN, extended, verbose) {
    if (verbose >= 1L) message("[readMethArray] Reading ", basename(g_file))
    gQ <- readIDAT(g_file)[["Quants"]]
    if (verbose >= 1L) message("[readMethArray] Reading ", basename(r_file))
    rQ <- readIDAT(r_file)[["Quants"]]
    n <- nrow(gQ)
    out <- list(
      n     = n,
      addr  = if (n == refN) NULL else rownames(gQ),
      gMean = gQ[, "Mean"],
      rMean = rQ[, "Mean"]
    )
    if (extended) {
      out$gSD <- gQ[, "SD"]
      out$rSD <- rQ[, "SD"]
      out$nBeads <- gQ[, "NBeads"]
    }
    out
  }

  stime <- system.time({
    sampleData <- .bpmapply_with_progress(
      FUN = read_one,
      iter_args = list(g_file = G.files, r_file = R.files),
      MoreArgs = list(
        refN = refN, extended = extended,
        verbose = verbose
      ),
      BPPARAM = BPPARAM,
      verbose = verbose,
      label = "readMethArray"
    )
  })[3]
  if (verbose >= 1L) {
    message(sprintf(
      "[readMethArray] Read idat files in %.1f seconds",
      stime
    ))
  }
  if (verbose >= 1L) {
    message("[readMethArray] Creating data matrices ... ",
      appendLF = FALSE
    )
  }
  ptime1 <- proc.time()
  allNProbes <- vapply(sampleData, function(s) s$n, integer(1L))
  arrayTypes <- cbind(do.call(rbind, lapply(allNProbes, .guessArrayTypes)),
    size = allNProbes
  )
  sameLength <- (length(unique(arrayTypes[, "size"])) == 1)
  sameArray <- (length(unique(arrayTypes[, "array"])) == 1)

  if (!sameLength && !sameArray) {
    cat("[readMethArray] Trying to parse IDAT files from different arrays.\n")
    cat("  Inferred Array sizes and types:\n")
    print(arrayTypes[, c("array", "size")])
    stop(
      "[readMethArray] Trying to parse different IDAT files, of ",
      "different size and type."
    )
  }
  if (!sameLength && verbose >= 1L) {
    message(
      sprintf(
        "[readMethArray] mixed probe counts auto-detected (%s); ",
        paste(sort(unique(arrayTypes[, "size"])), collapse = ", ")
      ),
      "output restricted to the address intersection"
    )
  }

  # IDATs from the same array with the same number of probes carry
  # identical addresses, so when sameLength holds we use the reference
  # addresses pre-read above and align every sample's vectors by position.
  # For mixed-size cohorts (different probe counts within the same array
  # type, e.g. early-access vs final-release EPIC) we *deduplicate by
  # probe count* before intersecting: a real-world mixed cohort almost
  # always has only 2-3 distinct probe counts, not N distinct address
  # sets. Running intersect over those few unique versions instead of N
  # copies turns an O(N*n) master-side cost into O(V*n), where V is
  # typically 2 or 3.
  aligned <- .alignProbeAddresses(sampleData, refAddr)
  commonAddresses <- aligned$commonAddresses
  sampleIdx <- aligned$sampleIdx

  n <- length(commonAddresses)
  k <- length(sampleData)
  sampleNames <- names(G.files)

  # Pre-allocate the output matrix and copy column-by-column, freeing each
  # sample's vector as we consume it so peak memory stays close to the
  # final matrix size rather than 2x it.
  pull <- function(field) {
    m <- matrix(integer(1L),
      nrow = n, ncol = k,
      dimnames = list(NULL, sampleNames)
    )
    storage.mode(m) <- storage.mode(sampleData[[1L]][[field]])
    if (is.null(sampleIdx)) {
      for (j in seq_len(k)) {
        m[, j] <- sampleData[[j]][[field]]
        sampleData[[j]][[field]] <<- NULL
      }
    } else {
      for (j in seq_len(k)) {
        m[, j] <- sampleData[[j]][[field]][sampleIdx[[j]]]
        sampleData[[j]][[field]] <<- NULL
      }
    }
    if (!is.null(BACKEND)) m <- realize(m, BACKEND = BACKEND)
    m
  }

  GreenMean <- pull("gMean")
  RedMean <- pull("rMean")
  if (extended) {
    GreenSD <- pull("gSD")
    RedSD <- pull("rSD")
    NBeads <- pull("nBeads")
  }
  rm(sampleData)

  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose >= 1L) {
    message(sprintf("done in %.1f seconds", stime))
  }
  if (verbose >= 2L) {
    message(sprintf(
      "[readMethArray] Array: %s / %s, %d x %d Red+Green matrices (~%s each)",
      arrayTypes[1L, 1L], arrayTypes[1L, 2L],
      nrow(GreenMean), ncol(GreenMean),
      .bytes_human(GreenMean)
    ))
  }
  if (verbose >= 1L) {
    message("[readMethArray] Instantiating final object ... ",
      appendLF = FALSE
    )
  }
  ptime1 <- proc.time()
  if (extended) {
    out <- RGChannelSetExtended(
      Red = RedMean,
      Green = GreenMean,
      RedSD = RedSD,
      GreenSD = GreenSD,
      NBeads = NBeads
    )
  } else {
    out <- RGChannelSet(Red = RedMean, Green = GreenMean)
  }
  rownames(out) <- commonAddresses
  out@annotation <- c(array = arrayTypes[1, 1], annotation = arrayTypes[1, 2])
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if (verbose >= 1L) {
    message(sprintf("done in %.1f seconds", stime))
  }
  out
}

readMethArraySheet <- function(base, pattern = "csv$", ignore.case = TRUE,
                               recursive = TRUE, verbose = TRUE) {
  readSheet <- function(file) {
    dataheader <- grep("^\\[DATA\\]", readLines(file), ignore.case = TRUE)
    if (length(dataheader) == 0) dataheader <- 0
    df <- read.csv(file, stringsAsFactor = FALSE, skip = dataheader)
    nam <- grep(
      pattern = "Sentrix_Position",
      x = names(df),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(nam) == 1) {
      df$Array <- as.character(df[, nam])
      df[, nam] <- NULL
    }
    nam <- grep(
      pattern = "Array[\\._]ID",
      x = names(df),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(nam) == 1) {
      df$Array <- as.character(df[, nam])
      df[, nam] <- NULL
    }
    if (!"Array" %in% names(df)) {
      warning(sprintf("Could not infer array name for file: %s", file))
    }
    nam <- grep("Sentrix_ID", names(df), ignore.case = TRUE, value = TRUE)
    if (length(nam) == 1) {
      df$Slide <- as.character(df[, nam])
      df[, nam] <- NULL
    }
    nam <- grep(
      pattern = "Slide[\\._]ID",
      x = names(df),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(nam) == 1) {
      df$Slide <- as.character(df[, nam])
      df[, nam] <- NULL
    }
    if (!"Slide" %in% names(df)) {
      warning(sprintf("Could not infer slide name for file: %s", file))
    } else {
      df[, "Slide"] <- as.character(df[, "Slide"])
    }
    nam <- grep(
      pattern = "Plate[\\._]ID",
      x = names(df),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(nam) == 1) {
      df$Plate <- as.character(df[, nam])
      df[, nam] <- NULL
    }
    for (nam in c("Pool_ID", "Sample_Plate", "Sample_Well")) {
      if (nam %in% names(df)) {
        df[[nam]] <- as.character(df[[nam]])
      }
    }

    if (!is.null(df$Array)) {
      patterns <- sprintf("%s_%s_Grn.idat", df$Slide, df$Array)
      allfiles <- list.files(
        path = dirname(file),
        recursive = recursive,
        full.names = TRUE
      )
      # TODO: Switch sapply() to vapply()
      basenames <- sapply(
        X = patterns,
        FUN = function(xx) grep(xx, allfiles, value = TRUE)
      )
      names(basenames) <- NULL
      basenames <- sub("_Grn\\.idat.*", "", basenames, ignore.case = TRUE)
      df$Basename <- basenames
    }
    df
  }
  if (!all(file.exists(base))) stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) {
    stop("'base needs to be either directories or files")
  }
  if (all(info$isdir)) {
    csvfiles <- list.files(
      path = base,
      recursive = recursive,
      pattern = pattern,
      ignore.case = ignore.case,
      full.names = TRUE
    )
    if (verbose) {
      message("[readMethArraySheet] Found the following CSV files:")
      print(csvfiles)
    }
  } else {
    csvfiles <- list.files(base, full.names = TRUE)
  }
  dfs <- lapply(csvfiles, readSheet)
  namesUnion <- Reduce(union, lapply(dfs, names))
  df <- do.call(
    what = rbind,
    args = lapply(dfs, function(df) {
      newnames <- setdiff(namesUnion, names(df))
      newdf <- matrix(
        data = NA,
        ncol = length(newnames),
        nrow = nrow(df),
        dimnames = list(NULL, newnames)
      )
      cbind(df, as.data.frame(newdf))
    })
  )
  df
}

readMethArrayExp <- function(base = NULL, targets = NULL, extended = FALSE,
                             recursive = FALSE, verbose = 0L,
                             BPPARAM = SerialParam()) {
  if (!is.null(targets)) {
    if (!"Basename" %in% names(targets)) {
      stop("Need 'Basename' amongst the column names of 'targets'")
    }
    if (!is.null(base)) {
      files <- file.path(base, basename(targets$Basename))
    } else {
      files <- targets$Basename
    }
    rgSet <- readMethArray(
      basenames = files,
      extended = extended,
      verbose = verbose,
      BPPARAM = BPPARAM
    )
    pD <- targets
    pD$filenames <- files
    rownames(pD) <- colnames(rgSet)
    colData(rgSet) <- as(pD, "DataFrame")
    return(rgSet)
  }
  # Now we just read all files in the directory
  Grn.files <- list.files(
    path = base,
    pattern = "_Grn.idat$",
    recursive = recursive,
    ignore.case = TRUE,
    full.names = TRUE
  )
  Red.files <- list.files(
    path = base,
    pattern = "_Red.idat$",
    recursive = recursive,
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(Grn.files) == 0 || length(Red.files) == 0) {
    stop("No IDAT files were found")
  }
  commonFiles <- intersect(
    x = sub("_Grn.idat$", "", Grn.files),
    y = sub("_Red.idat$", "", Red.files)
  )
  if (length(commonFiles) == 0) {
    stop("No IDAT files with both Red and Green channel were found")
  }
  commonFiles.Grn <- paste(commonFiles, "_Grn.idat", sep = "")
  commonFiles.Red <- paste(commonFiles, "_Red.idat", sep = "")
  if (!setequal(commonFiles.Grn, Grn.files)) {
    warning(
      sprintf(
        "the following files only exists for the green channel: %s",
        paste(
          setdiff(Grn.files, commonFiles.Grn),
          collapse = ", "
        )
      )
    )
  }
  if (!setequal(commonFiles.Red, Red.files)) {
    warning(
      sprintf(
        "the following files only exists for the red channel: %s",
        paste(
          setdiff(Red.files, commonFiles.Red),
          collapse = ", "
        )
      )
    )
  }
  readMethArray(
    basenames = commonFiles,
    extended = extended,
    verbose = verbose,
    BPPARAM = BPPARAM
  )
}
