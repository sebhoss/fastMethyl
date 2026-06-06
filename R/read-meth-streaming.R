# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Streaming IDAT reader: fuses readMethArray + detectionP + preprocessRaw +
# es2gds (plus the standard post-read sample/probe QC) into a single pass that
# writes a bigmelon-compatible GDS without ever assembling a full cohort-sized
# matrix in RAM.
#
# The function operates entirely on integer address positions in the IDAT's
# row order -- the master pre-reads the first file's rownames once, converts
# every manifest-derived address vector into integer offsets, and ships those
# (much cheaper to serialise and to index into) to the workers. The workers
# never see character addresses; their hot loop is pure positional indexing.
#
# Samples are processed in column blocks: workers compute each sample's
# methylated / unmethylated / detection-p columns in parallel, and the master
# appends the surviving columns straight into the extendable, LZ4_RA-compressed
# betas / methylated / unmethylated / pvals nodes. Peak RAM is therefore one
# block, not the whole cohort -- the property that lets large EPIC cohorts run
# without exhausting memory.
#
# Three optional QC arguments are applied during the stream:
#
#   drop_probes               character vector of CpG names to skip at
#                             compute time (sex/x-reactive probes, etc.)
#   sample_detP_threshold     drop samples with colMeans(detP) >= threshold
#   probe_detP_threshold      drop probes with any surviving sample's detP
#                             >= threshold
#
# Sample QC is per-column and decided as each sample streams past (dropped
# samples are simply not appended). Probe QC is a reduction over every surviving
# sample, so its mask is known only after the full pass; it is accumulated into
# a per-probe boolean and, when it drops any probe, the matrix nodes are
# compacted to their final probe set by a chunked row-subset rewrite. With all
# three QC args left at NULL the GDS values are identical to
# es2gds(preprocessRaw(readMethArray(...))).
#
# Mixed probe-count cohorts (same array type, e.g. early-access vs final-
# release EPIC) are auto-detected: the function always runs a quick parallel
# header pre-pass to learn each IDAT's probe count, then restricts the
# output to addresses present in every distinct version and dispatches the
# right per-version index table to each worker. For a uniform cohort the
# pre-pass adds no measurable cost (it amortises into the data pass via the
# OS file cache) and the multi-version path degenerates to a single entry.

# Internal helpers -------------------------------------------------------------

# Multi-version index builder. `version_addrs` is a list of per-version
# address vectors (typically one entry; multiple when the cohort has
# mixed probe counts -- e.g. early-access EPIC ~1.03M probes + final-
# release EPIC ~1.05M probes). The returned indices are arranged so a
# worker handed a `version_id` does pure positional indexing into its
# own IDAT's row order.
#
# The output locus set is the *intersection* of:
#   - probes whose addresses are present in EVERY version's IDAT address list,
#   - probes not in `drop_probes` (user filter).
#
# The address-presence filter applies even to a single-version cohort: a real
# manifest references addresses a given batch of IDATs may not all carry (the
# installed Bioconductor manifest is a fixed revision; production IDATs drift
# from it by a few hundred probes), and a worker indexes its IDAT positionally,
# so a probe whose address is absent cannot be computed and must be dropped
# rather than indexed as NA. On data where every manifest address is present
# (e.g. the minfiData/minfiDataEPIC fixtures) the filter keeps everything, so
# the output is byte-identical to es2gds(preprocessRaw(readMethArray(...))).
.streamingBuildIndices <- function(manifest, version_addrs, drop_probes) {
  stopifnot(is.list(version_addrs), length(version_addrs) >= 1L)

  TypeI.Red <- getProbeInfo(manifest, type = "I-Red")
  TypeI.Green <- getProbeInfo(manifest, type = "I-Green")
  TypeII <- getProbeInfo(manifest, type = "II")
  locusNames <- getManifestInfo(manifest, "locusNames")
  controlAddr <- getControlAddress(manifest, controlType = "NEGATIVE")

  # 1. "Available" means an address is present in every version's address list.
  # %in% hashes its table argument, so each test is O(length(addrs)) -- cheap
  # even for the single-version case, where it reduces to one membership test
  # against that version's IDAT addresses.
  addr_in_all <- function(addrs) {
    Reduce(`&`, lapply(version_addrs, function(v) addrs %in% v))
  }

  # 2. Filter probe tables to probes whose addresses are present in EVERY
  #    version's address list. AddressB only exists for TypeI probes.
  TR_addr_ok <- addr_in_all(TypeI.Red$AddressA) &
    addr_in_all(TypeI.Red$AddressB)
  TG_addr_ok <- addr_in_all(TypeI.Green$AddressA) &
    addr_in_all(TypeI.Green$AddressB)
  TII_addr_ok <- addr_in_all(TypeII$AddressA)
  TypeI.Red <- TypeI.Red[TR_addr_ok, ]
  TypeI.Green <- TypeI.Green[TG_addr_ok, ]
  TypeII <- TypeII[TII_addr_ok, ]

  # 3. Apply drop_probes filter (a user-supplied character vector of CpG
  #    names to skip at compute time).
  drop_set <- if (is.null(drop_probes)) character(0L) else drop_probes
  TypeI.Red <- TypeI.Red[!(TypeI.Red$Name %in% drop_set), ]
  TypeI.Green <- TypeI.Green[!(TypeI.Green$Name %in% drop_set), ]
  TypeII <- TypeII[!(TypeII$Name %in% drop_set), ]

  # 4. out_loci preserves the manifest locusNames order, minus dropped
  #    probes and probes whose addresses didn't survive the intersection.
  survivors <- c(TypeI.Red$Name, TypeI.Green$Name, TypeII$Name)
  out_loci <- locusNames[locusNames %in% survivors]

  # 5. Restrict control probes to addresses present in every version.
  controlAddr <- controlAddr[addr_in_all(controlAddr)]

  # 6. For each version, build the integer index arrays the worker uses.
  #    Every match() here MUST succeed (we already pruned to the
  #    intersection), so anyNA() is a hard assertion failure.
  build_one_version <- function(addr) {
    idx <- list(
      ctrl_pos        = match(controlAddr, addr),
      TI_Red_addrA    = match(TypeI.Red$AddressA, addr),
      TI_Red_addrB    = match(TypeI.Red$AddressB, addr),
      TI_Red_out      = match(TypeI.Red$Name, out_loci),
      TI_Green_addrA  = match(TypeI.Green$AddressA, addr),
      TI_Green_addrB  = match(TypeI.Green$AddressB, addr),
      TI_Green_out    = match(TypeI.Green$Name, out_loci),
      TII_addrA       = match(TypeII$AddressA, addr),
      TII_out         = match(TypeII$Name, out_loci),
      refN            = length(addr)
    )
    # Defensive invariant: out_loci is the address-intersection of
    # every version's IDAT, so every name in it must resolve in the
    # current version's address vector. Only fires if the addr_in_all
    # filtering above has been broken.
    # nocov start
    if (anyNA(unlist(idx[setdiff(names(idx), "refN")]))) {
      stop(
        "[processMethArray] internal: address present in ",
        "intersection but match failed -- this is a bug"
      )
    }
    # nocov end
    idx
  }
  per_version_indices <- lapply(version_addrs, build_one_version)

  list(
    per_version  = per_version_indices,
    out_loci     = out_loci,
    n_out        = length(out_loci),
    n_versions   = length(version_addrs)
  )
}

# Per-sample worker. Every address reference is a pre-resolved integer
# offset into THIS sample's IDAT row order; for multi-version cohorts the
# master picks the right index table per sample via `version_id`.
.streamingComputeColumns <- function(g_file, r_file, version_id,
                                     per_version_idx, n_out, verbose) {
  idx <- per_version_idx[[version_id]]
  refN <- idx$refN

  if (verbose >= 1L) message("[processMethArray] Reading ", basename(g_file))
  gQ <- readIDAT(g_file)[["Quants"]]
  if (verbose >= 1L) message("[processMethArray] Reading ", basename(r_file))
  rQ <- readIDAT(r_file)[["Quants"]]
  # Defensive invariant: master pre-scans every IDAT's header and
  # picks the version_id that matches its probe count, so worker-side
  # this should never see a mismatch. Only fires if the header scan
  # at line ~320 returned a different value from the data read here.
  # nocov start
  if (nrow(gQ) != refN || nrow(rQ) != refN) {
    stop(
      "[processMethArray] IDAT probe count for ", basename(g_file),
      " (", nrow(gQ), ") does not match version ", version_id,
      "'s expected ", refN, " probes"
    )
  }
  # nocov end
  gMean <- gQ[, "Mean"]
  rMean <- rQ[, "Mean"]
  rm(gQ, rQ)

  rBg <- rMean[idx$ctrl_pos]
  gBg <- gMean[idx$ctrl_pos]
  rMu <- median(rBg)
  rSd <- mad(rBg)
  gMu <- median(gBg)
  gSd <- mad(gBg)

  Meth <- integer(n_out)
  Unmeth <- integer(n_out)
  detP <- numeric(n_out)

  iR <- rMean[idx$TI_Red_addrA] + rMean[idx$TI_Red_addrB]
  detP[idx$TI_Red_out] <- pnorm(iR, mean = 2 * rMu, sd = 2 * rSd, lower.tail = FALSE)
  Meth[idx$TI_Red_out] <- rMean[idx$TI_Red_addrB]
  Unmeth[idx$TI_Red_out] <- rMean[idx$TI_Red_addrA]

  iG <- gMean[idx$TI_Green_addrA] + gMean[idx$TI_Green_addrB]
  detP[idx$TI_Green_out] <- pnorm(iG, mean = 2 * gMu, sd = 2 * gSd, lower.tail = FALSE)
  Meth[idx$TI_Green_out] <- gMean[idx$TI_Green_addrB]
  Unmeth[idx$TI_Green_out] <- gMean[idx$TI_Green_addrA]

  iII <- rMean[idx$TII_addrA] + gMean[idx$TII_addrA]
  detP[idx$TII_out] <- pnorm(iII, mean = rMu + gMu, sd = rSd + gSd, lower.tail = FALSE)
  Meth[idx$TII_out] <- gMean[idx$TII_addrA]
  Unmeth[idx$TII_out] <- rMean[idx$TII_addrA]

  list(Meth = Meth, Unmeth = Unmeth, detP = detP)
}

.streamingResolveManifest <- function(refN) {
  arrayAnn <- .guessArrayTypes(refN)
  if (arrayAnn["array"] == "Unknown") {
    stop(
      "[processMethArray] could not infer array type from ",
      refN, " probes"
    )
  }
  annTuple <- c(
    array = unname(arrayAnn["array"]),
    annotation = unname(arrayAnn["annotation"])
  )
  list(annTuple = annTuple, manifest = getManifest(annTuple))
}

# Per-block QC primitives, kept standalone so the mask logic is unit-testable on
# synthetic detection-p blocks without any IDAT or GDS machinery.
#
# Sample QC is decidable per column the moment a sample's detP column exists, so
# it streams trivially: a sample survives when its mean detection-p is below the
# threshold (every sample survives when no threshold is set).
.streamingChunkKeepSamples <- function(detP_chunk, sample_detP_threshold) {
  if (is.null(sample_detP_threshold)) {
    return(rep(TRUE, ncol(detP_chunk)))
  }
  colMeans(detP_chunk) < sample_detP_threshold
}

# Probe QC is a reduction over every surviving sample, so it cannot be decided
# until the whole cohort has streamed past. This returns one block's per-probe
# "fails in at least one of these (surviving) samples" contribution; the caller
# ORs it across blocks into a running probe_fails accumulator.
.streamingChunkProbeFails <- function(detP_chunk, probe_detP_threshold) {
  rowSums(detP_chunk >= probe_detP_threshold) > 0L
}

# Emit the verbose QC drop summary and stop with an actionable message if a
# threshold removed every sample or every probe. Operates on the final keep
# masks (sample QC accumulated per block, probe QC accumulated into
# probe_fails), so it is pure and unit-testable on synthetic mask vectors.
.streamingQCfinalize <- function(keepSamples, keepProbes, sampleNames,
                                 sample_detP_threshold, probe_detP_threshold,
                                 verbose) {
  k <- length(keepSamples)
  nLoci <- length(keepProbes)
  if (!all(keepSamples) || !all(keepProbes)) {
    if (verbose >= 1L) {
      message(sprintf(
        "[processMethArray] QC drops %d/%d samples, %d/%d probes",
        sum(!keepSamples), k, sum(!keepProbes), nLoci
      ))
    }
    if (verbose >= 2L && sum(!keepSamples) > 0L && sum(!keepSamples) <= 20L) {
      message(sprintf(
        "[processMethArray] dropped samples: %s",
        paste(sampleNames[!keepSamples], collapse = ", ")
      ))
    }
  }
  if (!any(keepSamples)) {
    stop(sprintf(
      paste0(
        "[processMethArray] all %d samples were dropped by ",
        "sample_detP_threshold = %g (every sample's mean detection-p is at ",
        "or above it). Raise the threshold or check IDAT quality."
      ),
      k, sample_detP_threshold
    ), call. = FALSE)
  }
  if (!any(keepProbes)) {
    stop(sprintf(
      paste0(
        "[processMethArray] all %d probes were dropped by ",
        "probe_detP_threshold = %g (every probe fails detection in at least ",
        "one surviving sample). Raise the threshold or check IDAT quality."
      ),
      nLoci, probe_detP_threshold
    ), call. = FALSE)
  }
  invisible(NULL)
}

# Restrict a probe-keep mask to probes the annotation actually covers. `out_loci`
# comes from the manifest, but the installed annotation package is a fixed
# revision that may not annotate every manifest locus (real EPIC manifests carry
# a few hundred probes the annotation lacks). A probe with no annotation cannot
# form a valid bigmelon fData row (no Name, no Type), so it is dropped here --
# folded into keepProbes so the same compaction removes it -- rather than
# written as an all-NA row, whose empty Name would duplicate across such probes
# and break feature-keyed steps like dasen.
.streamingAnnoKeep <- function(keepProbes, out_loci, annoNames) {
  keepProbes & (out_loci %in% annoNames)
}

# Drop the rows flagged by !keepProbes from an already-streamed matrix node, in
# place in the open GDS. Probe QC is only known after the full streaming pass,
# by which point every probe row has been written; this compacts a node to its
# final probe set by reading it back in column blocks (LZ4_RA supports random
# hyperslab reads), subsetting rows, and re-emitting into a fresh node that
# replaces the original. Chunked, so the rewrite never holds the full node in
# RAM. Only invoked when probe_detP_threshold actually dropped probes.
.streamingCompactRows <- function(bmln, node_name, keepProbes, compress,
                                  block_cols) {
  old <- index.gdsn(bmln, node_name)
  d <- objdesp.gdsn(old)$dim
  nrow_old <- d[1L]
  ncol_old <- d[2L]
  tmp_name <- paste0(node_name, "__compact")
  newnode <- NULL
  for (j0 in seq.int(1L, ncol_old, by = block_cols)) {
    ncols <- min(block_cols, ncol_old - j0 + 1L)
    blk <- read.gdsn(old, start = c(1L, j0), count = c(nrow_old, ncols))
    blk <- blk[keepProbes, , drop = FALSE]
    if (is.null(newnode)) {
      newnode <- add.gdsn(bmln, tmp_name, val = blk,
                          compress = compress, closezip = FALSE)
    } else {
      append.gdsn(newnode, blk)
    }
  }
  readmode.gdsn(newnode)
  delete.gdsn(old, force = TRUE)
  rename.gdsn(newnode, node_name)
  invisible(NULL)
}

# Compact ONE node into its own fresh GDS file and return that file's path. Reads
# the source GDS read-only (the caller must have closed its writable handle
# first, so forked workers can open it), so several of these can run in parallel
# -- one per matrix node -- without sharing a handle. The chunked column-block
# read keeps each worker's working set to one slab.
.streamingCompactNodeToFile <- function(gds_path, node_name, keepProbes,
                                        compress, block_cols) {
  src <- openfn.gds(gds_path, readonly = TRUE)
  on.exit(closefn.gds(src), add = TRUE)
  old <- index.gdsn(src, node_name)
  d <- objdesp.gdsn(old)$dim
  tmp <- tempfile(fileext = ".gds")
  dst <- createfn.gds(tmp)
  newnode <- NULL
  for (j0 in seq.int(1L, d[2L], by = block_cols)) {
    ncols <- min(block_cols, d[2L] - j0 + 1L)
    blk <- read.gdsn(old, start = c(1L, j0), count = c(d[1L], ncols))
    blk <- blk[keepProbes, , drop = FALSE]
    if (is.null(newnode)) {
      newnode <- add.gdsn(dst, node_name, val = blk,
                          compress = compress, closezip = FALSE)
    } else {
      append.gdsn(newnode, blk)
    }
  }
  readmode.gdsn(newnode)
  closefn.gds(dst)
  tmp
}

# Parallel probe-QC compaction: the four matrix nodes are independent, so each is
# compacted in its own worker into a temp file, then the master splices the
# compacted nodes back via copyto.gdsn (a cheap node copy, not a re-compress).
# The writable handle is closed first so workers can open the file read-only;
# this function reopens it and returns the new open handle (the caller must
# rebind it). The per-node block budget is divided by the worker count so the
# concurrent slabs sum to the serial path's single-slab budget -- the flat-memory
# invariant holds. Falls back to nothing here: the serial .streamingCompactRows
# is used directly when there is only one worker.
.streamingCompactParallel <- function(bmln, gds_path, node_order, keepProbes,
                                      compress, block_cap, BPPARAM) {
  nworkers <- max(1L, bpnworkers(BPPARAM))
  block_cols <- max(1L, block_cap %/% nworkers)
  closefn.gds(bmln)
  tmps <- bplapply(node_order, .streamingCompactNodeToFile,
                   gds_path = gds_path, keepProbes = keepProbes,
                   compress = compress, block_cols = block_cols,
                   BPPARAM = BPPARAM)
  bmln <- openfn.gds(gds_path, readonly = FALSE)
  # Delete every old node first so the copies can reuse the freed space rather
  # than appending (keeps the file from bloating to old + compacted).
  for (nm in node_order) delete.gdsn(index.gdsn(bmln, nm), force = TRUE)
  for (i in seq_along(node_order)) {
    s <- openfn.gds(tmps[[i]], readonly = TRUE)
    copyto.gdsn(bmln, index.gdsn(s, node_order[[i]]))
    closefn.gds(s)
    unlink(tmps[[i]])
  }
  bmln
}

# Write the non-matrix scaffolding (fData, pData, history, paths) into the open
# GDS once the matrix nodes have been streamed in. `annot_df` is the per-probe
# annotation for the final (compacted) probe set, already aligned to the matrix
# rows by the caller. fData is stored column by column into a folder node tagged
# as a data.frame: byte-equivalent to add.gdsn(val = <data.frame>) (gdsfmt
# represents a data.frame as exactly such a folder of one character child per
# column plus the R.class attribute, with no row names) but without
# materialising a second all-character copy of every column at once.
# check = FALSE keeps it silent -- annotation columns legitimately carry NA
# (gene maps, SNP fields), which String storage represents as "".
.streamingWriteMeta <- function(bmln, annot_df, sampleNames, pData_df) {
  fnode <- addfolder.gdsn(bmln, "fData")
  put.attr.gdsn(fnode, "R.class", "data.frame")
  for (col in names(annot_df)) {
    add.gdsn(fnode, col, val = as.character(annot_df[[col]]), check = FALSE)
  }
  fData_cols <- names(annot_df)
  rm(annot_df)

  if (is.null(pData_df)) {
    pData_df <- data.frame(barcode = sampleNames, stringsAsFactors = FALSE)
  } else {
    pData_df <- data.frame(
      barcode = sampleNames,
      lapply(pData_df, as.character),
      stringsAsFactors = FALSE
    )
  }
  add.gdsn(bmln, "pData", val = pData_df)

  history_df <- data.frame(
    submitted = as.character(Sys.time()),
    finished = as.character(Sys.time()),
    command = paste0(
      "Streaming processMethArray (minfi ",
      packageVersion("minfi"), ")"
    ),
    stringsAsFactors = FALSE
  )
  add.gdsn(bmln, "history", val = history_df)

  rowna <- fData_cols[head(
    c(
      grep("targetID", fData_cols, ignore.case = TRUE),
      grep("name", fData_cols, ignore.case = TRUE),
      grep("Probe_ID", fData_cols, ignore.case = TRUE)
    ),
    1L
  )]
  pData_cols <- names(pData_df)
  colna <- pData_cols[head(
    c(
      grep("sampleID", pData_cols, ignore.case = TRUE),
      grep("barcode", pData_cols, ignore.case = TRUE),
      grep("name", pData_cols, ignore.case = TRUE)
    ),
    1L
  )]
  add.gdsn(bmln, "paths",
    val = c(paste0("fData/", rowna), paste0("pData/", colna))
  )
  invisible(NULL)
}

# Exported function ------------------------------------------------------------

processMethArray <- function(basenames,
                             gds_path,
                             pData = NULL,
                             drop_probes = NULL,
                             sample_detP_threshold = NULL,
                             probe_detP_threshold = NULL,
                             compress = "LZ4_RA",
                             verbose = 0L,
                             BPPARAM = SerialParam()) {
  verbose <- .normalize_verbose(verbose)
  if (missing(gds_path) || is.null(gds_path) || !is.character(gds_path) ||
        length(gds_path) != 1L || is.na(gds_path) || !nzchar(gds_path)) {
    stop("`gds_path` must be a single non-empty string.")
  }
  if (!is.null(sample_detP_threshold)) {
    stopifnot(
      is.numeric(sample_detP_threshold),
      length(sample_detP_threshold) == 1L,
      sample_detP_threshold > 0, sample_detP_threshold < 1
    )
  }
  if (!is.null(probe_detP_threshold)) {
    stopifnot(
      is.numeric(probe_detP_threshold),
      length(probe_detP_threshold) == 1L,
      probe_detP_threshold > 0, probe_detP_threshold < 1
    )
  }

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

  # Pre-read first IDAT to anchor the reference probe count + address order.
  refQ <- readIDAT(G.files[[1L]])[["Quants"]]
  refN <- nrow(refQ)
  refAddr <- rownames(refQ)
  rm(refQ)

  mInfo <- .streamingResolveManifest(refN)
  manifest <- mInfo$manifest
  annTuple <- mInfo$annTuple

  # ---- Version auto-detection -----------------------------------------------
  # Mixed probe-count cohorts (same array type, different versions -- e.g.
  # early-access vs final-release EPIC) are detected automatically via a
  # quick parallel header-only scan. On a uniform cohort the scan amortises
  # into the OS file cache for the data pass that follows, so the overhead
  # is within run-to-run noise even at N=1000+.
  if (verbose >= 2L) {
    message("[processMethArray] scanning IDAT headers for probe counts")
  }
  nProbes_per_file <- unlist(bplapply(G.files, function(f) {
    as.integer(readIDAT(f, what = "nSNPsRead"))
  }, BPPARAM = BPPARAM), use.names = FALSE)
  distinct_n <- sort(unique(nProbes_per_file))

  # For each distinct probe count, capture its address order. We already
  # have refAddr for `refN`; read just the others.
  version_addrs <- list()
  version_addrs[[as.character(refN)]] <- refAddr
  for (n in distinct_n) {
    key <- as.character(n)
    if (!is.null(version_addrs[[key]])) next
    idx <- match(n, nProbes_per_file)
    version_addrs[[key]] <- rownames(readIDAT(G.files[[idx]])[["Quants"]])
  }
  # Order versions by probe count (ascending) for deterministic dispatch.
  version_addrs <- version_addrs[as.character(distinct_n)]

  sample_version_id <- match(nProbes_per_file, distinct_n)

  if (verbose >= 2L) {
    message(sprintf(
      "[processMethArray] %d IDAT pair(s) found; %d worker(s) via %s",
      length(G.files), bpnworkers(BPPARAM), class(BPPARAM)[1L]
    ))
    message(sprintf(
      "[processMethArray] Array resolved: %s / %s (refN=%d)",
      annTuple["array"], annTuple["annotation"], refN
    ))
    if (length(distinct_n) > 1L) {
      counts <- tabulate(sample_version_id, nbins = length(distinct_n))
      mix <- paste(sprintf(
        "%d probes x %d samples",
        distinct_n, counts
      ), collapse = "; ")
      message(sprintf(
        "[processMethArray] mixed-version cohort detected: %s", mix
      ))
    }
  }

  if (!is.null(drop_probes)) {
    drop_probes <- as.character(drop_probes)
  }
  indices <- .streamingBuildIndices(manifest, version_addrs, drop_probes)
  if (indices$n_out == 0L) {
    stop(
      "[processMethArray] no loci left after drop_probes + ",
      "multi-version intersection"
    )
  }
  if (verbose >= 1L) {
    total_loci <- length(getManifestInfo(manifest, "locusNames"))
    if (indices$n_out < total_loci) {
      message(sprintf(
        "[processMethArray] loci kept after drop_probes%s: %d -> %d (skipping %d)",
        if (length(distinct_n) > 1L) " + multi-version intersect" else "",
        total_loci, indices$n_out, total_loci - indices$n_out
      ))
    } else if (verbose >= 2L) {
      message(sprintf(
        "[processMethArray] no drop_probes: %d loci kept",
        indices$n_out
      ))
    }
  }

  k <- length(G.files)
  sampleNames <- names(G.files)
  nLoci <- indices$n_out
  out_loci <- indices$out_loci

  # Free the pre-streaming scaffolding before the parallel block loop: the
  # manifest and the per-version / reference address vectors are consumed by
  # .streamingBuildIndices and unused thereafter (workers index by the integer
  # tables in `indices`, not these). Holding them would bloat the master heap,
  # and every forked MulticoreParam worker copy-on-write-inflates against that
  # heap -- R's garbage collector in a worker touches the whole heap -- so a
  # large object retained here is effectively multiplied by the worker count.
  rm(manifest, mInfo, version_addrs, refAddr)

  # Stream samples to the GDS in column blocks so full cohort-sized matrices are
  # never assembled in RAM: workers compute each sample's columns in parallel,
  # the master appends the surviving columns to the (compressed, extendable)
  # betas/methylated/unmethylated/pvals nodes, and QC state is accumulated as it
  # goes -- peak RAM is one block, not the whole cohort. The GDS is opened here
  # and removed on any error so a partial file is never left to be mistaken for
  # a valid cache by the build-key reuse layer.
  bmln <- createfn.gds(gds_path)
  closed <- FALSE
  ok <- FALSE
  on.exit({
    if (!closed) try(closefn.gds(bmln), silent = TRUE)
    if (!ok && file.exists(gds_path)) unlink(gds_path)
  }, add = TRUE)

  nworkers <- max(1L, bpnworkers(BPPARAM))
  # The per-block forced GC (below) only earns its keep when the backend forks
  # multiple workers against the master heap. A non-forking backend (SerialParam,
  # SnowParam) has no shared heap to COW-inflate, and a single worker cannot
  # multiply anything, so the collection is pure overhead there and is skipped.
  forking <- methods::is(BPPARAM, "MulticoreParam") && nworkers > 1L
  # Bound a block's working set (worker results + assembled block + clip
  # temporaries) by capping columns per block at a fixed element budget, and at
  # nworkers * 4 so a block is never larger than the parallel fan-out needs.
  block_cap <- max(1L, as.integer(6e7 %/% max(nLoci, 1L)))
  block_size <- max(1L, min(k, nworkers * 4L, block_cap))
  blocks <- split(seq_len(k), ceiling(seq_len(k) / block_size))

  node_order <- c("betas", "methylated", "unmethylated", "pvals")
  nodes <- NULL
  keepSamples <- logical(k)
  probe_fails <- if (is.null(probe_detP_threshold)) NULL else logical(nLoci)

  if (verbose >= 1L) {
    message(sprintf(
      "[processMethArray] Streaming %d sample(s) to GDS in blocks of %d ...",
      k, block_size
    ))
  }
  # Split the streaming wall-clock into the parallel worker-compute (the
  # .bpmapply_try call) and the serial master-write (assemble + betas +
  # compress + append), reported at verbose >= 2. The master-write is dominated
  # by LZ4 compression and is the serial floor of the streaming pass.
  compute_secs <- 0
  write_secs <- 0
  start_t <- Sys.time()
  for (bi in blocks) {
    .t0 <- proc.time()[[3]]
    block_res <- .bpmapply_try(
      FUN = .streamingComputeColumns,
      iter_args = list(
        g_file = G.files[bi], r_file = R.files[bi],
        version_id = sample_version_id[bi]
      ),
      MoreArgs = list(
        per_version_idx = indices$per_version,
        n_out           = nLoci,
        verbose         = verbose
      ),
      BPPARAM = BPPARAM, label = "processMethArray"
    )
    .t1 <- proc.time()[[3]]
    Mblk <- vapply(block_res, `[[`, integer(nLoci), "Meth")
    Ublk <- vapply(block_res, `[[`, integer(nLoci), "Unmeth")
    Pblk <- vapply(block_res, `[[`, numeric(nLoci), "detP")
    rm(block_res)

    keep_cols <- .streamingChunkKeepSamples(Pblk, sample_detP_threshold)
    keepSamples[bi] <- keep_cols
    if (any(keep_cols)) {
      sel <- which(keep_cols)
      Mblk <- Mblk[, sel, drop = FALSE]
      Ublk <- Ublk[, sel, drop = FALSE]
      Pblk <- Pblk[, sel, drop = FALSE]
      # Clip intensities to compute Betas for this block only; the temporaries
      # are block-sized, never cohort-sized. 0/0 yields NaN for a probe with no
      # signal in a sample, matching the whole-matrix form.
      mclip <- pmax(Mblk, 0L)
      uclip <- pmax(Ublk, 0L)
      Bblk <- mclip / (mclip + uclip)
      rm(mclip, uclip)
      if (is.null(nodes)) {
        nodes <- list(
          betas        = add.gdsn(bmln, "betas", val = Bblk,
                                  compress = compress, closezip = FALSE),
          methylated   = add.gdsn(bmln, "methylated", val = Mblk,
                                  compress = compress, closezip = FALSE),
          unmethylated = add.gdsn(bmln, "unmethylated", val = Ublk,
                                  compress = compress, closezip = FALSE),
          pvals        = add.gdsn(bmln, "pvals", val = Pblk,
                                  compress = compress, closezip = FALSE)
        )
      } else {
        append.gdsn(nodes$betas, Bblk)
        append.gdsn(nodes$methylated, Mblk)
        append.gdsn(nodes$unmethylated, Ublk)
        append.gdsn(nodes$pvals, Pblk)
      }
      if (!is.null(probe_fails)) {
        probe_fails <- probe_fails |
          .streamingChunkProbeFails(Pblk, probe_detP_threshold)
      }
      rm(Bblk)
    }
    if (verbose >= 2L) {
      done <- max(bi)
      elapsed <- as.numeric(Sys.time() - start_t, units = "secs")
      message(sprintf(
        "[processMethArray] progress %d/%d (%5.1f%%) elapsed %.1fs (%.1f file/s)",
        done, k, 100 * done / k, elapsed, done / max(elapsed, 0.001)
      ))
    }
    # Drop this block's matrices and, when the next iteration will fork workers,
    # force a collection first. MulticoreParam workers copy-on-write against the
    # master heap (a worker's GC touches the whole heap), so leaving a block's
    # ~GB of now-collectable matrices live at fork time multiplies them by the
    # worker count -- the per-block growth that drove the peak up to OOM on a
    # large cohort. Collecting here keeps the fork base minimal; for a
    # non-forking backend it is skipped (see `forking`).
    rm(Mblk, Ublk, Pblk)
    if (forking) gc(FALSE)
    compute_secs <- compute_secs + (.t1 - .t0)
    write_secs <- write_secs + (proc.time()[[3]] - .t1)
  }
  if (verbose >= 1L) {
    elapsed <- as.numeric(Sys.time() - start_t, units = "secs")
    message(sprintf(
      "[processMethArray] Streamed in %.1f seconds (%.1f file/s)",
      elapsed, k / max(elapsed, 0.001)
    ))
  }
  if (verbose >= 2L) {
    message(sprintf(
      paste0("[processMethArray] streaming split: worker-compute %.1fs ",
             "(parallel) + master-write %.1fs (serial)"),
      compute_secs, write_secs
    ))
  }

  keepProbes <- if (is.null(probe_fails)) rep(TRUE, nLoci) else !probe_fails

  # Build the per-probe annotation now (matrices are on disk, so this ~1.8 GB
  # EPIC frame is the dominant live object, not an addition on top of them). It
  # serves two purposes: it is the fData source, and it defines which loci are
  # analysable -- a manifest locus the installed annotation does not cover has
  # no fData row, so it is folded out of keepProbes here and removed by the same
  # compaction below. getAnnotation() can be slow/memory-bound, so announce it
  # before the call: a stall here is then visible rather than silent.
  anno_str <- .getAnnotationString(annTuple)
  if (verbose >= 2L) {
    message(sprintf(
      "[processMethArray] building annotation '%s'; %s",
      anno_str, .mem_report()
    ))
  }
  atime <- system.time(annot_full <- getAnnotation(anno_str))[3]
  if (verbose >= 2L) {
    message(sprintf(
      "[processMethArray] annotation built: %d x %d in %.1fs; %s",
      nrow(annot_full), ncol(annot_full), atime, .mem_report()
    ))
  }
  n_before_anno <- sum(keepProbes)
  keepProbes <- .streamingAnnoKeep(keepProbes, out_loci, annot_full$Name)
  if (verbose >= 1L && sum(keepProbes) < n_before_anno) {
    message(sprintf(
      "[processMethArray] dropped %d probe(s) absent from the annotation",
      n_before_anno - sum(keepProbes)
    ))
  }

  # Stops (and triggers the on.exit cleanup) if a threshold removed every sample
  # or every probe; otherwise at least one sample survived so `nodes` below is
  # non-NULL.
  .streamingQCfinalize(keepSamples, keepProbes, sampleNames,
                       sample_detP_threshold, probe_detP_threshold, verbose)
  for (nm in node_order) readmode.gdsn(nodes[[nm]])

  if (!all(keepProbes)) {
    if (verbose >= 1L) {
      message(sprintf(
        "[processMethArray] Compacting %d -> %d probes ...",
        nLoci, sum(keepProbes)
      ))
    }
    if (max(1L, bpnworkers(BPPARAM)) > 1L) {
      # Compact the four nodes in parallel (one worker per node); reopens and
      # rebinds the GDS handle. block_cap (the memory-bound width), not the
      # streaming fan-out block_size, is the natural read width here.
      bmln <- .streamingCompactParallel(bmln, gds_path, node_order, keepProbes,
                                        compress, block_cap, BPPARAM)
    } else {
      for (nm in node_order) {
        .streamingCompactRows(bmln, nm, keepProbes, compress, block_cap)
      }
    }
  }

  final_loci <- out_loci[keepProbes]
  # Every final locus is in the annotation (keepProbes was filtered to it), so
  # this match has no NA -- fData carries a real Name/Type for every probe.
  final_annot <- annot_full[match(final_loci, annot_full$Name), , drop = FALSE]
  pData_kept <- if (is.null(pData)) NULL else pData[keepSamples, , drop = FALSE]
  .streamingWriteMeta(bmln, final_annot, sampleNames[keepSamples], pData_kept)

  closefn.gds(bmln)
  closed <- TRUE
  ok <- TRUE
  if (verbose >= 2L && file.exists(gds_path)) {
    message(sprintf(
      "[processMethArray] GDS file size: %s (%s)",
      .bytes_human_count(file.size(gds_path)), gds_path
    ))
  }
  invisible(list(
    gds_path = gds_path,
    keepSamples = keepSamples,
    keepProbes = keepProbes
  ))
}

# Experiment-level wrapper around processMethArray(). Bundles the
# orchestration that every fork-using pipeline repeats: resolve IDAT
# basenames from a samplesheet (+ optional base dir), load the
# annotation package, derive sex-chromosome probe names if asked,
# concatenate them with any caller-supplied drop_probes, and call
# processMethArray(). Returns the filtered
# samplesheet (post-QC sample drops applied) alongside the gds_path
# and keep masks, so the caller can plug them straight into a
# downstream analysis without re-deriving them.
#
# Mirrors readMethArrayExp() so users familiar with the upstream
# experiment-level reader find this naturally.
processMethArrayExp <- function(base = NULL,
                                targets,
                                gds_path,
                                annotation_package,
                                drop_probes = NULL,
                                drop_sex_chromosomes = TRUE,
                                sample_detP_threshold = NULL,
                                probe_detP_threshold = NULL,
                                compress = "LZ4_RA",
                                verbose = 0L,
                                BPPARAM = SerialParam()) {
  # ---- input validation -------------------------------------------------
  if (missing(targets) || !is.data.frame(targets)) {
    stop("`targets` must be a data.frame (or DataFrame) with a Basename column.")
  }
  if (!"Basename" %in% names(targets)) {
    stop("`targets` is missing the required `Basename` column.")
  }
  if (missing(gds_path) || !is.character(gds_path) || length(gds_path) != 1L ||
        is.na(gds_path) || !nzchar(gds_path)) {
    stop("`gds_path` must be a single non-empty string.")
  }
  if (missing(annotation_package) || !is.character(annotation_package) ||
        length(annotation_package) != 1L || is.na(annotation_package) ||
        !nzchar(annotation_package)) {
    stop(
      "`annotation_package` must be a single non-empty string ",
      "naming an installed Illumina*anno.* Bioconductor package."
    )
  }
  if (!is.logical(drop_sex_chromosomes) || length(drop_sex_chromosomes) != 1L ||
        is.na(drop_sex_chromosomes)) {
    stop("`drop_sex_chromosomes` must be a single TRUE or FALSE.")
  }
  verbose <- .normalize_verbose(verbose)

  # ---- resolve basenames ------------------------------------------------
  # Differs from readMethArrayExp() deliberately: file.path joins
  # base with targets$Basename without a basename() strip, so the
  # typical Illumina samplesheet convention of
  # `Basename = "<slide>/<array>"` (one slide-directory per chip) works
  # directly. Pass `base = NULL` and full paths in `targets$Basename`
  # if you have already pre-resolved.
  files <- if (!is.null(base)) {
    file.path(base, targets$Basename)
  } else {
    targets$Basename
  }

  # ---- annotation loading + sex-chromosome probe lookup -----------------
  sex_probes <- character(0L)
  if (drop_sex_chromosomes) {
    if (!requireNamespace(annotation_package, quietly = TRUE)) {
      stop(
        "annotation package not installed: ", annotation_package,
        ". Install with BiocManager::install(\"",
        annotation_package, "\")."
      )
    }
    # getAnnotation's internal dispatch looks the annotation object up
    # by `package:<name>` on the search path, so requireNamespace
    # alone is not enough -- we need library() to attach. The single
    # exported object has the same name as the package itself.
    #
    # Attaching the package loads its (large) annotation data into memory and,
    # like getAnnotation() below, can be slow or stall on a constrained box.
    # Announce it -- with the install path, its on-disk size, and current
    # memory -- before the call, so the last line a user sees during a stall
    # names exactly what was loading.
    pkg_dir <- tryCatch(find.package(annotation_package),
                        error = function(e) NA_character_)
    if (verbose >= 2L) {
      message(sprintf(
        "[processMethArrayExp] loading annotation package '%s' (%s on disk at %s); %s",
        annotation_package, .path_size(pkg_dir),
        if (is.na(pkg_dir)) "unknown" else pkg_dir, .mem_report()
      ))
    }
    atime <- system.time({
      suppressPackageStartupMessages(
        library(annotation_package, character.only = TRUE)
      )
      ann_obj <- get(annotation_package)
      ann_df <- getAnnotation(ann_obj)
    })[3]
    sex_probes <- ann_df$Name[ann_df$chr %in% c("chrX", "chrY")]
    if (verbose >= 2L) {
      message(sprintf(
        "[processMethArrayExp] annotation '%s' loaded in %.1fs: %d sex-chromosome probes of %d; %s",
        annotation_package, atime, length(sex_probes), nrow(ann_df),
        .mem_report()
      ))
    }
    # Free the (~1.8 GB for EPIC) annotation frame now -- only sex_probes is
    # needed downstream. Held, it would sit in the master heap through the
    # forked streaming pass in processMethArray, where every MulticoreParam
    # worker copy-on-write-inflates against the master heap (a worker's GC
    # touches the whole heap), multiplying it by the worker count.
    rm(ann_df, ann_obj)
  }

  all_drop <- unique(c(sex_probes, as.character(drop_probes)))
  if (length(all_drop) == 0L) all_drop <- NULL

  # ---- delegate to processMethArray ------------------------------------
  res <- processMethArray(
    basenames             = files,
    gds_path              = gds_path,
    pData                 = targets,
    drop_probes           = all_drop,
    sample_detP_threshold = sample_detP_threshold,
    probe_detP_threshold  = probe_detP_threshold,
    compress              = compress,
    verbose               = verbose,
    BPPARAM               = BPPARAM
  )

  # Filter the supplied samplesheet by the sample-QC mask so the caller
  # gets a one-to-one targets frame for the surviving GDS columns.
  res$targets <- targets[res$keepSamples, , drop = FALSE]
  invisible(res)
}
