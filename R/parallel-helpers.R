# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Parallel-reading helpers shared by readMethArray() and processMethArray():
# verbose-level normalisation and the fault-tolerant, progress-reporting
# bpmapply wrapper.

# Verbose level handling -------------------------------------------------------
# Functions take `verbose` as either a logical (back-compat with the historical
# TRUE/FALSE API) or an integer in {0L, 1L, 2L}.
#
#   0 - silent
#   1 - phase summaries + per-sample worker messages (the old TRUE behaviour)
#   2 - adds master-side periodic progress lines (every chunk completes), plus
#       extra context (array type, manifest, memory footprint, QC drops) so the
#       user always sees activity during long-running operations
#
# Backward compatibility: TRUE -> 1L, FALSE -> 0L.
.normalize_verbose <- function(verbose) {
  if (is.logical(verbose)) {
    if (length(verbose) != 1L || is.na(verbose)) {
      stop("`verbose` must be a single TRUE/FALSE or integer 0/1/2")
    }
    return(if (verbose) 1L else 0L)
  }
  if (!is.numeric(verbose) || length(verbose) != 1L || !is.finite(verbose)) {
    stop("`verbose` must be a single TRUE/FALSE or integer 0/1/2")
  }
  v <- as.integer(verbose)
  if (v < 0L || v > 2L) {
    stop("`verbose` must be 0, 1, or 2 (got ", verbose, ")")
  }
  v
}

# One bpmapply call wrapped for fault tolerance. A worker that dies on a
# single bad input (an unreadable IDAT, a transient NFS stall) would
# otherwise discard the work already done on every other input. bptry()
# captures per-element failures instead of aborting; we retry just the
# failed elements once via BPREDO (which re-runs only the not-ok ones),
# and only if a failure persists do we stop -- with a message naming the
# exact inputs that failed rather than an opaque worker traceback.
#
# On full success the return value is identical to a plain bpmapply call,
# so the hot path is unchanged for the common case.
.bpmapply_try <- function(FUN, iter_args, MoreArgs, BPPARAM, label) {
  call_bp <- function(redo) {
    args <- c(
      list(FUN = FUN), iter_args,
      list(
        MoreArgs = MoreArgs,
        SIMPLIFY = FALSE, USE.NAMES = TRUE,
        BPPARAM = BPPARAM
      )
    )
    if (!is.null(redo)) args$BPREDO <- redo
    bptry(do.call(bpmapply, args))
  }

  res <- call_bp(NULL)
  if (!all(bpok(res))) {
    failed <- which(!bpok(res))
    nm <- names(res)
    labels <- if (!is.null(nm)) nm[failed] else as.character(failed)
    message(sprintf(
      "[%s] %d/%d input(s) failed on first pass; retrying once: %s",
      label, length(failed), length(res),
      paste(utils::head(labels, 5L), collapse = ", ")
    ))
    res <- call_bp(res)
  }
  if (!all(bpok(res))) {
    failed <- which(!bpok(res))
    nm <- names(res)
    labels <- if (!is.null(nm)) nm[failed] else as.character(failed)
    n_show <- min(length(failed), 10L)
    # After bptry() a failed element is itself the captured condition.
    first_err <- tryCatch(conditionMessage(res[[failed[1L]]]),
      error = function(e) "<error message unavailable>"
    )
    stop(sprintf(
      "[%s] %d/%d input(s) failed after one retry:\n  %s%s\n  first error: %s",
      label, length(failed), length(res),
      paste(labels[seq_len(n_show)], collapse = "\n  "),
      if (length(failed) > n_show) {
        sprintf("\n  ... and %d more", length(failed) - n_show)
      } else {
        ""
      },
      first_err
    ), call. = FALSE)
  }
  res
}

# Splits a bpmapply call into chunks (so the master sees workers complete
# periodically and can emit a progress line per chunk). Used only when the
# caller asks for verbose >= 2; for verbose < 2 the helper degrades to a
# single bpmapply call so there is no chunking overhead on the hot path.
# Every underlying bpmapply runs through .bpmapply_try, so a single bad
# input is retried once and, if still failing, reported by name.
#
# Signature mirrors bpmapply: FUN, then iterating arguments by name, then
# MoreArgs / BPPARAM. Iterating arguments must be named and equal-length.
.bpmapply_with_progress <- function(FUN, iter_args, MoreArgs = list(),
                                    BPPARAM, verbose, label) {
  stopifnot(is.list(iter_args), length(iter_args) >= 1L)
  n <- length(iter_args[[1L]])
  stopifnot(all(vapply(iter_args, length, integer(1L)) == n))

  if (verbose < 2L) {
    return(.bpmapply_try(
      FUN = FUN, iter_args = iter_args,
      MoreArgs = MoreArgs, BPPARAM = BPPARAM,
      label = label
    ))
  }

  nworkers <- max(1L, bpnworkers(BPPARAM))
  chunk_size <- max(1L, min(n, nworkers * 4L))
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  result <- vector("list", n)
  result_names <- if (!is.null(names(iter_args[[1L]]))) {
    names(iter_args[[1L]])
  } else {
    character(n)
  }
  names(result) <- result_names
  start_t <- Sys.time()

  for (chunk_idx in chunks) {
    chunk_args <- lapply(iter_args, function(a) a[chunk_idx])
    chunk_res <- .bpmapply_try(
      FUN = FUN, iter_args = chunk_args,
      MoreArgs = MoreArgs, BPPARAM = BPPARAM,
      label = label
    )
    for (k in seq_along(chunk_idx)) result[[chunk_idx[k]]] <- chunk_res[[k]]
    done <- max(chunk_idx)
    elapsed <- as.numeric(Sys.time() - start_t, units = "secs")
    rate <- done / max(elapsed, 0.001)
    eta <- (n - done) / max(rate, 0.001)
    message(sprintf(
      "[%s] progress %d/%d (%5.1f%%) elapsed %.1fs, ETA %.1fs (%.1f file/s)",
      label, done, n, 100 * done / n, elapsed, eta, rate
    ))
  }
  result
}

# Friendly memory-footprint string for an in-memory R object (uses
# object.size() to walk the structure). For an on-disk byte count from
# file.size() use .bytes_human_count(), which skips object.size and
# just formats the number directly.
.bytes_human <- function(x) {
  bytes <- as.numeric(object.size(x))
  format(structure(bytes, class = "object_size"),
    units = "auto",
    standard = "SI"
  )
}
.bytes_human_count <- function(bytes) {
  format(structure(as.numeric(bytes), class = "object_size"),
    units = "auto", standard = "SI"
  )
}

# Process memory snapshot for verbose diagnostics: the resident set size and
# its high-water mark. When an external data load (annotation package,
# getAnnotation frame) balloons memory or stalls in swap, logging this before
# and after the load is what turns "it hung with no output" into a diagnosable
# event. Reads /proc/self/status on Linux (the figure that actually matters,
# including memory held by loaded packages); elsewhere falls back to the R heap
# reported by gc().
.mem_report <- function() {
  status <- "/proc/self/status"
  if (file.exists(status)) {
    lines <- readLines(status, warn = FALSE)
    grab <- function(key) {
      hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
      if (length(hit) == 0L) {
        return(NA_real_)
      }
      as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", hit[[1L]])) * 1024
    }
    rss <- grab("VmRSS")
    if (!is.na(rss)) {
      return(sprintf(
        "RSS %s (peak %s)", .bytes_human_count(rss),
        .bytes_human_count(grab("VmHWM"))
      ))
    }
  }
  g <- gc(verbose = FALSE)
  sprintf("R heap %.0f MB", sum(g[, 2L]))
}

# Human-readable on-disk size of a file or directory (directories summed
# recursively), for verbose load messages -- so a user sees how much data a
# load is about to pull in. Returns "missing" / "unknown" rather than erroring
# on an absent or unusable path, since this only ever feeds a log line.
.path_size <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    return("unknown")
  }
  if (!file.exists(path)) {
    return("missing")
  }
  bytes <- if (isTRUE(file.info(path)$isdir)) {
    sum(file.size(list.files(path, recursive = TRUE, full.names = TRUE)),
      na.rm = TRUE
    )
  } else {
    file.size(path)
  }
  .bytes_human_count(bytes)
}
