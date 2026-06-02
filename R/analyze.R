# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# analyze(): a loan-pattern entry point. It runs preprocessing, opens the
# resulting GDS, optionally normalises it, hands the open handle to a
# user-supplied FUN, and guarantees the GDS is closed afterwards -- on normal
# return AND on error -- then returns FUN's value.
#
# The closure (FUN) is how the caller supplies their analysis; the close
# guarantee comes from on.exit() in .withGds(), not from the closure itself.
#
# The heavy operations are injected (.preprocess / .open / .close /
# .normalize) with production defaults, so the control flow can be unit-tested
# without IDATs, minfi, gdsfmt, or bigmelon.
#
# At verbose = 2 (read from ..., the same level passed to preprocessing) each
# external-data touch is logged with the file's on-disk size and a memory
# snapshot before and after, so a slow or stalling load is visible while it
# happens instead of being a silent hang.

# Open `path`, run `fun(gds)` with the GDS guaranteed closed afterwards (even
# if `fun` errors), and return `fun`'s value. This is the reusable
# resource-management core; on.exit() is what makes the close unconditional.
.withGds <- function(path, fun,
                     .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
                     .close = gdsfmt::closefn.gds) {
  gds <- .open(path)
  on.exit(.close(gds), add = TRUE)
  fun(gds)
}

# Default normaliser: dasen into `node`. bigmelon is a Suggests, so guard on it
# and give an actionable error rather than an opaque one. `.have` is injectable
# so the guard branch is testable without manipulating the search path.
.dasenNormalize <- function(gds, node,
                            .have = function(p) {
                              requireNamespace(p, quietly = TRUE)
                            }) {
  if (!.have("bigmelon")) {
    stop("normalize = TRUE requires the 'bigmelon' package; install it ",
         "or call analyze(..., normalize = FALSE).", call. = FALSE)
  }
  bigmelon::dasen(gds, node = node)
  invisible(node)
}

# Exported -----------------------------------------------------------------

analyze <- function(..., FUN, normalize = TRUE, normNode = "normbetas",
                    .preprocess = runPreprocess,
                    .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
                    .close = gdsfmt::closefn.gds,
                    .normalize = .dasenNormalize) {
  if (!is.function(FUN)) {
    stop("`FUN` must be a function of the open GDS, e.g. ",
         "function(gds, res) { ... }.", call. = FALSE)
  }
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("`normalize` must be a single TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(normNode) || length(normNode) != 1L || is.na(normNode) ||
        !nzchar(normNode)) {
    stop("`normNode` must be a single non-empty string.", call. = FALSE)
  }

  # `verbose` is read out of ... (not a formal) so it still flows through to
  # .preprocess unchanged; here it only gates analyze's own diagnostics. At
  # verbose >= 2 every external-data touch -- preprocessing, opening the GDS,
  # normalising it, and the user's analysis -- is bracketed with the file's
  # on-disk size and a memory snapshot, so a long or stalling step announces
  # itself (and its memory footprint) before it blocks rather than after.
  dots <- list(...)
  verbose <- .normalize_verbose(if (is.null(dots$verbose)) 0L else dots$verbose)
  vlog <- function(fmt, ...) {
    if (verbose >= 2L) message(sprintf(paste0("[analyze] ", fmt), ...))
  }

  vlog("starting; %s", .mem_report())
  res <- .preprocess(...)
  vlog("preprocessing complete: GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())

  vlog("opening GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())
  .withGds(
    res$gds_path,
    function(gds) {
      vlog("GDS open; %s", .mem_report())
      if (normalize) {
        vlog("normalising node '%s' (dasen); %s", normNode, .mem_report())
        .normalize(gds, normNode)
        vlog("normalisation done: GDS now %s; %s",
             .path_size(res$gds_path), .mem_report())
      }
      vlog("running analysis function; %s", .mem_report())
      out <- FUN(gds, res)
      vlog("analysis function returned; %s", .mem_report())
      out
    },
    .open = .open, .close = .close
  )
}
