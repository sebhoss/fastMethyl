# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# analyze(): the package's sole entry point and a loan-pattern wrapper. One call
# runs the preprocessing pipeline (config in, QC'd GDS out), optionally hands
# the open GDS handle to a user-supplied FUN, and *guarantees the GDS is closed
# afterwards* -- on normal return AND on error.
#
# It is deliberately unopinionated about what to do with the GDS: normalisation
# (dasen) and outlier detection (outlyx) both live in FUN, where the caller
# controls their order -- and the statistically correct order is outlyx on the
# raw betas first, then dasen (dasen can mask the technical artefacts outlyx
# keys on; removing flagged samples first also keeps them out of dasen's
# quantile reference). With FUN omitted, analyze() just returns the freshly
# built QC'd GDS's path and metadata -- the "build the GDS" path.
#
# The close guarantee comes from on.exit() in .withGds(), not from the closure,
# so a leftover open handle can never block the next run even if FUN errors.
#
# The heavy operations are injected (.preprocess / .open / .close) with
# production defaults, so the control flow can be unit-tested without IDATs,
# minfi, gdsfmt, or bigmelon.
#
# At verbose = 2 (read from ..., the same level passed to preprocessing) each
# external-data touch -- preprocessing, opening the GDS, the analysis function
# -- is logged with the file's on-disk size and a memory snapshot, so a slow or
# stalling step is visible while it happens instead of being a silent hang.

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

# Exported -----------------------------------------------------------------

analyze <- function(..., FUN = NULL,
                    .preprocess = runPreprocess,
                    .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
                    .close = gdsfmt::closefn.gds) {
  # FUN is optional. When omitted, analyze() just builds the GDS and returns the
  # preprocessing result; when supplied it must be a function of the open GDS.
  if (!is.null(FUN) && !is.function(FUN)) {
    stop("`FUN`, when supplied, must be a function of the open GDS, e.g. ",
         "function(gds, res) { ... }; omit it to just build the GDS.",
         call. = FALSE)
  }

  # `verbose` is read out of ... (not a formal) so it still flows through to
  # .preprocess unchanged; here it only gates analyze's own diagnostics.
  dots <- list(...)
  verbose <- .normalize_verbose(if (is.null(dots$verbose)) 0L else dots$verbose)
  vlog <- function(fmt, ...) {
    if (verbose >= 2L) message(sprintf(paste0("[analyze] ", fmt), ...))
  }

  vlog("starting; %s", .mem_report())
  res <- .preprocess(...)
  vlog("preprocessing complete: GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())

  # No FUN: preprocessing already built and closed the GDS, so just hand back
  # its path and metadata.
  if (is.null(FUN)) {
    return(invisible(res))
  }

  vlog("opening GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())
  .withGds(
    res$gds_path,
    function(gds) {
      vlog("GDS open; running analysis function; %s", .mem_report())
      out <- FUN(gds, res)
      vlog("analysis function returned; %s", .mem_report())
      out
    },
    .open = .open, .close = .close
  )
}
