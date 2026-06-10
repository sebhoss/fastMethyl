# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# analyze(): the package's primary entry point and a loan-pattern orchestrator.
# One call runs the preprocessing pipeline (config in, QC'd GDS out), then drives
# the analysis phases in the statistically correct order and *guarantees the GDS
# is closed afterwards* -- on normal return AND on error.
#
# The order is the whole point. Outlier detection must precede normalisation
# (dasen's between-array quantile reference must be built from outlier-free
# samples; an outlier left in pollutes the reference every sample is mapped onto),
# so analyze() runs them in a fixed sequence -- build, then outlierRemoval, then
# normalize over the survivors, then FUN -- and the caller cannot get it wrong. `outlierRemoval` and `normalize` take
# the SAME shape -- a built-in method name ("outlyx" / "dasen"), a user function,
# or "none" -- so a study swaps in its own step or accepts the fast default
# without restructuring anything. Defaults are "none" for both, so an analyze()
# call with neither set is unopinionated and behaves as a plain build + FUN.
#
#   outlierRemoval(gds, res) -> a logical keep mask over the GDS samples (TRUE =
#     survivor). Built-in "outlyx" returns !outlyx(gds)$outliers.
#   normalize(gds, res, keep) -> writes a normalised node (default "normbetas").
#     Built-in "dasen" is gdsDasen(gds, keep = keep): the reference is built from
#     the survivors, but normbetas is written for every sample, so the GDS stays
#     rectangular and a flagged-but-kept outlier still gets a column.
#
# The keep mask is exposed to FUN as res$outlierKeep. With FUN and both hooks
# omitted, analyze() just returns the freshly built QC'd GDS's path and metadata.
#
# The close guarantee comes from on.exit() in .withGds(), not from the closure,
# so a leftover open handle can never block the next run even if a hook or FUN
# errors. The heavy operations are injected (.preprocess / .open / .close) with
# production defaults, so the control flow can be unit-tested without IDATs,
# minfi, gdsfmt, or bigmelon.
#
# `verbose` defaults to 0 (quiet) and analyze is its single authority: the
# resolved level is forwarded explicitly to preprocessing, so an unspecified
# verbose is silent everywhere. At verbose = 2 each external-data touch --
# preprocessing, opening the GDS, each hook, FUN -- is logged with a memory
# snapshot, so a slow or stalling step is visible while it happens.

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

# Resolve an `outlierRemoval` spec to a function(gds, res) -> keep mask, or NULL
# to skip. A function is used as-is; a method name selects a built-in; "none"
# (the default) skips. Built-in "outlyx" needs bigmelon, so the dependency error
# is raised only when that method is actually requested.
.resolveOutlierRemoval <- function(spec) {
  if (is.null(spec) || is.function(spec)) {
    return(spec)
  }
  if (is.character(spec) && length(spec) == 1L && !is.na(spec)) {
    return(switch(spec,
      none = NULL,
      outlyx = function(gds, res) {
        # The outlyx generic lives in wateRmelon; bigmelon registers the
        # gds.class method. Loading bigmelon's namespace registers that method,
        # so the wateRmelon generic then dispatches on the open GDS handle.
        if (!requireNamespace("bigmelon", quietly = TRUE)) {
          stop("outlierRemoval = \"outlyx\" needs the bigmelon package.",
               call. = FALSE)
        }
        out <- wateRmelon::outlyx(gds, plot = FALSE)
        !out$outliers
      },
      stop(sprintf("unknown outlierRemoval method '%s'; use 'outlyx', 'none', or a function.", spec), call. = FALSE)
    ))
  }
  stop("`outlierRemoval` must be a method name, a function(gds, res), or \"none\".",
       call. = FALSE)
}

# Resolve a `normalize` spec to a function(gds, res, keep) that writes a node, or
# NULL to skip. Built-in "dasen" is the streaming gdsDasen with the survivor mask.
.resolveNormalize <- function(spec) {
  if (is.null(spec) || is.function(spec)) {
    return(spec)
  }
  if (is.character(spec) && length(spec) == 1L && !is.na(spec)) {
    return(switch(spec,
      none = NULL,
      # Cache-aware: skip when a normbetas node already exists and nothing
      # upstream rebuilt this run (res$rebuildDownstream), so a re-run that did
      # not touch the build/normalize stages reuses it. analyze() sets
      # res$rebuildDownstream from the resolved `rebuild` request.
      dasen = function(gds, res, keep) {
        if (isTRUE(res$rebuildDownstream) ||
              !("normbetas" %in% gdsfmt::ls.gdsn(gds))) {
          gdsDasen(gds, node = "normbetas", keep = keep)
        }
        invisible(NULL)
      },
      stop(sprintf("unknown normalization method '%s'; use 'dasen', 'none', or a function.", spec), call. = FALSE)
    ))
  }
  stop("`normalize` must be a method name, a function(gds, res, keep), or \"none\".",
       call. = FALSE)
}

# Resolve the staged `rebuild` request to per-stage flags. The pipeline is a
# linear dependency chain (build -> normalize -> the FUN's own outputs), so
# naming a stage rebuilds it AND every stage after it. `rebuild` accepts a stage
# name ("build" / "normalize"), "none"/"all", or a logical (FALSE = "none",
# TRUE = "all"); `forceRebuild` is the deprecated logical alias. Returns a list
# of booleans, one per analyze-owned stage.
.resolveRebuild <- function(rebuild = NULL, forceRebuild = NULL) {
  spec <- if (!is.null(rebuild)) {
    rebuild
  } else if (!is.null(forceRebuild)) {
    warning("`forceRebuild` is deprecated; use `rebuild` ",
            "(\"none\"/\"build\"/\"normalize\"/\"all\", or a logical).",
            call. = FALSE)
    forceRebuild
  } else {
    "none"
  }
  if (is.logical(spec) && length(spec) == 1L && !is.na(spec)) {
    spec <- if (spec) "all" else "none"
  }
  if (!(is.character(spec) && length(spec) == 1L && !is.na(spec))) {
    stop("`rebuild` must be a stage name (\"build\" / \"normalize\"), ",
         "\"none\" / \"all\", or a single logical.", call. = FALSE)
  }
  # Ordered analyze-owned stages; a named stage selects it and all following.
  stages <- c("build", "normalize")
  set <- if (spec == "none") {
    character(0L)
  } else if (spec %in% c("all", "build")) {
    stages
  } else if (spec == "normalize") {
    "normalize"
  } else {
    msg <- sprintf("unknown rebuild stage '%s'; use 'none'/'build'/'normalize'/'all' or a logical.", spec)
    stop(msg, call. = FALSE)
  }
  list(build = "build" %in% set, normalize = "normalize" %in% set)
}

# Exported -----------------------------------------------------------------

analyze <- function(..., analysis = NULL,
                    outlierRemoval = "none",
                    normalization = "none",
                    .preprocess = runPreprocess,
                    .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
                    .close = gdsfmt::closefn.gds) {
  dots <- list(...)

  # Deprecated argument aliases -> remap to the current names with a warning.
  # The two analysis hooks (FUN -> analysis, normalize -> normalization) are
  # formals, so their old names arrive via `...`; the preprocessing renames are
  # forwarded onward through `dots`. An explicit new-name value always wins.
  .warnDep <- function(old, new) {
    warning(sprintf("`%s` is deprecated; use `%s`.", old, new), call. = FALSE)
  }
  if (!is.null(dots$FUN)) {
    .warnDep("FUN", "analysis")
    if (is.null(analysis)) analysis <- dots$FUN
    dots$FUN <- NULL
  }
  if (!is.null(dots$normalize)) {
    .warnDep("normalize", "normalization")
    if (identical(normalization, "none")) normalization <- dots$normalize
    dots$normalize <- NULL
  }
  # nonSpecificProbesPath -> crossReactiveProbes is a pure rename.
  if (!is.null(dots$nonSpecificProbesPath)) {
    .warnDep("nonSpecificProbesPath", "crossReactiveProbes")
    if (is.null(dots$crossReactiveProbes)) {
      dots$crossReactiveProbes <- dots$nonSpecificProbesPath
    }
    dots$nonSpecificProbesPath <- NULL
  }
  # targetPattern named a filename fragment under dataDirectory; the current
  # `samplesheet` is the path itself, so reconstruct the old convention to keep
  # legacy calls working.
  if (!is.null(dots$targetPattern)) {
    .warnDep("targetPattern", "samplesheet")
    if (is.null(dots$samplesheet)) {
      sheet <- paste0("samplesheet_", dots$targetPattern, ".csv")
      dots$samplesheet <- file.path(dots$dataDirectory, sheet)
    }
    dots$targetPattern <- NULL
  }
  # datasetClass was a prefix (the GDS was "<datasetClass>.gds"); the current
  # `gdsOutput` is the full path, so append the extension for legacy calls.
  if (!is.null(dots$datasetClass)) {
    .warnDep("datasetClass", "gdsOutput")
    if (is.null(dots$gdsOutput)) dots$gdsOutput <- paste0(dots$datasetClass, ".gds")
    dots$datasetClass <- NULL
  }

  if (!is.null(analysis) && !is.function(analysis)) {
    stop("`analysis`, when supplied, must be a function of the open GDS, e.g. ",
         "function(gds, res) { ... }; omit it to just build the GDS.",
         call. = FALSE)
  }
  outlierFn <- .resolveOutlierRemoval(outlierRemoval)
  normalizeFn <- .resolveNormalize(normalization)

  # Resolve `verbose` once and forward the normalised value explicitly to
  # .preprocess, so analyze is the single authority on verbosity. Writing it back
  # into `dots` keeps the raw value from also reaching .preprocess via `...`.
  verbose <- .normalize_verbose(if (is.null(dots$verbose)) 0L else dots$verbose)
  dots$verbose <- verbose

  # Resolve the staged `rebuild` request (with the deprecated `forceRebuild`
  # alias) to per-stage flags. The build stage is realised through .preprocess's
  # forceRebuild, so translate to that boolean and drop `rebuild` from the
  # forwarded args; the normalize stage is enforced below via res$rebuildDownstream.
  rb <- .resolveRebuild(dots$rebuild, dots$forceRebuild)
  dots$rebuild <- NULL
  dots$forceRebuild <- rb$build
  vlog <- function(fmt, ...) {
    if (verbose >= 2L) message(sprintf(paste0("[analyze] ", fmt), ...))
  }

  vlog("starting; %s", .mem_report())
  res <- do.call(.preprocess, dots)
  vlog("preprocessing complete: GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())

  # Nothing to run on the open GDS: preprocessing already built and closed it, so
  # hand back its path and metadata.
  if (is.null(analysis) && is.null(outlierFn) && is.null(normalizeFn)) {
    return(invisible(res))
  }

  vlog("opening GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())
  .withGds(
    res$gds_path,
    function(gds) {
      # Cascade rebuild to the downstream stages: normbetas (and the FUN's own
      # outputs) are stale if the normalize stage was named OR the GDS was
      # actually rebuilt (a forced build or a cache miss). The cache-aware
      # "dasen" hook and the user's FUN both read res$rebuildDownstream.
      res$rebuildDownstream <- rb$normalize || isTRUE(res$rebuilt)
      keep <- NULL
      if (!is.null(outlierFn)) {
        vlog("running outlier removal; %s", .mem_report())
        keep <- outlierFn(gds, res)
        if (is.logical(keep)) {
          vlog("outlier removal kept %d/%d samples", sum(keep), length(keep))
        }
      }
      # Expose the survivor mask to the normalization and analysis steps without
      # changing their arity (steps that ignore it see the extra field harmlessly).
      res$outlierKeep <- keep
      if (!is.null(normalizeFn)) {
        vlog("running normalization; %s", .mem_report())
        normalizeFn(gds, res, keep)
        vlog("normalization complete; %s", .mem_report())
      }
      if (!is.null(analysis)) {
        vlog("running analysis function; %s", .mem_report())
        out <- analysis(gds, res)
        vlog("analysis function returned; %s", .mem_report())
        out
      } else {
        invisible(res)
      }
    },
    .open = .open, .close = .close
  )
}
