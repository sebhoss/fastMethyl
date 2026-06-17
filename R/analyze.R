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
          .userStop("outlierRemoval = \"outlyx\" needs the bigmelon package. ",
                    "Install it with BiocManager::install(\"bigmelon\"), or pass ",
                    "outlierRemoval = \"none\" / your own function instead.")
        }
        out <- wateRmelon::outlyx(gds, plot = FALSE)
        !out$outliers
      },
      .userStopf(
        "unknown outlierRemoval method '%s'; use 'outlyx', 'none', or a function(gds, res).",
        spec
      )
    ))
  }
  .userStop("`outlierRemoval` must be a method name (\"outlyx\" / \"none\"), ",
            "a function(gds, res), or \"none\".")
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
      .userStopf(
        "unknown normalization method '%s'; use 'dasen', 'none', or a function(gds, res, keep).",
        spec
      )
    ))
  }
  .userStop("`normalization` must be a method name (\"dasen\" / \"none\"), ",
            "a function(gds, res, keep), or \"none\".")
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
    .userStop("`rebuild` must be a stage name (\"build\" / \"normalize\"), ",
              "\"none\" / \"all\", or a single logical.")
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
    .userStopf(
      "unknown rebuild stage '%s'; use 'none' / 'build' / 'normalize' / 'all', or a logical.",
      spec
    )
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
  #
  # Detection is by NAME PRESENCE (`has()`), not by value. A user who passes a
  # deprecated name explicitly as NULL -- e.g. the legacy `FUN = NULL` "just
  # build the GDS" form -- must still have it stripped from `dots`; checking
  # `!is.null(dots$FUN)` would skip the strip (the value *is* NULL), leaving FUN
  # in `dots` to be forwarded to .preprocess, which has no such formal and would
  # abort with an opaque "unused argument (FUN = NULL)".
  .warnDep <- function(old, new) {
    warning(sprintf("`%s` is deprecated; use `%s`.", old, new), call. = FALSE)
  }
  has <- function(k) k %in% names(dots)
  if (has("FUN")) {
    .warnDep("FUN", "analysis")
    if (is.null(analysis)) analysis <- dots$FUN
    dots$FUN <- NULL
  }
  if (has("normalize")) {
    .warnDep("normalize", "normalization")
    if (identical(normalization, "none")) normalization <- dots$normalize
    dots$normalize <- NULL
  }
  # nonSpecificProbesPath -> crossReactiveProbes is a pure rename.
  if (has("nonSpecificProbesPath")) {
    .warnDep("nonSpecificProbesPath", "crossReactiveProbes")
    if (!has("crossReactiveProbes")) {
      dots$crossReactiveProbes <- dots$nonSpecificProbesPath
    }
    dots$nonSpecificProbesPath <- NULL
  }
  # targetPattern named a filename fragment under dataDirectory; the current
  # `samplesheet` is the path itself, so reconstruct the old convention to keep
  # legacy calls working.
  if (has("targetPattern")) {
    .warnDep("targetPattern", "samplesheet")
    if (!has("samplesheet") && !is.null(dots$targetPattern)) {
      sheet <- paste0("samplesheet_", dots$targetPattern, ".csv")
      dots$samplesheet <- file.path(dots$dataDirectory, sheet)
    }
    dots$targetPattern <- NULL
  }
  # datasetClass was a prefix (the GDS was "<datasetClass>.gds"); the current
  # `gdsOutput` is the full path, so append the extension for legacy calls.
  if (has("datasetClass")) {
    .warnDep("datasetClass", "gdsOutput")
    if (!has("gdsOutput") && !is.null(dots$datasetClass)) {
      dots$gdsOutput <- paste0(dots$datasetClass, ".gds")
    }
    dots$datasetClass <- NULL
  }

  if (!is.null(analysis) && !is.function(analysis)) {
    .userStop("`analysis`, when supplied, must be a function of the open GDS, ",
              "e.g. function(gds, res) { ... }; omit it to just build the GDS.")
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

  # Guard the forwarded arguments. The real .preprocess (runPreprocess) has fixed
  # formals and no `...`, so a stray or misspelled name in `dots` would otherwise
  # surface as R's opaque "unused argument" error from inside do.call. Catch it
  # here and name both the offender and the valid arguments. A mock .preprocess
  # that declares `...` accepts anything, so the check is skipped for it.
  preFormals <- names(formals(.preprocess))
  if (!("..." %in% preFormals)) {
    named <- names(dots)
    unknown <- setdiff(named[nzchar(named)], preFormals)
    if (length(unknown) > 0L) {
      # Hide internal (dotted) formals and the deprecated forceRebuild from the
      # advertised list -- forceRebuild is named below under the aliases instead.
      hide <- c(grep("^\\.", preFormals, value = TRUE), "forceRebuild")
      valid <- setdiff(preFormals, hide)
      .userStopf(
        paste0(
          "unknown argument(s) to analyze(): %s.\n",
          "  Valid preprocessing arguments are: %s.\n",
          "  analyze() also takes: analysis, outlierRemoval, normalization, ",
          "rebuild, verbose.\n",
          "  Using an old name? Deprecated aliases: FUN, normalize, ",
          "nonSpecificProbesPath, targetPattern, datasetClass, forceRebuild. ",
          "See ?analyze."
        ),
        paste(sprintf("`%s`", unknown), collapse = ", "),
        paste(valid, collapse = ", ")
      )
    }
  }
  vlog <- function(fmt, ...) {
    if (verbose >= 2L) message(sprintf(paste0("[analyze] ", fmt), ...))
  }

  vlog("starting; %s", .mem_report())
  # Every heavy step runs inside .withPhase: an unexpected failure anywhere in
  # the call graph (minfi/gdsfmt/BiocParallel/base R) is caught and re-raised
  # with what analyze() was doing, a remedy to try, and how to report a genuine
  # bug. Errors we raise ourselves as fastMethyl_user_error already carry an
  # actionable message and pass through unchanged.
  res <- .withPhase(
    "reading the IDATs and building the QC'd GDS",
    do.call(.preprocess, dots),
    hint = paste0(
      "Re-run with verbose = 2 to see the last step before the failure. ",
      "Confirm dataDirectory, samplesheet and crossReactiveProbes point to the ",
      "intended files, that annotationPackage matches your array, and that the ",
      "machine has enough free memory (lower readerWorkers if it ran out)."
    )
  )
  vlog("preprocessing complete: GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())

  # Nothing to run on the open GDS: preprocessing already built and closed it, so
  # hand back its path and metadata.
  if (is.null(analysis) && is.null(outlierFn) && is.null(normalizeFn)) {
    return(invisible(res))
  }

  vlog("opening GDS %s (%s); %s",
       res$gds_path, .path_size(res$gds_path), .mem_report())
  # Wrap the opener so a corrupt/locked/incomplete GDS reports with context; the
  # close still happens via .withGds's on.exit, so the wrapper cannot leak a handle.
  openWithContext <- function(p) {
    .withPhase(
      sprintf("opening the QC'd GDS file at %s", p),
      .open(p),
      hint = paste0("The GDS may be incomplete or corrupt (for example from an ",
                    "interrupted earlier run), or open in another process. ",
                    "Re-run with rebuild = \"build\" to rebuild it from the IDATs.")
    )
  }
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
        keep <- .withPhase(
          "removing outliers",
          outlierFn(gds, res),
          hint = paste0(
            "Built-in outlierRemoval = \"outlyx\" needs the bigmelon package ",
            "and a GDS with methylated/unmethylated nodes. If you passed your ",
            "own function, check it takes (gds, res) and returns a logical keep ",
            "mask, one entry per sample."
          )
        )
        if (is.logical(keep)) {
          vlog("outlier removal kept %d/%d samples", sum(keep), length(keep))
        }
      }
      # Expose the survivor mask to the normalization and analysis steps without
      # changing their arity (steps that ignore it see the extra field harmlessly).
      res$outlierKeep <- keep
      if (!is.null(normalizeFn)) {
        vlog("running normalization; %s", .mem_report())
        .withPhase(
          "normalising the methylation values",
          normalizeFn(gds, res, keep),
          hint = paste0(
            "Built-in normalization = \"dasen\" needs methylated/unmethylated ",
            "intensity nodes and an fData$Type column of \"I\"/\"II\". If you ",
            "passed your own function, check it takes (gds, res, keep)."
          )
        )
        vlog("normalization complete; %s", .mem_report())
      }
      if (!is.null(analysis)) {
        vlog("running analysis function; %s", .mem_report())
        out <- .withPhase(
          "running your analysis function",
          analysis(gds, res),
          hint = paste0(
            "This error came from the analysis function you supplied. Check it ",
            "takes (gds, res) and handles this GDS; re-run with verbose = 2 to ",
            "see how far it got."
          )
        )
        vlog("analysis function returned; %s", .mem_report())
        out
      } else {
        invisible(res)
      }
    },
    .open = openWithContext, .close = .close
  )
}
