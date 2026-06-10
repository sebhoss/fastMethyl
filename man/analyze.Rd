\name{analyze}
\alias{analyze}
\title{
  Preprocess an Illumina methylation cohort and (optionally) analyse it, with
  guaranteed GDS cleanup
}
\description{
  The package's primary entry point. One call validates the inputs and builds a
  QC'd, bigmelon-compatible GDS from raw IDATs, then drives the analysis phases
  in a fixed, statistically correct order and \strong{guarantees the GDS is
  closed afterwards} -- on normal return and on error.

  \code{analyze} orchestrates the phases

  \preformatted{  build  ->  outlierRemoval  ->  normalization(keep)  ->  analysis}

  enforcing the order so outliers are excluded from the normalisation reference
  (an outlier left in pollutes the shared quantile reference). \code{outlierRemoval}
  and \code{normalization} are uniform hooks: each takes a built-in method name,
  a user function, or \code{"none"} (the default for both, so a bare call is
  unopinionated). The close comes from \code{\link[base]{on.exit}}, so a leftover
  open handle can never block the next run even if a hook or your analysis errors.
}
\usage{
analyze(..., analysis = NULL,
        outlierRemoval = "none",
        normalization = "none",
        .preprocess = runPreprocess,
        .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
        .close = gdsfmt::closefn.gds)
}
\arguments{
  \item{\dots}{The preprocessing inputs: \code{dataDirectory} (folder the IDATs
    live under), \code{samplesheet} (path to the samplesheet CSV; needs a
    \code{Basename} column), \code{crossReactiveProbes} (path to a CSV of probes
    to drop; needs a \code{TargetID} column), \code{gdsOutput} (path to the
    output GDS file, e.g. \code{".../study.gds"}), and \code{annotationPackage};
    plus optionally \code{sampleDetPThreshold}, \code{probeDetPThreshold},
    \code{rebuild} (\code{"none"}/\code{"normalize"}/\code{"build"}/\code{"all"},
    or a logical), \code{readerWorkers}, \code{BPPARAM}, \code{compress}, and
    \code{verbose} (\code{0L}/\code{1L}/\code{2L}).}
  \item{analysis}{Optional. A function \code{function(gds, res)} run with the
    open GDS handle and the preprocessing result list (\code{res}, see Value),
    after the outlier and normalisation hooks. Its return value becomes
    \code{analyze}'s return value, and must be \emph{materialised} (matrices,
    data.frames, model objects), not the open \code{gds} handle. If omitted
    (\code{NULL}) and both hooks are \code{"none"}, \code{analyze} returns
    \code{res}.}
  \item{outlierRemoval}{Outlier step. A method name (\code{"outlyx"}), a function
    \code{function(gds, res)} returning a logical keep mask over the samples, or
    \code{"none"}. The mask is fed to \code{normalization} and exposed to
    \code{analysis} as \code{res$outlierKeep}.}
  \item{normalization}{Normalisation step. A method name (\code{"dasen"} -- the
    streaming \code{gdsDasen}), a function \code{function(gds, res, keep)}
    that writes a normalised node, or \code{"none"}.}
  \item{.preprocess, .open, .close}{Injection points for the heavy operations,
    defaulting to the production implementations. Exposed so the orchestration
    can be unit-tested without IDATs or Bioconductor packages; production callers
    should not pass these.}
}
\value{
  When \code{analysis} is supplied, whatever it returns. Otherwise, invisibly the
  preprocessing result: a list with \code{gds_path} (the QC'd GDS), \code{targets}
  (the post-QC samplesheet), \code{keepSamples}, \code{keepProbes}, and
  \code{rebuildDownstream}.
}
\details{
  \code{forceRebuild} (logical) is a deprecated alias for \code{rebuild}; the
  preprocessing argument names \code{targetPattern}, \code{datasetClass}, and
  \code{nonSpecificProbesPath} are deprecated aliases for \code{samplesheet},
  \code{gdsOutput}, and \code{crossReactiveProbes}. A worked example pipeline
  ships with the package; retrieve it with
  \code{system.file("scripts", "pipeline.R", package = "fastMethyl")}.
}
\examples{
\dontrun{
library(bigmelon)   # prcomp / estimateCellCounts.gds on the open GDS

# Full analysis: built-in outlier removal + streaming normalisation run before
# your analysis, in the correct order. The GDS is closed for you afterwards.
res <- analyze(
  dataDirectory       = "/path/to/idats",
  samplesheet         = "/path/to/samplesheet.csv",
  crossReactiveProbes = "/path/to/cross-reactive.csv",
  gdsOutput           = "/path/to/my_study.gds",
  annotationPackage   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  outlierRemoval      = "outlyx",
  normalization       = "dasen",
  analysis = function(gds, res) {
    prcomp(gds, node.name = "normbetas", method = "quick")
  })

# Just build the QC'd GDS (no hooks, no analysis):
res <- analyze(
  dataDirectory       = "/path/to/idats",
  samplesheet         = "/path/to/samplesheet.csv",
  crossReactiveProbes = "/path/to/cross-reactive.csv",
  gdsOutput           = "/path/to/my_study.gds",
  annotationPackage   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
res$gds_path   # path to the QC'd GDS
}
}
\author{
  Sebastian Hoss.
}
