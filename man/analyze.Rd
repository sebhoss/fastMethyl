\name{analyze}
\alias{analyze}
\title{
  Preprocess an Illumina methylation cohort and (optionally) analyse it, with
  guaranteed GDS cleanup.
}
\description{
  The package's single entry point. One call validates the inputs and builds a
  QC'd, LZ4_RA-compressed, bigmelon-compatible GDS from raw IDATs, optionally
  hands the open GDS handle to a user-supplied function \code{FUN}, and
  \strong{guarantees the GDS is closed afterwards} -- on normal return and on
  error.

  With \code{FUN} supplied this is the "loan pattern": you provide a closure
  describing the analysis and \code{analyze} owns the resource lifecycle, so a
  leftover open handle can never block the next run even if your analysis
  errors (the close comes from \code{\link[base]{on.exit}}). With \code{FUN}
  omitted, \code{analyze} simply builds the GDS and returns its path and
  metadata.

  \code{analyze} is deliberately unopinionated about what to do with the GDS:
  normalisation (\code{dasen}) and outlier detection (\code{outlyx}) both live
  in \code{FUN}, where you control their order. The statistically correct order
  is \code{outlyx} on the raw betas \emph{first}, then \code{dasen} -- see
  Examples.
}
\usage{
analyze(..., FUN = NULL,
        .preprocess = runPreprocess,
        .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
        .close = gdsfmt::closefn.gds)
}
\arguments{
  \item{\dots}{The preprocessing inputs: \code{dataDirectory},
    \code{nonSpecificProbesPath}, \code{targetPattern}, \code{datasetClass},
    \code{annotationPackage}, and optionally \code{sampleDetPThreshold},
    \code{probeDetPThreshold}, \code{forceRebuild}, \code{readerWorkers},
    \code{BPPARAM}, and \code{verbose} (\code{0L}/\code{1L}/\code{2L};
    \code{2L} logs each external-data step with on-disk size and a memory
    snapshot).}
  \item{FUN}{Optional. A function \code{function(gds, res)} run with the open
    GDS handle (\code{gds}) and the preprocessing result list (\code{res}, see
    Value). Its return value becomes \code{analyze}'s return value. It must
    return \emph{materialised} results (matrices, data.frames, model objects),
    not the open \code{gds} handle -- the handle is closed by the time
    \code{analyze} returns. If omitted (\code{NULL}), \code{analyze} does no
    analysis and returns \code{res}.}
  \item{.preprocess, .open, .close}{Injection points for the heavy operations,
    defaulting to the production implementations (the internal preprocessing
    pipeline; \code{openfn.gds} read-write; \code{closefn.gds}). Exposed so the
    control flow can be unit-tested without IDATs or Bioconductor packages;
    production callers should not pass these.}
}
\value{
  When \code{FUN} is supplied, whatever \code{FUN} returns. When \code{FUN} is
  omitted, invisibly a list with \code{gds_path} (the QC'd GDS), \code{targets}
  (the post-QC samplesheet), \code{keepSamples}, and \code{keepProbes}.
}
\details{
  A worked example pipeline ships with the package; retrieve it with
  \code{system.file("scripts", "pipeline.R", package = "fastMethyl")}.
}
\examples{
\dontrun{
library(bigmelon)   # dasen / outlyx / prcomp on the open GDS

# Full analysis: outlier detection on the RAW betas first, then dasen, then
# PCA. The GDS is closed for you afterwards, even if FUN errors.
res <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",
  targetPattern         = "batch1",
  datasetClass          = "my_study",
  annotationPackage     = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  FUN = function(gds, res) {
    outliers <- outlyx(gds, plot = FALSE)              # 1. detect on raw betas
    dasen(gds, node = "normbetas")                     # 2. normalise
    pca <- prcomp(gds, node.name = "normbetas", method = "quick")  # 3. analyse
    list(outliers = outliers, pca = pca)
  })

# Just build the QC'd GDS (no analysis):
res <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",
  targetPattern         = "batch1",
  datasetClass          = "my_study",
  annotationPackage     = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
res$gds_path   # path to the QC'd GDS
}
}
\author{
  Sebastian Hoss.
}
