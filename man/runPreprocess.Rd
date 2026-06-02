\name{runPreprocess}
\alias{runPreprocess}
\title{
  End-to-end pre-processing from per-study arguments to a QC'd GDS.
}
\description{
  Takes every input a methylation pre-processing pipeline needs as
  ordinary named arguments, validates each one (types, ranges, file
  existence, annotation-package availability, samplesheet structure,
  IDAT presence, cross-reactive-probes CSV structure), and then calls
  \code{\link{processMethArrayExp}} to produce a bigmelon-compatible,
  LZ4_RA-compressed GDS file.

  The shipped example (\code{system.file("scripts", "pipeline.R", package = "fastMethyl")}) demonstrates the intended usage: declare
  each setting as an R variable at the top of the file, then call
  \code{runPreprocess(arg = arg, ...)} once. The function returns a
  result object the caller's analysis code plugs into.
}
\usage{
runPreprocess(dataDirectory,
                 nonSpecificProbesPath,
                 targetPattern,
                 datasetClass,
                 annotationPackage,
                 sampleDetPThreshold = 0.01,
                 probeDetPThreshold  = 0.01,
                 forceRebuild        = FALSE,
                 readerWorkers       = 1L,
                 verbose             = 2L,
                 BPPARAM             = NULL,
                 .check_annotation   = function(pkg)
                                       requireNamespace(pkg, quietly = TRUE))
}
\arguments{
  \item{dataDirectory}{Path to the directory containing IDAT files and
    the samplesheet CSV.}
  \item{nonSpecificProbesPath}{Path to a CSV listing cross-reactive
    probes. Must contain a \code{TargetID} column.}
  \item{targetPattern}{The samplesheet under \code{dataDirectory} is
    expected at \code{samplesheet_<targetPattern>.csv}.}
  \item{datasetClass}{Prefix for every output file the pipeline writes
    (GDS, txt). The GDS is written to
    \code{paste0(datasetClass, ".gds")} in the current working dir.}
  \item{annotationPackage}{Name of an installed Illumina annotation
    package, e.g.
    \code{"IlluminaHumanMethylation450kanno.ilmn12.hg19"}.}
  \item{sampleDetPThreshold}{A sample is dropped if its mean
    detection-p \eqn{\geq}{>=} this threshold. Default \code{0.01}.}
  \item{probeDetPThreshold}{A probe is dropped if any surviving
    sample's detection-p \eqn{\geq}{>=} this threshold. Default
    \code{0.01}.}
  \item{forceRebuild}{If \code{FALSE} (default) and the GDS already
    exists, the read pass is skipped and the existing file is reused.
    Set to \code{TRUE} to overwrite.}
  \item{readerWorkers}{IDAT-reader worker count for the default
    \code{MulticoreParam} backend. Default \code{1L} (serial). Ignored
    when \code{BPPARAM} is supplied.}
  \item{verbose}{Verbosity for the underlying reader; passed straight
    through to \code{\link{processMethArrayExp}}. \code{0L} silent,
    \code{1L} phase summaries + per-sample messages, \code{2L} adds
    master-side progress and footprint lines.}
  \item{BPPARAM}{A \pkg{BiocParallel} \code{BiocParallelParam} object
    defining the parallel backend for IDAT reading. \code{NULL}
    (default) derives a \code{MulticoreParam} from \code{readerWorkers},
    which suits most users. Pass an explicit backend
    (\code{\link[BiocParallel]{SnowParam}} on Windows or for a socket
    cluster, \code{BatchtoolsParam} under an HPC scheduler) when
    \code{MulticoreParam} is unavailable or unsuitable; set the worker
    count on that object rather than via \code{readerWorkers}.}
  \item{.check_annotation}{Injectable function used to check whether
    \code{annotationPackage} is installed. Default is
    \code{requireNamespace}. Exists for test stubbing; production
    callers should not pass this.}
}
\value{
  Invisibly returns a list with components:
  \describe{
    \item{\code{gds_path}}{Path to the QC'd GDS file.}
    \item{\code{targets}}{The samplesheet re-indexed by the post-QC
      sample-keep mask (matches the GDS columns one-to-one).}
    \item{\code{keepSamples}}{Logical vector of sample-detection-p
      survivors, \code{TRUE} for samples that survived. All
      \code{TRUE} when the existing GDS was reused.}
    \item{\code{keepProbes}}{Logical vector of probe-detection-p
      survivors. \code{NULL} when the existing GDS was reused.}
  }
}
\author{
  Sebastian Hoss, with assistance from Claude (Anthropic AI).
}
\seealso{
  \code{\link{processMethArrayExp}} for the lower-level reader
  \code{runPreprocess} wraps;
  \code{\link{processMethArray}} for the per-sample fused reader at the
  bottom of the stack.
}
