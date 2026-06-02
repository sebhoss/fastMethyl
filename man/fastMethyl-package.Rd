\name{fastMethyl-package}
\alias{fastMethyl-package}
\alias{fastMethyl}
\docType{package}
\title{
  Fast, parallel preprocessing of Illumina methylation arrays
}
\description{
  A performance layer over \pkg{minfi} and \pkg{bigmelon} for large Illumina
  methylation cohorts. It imports those packages (rather than reimplementing
  them) and adds faster paths for the steps that dominate runtime and memory
  on big cohorts.
}
\details{
  Two tiers of functionality:

  \strong{Preprocessing.}
  \itemize{
    \item \code{\link{readMethArray}} / \code{\link{readMethArrayExp}} --
      minfi-compatible IDAT readers with a \code{BPPARAM=} argument for
      parallel reading and an integer-indexed worker assembly hot path.
    \item \code{\link{processMethArray}} -- a fused per-sample pass (read +
      detection p-value + raw preprocessing + sample/probe QC) written once to
      an LZ4_RA-compressed, bigmelon-compatible GDS file.
    \item \code{\link{processMethArrayExp}} -- the experiment-level wrapper
      (samplesheet + base dir + annotation package).
    \item \code{\link{runPreprocess}} -- one validated entry point:
      config in, QC'd GDS out, with a build-key cache so re-runs that change
      thresholds or inputs rebuild instead of silently reusing a stale GDS.
  }

  The downstream analysis (normalisation, outlier detection, PCA, cell
  composition) runs on the resulting GDS via \pkg{bigmelon} directly; see the
  example pipeline.

  A worked example pipeline ships with the package; retrieve it with
  \code{system.file("scripts", "pipeline.R", package = "fastMethyl")}.
}
\author{
  Sebastian Hoss. Portions derived from \pkg{minfi} (Kasper D. Hansen and
  contributors).
}
\seealso{
  \code{\link{runPreprocess}}, \code{\link{processMethArray}}.
}
\keyword{package}
