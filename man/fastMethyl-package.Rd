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
  The single public entry point is \code{\link{analyze}}. One call validates
  the inputs and builds a QC'd, LZ4_RA-compressed, bigmelon-compatible GDS from
  raw IDATs -- a fused, column-streaming pass (read + detection p-value + raw
  preprocessing + sample/probe QC) whose peak memory stays bounded regardless
  of cohort size -- then optionally normalises it (\code{dasen}), optionally
  runs a user-supplied analysis function on the open GDS, and closes the GDS
  for you. Omitting that function makes \code{analyze} simply build (and
  optionally normalise) the GDS and return its path.

  The downstream analysis (outlier detection, PCA, cell composition) runs on
  the resulting GDS via \pkg{bigmelon} directly; see the example pipeline,
  retrievable with
  \code{system.file("scripts", "pipeline.R", package = "fastMethyl")}.
}
\author{
  Sebastian Hoss. Portions derived from \pkg{minfi} (Kasper D. Hansen and
  contributors).
}
\seealso{
  \code{\link{analyze}}.
}
\keyword{package}
