\name{analyze}
\alias{analyze}
\title{
  Preprocess, open, normalise, and run an analysis on a GDS with guaranteed cleanup.
}
\description{
  A one-call entry point that runs \code{\link{runPreprocess}}, opens the
  resulting GDS, optionally normalises it (\code{dasen}), hands the open handle
  to a user-supplied function \code{FUN}, and \strong{guarantees the GDS is
  closed afterwards} -- on normal return and on error -- before returning
  \code{FUN}'s value.

  This is the "loan pattern": you supply a closure describing what to do with
  the data; \code{analyze} owns the lifecycle. The close guarantee comes from
  \code{\link[base]{on.exit}}, so a leftover open handle can never block the
  next run, even if your analysis errors.
}
\usage{
analyze(..., FUN, normalize = TRUE, normNode = "normbetas",
        .preprocess = runPreprocess,
        .open = function(p) gdsfmt::openfn.gds(p, readonly = FALSE),
        .close = gdsfmt::closefn.gds,
        .normalize = .dasenNormalize)
}
\arguments{
  \item{\dots}{Arguments forwarded to \code{\link{runPreprocess}}
    (\code{dataDirectory}, \code{nonSpecificProbesPath}, \code{targetPattern},
    \code{datasetClass}, \code{annotationPackage}, ...).}
  \item{FUN}{A function \code{function(gds, res)} run with the open GDS handle
    (\code{gds}) and the \code{runPreprocess} result list (\code{res}: see its
    Value). Its return value becomes \code{analyze}'s return value. It must
    return \emph{materialised} results (matrices, data.frames, model objects),
    not the open \code{gds} handle -- the handle is closed by the time
    \code{analyze} returns.}
  \item{normalize}{If \code{TRUE} (default), run the normaliser into
    \code{normNode} before calling \code{FUN}.}
  \item{normNode}{GDS node the normalised betas are written to. Default
    \code{"normbetas"}.}
  \item{.preprocess, .open, .close, .normalize}{Injection points for the heavy
    operations, defaulting to the production implementations
    (\code{runPreprocess}; \code{openfn.gds} read-write; \code{closefn.gds}; an
    internal \code{dasen} wrapper that requires \pkg{bigmelon}). Exposed so the
    control flow can be unit-tested without IDATs or Bioconductor packages;
    production callers should not pass these.}
}
\value{
  Whatever \code{FUN} returns.
}
\details{
  Equivalent in spirit to the shipped example pipeline, but reusable: the
  fixed preamble (preprocess, open, normalise, guaranteed close) is owned by
  \code{analyze}, and the study-specific analysis lives in \code{FUN}.

  \code{normalize = TRUE} requires the \pkg{bigmelon} package; without it,
  either install bigmelon or pass \code{normalize = FALSE} and normalise inside
  \code{FUN}.
}
\examples{
\dontrun{
pcs <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",
  targetPattern         = "batch1",
  datasetClass          = "my_study",
  annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  FUN = function(gds, res) {
    prcomp(gds, node.name = "normbetas", method = "quick")
  })
# `pcs` is the PCA result; the GDS is already closed.
}
}
\author{
  Sebastian Hoss.
}
\seealso{
  \code{\link{runPreprocess}} for the preprocessing it wraps.
}
