\name{gdsDasen}
\alias{gdsDasen}
\title{
  Streaming between-array (dasen) normalisation of a GDS
}
\description{
  A streaming, two-pass between-array quantile normalisation that reproduces
  \code{wateRmelon::dasen} bit-for-bit while reading the GDS one column block at
  a time -- it never assembles the full samples-x-probes matrix and writes no
  scratch GDS. It is the default normaliser the \code{\link{analyze}}
  \code{normalization = "dasen"} hook uses, and depends only on \pkg{stats} +
  \pkg{gdsfmt} (it does not call \pkg{wateRmelon}).

  The quantile-normalisation reference is the mean of the sorted columns, which
  accumulates one sample at a time: pass 1 builds the four references
  (methylated / unmethylated x Type I / II) over the \code{keep} samples; pass 2
  re-reads each column, dye/background-corrects it, maps its ranks onto the
  reference, forms \eqn{beta = M / (M + U + fudge)}, and appends.
}
\usage{
gdsDasen(gds, node = "normbetas", keep = NULL, fudge = 100, block = 64L,
         BPPARAM = NULL)
}
\arguments{
  \item{gds}{The open read-write GDS handle that \code{\link{analyze}}'s analysis
    function receives, or a path to a GDS file (opened read-write, closed on
    return). It must contain \code{methylated}, \code{unmethylated}, and an
    \code{fData$Type} column of \code{"I"} / \code{"II"} probe designs.}
  \item{node}{Name of the beta node to write (deleted first if it exists).}
  \item{keep}{Optional logical mask over the samples selecting which build the
    quantile reference -- pass the outlier survivors so the reference excludes
    outliers. The default (\code{NULL}) uses every sample. The normalised
    \code{node} is written for \emph{every} sample regardless, so the GDS stays
    rectangular.}
  \item{fudge}{Offset added to the denominator when forming betas (the
    \code{dasen} default, 100).}
  \item{block}{Column-block size that bounds peak memory: larger reads fewer
    times, smaller holds less at once.}
  \item{BPPARAM}{Optional, opt-in parallel backend. The default (\code{NULL} ->
    \code{SerialParam}) runs serial and lean. Pass a BiocParallel backend (e.g.
    \code{MulticoreParam(4)}) to fan the per-column compute (the ~90% that is
    dye-offset densities + the rank/map) across workers -- value-identical to
    serial up to floating-point reassociation, roughly twice as fast, at higher
    but block-bound (not cohort-bound) peak memory. It is opt-in precisely so the
    default keeps the streaming memory guarantee.}
}
\value{
  Invisibly, \code{node}. Called for its side effect: the normalised node is
  written into \code{gds}.
}
\details{
  This is the per-sample \code{dfsfit} case (\code{roco = NULL}); sentrix-position
  smoothing is not applied. Peak memory is bounded and roughly flat in cohort
  size, where an in-RAM \code{dasen} grows with it.
}
\examples{
\dontrun{
analyze(..., outlierRemoval = "outlyx",
        normalization = function(gds, res, keep) {
          gdsDasen(gds, node = "normbetas", keep = keep)
        })
}
}
\author{
  Sebastian Hoss.
}
