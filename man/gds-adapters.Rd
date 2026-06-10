\name{gdsBetaMatrix}
\alias{gdsBetaMatrix}
\alias{gdsSummarizedExperiment}
\title{
  Project a QC'd GDS into a matrix or a SummarizedExperiment
}
\description{
  Read-side adapters that materialise a finished GDS into the container a
  non-GDS-aware downstream package expects. The GDS that \code{\link{analyze}}
  builds is the canonical, low-memory store; these adapters are the opposite
  direction -- they pull a representation into RAM on demand, sized to the
  caller, for tools that cannot read a GDS directly. GDS-native analysis
  (\code{bigmelon}'s \code{dasen}, \code{outlyx}, \code{prcomp},
  \code{estimateCellCounts.gds}) needs neither.

  \code{gdsBetaMatrix} returns a dense matrix; \code{gdsSummarizedExperiment}
  returns a \code{SummarizedExperiment}.
}
\usage{
gdsBetaMatrix(gds, node = "betas", transpose = TRUE)

gdsSummarizedExperiment(gds,
  assays = c("betas", "methylated", "unmethylated", "pvals"))
}
\arguments{
  \item{gds}{Either the open GDS handle that \code{\link{analyze}}'s analysis
    function receives (a gdsfmt \code{gds.class}) or a path to a GDS file on disk.
    A path is opened read-only and closed before the call returns.}
  \item{node}{The 2-D node to extract -- a \code{"betas"} / \code{"methylated"} /
    \code{"unmethylated"} / \code{"pvals"} matrix, or a node you added such as
    \code{"normbetas"}.}
  \item{transpose}{If \code{TRUE} (the default) the GDS's probes-x-samples node
    is returned as \code{samples x probes} -- the individuals-x-variables
    orientation \code{FactoMineR::PCA}, \code{stats::prcomp}, \code{irlba},
    \code{umap} and \code{Rtsne} expect. \code{FALSE} keeps the on-disk
    probes-x-samples layout.}
  \item{assays}{Names of the matrix nodes to load as assays. Any not present in
    the GDS are dropped; at least one must be present.}
}
\value{
  \code{gdsBetaMatrix}: a base \code{matrix}, with probe and sample names taken
  from the GDS when available.

  \code{gdsSummarizedExperiment}: a \code{SummarizedExperiment} whose assays keep
  the GDS probes-x-samples orientation, with \code{rowData} from the GDS
  \code{fData} and \code{colData} from \code{pData}.
}
\details{
  An adapter materialises the matrix in RAM, so it is opt-in and sized to the
  cohort -- reach for it when the destination package needs it, and stay on the
  GDS-native route otherwise.
}
\examples{
\dontrun{
analyze(..., normalization = "dasen", analysis = function(gds, res) {
  X <- gdsBetaMatrix(gds, node = "normbetas")   # samples x probes
  FactoMineR::PCA(X, graph = FALSE)
})

se <- gdsSummarizedExperiment("/path/to/my_study.gds")   # for limma / sva
}
}
\author{
  Sebastian Hoss.
}
