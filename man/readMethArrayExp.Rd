\name{readMethArrayExp}
\alias{readMethArrayExp}
\title{
  Reads an entire metharray experiment using a sample sheet
}
\description{
  Reads an entire methylation array experiment using a sample sheet or (optionally)
  a target like data.frame.
}
\usage{
readMethArrayExp(base = NULL, targets = NULL, extended = FALSE,
    recursive = FALSE, verbose = 0L, BPPARAM = SerialParam())
}
\arguments{
  \item{base}{The base directory.}
  \item{targets}{A targets \code{data.frame}, see details}
  \item{extended}{Should the output of the function be a
    \code{"RGChannelSetExtended"} (default is \code{"RGChannelSet"}).}
  \item{recursive}{Should the search be recursive (see details)}
  \item{verbose}{Verbosity level (\code{0L} silent, \code{1L} phase
    summaries + per-sample worker messages, \code{2L} adds master-side
    periodic progress). \code{TRUE}/\code{FALSE} also accepted for back
    compat; see \code{\link{readMethArray}}.}
  \item{BPPARAM}{A \code{\link[BiocParallel]{BiocParallelParam}} instance
    passed to \code{\link{readMethArray}} to control how IDAT files are
    read.  See \code{\link{readMethArray}} for details.}
}
\details{
  If the \code{targets} argument is \code{NULL}, the function finds all
  two-color IDAT files in the directory given by \code{base}.  If
  \code{recursive} is \code{TRUE}, the function searches \code{base} and
  all subdirectories.  A two-color IDAT files are pair of files with
  names ending in \code{_Red.idat} or \code{_Grn.idat}.

  If the \code{targets} argument is not \code{NULL} it is assumed it has
  a columned named \code{Basename}, and this is assumed to be pointing
  to the base name of a two color IDAT file, ie. a name that can be made
  into a real IDAT file by appending either \code{_Red.idat} or
  \code{_Grn.idat}.

  The type of methylation array is guess by looking at the number of
  probes in the IDAT files.
}
\value{
  An object of class \code{"RGChannelSet"} or
  \code{"RGChannelSetExtended"}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{readMethArray}} for the workhorse function,
  \code{\link{readMethArraySheet}} for reading a sample sheet and
  \code{\link[minfi:RGChannelSet-class]{RGChannelSet}} for the output class.
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGset <- readMethArrayExp(file.path(baseDir, "5723646052"))

}
}
