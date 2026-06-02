\name{readMethArray}
\alias{readMethArray}
\title{
  Parsing IDAT files from Illumina methylation arrays.
}
\description{
  Parsing IDAT files from Illumina methylation arrays.
}
\usage{
readMethArray(basenames, extended = FALSE, verbose = 0L,
               BPPARAM = SerialParam())
}
\arguments{
  \item{basenames}{The basenames or filenames of the IDAT files.  By
    basenames we mean the filename without the ending \code{_Grn.idat}
    or \code{_Red.idat} (such that each sample occur once).  By
    filenames we mean filenames including \code{_Grn.idat}
    or \code{_Red.idat} (but only one of the colors)}
  \item{extended}{Should a \code{RGChannelSet} or a
    \code{RGChannelSetExtended} be returned. }
  \item{verbose}{Verbosity level. Either an integer \code{0L}, \code{1L},
    or \code{2L}, or a logical (\code{FALSE} maps to \code{0L},
    \code{TRUE} to \code{1L}, for back-compat). \code{0L} is silent;
    \code{1L} prints phase summaries and per-sample worker messages;
    \code{2L} additionally prints master-side progress updates every chunk
    of samples (so long-running parallel runs always show activity), plus
    extra context (array type, output matrix size, memory footprint).}
  \item{BPPARAM}{A \code{\link[BiocParallel]{BiocParallelParam}} instance
    controlling how IDAT files are read.  Defaults to
    \code{\link[BiocParallel]{SerialParam}}; pass e.g.
    \code{\link[BiocParallel]{MulticoreParam}(4)} to read files in parallel.}
}
\details{
  The type of methylation array is inferred from the number of probes in
  the IDAT files.

  IDAT files from the same array type but with different probe counts
  (e.g. early-access vs final-release EPIC) are supported automatically:
  the output is restricted to the intersection of probe addresses found
  across every distinct version. The intersect is computed per distinct
  probe count rather than per sample, so a cohort of N samples with V
  distinct versions costs O(V*n) instead of O(N*n) master-side work.
  Mixing different array types entirely (e.g. 450k + EPIC) is still a
  hard error -- use \code{combineArrays()} for that.
}
\value{
  An object of class \code{RGChannelSet} or
  \code{RGChannelSetExtended}. 
}
\author{
  Kasper Daniel Hansen\email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{readMethArrayExp}} for a convenience function for
  reading an experiment, \code{\link{readMethArraySheet}} for
  reading a sample sheet and \code{\link[minfi:RGChannelSet-class]{RGChannelSet}} for the
  output class. 
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGset1 <- readMethArray(file.path(baseDir, "5723646052", "5723646052_R02C02"))

}
}
