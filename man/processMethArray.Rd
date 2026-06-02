\name{processMethArray}
\alias{processMethArray}
\title{
  Read IDAT files and stream preprocessed methylation data to a GDS in one
  fused pass.
}
\description{
  Fuses \code{\link{readMethArray}}, \code{\link[minfi]{detectionP}},
  \code{\link[minfi]{preprocessRaw}} and \code{\link[bigmelon]{es2gds}} into a
  single streaming pass. The full \code{RGChannelSet}, the intermediate
  \code{MethylSet}, and the cohort-sized output matrices are never assembled in
  master memory: each worker reads one IDAT pair and computes its detection
  p-values and methylated / unmethylated intensity columns, and the master
  appends the surviving columns straight into a bigmelon-compatible GDS file.
  Peak memory is one column block, not the whole cohort.
}
\usage{
processMethArray(basenames, gds_path, pData = NULL,
                  drop_probes = NULL,
                  sample_detP_threshold = NULL,
                  probe_detP_threshold = NULL,
                  compress = "LZ4_RA",
                  verbose = 0L,
                  BPPARAM = SerialParam())
}
\arguments{
  \item{basenames}{The basenames of the IDAT files (no \code{_Grn.idat} /
    \code{_Red.idat} suffix).}
  \item{gds_path}{Output path for the GDS file (required). A partial file is
    removed if the run errors, so it is never mistaken for a valid cache.}
  \item{pData}{Optional \code{data.frame} / \code{DataFrame} of sample
    metadata (one row per basename, in the same order). Written into the
    GDS \code{pData} node and re-indexed by \code{keepSamples} if sample
    QC drops anything.}
  \item{drop_probes}{Character vector of CpG names to skip at compute time
    (e.g. sex-chromosome and cross-reactive probes). Workers do not waste
    cycles on dropped probes; the GDS contains only kept probes.}
  \item{sample_detP_threshold}{Optional numeric in \code{(0, 1)}. Drops
    any sample whose \code{colMeans(detP) >= threshold} as it streams past
    (its columns are simply not written). The result list's
    \code{keepSamples} component reports which samples survived.}
  \item{probe_detP_threshold}{Optional numeric in \code{(0, 1)}. Drops any
    probe where at least one surviving sample has \code{detP >= threshold}.
    Because this is known only after the full pass, the matrix nodes are
    compacted to their final probe set by a chunked row-subset rewrite once
    the stream completes. \code{keepProbes} reports the mask.}
  \item{compress}{Compression mode for the GDS bulk-numeric nodes
    (\code{betas}, \code{methylated}, \code{unmethylated}, \code{pvals}).
    Default \code{"LZ4_RA"} (LZ4 with random access) trades minimal CPU
    for ~2-4x disk savings on typical methylation matrices, and
    bigmelon's downstream tooling reads compressed nodes transparently.
    Pass \code{""} for no compression, or any other valid \code{gdsfmt}
    code (\code{"ZIP_RA"}, \code{"LZMA_RA"}, ...) as needed.}
  \item{verbose}{Verbosity level. Integer \code{0L}/\code{1L}/\code{2L} or
    logical (\code{FALSE} -> \code{0L}, \code{TRUE} -> \code{1L}).
    \code{0L} silent; \code{1L} prints phase summaries and per-sample
    worker messages; \code{2L} adds master-side per-block progress lines
    so long parallel runs always show activity.}
  \item{BPPARAM}{A \code{\link[BiocParallel]{BiocParallelParam}} instance.
    Defaults to \code{SerialParam()}; pass e.g.
    \code{\link[BiocParallel]{MulticoreParam}(4)} to parallelise.}
}
\details{
  \strong{Mixed probe-count cohorts.} The function always does a quick
  parallel header-only pre-pass of every IDAT to learn its probe count.
  When the cohort contains more than one distinct probe count (e.g.
  early-access and final-release EPIC IDATs side-by-side) the output is
  restricted to addresses present in every version's IDAT, and each worker
  dispatches to the right per-version index table. Uniform cohorts incur
  no measurable overhead from the pre-pass.

  \strong{Streaming GDS write.} Samples are processed in column blocks sized
  to bound a block's working set. Workers compute each sample's columns in
  parallel and the master appends the surviving columns to the extendable,
  compressed \code{betas} / \code{methylated} / \code{unmethylated} /
  \code{pvals} nodes, so a full cohort-sized matrix is never held in RAM. The
  resulting node layout matches
  \code{bigmelon::es2gds(preprocessRaw(readMethArray(...)))}, so downstream
  bigmelon functions (\code{outlyx}, \code{dasen}, \code{prcomp.gds},
  \code{estimateCellCounts.gds}) consume it unchanged. Nodes: \code{betas},
  \code{methylated}, \code{unmethylated}, \code{pvals}, \code{fData},
  \code{pData}, \code{history}, \code{paths}. With all three QC arguments left
  at \code{NULL} the GDS values are identical to
  \code{es2gds(preprocessRaw(readMethArray(...)))}.
}
\value{
  Invisibly returns a list with \code{gds_path}, \code{keepSamples},
  \code{keepProbes}.
}
\author{
  Sebastian Hoss, with assistance from Claude (Anthropic AI).
}
\seealso{
  \code{\link{readMethArray}} for the raw IDAT reader (returns an
  \code{RGChannelSet} without preprocessing or QC),
  \code{\link[bigmelon]{es2gds}} for the upstream GDS builder.
}
