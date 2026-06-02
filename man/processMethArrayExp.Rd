\name{processMethArrayExp}
\alias{processMethArrayExp}
\title{
  Experiment-level pre-processing of an Illumina methylation array
  cohort into a bigmelon-compatible GDS file.
}
\description{
  Wraps \code{\link{processMethArray}} with the orchestration that
  every fork-using pipeline repeats: resolve IDAT basenames from a
  samplesheet (and an optional base directory), load the annotation
  package, derive the sex-chromosome probe names from the annotation,
  concatenate them with any caller-supplied \code{drop_probes}, and call
  \code{processMethArray}. Returns the
  filtered samplesheet (with post-QC sample drops applied) alongside the
  GDS path and the keep-masks so the caller can plug straight into a
  downstream analysis without re-deriving them.

  Mirrors \code{\link{readMethArrayExp}} so users familiar with the
  upstream experiment-level reader find this naturally.
}
\usage{
processMethArrayExp(base = NULL, targets, gds_path,
                      annotation_package,
                      drop_probes = NULL,
                      drop_sex_chromosomes = TRUE,
                      sample_detP_threshold = NULL,
                      probe_detP_threshold = NULL,
                      compress = "LZ4_RA",
                      verbose = 0L,
                      BPPARAM = SerialParam())
}
\arguments{
  \item{base}{Optional directory path; if supplied, IDAT files are
    looked up at \code{file.path(base, targets$Basename)} (the
    \code{Basename} is preserved as-is, including any slide-subdirectory
    prefix). Otherwise \code{targets$Basename} is used as-is. Differs
    from \code{\link{readMethArrayExp}}, which calls
    \code{basename(targets$Basename)} and so requires a flat layout; here
    the typical Illumina \code{<slide>/<array>} layout works directly.}
  \item{targets}{A \code{data.frame} (or \code{DataFrame}) describing
    the experiment, one row per IDAT pair. Must contain a
    \code{Basename} column. The full frame is written into the GDS
    \code{pData} node; on sample-QC drops it is re-indexed so the
    returned \code{targets} matches the surviving GDS columns one to
    one.}
  \item{gds_path}{Output path for the GDS file. Required.}
  \item{annotation_package}{Name of an installed Illumina annotation
    package, e.g. \code{"IlluminaHumanMethylation450kanno.ilmn12.hg19"}
    or \code{"IlluminaHumanMethylationEPICanno.ilm10b4.hg19"}. Loaded
    via \code{requireNamespace} (no attach).}
  \item{drop_probes}{Optional character vector of CpG names to skip at
    compute time. Concatenated with the auto-derived sex-chromosome set
    when \code{drop_sex_chromosomes = TRUE}. Typically the
    cross-reactive / non-specific probe list (e.g. Chen et al. 2013 for
    450k, McCartney et al. 2016 for EPIC).}
  \item{drop_sex_chromosomes}{If \code{TRUE} (default), looks up the
    \code{chr} column in the annotation and drops every probe on
    \code{chrX} or \code{chrY}.}
  \item{sample_detP_threshold}{Optional numeric in \code{(0, 1)}, passed
    straight through to \code{\link{processMethArray}}.}
  \item{probe_detP_threshold}{Optional numeric in \code{(0, 1)}, passed
    straight through to \code{\link{processMethArray}}.}
  \item{compress}{GDS compression mode for the bulk-numeric nodes;
    forwarded to \code{\link{processMethArray}}. Default
    \code{"LZ4_RA"}.}
  \item{verbose}{Verbosity level; forwarded to
    \code{\link{processMethArray}}.}
  \item{BPPARAM}{A \code{BiocParallel} backend; forwarded to
    \code{\link{processMethArray}}.}
}
\value{
  Invisibly returns a list with components \code{gds_path},
  \code{keepSamples}, \code{keepProbes}, and \code{targets} (the input
  samplesheet re-indexed by \code{keepSamples}).
}
\author{
  Sebastian Hoss, with assistance from Claude (Anthropic AI).
}
\seealso{
  \code{\link{processMethArray}} for the low-level fused reader;
  \code{\link{readMethArrayExp}} for the raw-RGChannelSet
  experiment-level reader this function mirrors.
}
