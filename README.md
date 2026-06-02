<!--
SPDX-FileCopyrightText: 2026 The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# fastMethyl

A performance layer over [minfi](https://bioconductor.org/packages/minfi) and
[bigmelon](https://bioconductor.org/packages/bigmelon) for large Illumina
methylation cohorts. It does **not** reimplement minfi/bigmelon — it imports
them and adds faster paths for the steps that dominate runtime and memory on
big cohorts.

## What it provides

**Preprocessing (the shipped path):**

- `readMethArray()` / `readMethArrayExp()` — minfi-compatible IDAT readers
  with a `BPPARAM=` argument for parallel reading and an integer-indexed worker
  assembly hot path (~11× faster than the serial upstream reader on large
  cohorts).
- `processMethArray()` — a fused, **column-streaming** pass: read + detection
  p-value + raw preprocessing + sample/probe QC, with each sample's columns
  appended straight into an LZ4_RA-compressed, bigmelon-compatible GDS. Full
  cohort-sized matrices are never assembled in RAM, so **peak memory stays flat
  in cohort size** — it does not grow with the number of samples (a 1000-sample
  EPIC run peaks at roughly the same memory as a 100-sample one). Sample QC
  drops columns as they stream; probe QC, known only after the full pass,
  triggers a chunked row-compaction so the GDS still ends at its final size.
- `processMethArrayExp()` — the experiment-level wrapper (samplesheet + base
  dir + annotation package).
- `runPreprocess()` — one validated entry point: config in, QC'd GDS out,
  with a build-key cache so re-runs that change thresholds/inputs rebuild
  instead of silently reusing a stale GDS.

**High-level driver:**

- `analyze()` — the recommended way to run a full analysis. One call runs
  `runPreprocess`, opens the GDS, (optionally) normalises it, runs *your*
  analysis function on the open handle, and **guarantees the GDS is closed
  afterwards** — on normal return and on error. You supply a closure; it owns
  the resource lifecycle. See Usage.

**Diagnostics.** Pass `verbose = 2L` to log every external-data step —
preprocessing, opening the GDS, normalisation, loading the annotation package,
and your analysis function — each annotated with the file's on-disk size and a
memory snapshot before and after. A slow or stalling load (e.g. materialising
the EPIC annotation) announces itself *before* it blocks, so an apparent hang
becomes a diagnosable event.

**Analysis phase.** The downstream analysis (normalisation, outlier
detection, PCA, cell composition) runs on the GDS via `bigmelon` directly —
see the example pipeline. Parallel/Gram-matrix reimplementations of `dasen` and
`outlyx` were prototyped and benchmarked but **dropped**: on a GDS, bigmelon's
paths are already well-optimised (outlier detection via gdsfmt's C-level
`apply.gdsn`; `dasen`'s cost is quantile normalisation + I/O, not the
parallelisable per-sample fits), so the reimplementations were net-negative.

## Usage

The high-level entry point is `analyze()`: give it the same inputs
`runPreprocess()` takes plus a function describing your analysis. It
preprocesses, opens and normalises the GDS, runs your function on the open
handle, and closes the GDS for you (even if your code errors).

```r
library(fastMethyl)
library(bigmelon)   # the normalisation + analysis below use bigmelon

result <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",
  targetPattern         = "batch1",
  datasetClass          = "my_study",
  annotationPackage     = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  FUN = function(gds, res) {
    # gds is open and normalised into the `normbetas` node; res carries
    # gds_path, the post-QC targets samplesheet, and the keep masks.
    prcomp(gds, node.name = "normbetas", method = "quick")
  })
# `result` is whatever FUN returned (here, a prcomp object).
# The GDS is already closed.
```

Pass `normalize = FALSE` to skip the built-in `dasen` and normalise inside
`FUN` yourself.

For a full worked pipeline — outlier detection, cache-aware normalisation, PCA,
cell composition, and optional outlier removal — copy the shipped template and
edit its CONFIG block:

```r
file.copy(system.file("scripts", "pipeline.R", package = "fastMethyl"),
          "pipeline.R")
# edit the CONFIG block, then:  Rscript pipeline.R
```

## Provenance & license

Derived from minfi (Hansen et al.) and interoperates with the bigmelon /
wateRmelon GDS ecosystem. Distributed under Artistic-2.0, matching minfi.
