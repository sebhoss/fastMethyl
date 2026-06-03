<!--
SPDX-FileCopyrightText: The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# fastMethyl

A performance layer over [minfi](https://bioconductor.org/packages/minfi) and
[bigmelon](https://bioconductor.org/packages/bigmelon) for large Illumina
methylation cohorts. It does **not** reimplement them — it imports them and
adds a single, fast, memory-bounded entry point: **`analyze()`**.

> **Notice.** This package — the streaming preprocessing engine, the `analyze()`
> driver, and this documentation — was prepared with the help of
> [Claude](https://claude.ai/) (Anthropic's AI assistant) acting under the
> maintainer's direction. Outputs are equivalence-checked against
> `minfi`/`bigmelon` (the QC'd GDS is value-identical to
> `bigmelon::es2gds(minfi::preprocessRaw(minfi::read.metharray(...)))`), and the
> full pipeline has been run end-to-end on a real 582-sample EPIC cohort, but it
> has not been independently peer-reviewed. Please report anything you observe
> at <https://github.com/sebhoss/fastMethyl/issues>.

## What it does

`analyze()` is the one public function. A single call:

1. **validates** every input, then builds a QC'd, LZ4_RA-compressed,
   bigmelon-compatible **GDS** from raw IDATs in a fused, **column-streaming**
   pass (read + detection p-value + raw preprocessing + sample/probe QC). Full
   cohort-sized matrices are never assembled in RAM, so **peak memory stays
   bounded regardless of cohort size** — a 1000-sample EPIC run peaks at roughly
   the same memory as a 100-sample one;
2. optionally runs **your analysis function** on the open GDS handle; and
3. **closes the GDS for you** — on normal return *and* on error.

`analyze()` is deliberately unopinionated about the analysis: normalisation
(`dasen`) and outlier detection (`outlyx`) live in **your** function, so you
control their order — the statistically correct one being **outliers first,
then normalise**. Omit the function to just get the QC'd GDS.

## Benchmark

The streaming preprocessing inside `analyze()` is where the speed comes from.
All timings on a 16-core Linux box with `MulticoreParam(workers = 4)`.
*Upstream* is the released Bioconductor minfi pipeline; *fastMethyl* is the same
work done in `analyze()`'s fused, column-streaming pass. Outputs are
equivalence-verified (the GDS is value-identical) before timings are reported.

Full IDAT → fully-QC'd GDS: `read.metharray()` + `detectionP()` +
`preprocessRaw()` + sample/probe QC + `bigmelon::es2gds()` (upstream) vs the
single streaming pass (fastMethyl).

| Sample count | Upstream pipeline | fastMethyl | Speedup   |
|--------------|-------------------|------------|-----------|
| 30           | 87.9 s            | 15.1 s     | **5.8×**  |
| 100          | 267.5 s           | 24.8 s     | **10.8×** |
| 200          | 537.3 s           | 42.1 s     | **12.8×** |

The speedup *grows* with cohort size because the upstream pipeline's dominant
phase (`detectionP()`) scales linearly in the master process while fastMethyl
keeps it inside the parallel per-sample loop. Extrapolated to 1000 samples:
~45 min upstream → ~4 min with fastMethyl — and, unlike the upstream pipeline,
peak memory stays bounded regardless of cohort size (a real 582-sample EPIC
cohort ran the full pipeline in ~14.6 GB; see the Notice above).

## Installation

### Requirements

- **R 4.4 or later** (<https://cran.r-project.org/>).
- **An internet connection** during install — R downloads the package and its
  ~30 Bioconductor dependencies.
- **`bigmelon`** for the GDS output and the downstream analysis (`dasen`,
  `outlyx`, `prcomp`, `estimateCellCounts.gds`). The install command below
  installs it for you.

### Install command

Paste the whole block into the R console:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sebhoss/fastMethyl")
BiocManager::install("bigmelon")
```

The first install takes **5–20 minutes** because of the Bioconductor
dependencies. R may interrupt with two questions; the answer to both is `n`:

- `Update all/some/none? [a/s/n]:` — type `n`.
- `Do you want to install from sources the packages which need compilation? [y/n]:` — type `n` (the pre-built binaries are fine, so you don't need Rtools/Xcode).

You'll also need the annotation package for your array, e.g.
`BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")` (EPIC)
or `...450kanno.ilmn12.hg19` (450k).

### Verify

```r
library(fastMethyl)
#> fastMethyl <version>
#>   analyze(): one call -- preprocess ... -> your analysis function ...
```

If something fails to compile on Windows/macOS, install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows) or the Xcode
command-line tools (`xcode-select --install`, macOS) and re-run; for anything
else, open an issue at <https://github.com/sebhoss/fastMethyl/issues>.

## Usage

### Run a full analysis (supply a closure)

`analyze()` builds the GDS, runs your function on the open handle, and closes
the GDS even if your code errors. Your function receives the open handle `gds`
and the preprocessing result `res` (`gds_path`, the post-QC `targets`
samplesheet, and the keep masks). Do the normalisation and outlier detection
**inside** your function, in the right order — **`outlyx` on the raw betas
first, then `dasen`** (dasen can mask the technical artefacts `outlyx` keys on,
and removing flagged samples first keeps them out of dasen's quantile
reference):

```r
library(fastMethyl)
library(bigmelon)   # dasen / outlyx / prcomp on the open GDS

res <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",  # CSV with a TargetID column
  targetPattern         = "batch1",   # expects samplesheet_batch1.csv in dataDirectory
  datasetClass          = "my_study", # output GDS is my_study.gds
  annotationPackage     = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  FUN = function(gds, res) {
    outliers <- outlyx(gds, plot = FALSE)             # 1. detect on RAW betas
    dasen(gds, node = "normbetas")                    # 2. then normalise
    pca <- prcomp(gds, node.name = "normbetas", method = "quick")  # 3. analyse
    list(outliers = outliers, pca = pca)
  })
# `res` is whatever FUN returned; the GDS is already closed.
```

Pass `verbose = 2L` to log every step with its on-disk size and a memory
snapshot — so a slow or stalling load (e.g. materialising the EPIC annotation)
announces itself before it blocks instead of looking like a hang.

### Just build a GDS (no closure)

Omit `FUN` to stop after preprocessing. `analyze()` returns the result list
(`gds_path`, `targets`, `keepSamples`, `keepProbes`) and closes the GDS, so you
can open it later and normalise/analyse however you like:

```r
res <- analyze(
  dataDirectory         = "/path/to/idats",
  nonSpecificProbesPath = "/path/to/cross-reactive.csv",
  targetPattern         = "batch1",
  datasetClass          = "my_study",
  annotationPackage     = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
res$gds_path             # the QC'd GDS; open it with bigmelon/gdsfmt for any analysis
```

### Worked pipeline template

For a full end-to-end pipeline — outlier detection, cache-aware normalisation,
PCA, cell composition, optional outlier removal — copy the shipped template and
edit its CONFIG block:

```r
file.copy(system.file("scripts", "pipeline.R", package = "fastMethyl"), "pipeline.R")
# edit the CONFIG block, then:  Rscript pipeline.R
```

## Provenance & license

Derived from minfi (Hansen et al.) and interoperates with the
bigmelon / wateRmelon GDS ecosystem. Distributed under Artistic-2.0, matching
minfi.
