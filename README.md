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
The benchmark compares the IDAT → fully-QC'd-data step both ecosystems perform:

- **Upstream** — `read.metharray.exp()` + `detectionP()` + `preprocessRaw()` +
  sample / probe / sex-chromosome / cross-reactive QC + writing the QC'd matrices
  to an LZ4_RA-compressed GDS.
- **fastMethyl** — the same end result in `analyze()`'s single fused,
  column-streaming pass.

Both pipelines produce the *same* artifact — a QC'd, LZ4_RA-compressed,
bigmelon-compatible GDS — and both pay the *same* compression cost. (bigmelon's
`es2gds` convenience writer is deliberately not used on the upstream side: it
stores the matrices **uncompressed**, which would let it skip the compression
work fastMethyl performs; the upstream column instead writes the same four nodes
with the same `LZ4_RA` codec.) What *is* excluded from both columns is the
downstream bigmelon **analysis** — `dasen`, `outlyx`, `prcomp`,
`estimateCellCounts.gds` — which is the same code regardless of which reader
produced the GDS. Outputs are equivalence-verified (the QC'd betas are
value-identical) before timings are reported.

**Data.** The benchmark runs on *real* Illumina 450k IDAT files — the six-sample
dataset shipped in the Bioconductor [`minfiData`](https://bioconductor.org/packages/minfiData)
package — not synthetic intensities. To reach a target cohort size *N*, those six
real Grn/Red IDAT pairs are symlinked in a cycle into an *N*-sample directory
with a matching samplesheet, so both pipelines do genuine on-disk IDAT reads and
the full read → detection-p → preprocess → QC → compressed-GDS-write path runs on
real array data at the chosen scale. (Repeating the six samples exercises the
per-sample compute and I/O faithfully; it does not bias the comparison, which is
about pipeline *structure* — fused streaming vs. read-everything-then-write — not
biological variety.)

**How the benchmark is executed.** To make the timings reproducible and to
protect the host, the whole run is launched inside a `systemd-run --user`
cgroup v2 scope with a **single fixed resource envelope used for every sample
size**: CPU via `CPUQuota` (**4 cores**), memory via `MemoryHigh` + `MemoryMax`
(sized once from a 200-sample reference, ~14 GiB / 17 GiB, so it never OOMs at
the largest cohort), and disk-write bandwidth via `IOWriteBandwidthMax`. The
memory/IO envelope is deliberately *not* re-sized per cohort — every run uses the
identical ~14 GiB / 17 GiB cap — so the across-size *and* across-core comparisons
reflect the workload, not changing limits (the core-scaling section below varies
only the CPU quota). Every process the benchmark spawns — the R
master *and* its forked reader workers — nests under that scope, so the caps
apply to the whole tree rather than to one process. fastMethyl's reader is
pinned to a fixed worker count (`MulticoreParam(workers = 4)`, matching the CPU
quota) rather than auto-detected, because under a quota the host's core count is
not the budget; upstream minfi reads serially, as it has no parallel read path.
Peak RAM is the **absolute peak of each pipeline run in a fresh process**,
measured from the cgroup so it spans the whole process tree (master + workers).
Running each pipeline fresh — rather than measuring fastMethyl's increment over a
still-warm minfi heap — is what makes the two columns honest and comparable:
each number is what that pipeline would peak at if you ran it on its own.

| Samples | Upstream time | fastMethyl time | Speedup  | Upstream peak RAM | fastMethyl peak RAM | Reduction |
|--------:|--------------:|----------------:|:--------:|------------------:|--------------------:|:---------:|
|      50 |     156.1 s   |        56.0 s   | **2.8×** | 4.4 GiB           | 4.0 GiB             | **8%**    |
|     100 |     282.8 s   |        99.3 s   | **2.8×** | 6.7 GiB           | 4.1 GiB             | **40%**   |
|     200 |     657.2 s   |       187.0 s   | **3.5×** | 12.7 GiB          | 4.2 GiB             | **67%**   |

**Time** is a clean, growing win — **2.8× → 3.5×** — because upstream's dominant
phases (`detectionP`, `preprocessRaw`) run serially in the master while
fastMethyl keeps them in the parallel per-sample loop, so the lead widens with
the cohort.

**Memory** is where the streaming design shows. Upstream holds the full
Red/Green, methylated/unmethylated and detection-p matrices in RAM at once, so its
peak grows roughly linearly (**4.4 → 12.7 GiB**). fastMethyl never assembles a
cohort-sized matrix, so its peak stays **essentially flat (~4 GiB)** — at these
sizes it is dominated by the one-off 450k annotation frame, not the cohort. The
consequence is honest about where the win is: at small N both pipelines are
annotation-bound and roughly even (8% less), but the gap opens to **67% less at
N=200** and keeps widening — exactly the property that keeps large EPIC cohorts
bounded where upstream OOMs (a real 582-sample EPIC cohort ran in ~14.6 GB; see
the Notice above).

### Scaling the reader across cores

fastMethyl's reader is parallel; upstream minfi reads serially and is therefore
**core-independent** (its column is unchanged below). Re-running only fastMethyl
at 8 and 12 worker cores against the same minfi baseline:

| Cores | N=50 time / peak | N=100 time / peak | N=200 time / peak |
|------:|------------------|-------------------|-------------------|
|   **4** | 56.0 s / 4.0 GiB | 99.3 s / 4.1 GiB | 187.0 s / 4.2 GiB |
|   **8** | 51.3 s / 6.1 GiB | 93.4 s / 6.1 GiB | 176.7 s / 4.9 GiB |
|  **12** | 56.8 s / 7.2 GiB | 94.4 s / 7.5 GiB | 177.7 s / 6.0 GiB |

The takeaway: **more cores buy almost no speed but cost real memory.** The host
has 8 physical cores, so wall-clock barely moves from 4 → 8 → 12 workers (at these
cohort sizes fastMethyl is not read-bound — the annotation load and GDS write
dominate), while peak RAM climbs steadily with the worker count. That climb is
**copy-on-write inflation**: each forked `MulticoreParam` worker's garbage
collector touches, and so privately copies, the shared master heap, so peak grows
roughly as *master-heap × workers*. (Its exact size is GC-timing dependent, so
the per-cohort peaks above are noisy to ±~1 GiB and not strictly monotonic in N;
the robust signal is that peak rises with worker count.) At 8–12 workers this can
even push fastMethyl's small-cohort peak above upstream's, though fastMethyl still
wins decisively at N=200 where upstream's matrices dominate. **Recommendation:
size workers to roughly the physical-core count (≈4–8 here); beyond that you pay
memory for no speed.**

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

BiocManager::install("sebhoss/fastMethyl", update = FALSE)
BiocManager::install("bigmelon")
```

`update = FALSE` is important: it tells the installer to leave your existing
dependencies in place and only install fastMethyl. Without it the installer may
try to re-install `minfi` (a dependency); if your `minfi` lives in a system
library you cannot write to, or is loaded in another R session, that re-install
fails with a permission error (see Troubleshooting below).

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

### Troubleshooting: a permission error mentioning `minfi`

```
Warning: cannot remove prior installation of package 'minfi'
Error: ... Permission denied
```

This means the installer tried to re-install `minfi` and could not overwrite the
copy you already have. fastMethyl depends only on `minfi`'s public API and does
**not** need a special build of it, so the re-install is unnecessary — stop it
from happening:

1. **Always pass `update = FALSE`** (as in the install command above). That alone
   prevents the installer from touching your existing `minfi`.
2. **Install from a fresh R session** (no `library(minfi)` loaded), so the
   package files are not locked.

If you previously installed the **`sebhoss/minfi` development fork** (it reports a
version like `1.59.0.9000`, higher than the Bioconductor release), the installer
sees a `minfi` "newer than Bioconductor" and is especially prone to trying to
replace it — which is what triggers the error when that copy is read-only or
loaded. fastMethyl works with both the fork and the stock release, so the
simplest fix is to return to the official `minfi` once:

```r
remove.packages("minfi")
BiocManager::install("minfi")        # the official Bioconductor build
BiocManager::install("sebhoss/fastMethyl", update = FALSE)
```

## Usage

### Run a full analysis (supply a closure)

`analyze()` builds the GDS, runs your function on the open handle, and closes
the GDS even if your code errors.

**Your function is `function(gds, res)`** — it receives two arguments:

- **`gds`** — the **open GDS handle** (a gdsfmt `gds.class` object), opened
  read-write. This is the bigmelon-compatible GDS, so you operate on it with
  `bigmelon`/`gdsfmt`: `dasen()`, `outlyx()`, `prcomp()`,
  `estimateCellCounts.gds()`, `index.gdsn()` / `read.gdsn()`, etc. Its nodes are
  `betas`, `methylated`, `unmethylated`, `pvals` (the raw QC'd data) plus
  `fData` / `pData` / `history` / `paths`; anything you add (e.g. a `normbetas`
  node from `dasen`) persists in the file. **Do not return `gds`** — it is
  closed by the time `analyze()` returns; return materialised results instead
  (matrices, data.frames, model objects).
- **`res`** — the preprocessing result, a list with:
  - `res$gds_path` — path to the GDS file on disk;
  - `res$targets` — the post-QC samplesheet (a data.frame, one row per
    surviving sample, in GDS column order);
  - `res$keepSamples` — logical vector marking which input samples passed
    sample-level QC;
  - `res$keepProbes` — logical vector marking which probes passed probe-level
    QC (or `NULL` when no probe threshold was applied).

Do the normalisation and outlier detection **inside** your function, in the
right order — **`outlyx` on the raw betas first, then `dasen`** (dasen can mask
the technical artefacts `outlyx` keys on, and removing flagged samples first
keeps them out of dasen's quantile reference):

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
