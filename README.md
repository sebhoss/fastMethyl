<!--
SPDX-FileCopyrightText: The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# fastMethyl

A performance layer over [minfi](https://bioconductor.org/packages/minfi) and
[bigmelon](https://bioconductor.org/packages/bigmelon) for large Illumina
methylation cohorts. It does **not** reimplement them — it imports them and
adds a single, fast, memory-bounded entry point: **`analyze()`**.

> **Based on the work of Marisol Herrera-Rivero.** The methylation-preprocessing
> approach this package implements and accelerates is from her work; **please
> cite the paper** if you use fastMethyl:
>
> Herrera-Rivero, M., Nauck, M., Berger, K., & Baune, B. T. (2025). Immune DNA
> methylation in depression: cross-sectional and longitudinal study. *BJPsych
> open*, 11(4), e129. <https://doi.org/10.1192/bjo.2025.10065>

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

1. **validates** every input, then builds a QC'd,
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
  to a GDS.
- **fastMethyl** — the same end result in `analyze()`'s single fused,
  column-streaming pass.

Both pipelines produce the *same* artifact — a QC'd, bigmelon-compatible GDS —
written with the *same* codec. The benchmark runs both **uncompressed**
(`compress = ""`, fastMethyl's default), for two reasons: it is the configuration
the package recommends for a working GDS, and it isolates the *structural* win.
GDS compression is a cost both sides pay almost equally — it is the same four
matrices either way — so turning it on mostly adds the same seconds to both
columns and **shrinks** the ratio rather than changing who wins. (With `LZ4_RA`
on both, the speed-up narrows from ≈9× to ≈4–5× while fastMethyl still leads at
every size; that compressed comparison lives in
[Speed vs. disk size](#speed-vs-disk-size-the-compress-knob) below.) bigmelon's
`es2gds` writer is not used on the upstream side — it would impose its own node
layout; the upstream column writes the same four nodes through gdsfmt so the two
artifacts match. What *is* excluded from both columns is the downstream bigmelon
**analysis** — `dasen`, `outlyx`, `prcomp`, `estimateCellCounts.gds` — which is
the same code regardless of which reader produced the GDS. Outputs are
equivalence-verified (the QC'd betas are value-identical) before timings are
reported.

**Data.** The benchmark runs on *real* Illumina 450k IDAT files — the six-sample
dataset shipped in the Bioconductor [`minfiData`](https://bioconductor.org/packages/minfiData)
package — not synthetic intensities. To reach a target cohort size *N*, those six
real Grn/Red IDAT pairs are symlinked in a cycle into an *N*-sample directory
with a matching samplesheet, so both pipelines do genuine on-disk IDAT reads and
the full read → detection-p → preprocess → QC → GDS-write path runs on
real array data at the chosen scale. (Repeating the six samples exercises the
per-sample compute and I/O faithfully; it does not bias the comparison, which is
about pipeline *structure* — fused streaming vs. read-everything-then-write — not
biological variety.)

**How the benchmark is executed.** To make the timings reproducible and to
protect the host, the whole run is launched inside a `systemd-run --user`
cgroup v2 scope with a **single fixed resource envelope used for every sample
size**: CPU via `CPUQuota` (**4 cores**), memory via `MemoryHigh` + `MemoryMax`
(sized once from a 200-sample reference, ≈14 GiB / 17 GiB, so it never OOMs at
the largest cohort), and disk-write bandwidth via `IOWriteBandwidthMax`. The
memory/IO envelope is deliberately *not* re-sized per cohort — every run uses the
identical ≈14 GiB / 17 GiB cap — so the across-size *and* across-core comparisons
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
|      50 |     129.0 s   |        23.1 s   | **5.6×** | 4.6 GiB           | 4.0 GiB             | **13%**   |
|     100 |     258.4 s   |        34.8 s   | **7.4×** | 7.1 GiB           | 4.0 GiB             | **44%**   |
|     200 |     612.8 s   |        69.2 s   | **8.9×** | 12.5 GiB          | 6.9 GiB             | **45%**   |

**Time** is a large, growing win — **5.6× → 8.9×** — because upstream's dominant
phases (`detectionP`, `preprocessRaw`) run serially in the master while
fastMethyl keeps them in the parallel per-sample loop, so the lead widens with
the cohort.

**Memory** is where the streaming design shows. Upstream holds the full
Red/Green, methylated/unmethylated and detection-p matrices in RAM at once, so its
peak grows roughly linearly (**4.6 → 12.5 GiB**). fastMethyl never assembles a
cohort-sized matrix, so its peak stays low and grows far slower (**4.0 → 6.9
GiB**) — at small N it is dominated by the one-off 450k annotation frame, not the
cohort, and the rise at N=200 is the parallel probe-QC compaction briefly holding
the four node buffers concurrently. The consequence is honest about where the win
is: at small N both pipelines are annotation-bound and closer (13% less), but the
gap opens to **45% less at N=200** and keeps widening — exactly the property that
keeps large EPIC cohorts bounded where upstream OOMs (a real 582-sample EPIC
cohort ran in ≈14.6 GB; see the Notice above).

### Past where minfi fits: N = 500 and 1000

The comparison table stops at N=200 because that is where a like-for-like
comparison stops being possible — **upstream minfi runs out of RAM first.** Its
peak grows at ≈0.05 GiB/sample, so on this 32 GB host it needs ≈29 GiB at N=500
and ≈57 GiB at N=1000, i.e. it OOMs past **N ≈ 600**. fastMethyl never
materialises a cohort-sized matrix, so it just keeps going:

| Samples | fastMethyl time | Peak RAM | Working set (avg) | GDS size |
|--------:|----------------:|---------:|------------------:|---------:|
|     500 |       126.8 s    | 12.6 GiB |      5.2 GiB       |  5.5 GB  |
|    1000 |       247.0 s    | 16.3 GiB |      8.4 GiB       | 10.8 GB  |

4 cores, 450k, uncompressed; throughput is a flat **≈0.25 s/sample**. The real
footprint is the **working set** (≈5–8 GiB) — the higher *peak* is mostly the
uncompressed GDS's reclaimable page cache (a 10.8 GB file at N=1000), not live
data, which is why the host never came under memory pressure. So fastMethyl's
ceiling here is **disk, not RAM**: a 450k 1000-sample uncompressed GDS is ≈11 GB
(EPIC ≈2×), so the scratch disk would hold several thousand samples — while minfi
cannot produce these at all on a 32 GB machine. (These two rows use a relaxed
I/O cap, so their per-sample time is not directly comparable to the throttled
50–200 rows above; the point is feasibility and the bounded footprint.)

### Scaling the reader across cores

fastMethyl's reader is parallel; upstream minfi reads serially and is therefore
**core-independent** (its column is unchanged below). Re-running only fastMethyl
at 8 and 12 worker cores against the same minfi baseline:

| Cores | N=50 time / peak | N=100 time / peak | N=200 time / peak |
|------:|------------------|-------------------|-------------------|
|   **4** | 23.1 s / 4.0 GiB | 34.8 s / 4.0 GiB | 69.2 s / 6.9 GiB |
|   **8** | 21.6 s / 6.2 GiB | 31.1 s / 6.1 GiB | 64.8 s / 7.1 GiB |
|  **12** | 24.1 s / 7.5 GiB | 38.1 s / 7.5 GiB | **50.4 s** / 8.5 GiB |

The takeaway: **with compression off, more cores genuinely speed up large
cohorts.** At N=200 the reads are the bottleneck again — once the serial
master-side compression is gone — so wall-clock drops with worker count (69 → 65
→ **50 s**, ≈27% from 4 → 12). At small N the per-worker overhead dominates
instead, so extra cores do nothing or slightly *regress* (12 workers is the
slowest at N=50). Peak RAM, meanwhile, always climbs with the worker count: that
is **copy-on-write inflation** — each forked `MulticoreParam` worker's garbage
collector touches, and so privately copies, the shared master heap, so peak grows
roughly as *master-heap × workers*. (Its exact size is GC-timing dependent, so
the per-cohort peaks above are noisy to ±1 GiB; the robust signal is that peak
rises with worker count.) **Recommendation: scale workers toward the core count
for large cohorts where the read parallelises, but expect diminishing returns and
rising memory; for small cohorts a handful of workers is plenty.** (Under the
non-default `compress = "LZ4_RA"`, the serial compression caps the gain and more
cores buy almost nothing — see [Speed vs. disk size](#speed-vs-disk-size-the-compress-knob).)

## Installation

### Requirements

- **R 4.4 or later** (<https://cran.r-project.org/>).
- **An internet connection** during install — R downloads the package and its
  ≈30 Bioconductor dependencies.
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
  readerWorkers         = parallel::detectCores(),  # parallel reads; see below
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

By default `analyze()` writes the GDS **uncompressed** (`compress = ""`) — the
fastest option and the right one for a working file. If you intend to archive or
share the GDS, pass `compress = "LZ4_RA"` to shrink it; see
[Speed vs. disk size](#speed-vs-disk-size-the-compress-knob).

### Choosing the worker count

`readerWorkers` sets how many IDATs are read in parallel (it defaults to `1L` —
serial). With compression off the reads are the bottleneck, so more workers
genuinely speed up large cohorts; the trade-off is that peak RAM rises with the
worker count. Pick by your situation:

| Goal | Set | Use when |
|---|---|---|
| Serial / reproducible / lowest memory | `readerWorkers = 1L` *(default)* | debugging, tiny cohorts, memory-tight machines, non-fork OS (Windows) |
| Fixed budget | `readerWorkers = 8L` | shared machine or a cgroup quota; bounds peak RAM |
| Use the whole machine | `readerWorkers = parallel::detectCores()` | a dedicated box with a large cohort — the biggest wins |

The speed-up concentrates at large cohorts (small ones are overhead-bound, where
extra workers can even regress); peak RAM climbs with workers; and inside a
container `detectCores()` over-reports the host's cores, so under a CPU quota pass
your real allotment instead. For a non-fork backend (Windows, or a socket
cluster) supply `BPPARAM = SnowParam(...)` rather than `readerWorkers`.

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
  annotationPackage     = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  readerWorkers         = parallel::detectCores())
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

### Speed vs. disk size: the `compress` knob

> **GDS compression is the single biggest cost in `analyze()`, so it is off by
> default.** `analyze()` writes the matrices **uncompressed** (`compress = ""`),
> the fastest option. Passing **`compress = "LZ4_RA"`** shrinks the GDS but is
> **≈2× slower**, because the compression CPU dominates both the streaming write
> *and* the probe-QC compaction. Turn it on only when you want a smaller file:

```r
analyze(...)                      # default compress = "": fastest, larger GDS
analyze(..., compress = "LZ4_RA") # ~2x slower; smaller GDS (~1.3x) for archiving
```

The trade is asymmetric: compression roughly **doubles the runtime** but shrinks
the file by only about a quarter (uncompressed is ≈1.3× larger), because
methylation intensity/beta matrices do not compress dramatically. Same 450k
cohorts, 4 cores, identical envelope:

| Samples | `""` time *(default)* | `LZ4_RA` time | Slow-down | `""` size | `LZ4_RA` size | Size saved |
|--------:|----------------------:|--------------:|:---------:|----------:|--------------:|:----------:|
|      50 |        23.1 s          |     44.7 s     | **1.9×** |  771 MiB  |    629 MiB     | **−18%**   |
|     100 |        34.8 s          |     76.4 s     | **2.2×** | 1314 MiB  |    972 MiB     | **−26%**   |
|     200 |        69.2 s          |    144.3 s     | **2.1×** | 2399 MiB  |   1765 MiB     | **−26%**   |

It changes only *how the data is stored*, never the values — the GDS is
**bigmelon-compatible either way** (uncompressed nodes read back identically and
still support the row-subset reads `dasen`/`outlyx` need), and it is not part of
the build-key cache. Recommendations:

| Situation | Use | Why |
|---|---|---|
| **You re-run often** (development, threshold tuning) | `compress = ""` *(default)* | the time saved each run dwarfs the disk cost; the file is transient |
| **Large cohort + ample/fast disk** | `compress = ""` *(default)* | the absolute time saving grows with cohort size (minutes on 1000-sample EPIC) |
| **Archiving, sharing, or disk-constrained** | `"LZ4_RA"` | a 1000-sample EPIC GDS can be tens of GB; ≈26% off is real |
| **One-off run you keep long-term** | `"LZ4_RA"` | you pay the compression once, save disk forever |
| **Small cohort** (≤ a few hundred) | either | the time difference is seconds; pick by whether you keep the file |

In short: **keep the default `compress = ""` whenever the GDS is a working
intermediate; pass `"LZ4_RA"` when the GDS is an artifact you store.**

## Provenance & license

Derived from minfi (Hansen et al.) and interoperates with the
bigmelon / wateRmelon GDS ecosystem. Distributed under Artistic-2.0, matching
minfi.
