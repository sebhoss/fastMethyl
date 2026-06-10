<!--
SPDX-FileCopyrightText: The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# fastMethyl

A performance layer over [minfi](https://bioconductor.org/packages/minfi) and
[bigmelon](https://bioconductor.org/packages/bigmelon) for large Illumina
methylation cohorts. It does **not** reimplement them — it imports them and
adds a single, fast, memory-lean entry point: **`analyze()`**.

> **Based on the work of Marisol Herrera-Rivero.** The methylation-preprocessing
> approach this package implements and accelerates is from her work:
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

`analyze()` is the main public function. A single call:

1. **validates** every input, then builds a QC'd,
   bigmelon-compatible **GDS** from raw IDATs in a fused, **column-streaming**
   pass (read + detection p-value + raw preprocessing + sample/probe QC). Full
   cohort-sized matrices are never assembled in RAM, so **peak memory grows far
   more slowly than the cohort** — large cohorts that make upstream minfi run out
   of RAM stay feasible (see the benchmark below);
2. **detects outliers** and **normalises**, via pluggable hooks, in the
   statistically correct order; then
3. runs **your analysis function** on the result and **closes the GDS for you** —
   on normal return *and* on error.

`analyze()` *orchestrates* the analysis phases — **build → outlier removal →
normalisation → your function** — and **enforces the order**, so outliers are
excluded from the normalisation reference (the one ordering that's easy to get
wrong, and statistically important). The `outlierRemoval` and `normalization` steps
each take a **built-in method name** (`"outlyx"` / `"dasen"`), **your own
function**, or `"none"` — all default to `"none"`, so a bare call stays
unopinionated. Omit everything to just get the QC'd GDS.

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

### The whole pipeline, end to end

The table above isolates the **preprocessing** step (IDAT → QC'd GDS), which is
where fastMethyl actually differs. But you run a *whole* pipeline — outlier
detection, normalisation, PCA, cell composition — so here is the complete
workflow both ways: upstream (`read.metharray.exp` + `preprocessRaw` + QC +
`es2gds`, then `outlyx` / `dasen` / `prcomp` / `estimateCellCounts.gds`) vs
fastMethyl (`analyze()`'s streaming build + `gdsDasen`, same analysis). 450k, the
same fixed envelope.

| Samples | Upstream total | fastMethyl total | Speedup | of which build (up → fm) | Upstream peak | fastMethyl peak |
|--------:|---------------:|-----------------:|:-------:|:------------------------:|--------------:|----------------:|
|      50 |     307.5 s    |      156.7 s     | **2.0×** | 184.8 → 29.6 s (6.2×)    | 3.9 GiB       | 3.2 GiB         |
|     100 |     573.1 s    |      218.0 s     | **2.6×** | 396.1 → 43.5 s (9.1×)    | 7.5 GiB       | 3.0 GiB         |
|     200 |    1135.9 s    |      347.4 s     | **3.3×** | 818.7 → 98.7 s (8.3×)    | 12.5 GiB      | 4.0 GiB         |

The full-pipeline speedup (**2–3.3×**) is *smaller* than the preprocessing speedup
above — and that's the honest framing. The analysis half (`outlyx`, `prcomp`,
`estimateCellCounts.gds`) is the **same bigmelon code** whichever reader built the
GDS, so it adds the same time to both columns and dilutes the win; the only
analysis difference is fastMethyl's streaming `gdsDasen`, which edges bigmelon's
`dasen.gds` at N=200 (249 s vs 317 s). The lead lives in the **build — 6–9×**.

**Memory is again the structural story.** Upstream's peak climbs with the cohort
(**3.9 → 12.5 GiB** — it holds the full matrices during the read); fastMethyl
stays **bounded at ~3–4 GiB**. fastMethyl's *build* peak is near-nothing (110 MiB
at N=50), so its full-pipeline peak is set by the **shared** cell-composition
reference load, **not** by anything fastMethyl does — which is why it flattens
while upstream grows without bound. At a 1000-sample cohort that's the difference
between a run that fits and one that OOMs.

> The fastMethyl peak here (≈3–4 GiB) reads a little *lower* than the
> preprocessing-only table above (≈4–6.9 GiB) — same operation, separate runs.
> The difference is run-to-run variance in the build's copy-on-write inflation
> across forked reader workers (sensitive to GC timing; the isolated
> preprocessing harness, which forks against a pre-warmed ~2 GiB annotation
> parent, tends to *overstate* it). Both stay in the few-GiB, bounded regime —
> the contrast with upstream's growth is the point, not the second decimal.

### Past where minfi fits: N = 500 to 2000

The comparison table stops at N=200 because that is where a like-for-like
comparison stops being possible — **upstream minfi runs out of RAM first.** Its
peak grows at ≈0.05 GiB/sample, so on this 32 GB host it needs ≈29 GiB at N=500,
≈57 GiB at N=1000 and ≈110 GiB at N=2000 — it OOMs past **N ≈ 600**. fastMethyl
never materialises a cohort-sized matrix, so it just keeps going:

| Samples | fastMethyl time | Peak RAM | Working set (avg) | GDS size |
|--------:|----------------:|---------:|------------------:|---------:|
|     500 |       126.8 s    | 12.6 GiB |      5.2 GiB       |  5.5 GB  |
|    1000 |       247.0 s    | 16.3 GiB |      8.4 GiB       | 10.8 GB  |
|    2000 |       509.6 s    | 14.9 GiB |     10.5 GiB       | 21.3 GB  |

4 cores, 450k, uncompressed. Throughput is a flat **≈0.25 s/sample** all the way
to N=2000 (≈8.5 min) and the GDS grows linearly. The real footprint is the
**working set** (≈5–11 GiB); the *peak* is mostly the uncompressed GDS's
reclaimable page cache, so it is writeback-/reclaim-noisy and not even monotonic
in N (N=2000 peaks below N=1000), and the host never came under memory pressure.
fastMethyl's ceiling here is therefore **disk, not RAM**: a 450k GDS is ≈11 GB per
1000 samples uncompressed (EPIC ≈2×), so the scratch disk would hold several
thousand samples — while minfi cannot produce any of these on a 32 GB machine.
(These rows use a relaxed I/O cap, so their per-sample time is not directly
comparable to the throttled 50–200 rows above; the point is feasibility and the
bounded footprint.)

At **8 cores** the same cohorts run **≈25–30% faster** — 94.9 / 186.5 / 367.3 s
for N=500 / 1000 / 2000 (≈0.185 s/sample, vs ≈0.25 s/sample at 4 cores) — and the
host stayed at 20–25 GiB free throughout, so 8 cores is well within this 32 GB
machine. This is the same read-bound scaling the core sweep below shows, now at
cohort sizes minfi cannot reach. (Peak RAM at this scale is dominated by
reclaimable page-cache/writeback timing and is too noisy to tabulate; the
throughput and the comfortable host headroom are the reliable signals.)

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

### Run a full analysis (the hooks + your closure)

`analyze()` builds the GDS, then runs **outlier removal → normalisation → your
function** on it, and closes the GDS even if your code errors. The two analysis
hooks share one shape — a built-in method name, your own function, or `"none"`:

```r
library(fastMethyl)
library(bigmelon)   # outlyx / prcomp / estimateCellCounts.gds on the open GDS

res <- analyze(
  dataDirectory       = "/path/to/idats",                  # folder the IDATs live under
  samplesheet         = "/path/to/samplesheet.csv",        # CSV with a Basename column
  crossReactiveProbes = "/path/to/cross-reactive.csv",     # CSV with a TargetID column
  gdsOutput           = "my_study.gds",                    # path to the output GDS
  annotationPackage   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  readerWorkers       = min(parallel::detectCores(), 16L), # parallel reads; cap so containers/huge boxes don't oversubscribe -- see below

  outlierRemoval = "outlyx",   # flag outliers on the RAW betas, exclude them from the reference
  normalization  = "dasen",    # streaming dasen -> a `normbetas` node (exact, memory-bounded)

  analysis = function(gds, res) {   # your analysis, on the normalised GDS
    pca <- prcomp(gds, node.name = "normbetas", method = "quick")
    list(pca = pca, keep = res$outlierKeep)
  })
# `res` is whatever `analysis` returned; the GDS is already closed.
```

`analyze()` runs the steps **in order** and feeds the outlier **keep mask** from
`outlierRemoval` to `normalization`, so `dasen`'s quantile reference is built from the
survivors (the mask is also exposed to `analysis` as `res$outlierKeep`). You can't get
outliers-before-normalisation wrong, because `analyze()` owns the order.

**The hooks** (`outlierRemoval`, `normalization`) each accept:

- a **method name** — `"outlyx"` / `"dasen"` (the fast, known-good defaults), or
  `"none"` to skip the step (the default for both);
- **your own function** — `outlierRemoval = function(gds, res)` must return a
  logical **keep mask** over the samples (`TRUE` = survivor); `normalization =
  function(gds, res, keep)` writes a normalised node, using `keep` to choose the
  reference samples. This is how you swap in a different outlier rule or
  normalisation (e.g. `funnorm`, a custom QC) while keeping the orchestration.

Built-in `normalization = "dasen"` is fastMethyl's own `gdsDasen()` — a streaming
quantile normalisation that reproduces `wateRmelon::dasen` **bit-for-bit** while
reading the GDS one block at a time, so peak memory stays **bounded and flat**
even on cohorts where the in-RAM `dasen` would not fit. It writes `normbetas` for
**every** sample (the GDS stays rectangular); a flagged-but-kept outlier just gets
a column normalised against an outlier-free reference.

#### Parallel normalisation (opt-in `BPPARAM`)

`dasen` is the heaviest analysis phase, and it's ~90% per-column compute (the
dye-offset densities and the rank-mapping), so it parallelises well. `gdsDasen()`
takes an opt-in **`BPPARAM`**: the default (serial) is lean — the master holds one
block at a time and forks nothing; pass a BiocParallel backend to fan the
per-column work across workers for roughly **2× speed**. The result is
**value-identical** either way (the reference is summed per worker, then combined,
so it matches the serial result to floating-point precision).

It's opt-in — not the default — for one reason: forking workers raises peak memory
(copy-on-write), so you trade RAM for speed. The cost is **block-bound, not
cohort-bound** (it tracks the fixed column-block, not N), so it stays within a sane
envelope at any cohort size, but it is real (~2–3× the serial peak). The lean
serial default keeps fastMethyl's memory guarantee intact; reach for the parallel
path when CPU, not RAM, is your constraint. (Measured at 450k, N=200, 4 workers:
normalisation **~180 s → ~95 s**, peak **~2.5 GiB → ~5.8 GiB**.)

The `"dasen"` string runs serial. To parallelise, pass `gdsDasen()` yourself from a
**custom normalisation hook**:

```r
library(BiocParallel)

analyze(
  ...,
  outlierRemoval = "outlyx",
  normalization  = function(gds, res, keep) {       # parallel dasen, ~2x faster
    gdsDasen(gds, node = "normbetas", keep = keep,
             BPPARAM = MulticoreParam(4))
  },
  analysis = function(gds, res) prcomp(gds, node.name = "normbetas", method = "quick"))
```

`gdsDasen()` is also callable standalone (on the open handle your `analysis`
receives, or on a GDS path) — handy if you build the GDS with `analyze(FUN = NULL)`
and normalise it separately.

**Your function is `function(gds, res)`** — it receives two arguments:

- **`gds`** — the **open GDS handle** (a gdsfmt `gds.class` object), opened
  read-write. By the time `analysis` runs, the hooks have already added the
  `normbetas` node (if `normalization` was set). It is the bigmelon-compatible GDS,
  so you operate on it with `bigmelon`/`gdsfmt`: `prcomp()`,
  `estimateCellCounts.gds()`, `index.gdsn()` / `read.gdsn()`, etc. Its nodes are
  `betas`, `methylated`, `unmethylated`, `pvals` (the raw QC'd data) plus
  `fData` / `pData` / `history` / `paths`, plus `normbetas`; anything you add
  persists in the file. **Do not return `gds`** — it is closed by the time
  `analyze()` returns; return materialised results instead (matrices,
  data.frames, model objects).
- **`res`** — the preprocessing result, a list with:
  - `res$gds_path` — path to the GDS file on disk;
  - `res$targets` — the post-QC samplesheet (a data.frame, one row per
    surviving sample, in GDS column order);
  - `res$keepSamples` — logical vector marking which input samples passed
    sample-level QC;
  - `res$keepProbes` — logical vector marking which probes passed probe-level
    QC (or `NULL` when no probe threshold was applied);
  - `res$outlierKeep` — the keep mask from the `outlierRemoval` hook (`NULL` when
    that hook was `"none"`).

**Prefer the hooks over doing outlier removal / normalisation inside `analysis`.**
Putting them in the hooks lets `analyze()` enforce the order and exclude outliers
from the normalisation reference for you. You *can* still do everything in `analysis`
with both hooks left `"none"` — that's the unopinionated mode — but then the
order is yours to get right (outliers on the raw betas first, then `dasen`,
because `dasen` can mask the artefacts `outlyx` keys on and outliers skew its
quantile reference). The shippable
[`pipeline.R`](inst/scripts/pipeline.R) template wires the hooks up end to end.

Pass `verbose = 2L` to log every step with its on-disk size and a memory
snapshot — so a slow or stalling load (e.g. materialising the EPIC annotation)
announces itself before it blocks instead of looking like a hang.

By default `analyze()` writes the GDS **uncompressed** (`compress = ""`) — the
fastest option and the right one for a working file. If you intend to archive or
share the GDS, pass `compress = "LZ4_RA"` to shrink it; see
[Speed vs. disk size](#speed-vs-disk-size-the-compress-knob).

### Feeding other packages — the adapters

GDS-native tools (`bigmelon`'s `dasen`, `outlyx`, `prcomp`,
`estimateCellCounts.gds`) work on the open `gds` handle directly — no conversion
needed, and that route never materialises a cohort-sized matrix. For a package
that **cannot** read a GDS, two adapters project the QC'd data into the container
it expects:

- **`gdsBetaMatrix(gds, node = "betas", transpose = TRUE)`** → a dense
  `samples × probes` matrix (the *individuals × variables* orientation
  `FactoMineR::PCA`, `stats::prcomp`, `irlba`, `umap` and `Rtsne` expect). Pass
  `transpose = FALSE` for `probes × samples`, or a different `node` (e.g.
  `"normbetas"`) to read a node you added.
- **`gdsSummarizedExperiment(gds)`** → a `SummarizedExperiment` whose assays are
  the matrix nodes, with `rowData` from `fData` and `colData` from `pData` — ready
  for `limma`, `sva`, and the rest of the `SummarizedExperiment` ecosystem.

Both take the open handle your `analysis` receives *or* a path on disk (a path is
opened read-only and closed for you), so they work inside `analysis` or standalone on
`analyze(analysis = NULL)$gds_path`:

```r
analyze(..., normalization = "dasen",                # writes normbetas before analysis runs
        analysis = function(gds, res) {
  X   <- gdsBetaMatrix(gds, node = "normbetas")  # n × p matrix, materialised here
  pca <- FactoMineR::PCA(X, graph = FALSE)       # graph = FALSE: no auto-plot
  list(pca = pca)
})
```

An adapter materialises the matrix in RAM, so it is **opt-in and sized to your
cohort** — reach for it when the destination package needs it, and stay on the
GDS-native route otherwise. (For PCA specifically, the GDS-native
`prcomp(node.name = "normbetas", method = "quick")` is both the fastest and the
leanest option for this `p ≫ n` data — a materialise-then-decompose route such as
`FactoMineR::PCA` or `irlba` reads the whole matrix into RAM first and does not
win the time back, so `pipeline.R` stays on `prcomp`.)

### Parameters

`analysis`, `outlierRemoval` and `normalization` are `analyze()`'s **own** arguments (the
orchestration); **every other argument is forwarded to the preprocessing step**,
so the full set is below. The ones whose value takes some thought link to a
dedicated guidance section.

**Required** — every run must supply these (they have no default):

| Parameter | What it is | Picking a value |
|---|---|---|
| `dataDirectory` | Directory the raw `*_Grn.idat` / `*_Red.idat` files live under. | Point it at the folder your IDATs live in; both `.idat` and `.idat.gz` are accepted. The samplesheet `Basename` entries resolve relative to this. |
| `samplesheet` | Path to the samplesheet CSV (it may live anywhere). | The sheet needs a `Basename` column whose entries resolve to the IDAT basenames under `dataDirectory`; the slide subdirectory in a `Basename` is preserved. |
| `gdsOutput` | Path to the output GDS file (written verbatim, e.g. `/data/my_study.gds`). | A full path including the `.gds` name; the [build-cache](#the-build-cache-and-rebuild) sidecar (`<gdsOutput>.buildkey.rds`) is written beside it. |
| `annotationPackage` | The Illumina annotation package for **your array type** (probe annotation, incl. the sex-chromosome map). | Must match the array **and** be installed. 450k: `IlluminaHumanMethylation450kanno.ilmn12.hg19`; EPIC v1: `IlluminaHumanMethylationEPICanno.ilm10b4.hg19`. A wrong/missing package fails validation before any IDAT is read. |
| `crossReactiveProbes` | CSV of cross-reactive / non-specific probes to drop; must contain a `TargetID` column (one probe ID per row). | Use a published cross-reactive list for your array (e.g. Chen et al. 2013 for 450k, Pidsley et al. 2016 for EPIC). The column **must** be named `TargetID`. |

**Optional** — sensible defaults; tune as needed:

| Parameter | Default | What it is | Picking a value |
|---|---|---|---|
| `analysis` | `NULL` | Your analysis closure `function(gds, res)`, run on the (outlier-aware, normalised) GDS. | Omit to stop after the hooks. See [Run a full analysis](#run-a-full-analysis-the-hooks--your-closure) / [Just build a GDS](#just-build-a-gds-no-closure). |
| `outlierRemoval` | `"none"` | Outlier step: a method name (`"outlyx"`), your own `function(gds, res)` returning a keep mask, or `"none"`. | `"outlyx"` is the fast default; runs on the raw betas and excludes flagged samples from the normalisation reference. |
| `normalization` | `"none"` | Normalisation step: a method name (`"dasen"`), your own `function(gds, res, keep)` that writes a node, or `"none"`. | `"dasen"` is the streaming, memory-bounded `gdsDasen()`; writes `normbetas` for every sample using the survivors as the reference. |
| `sampleDetPThreshold` | `0.01` | Sample QC: drop a sample whose **mean** detection p-value is at or above this. | See [Detection p-value QC thresholds](#detection-p-value-qc-thresholds). |
| `probeDetPThreshold` | `0.01` | Probe QC: drop a probe that fails detection in **any** surviving sample. | See [Detection p-value QC thresholds](#detection-p-value-qc-thresholds). |
| `readerWorkers` | `1L` | How many IDATs are read in parallel. | See [Choosing the worker count](#choosing-the-worker-count). |
| `BPPARAM` | `NULL` | An explicit BiocParallel backend, overriding `readerWorkers`. | Only for non-fork backends (Windows / clusters). See [Choosing the worker count](#choosing-the-worker-count). |
| `compress` | `""` | gdsfmt codec for the matrix nodes (`""` = uncompressed). | See [Speed vs. disk size](#speed-vs-disk-size-the-compress-knob). |
| `rebuild` | `"none"` | Which pipeline stage to recompute (it + everything downstream): `"none"`, `"normalize"`, `"build"`/`"all"`, or a logical. | See [The build cache and `rebuild`](#the-build-cache-and-rebuild). |
| `verbose` | `0L` | Logging level: `0L` silent, `1L` phase + per-sample lines, `2L` adds a size + memory snapshot around every external-data load. | Use `2L` when a run seems to hang — it shows *which* load is slow before it blocks. Legacy `TRUE`/`FALSE` is accepted. |

### Detection p-value QC thresholds

minfi assigns every (probe, sample) a **detection p-value** — the chance the
measured signal is indistinguishable from the array's background (negative-control)
probes. A *small* p is a confident measurement; a *large* p is a failed one.
fastMethyl uses these for two independent filters, both on the same 0–1 scale and
both defaulting to `0.01` (the common literature default). Each must be a single
number strictly between 0 and 1.

**`sampleDetPThreshold` — drops whole samples.** A sample is removed when the
**mean** of its detection p-values across all probes is at or above the threshold.
This catches *globally* failed arrays (poor bisulfite conversion, too little input
DNA, scan failure), not the odd bad probe.

- A healthy array's mean detection p is tiny (often ≈1e-3 or smaller), so the
  `0.01` default removes only clearly broken samples and rarely fires on good data.
- Lower it (e.g. `0.005`) to be stricter; raise it (`0.05` is sometimes used) to
  be more permissive. Because it is a *mean*, a sample with a handful of failed
  probes still passes — culling those is the probe filter's job.

**`probeDetPThreshold` — drops probes cohort-wide.** A probe is removed when it
fails detection (its detection p is at or above the threshold) in **at least one**
surviving sample. This is the strictest per-probe policy: every retained probe is
reliably measured in *every* sample.

- Because a single bad sample can knock out a probe for the whole cohort, the
  number of probes dropped grows with cohort size and heterogeneity — on large or
  noisy cohorts it can be substantial.
- If you want a softer rule (drop a probe only when it fails in some *fraction* of
  samples), do it yourself in `analysis`: set `probeDetPThreshold` close to `1` (e.g.
  `0.999999`) to effectively disable the built-in all-or-nothing filter, then apply
  your own per-probe mask to the GDS. The built-in rule has no fraction knob.

### The build cache and `rebuild`

To avoid re-reading IDATs you have already processed, `analyze()` caches its work.
Alongside the `gdsOutput` GDS it writes a small **build-key sidecar** recording
the inputs that determine the GDS's *content*: the annotation package, both
detection-p thresholds, and content hashes of the samplesheet and the
cross-reactive CSV. On a later run with the same `gdsOutput`:

- **Key matches** ⇒ the existing GDS is reused silently.
- **Key differs** (you changed a threshold, the annotation package, the
  samplesheet, or the probe list) ⇒ it is rebuilt automatically.
- **No sidecar** (an older or hand-built GDS) ⇒ it is reused with a warning that
  the configuration cannot be verified.

#### Staged rebuilds

The pipeline is a linear chain — **build → normalize → your analysis** — so the
`rebuild` argument names *where to start over*: the named stage **and every stage
after it** are recomputed, everything before is reused.

| `rebuild` | rebuilds | reuses |
|---|---|---|
| `"none"` (default) / `FALSE` | nothing | the build, `normbetas`, your cached outputs |
| `"normalize"` | normalisation + everything downstream | the built GDS |
| `"build"` / `"all"` / `TRUE` | everything, from the IDAT read down | — |

This replaces a pile of independent toggles with one coherent knob: you can't ask
to reuse a stage whose input you just rebuilt. The headline case is *"I changed
`outlyxPerc` — redo from normalisation, don't re-read 200 IDATs"* (`rebuild =
"normalize"`). A **build-cache miss** (changed thresholds / annotation /
samplesheet) rebuilds the GDS *and* cascades downstream automatically, so you only
name a stage to force a redo whose inputs the key doesn't track.

`analyze()` exposes the resolved decision to your `analysis` as
**`res$rebuildDownstream`** (`TRUE` if anything upstream was rebuilt) — use it to
gate your own cached steps. The pattern: recompute when the input is fresh,
otherwise load the saved result. For example, cache a PCA and only redo it when
`normbetas` actually changed:

```r
res <- analyze(
  dataDirectory       = "/path/to/idats",
  samplesheet         = "/path/to/samplesheet.csv",
  crossReactiveProbes = "/path/to/cross-reactive.csv",
  gdsOutput           = "/path/to/my_study.gds",
  annotationPackage   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  outlierRemoval = "outlyx",
  normalization  = "dasen",          # writes/refreshes normbetas (or reuses it)
  analysis = function(gds, res) {
    pcaCache <- "my_study_pca.rds"
    if (res$rebuildDownstream || !file.exists(pcaCache)) {
      message("Computing PCA")        # normbetas is fresh -> (re)compute
      pca <- prcomp(gds, node.name = "normbetas", method = "quick")
      saveRDS(pca, pcaCache)          # cache it for next time
    } else {
      message("Reusing cached PCA")   # nothing upstream changed -> reuse
      pca <- readRDS(pcaCache)
    }
    plot(pca)                         # plot either way
    invisible(pca)
  })
```

On the first run (or after `rebuild = "normalize"`, or any build-cache miss)
`res$rebuildDownstream` is `TRUE`, so the PCA is computed and saved; a plain re-run
that touched nothing upstream leaves it `FALSE`, so the saved PCA is read back and
plotted without recomputing. The built-in `normalization = "dasen"` hook honours
the same flag: it skips when a `normbetas` node exists and nothing upstream
changed. (`pipeline.R` uses this exact pattern to cache its cell-count estimates.)

`compress` is deliberately *not* part of the build key: it changes only how the
data is stored, not its values, so switching codecs alone does not trigger a
rebuild. (`forceRebuild` still works as a deprecated logical alias —
`TRUE → "all"`, `FALSE → "none"` — but emits a warning; prefer `rebuild`.)

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

Omit `analysis` **and** leave both hooks `"none"` (the defaults) to stop after
preprocessing. `analyze()` returns the result list (`gds_path`, `targets`,
`keepSamples`, `keepProbes`) and closes the GDS, so you can open it later and
normalise/analyse however you like. (Setting a hook but no `analysis` is also valid —
e.g. `normalization = "dasen"` alone builds a *normalised* GDS and returns its path.)

```r
res <- analyze(
  dataDirectory       = "/path/to/idats",
  samplesheet         = "/path/to/samplesheet.csv",
  crossReactiveProbes = "/path/to/cross-reactive.csv",
  gdsOutput           = "my_study.gds",
  annotationPackage   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  readerWorkers       = min(parallel::detectCores(), 16L))
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
