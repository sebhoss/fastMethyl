<!--
SPDX-FileCopyrightText: The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## What this repo is

`fastMethyl` is a standalone R/Bioconductor package — a **performance layer over
minfi and bigmelon** for large Illumina methylation cohorts. It does not
reimplement or vendor either dependency: it imports them and adds faster paths
for the steps that dominate runtime and memory.

**`analyze()` (`R/analyze.R`) is the primary exported function** (`NAMESPACE` is
hand-maintained, not roxygen-generated). It is a loan-pattern **orchestrator**:
one call runs the preprocessing pipeline, then drives the analysis phases **in a
fixed, statistically correct order** and **guarantees the GDS is closed** on
normal return and on error (via `on.exit`, not the closure). The order is:

    build  →  outlierRemoval  →  normalize(keep)  →  FUN

`outlierRemoval` and `normalize` are **uniform hooks**: each takes a built-in
method name (`"outlyx"` / `"dasen"`), a user function, or `"none"` (the default).
`outlierRemoval(gds, res)` returns a logical **keep mask** over the GDS samples;
`normalize(gds, res, keep)` writes a normalised node. Built-in `"dasen"` is
`gdsDasen(gds, keep = keep)` — the reference is built from the survivors but
`normbetas` is written for **every** sample, so the GDS stays rectangular and a
flagged-but-kept outlier still gets a column. The keep mask is exposed to `FUN`
as `res$outlierKeep`. **The ordering is the point** — analyze() enforces
outliers-before-normalisation (dasen's between-array quantile reference must
exclude outliers) so a caller cannot get it wrong. Both hooks default to `"none"`,
so an `analyze()` call with neither set is unopinionated and equivalent to a plain
build + `FUN` (the legacy "do outlyx/dasen inside `FUN`" pattern still works).
`analyze(..., FUN = NULL)` with both hooks `"none"` just builds the QC'd GDS and
returns its path + metadata, subsuming the old `runPreprocess()` public role. The
heavy operations are dependency-injected (`.preprocess`/`.open`/`.close`) so the
orchestration is unit-testable without IDATs/minfi/gdsfmt/bigmelon.

`gdsDasen()` (`R/streaming-normalize.R`) is the streaming between-array
normaliser the `"dasen"` hook uses: a two-pass quantile normalisation that
reproduces `wateRmelon::dasen` **bit-for-bit** (pinned by `test-unit-normalize.R`)
while reading the GDS one column block at a time — **no full matrix, no scratch
GDS, bounded/flat peak memory** (see Performance). It depends only on
stats+gdsfmt, so it adds no package dependency. `keep` selects the reference
samples; `normbetas` is written for all.

The read-side **adapters** `gdsBetaMatrix()` and `gdsSummarizedExperiment()`
(`R/adapters.R`) project a finished GDS into the container a non-GDS-aware
downstream package wants — a dense `individuals × variables` matrix for
FactoMineR / `stats::prcomp` / irlba / umap, or a `SummarizedExperiment` with
rowData/colData from fData/pData for limma / sva. They accept the open handle
`FUN` receives or a GDS path (opened read-only, closed on return). The GDS stays
the canonical low-memory store; an adapter materialises a representation **on
demand, sized to the caller** — so it is an explicit call, not an `analyze()`
output mode (that would reintroduce whole-matrix assembly). GDS-native analysis
(bigmelon `outlyx`/`prcomp`/`estimateCellCounts.gds`) needs no adapter, and
bigmelon's `prcomp(node.name=, method="quick")` is the fastest PCA route here
(see Performance).

Everything `analyze()` composes is **internal** (defined in `R/`, reached from
tests via `fastMethyl:::` or the bare-name bindings in `helper-internals.R`):

- `readMethArray()` / `readMethArrayExp()` / `readMethArraySheet()`
  (`R/read-meth.R`) — minfi-compatible IDAT readers with a `BPPARAM=` argument
  for parallel reading and an integer-indexed worker assembly hot path
  (~11× faster than the serial upstream reader). Return minfi `RGChannelSet`s.
  Internal API is camelCase (Bioconductor convention).
- `processMethArray()` / `processMethArrayExp()` (`R/read-meth-streaming.R`) —
  the fused, column-streaming pipeline: read + detectionP + preprocessRaw +
  es2gds + sample/probe QC in one pass that writes a bigmelon-compatible GDS
  **without ever assembling a full cohort-sized matrix in RAM**.
  `processMethArrayExp()` is the experiment-level wrapper.
- `runPreprocess()` (`R/pipeline-preprocess.R`) — the validated config-in /
  QC'd-GDS-out step `analyze()`'s default `.preprocess` calls. Its user-facing
  inputs are **paths**: `dataDirectory` (IDAT root), `samplesheet` (the CSV path
  itself — *not* a `samplesheet_<x>.csv` convention), `crossReactiveProbes` (the
  drop-list CSV), `gdsOutput` (the literal output `.gds` path, used verbatim),
  and `annotationPackage`. `analyze()` resolves **deprecated aliases** (with a
  warning) for the old names: `targetPattern`→`samplesheet` (reconstructing the
  old `dataDirectory/samplesheet_<pattern>.csv` path), `datasetClass`→`gdsOutput`
  (appending `.gds`), `nonSpecificProbesPath`→`crossReactiveProbes`, `FUN`→
  `analysis`, `normalize`→`normalization`, `forceRebuild`→`rebuild`.
  `.validateArgs()` checks every input up front (types, ranges, file existence,
  annotation-package availability, samplesheet structure, IDAT presence,
  cross-reactive CSV structure) before any IDAT byte is read, and a build-key
  sidecar caches results so re-runs that change thresholds/inputs rebuild instead
  of silently reusing a stale GDS.

Supporting internals: `R/minfi-internals.R` (small pure minfi helpers
reproduced so the package needs only minfi's public API), `R/parallel-helpers.R`
(verbose normalisation, the fault-tolerant `.bpmapply_*` wrappers, memory/size
diagnostics), `R/zzz.R` (`.onAttach` banner).

The README is the user-facing entry point. `inst/scripts/pipeline.R` is the
shippable worked example and is built around `analyze()`.

## Running anything

This machine has **no local R**. Every R / R CMD / Rscript command must run
inside the `ilo` container defined by `dev/Containerfile`. Invoke as:

```
ilo bash -c '<COMMAND>'
```

`.ilo.rc` (not committed; local paths) selects the image:

```
shell
--containerfile dev/Containerfile
dev/fastmethyl:latest
```

Docker tags must be lowercase (`dev/fastmethyl`, not `dev/fastMethyl`). `ilo`
expands `$()`/`$VAR` itself — avoid shell substitution in commands; use fixed
paths. The project is mounted at the same absolute path inside the container;
`~` is not writable, so use `/tmp` or `dev/scratch` for scratch output.

For long-running tests/benchmarks wrap with `nice -n 19 ionice -c 3` (the
maintainer has hit OOM/crash issues on large IDAT cohorts). Note `ionice -c 3`
is *idle* I/O priority — if the box has other disk activity it can starve a
heavy GDS write; drop to `ionice -c 2 -n 7` if a run stalls on I/O rather than
memory.

## Tests

The testthat suite under `tests/testthat/` is split into three buckets by
speed, distinguished by filename:

| Bucket | Pattern | Loads minfi? | Covers |
|---|---|---|---|
| **Pure** | `test-pure-*.R` | no | base-R helpers: verbose normalisation, array-type guessing, probe-address alignment, QC mask logic, arg validation, build-key inspection, `analyze()` control flow, memory/size helpers |
| **Unit** | `test-unit-*.R` | yes | mocked mixed-cohort dispatch (needs real minfi classes); the install contract (`test-unit-install.R`: the resolved minfi exports every symbol/class fastMethyl imports, every NAMESPACE-imported package is in DESCRIPTION, `analyze()`'s signature) |
| **Integration** | `test-integration-*.R` | yes | real minfiData / minfiDataEPIC IDAT reads, GDS writes, byte-identity equivalence, QC + compaction, `analyze()` end-to-end. Gated by `FASTMETHYL_RUN_INTEGRATION=1`. |

Integration tests `skip_if_slow()` (which skips unless
`FASTMETHYL_RUN_INTEGRATION=1`) and `skip_if_no_gds()` (needs
bigmelon + gdsfmt). Pure helpers are exposed to tests by bare name via
`helper-internals.R`; shared fixtures live in `helper-fixtures.R`.

Iteration loop (host-persisted install lib so repeated runs skip reinstall):

```
# install once into dev/scratch/rlib (host-mounted, persists across runs)
ilo bash -c 'R CMD INSTALL --library=dev/scratch/rlib .'

# fast buckets only (pure + unit; integration skipped):
ilo bash -c 'Rscript -e ".libPaths(\"dev/scratch/rlib\"); library(testthat); library(fastMethyl); test_dir(\"tests/testthat\")"'

# everything, including the real-IDAT tier:
ilo bash -c 'Rscript -e ".libPaths(\"dev/scratch/rlib\"); Sys.setenv(FASTMETHYL_RUN_INTEGRATION=\"1\"); library(testthat); library(fastMethyl); test_dir(\"tests/testthat\")"'
```

**Re-install after editing any `R/` file** (`R CMD INSTALL --library=dev/scratch/rlib .`);
testthat runs against the installed namespace, not the source tree. For heavy
iteration, `devtools::load_all(".")` inside a persistent `ilo shell` reloads
`R/` without a full reinstall.

The **equivalence oracle**: the unfiltered GDS run must be byte-identical
(read-back values, compared via `read.gdsn`) to
`preprocessRaw(readMethArray(...))` + `detectionP(...)`, and every QC-arg path is
checked against a post-hoc subset of that unfiltered run. Any change to
`R/read-meth.R` or `R/read-meth-streaming.R` must keep those tests green.

## Linting

CI runs `lintr::lint_dir("R")` over the whole `R/` tree. Config is `.lintr`
(read automatically): defaults plus `line_length_linter(120L)`, with
`object_usage_linter` and `object_name_linter` disabled (the latter because the
public API is camelCase by design). `lintr` is preinstalled in the container:

```
ilo bash -c 'Rscript -e "print(lintr::lint_dir(\"R\"))"'
```

Keep `R/` lint-clean — CI fails on any lint. Multi-line function signatures use
paren-aligned hanging indent; the linter is picky about the exact column.

## CI (`.github/workflows/pr-checks.yml`)

- **lint** — `lintr::lint_dir("R")`, fails on any lint.
- **R CMD check** — runs `tests/testthat.R` (`test_check`); the integration tier
  self-skips because `FASTMETHYL_RUN_INTEGRATION` is unset there.
- **install** — `dev/install_test.sh`: builds the tarball, installs fastMethyl
  alone into a fresh empty library, and loads it in a clean process (the closest
  reproduction of `BiocManager::install("sebhoss/fastMethyl")`). Run it locally
  with `ilo bash dev/install_test.sh`.
- **test-integration** — `R CMD INSTALL .` then
  `testthat::test_dir("tests/testthat", filter = "integration")` with
  `FASTMETHYL_RUN_INTEGRATION=1`; runs only after the fast checks pass.
- **bioc-check** — `BiocCheck::BiocCheck(".")`, advisory.
- **reuse** (`reuse.yml`) — REUSE/SPDX compliance over every file.

## Architecture notes worth knowing up front

- **Column-streaming GDS write (the core of `processMethArray`).** Samples are
  processed in column blocks: workers compute each sample's
  methylated/unmethylated/detection-p columns in parallel, and the master
  **appends** the surviving columns straight into the extendable
  `betas`/`methylated`/`unmethylated`/`pvals` nodes, written with the `compress`
  codec (**default `compress = ""`, uncompressed**; pass `LZ4_RA` to shrink — see
  the Performance section). Full matrices are never held in RAM, so **peak memory
  grows far slower than cohort size** — a 1000-sample EPIC run peaks far below a
  whole-matrix build. The thing being avoided: holding `Meth`+`Unmeth`+`detP`+
  `Betas` at once is **24 bytes/cell** (`int+int+double+double`), a ~21 GB floor
  for a 1000-sample EPIC cohort. `append.gdsn` extends a node along the sample
  (column) dimension and works on a still-in-write-mode compressed node; once a
  node hits readmode it cannot be modified, which is why probe-row drops need the
  compaction rewrite rather than an in-place edit. The streaming write is the
  **only** path: `processMethArray` always streams into a GDS — there is no
  in-memory / MethylSet output mode, and the read-side `gdsBetaMatrix()` /
  `gdsSummarizedExperiment()` adapters are how an in-memory representation is
  obtained, sized to the caller, *after* the GDS exists. Do not reintroduce
  whole-matrix assembly on the write side.
- **QC interacts with streaming asymmetrically.** Sample QC
  (`sample_detP_threshold`) is per-column and decided as each sample streams
  past (dropped samples are simply not appended). Probe QC
  (`probe_detP_threshold`) is a reduction over all surviving samples, so its
  mask is known only after the full pass; it is accumulated into a per-probe
  boolean and, when it drops any probe, the matrix nodes are compacted to their
  final probe set by a **chunked row-subset rewrite**, run in parallel across the
  four independent nodes (`.streamingCompactParallel`, each spliced back via
  `copyto.gdsn`; `.streamingCompactRows` is the 1-worker fallback).
  With all QC args `NULL` the GDS values equal
  `es2gds(preprocessRaw(readMethArray(...)))`.
- **Integer-indexed worker hot path.** Both readers pre-compute integer
  position arrays (`match(address, refAddr)`) once in master and ship them to
  workers via `MoreArgs`; workers never do character row lookups. This is the
  dominant perf win versus upstream. Don't reintroduce
  `Quants[address, "Mean"]`-style character indexing.
- **Multi-version dispatch is always on.** A parallel header-only pre-pass
  (`readIDAT(f, what = "nSNPsRead")`) detects distinct probe counts; for uniform
  cohorts it degenerates to a list-of-one at negligible cost, for mixed cohorts
  (e.g. early-access vs final-release EPIC) each worker gets a per-version index
  table. No `force =` argument — auto-detection makes it unnecessary.
- **`verbose` is an integer** (`0L` silent, `1L` phase + per-sample lines, `2L`
  adds master-side progress and detailed diagnostics). Legacy `TRUE`/`FALSE` is
  accepted via `.normalize_verbose()`. At **`verbose = 2L`** every external-data
  load — preprocessing, opening the GDS, normalisation, loading the annotation
  package, `getAnnotation()` — is logged with the file's on-disk size and a
  memory snapshot (`.mem_report()`, RSS+peak on Linux) before and after, so a
  slow or stalling load announces itself before it blocks. Preserve this when
  adding code that loads external data.
- **`analyze()` reads `verbose` out of `...`** (not a formal), so it still flows
  unchanged to `runPreprocess`; it only gates `analyze`'s own diagnostics.
- **Build-key cache.** `runPreprocess` writes a sidecar keyed on annotation
  package + thresholds + samplesheet + cross-reactive CSV next to the GDS; a
  matching key reuses the GDS silently, a missing one reuses with a warning, a
  mismatched one rebuilds. A failed/partial GDS write is unlinked so it is never
  mistaken for a valid cache. `runPreprocess` returns `res$rebuilt` (TRUE when it
  actually (re)wrote the GDS, i.e. a forced build OR a cache miss/mismatch).
- **Staged `rebuild`.** The pipeline is a linear chain
  (`build → normalize → FUN's outputs`), so `analyze()`'s `rebuild` arg names a
  stage and rebuilds it **and every stage after it**; earlier stages are reused.
  `.resolveRebuild()` (`R/analyze.R`) maps `"none"`/`"build"`/`"normalize"`/`"all"`
  (or a logical: `FALSE`→none, `TRUE`→all) to per-stage booleans for the two
  analyze-owned stages. The build stage is realised through `runPreprocess`'s
  `forceRebuild` boolean (analyze translates `rebuild` → `forceRebuild` in `...`
  and drops `rebuild`); the normalize stage is enforced via
  **`res$rebuildDownstream`** = `rebuild names normalize` OR `res$rebuilt` — so a
  build-cache miss cascades downstream automatically. The cache-aware built-in
  `"dasen"` hook skips when a `normbetas` node exists and `!res$rebuildDownstream`;
  `FUN` reads `res$rebuildDownstream` to gate its own cached steps (cell counts in
  `pipeline.R`). `forceRebuild` is a **deprecated logical alias** (warns) for
  back-compat; `runPreprocess`/`.validateArgs` keep `forceRebuild` as their
  internal formal (the build mechanism), not deprecated.

## Performance: what moves the needle, and the dead ends

Recorded so the same experiments are not repeated; this is the operative record
(distilled from `dev/` scratch design docs that have since been folded in here).
The harness is `dev/benchmark_analysis.R` (preprocessing: IDAT → QC'd GDS,
fastMethyl vs upstream minfi), `dev/benchmark_pca.R` (analysis-phase PCA routes on
a normalised GDS), `dev/benchmark_pipeline.R` (per-phase breakdown of
`pipeline.R`'s `FUN`: build / outlyx / dasen / prcomp / estimateCellCounts.gds),
`dev/benchmark_streamdasen.R` (streaming `gdsDasen` vs bigmelon `dasen.gds`), and
`dev/benchmark_fullpipeline.R` (end-to-end upstream minfi+bigmelon vs fastMethyl,
both ways), all run **only** via `dev/run-benchmark.sh` (a `systemd-run --user
--scope` envelope — no resource caps belong in `.ilo.rc`); select the script with
`BENCH_SCRIPT=dev/benchmark_pca.R dev/run-benchmark.sh`. Figures are 450k, the
fixed ~14/17 GiB envelope, N = 50/100/200. **Run a timing benchmark with NOTHING
else touching the box** — concurrent installs/tests/lint starve the throttled GDS
writes and inflate timings (a contaminated run logged a N=50 upstream build at
1851 s vs 185 s clean); the benchmarks now print each cohort live (a `-> N=...`
line via stderr) so a long sweep shows progress.

- **Full-pipeline apples-to-apples (`dev/benchmark_fullpipeline.R`, clean run,
  450k, N = 50/100/200).** The complete workflow both ecosystems run — upstream
  `read.metharray.exp`+`preprocessRaw`+QC+`es2gds` then
  `outlyx`/`dasen`/`prcomp`/`estimateCellCounts.gds`, vs `analyze(FUN=NULL)` then
  the same analysis with `gdsDasen`. The **analysis phase is identical code both
  sides except dasen**, so it dilutes the win; the difference is the build + the
  streaming dasen. Total time **307.5/573.1/1135.9 s upstream vs
  156.7/218.0/347.4 s fastMethyl** = **1.96×/2.63×/3.27×** (growing). The
  **build** alone (the part fastMethyl reimplements) is **6–9×**
  (185→30, 396→44, 819→99 s); the analysis is a wash-to-slight-win
  (gdsDasen edges dasen.gds at N=200: 249 vs 317 s). **Memory is the headline:**
  upstream peak **3.9/7.5/12.5 GiB** (climbs — holds the full matrices), fastMethyl
  **3.2/3.0/4.0 GiB** (bounded). fastMethyl's *build* peak is near-nothing
  (110 MiB at N=50, 3.5 GiB at N=200), so its full-pipeline peak is set by the
  **shared `estimateCellCounts` reference load**, not by anything fastMethyl does
  — which is why it stays ~4 GiB while upstream grows unbounded (a N=1000+ cohort
  OOMs upstream but stays feasible on fastMethyl). ~3.1× less memory at N=200,
  widening.

- **GDS compression is the single biggest runtime cost — which is why `compress`
  defaults to `""` (uncompressed).** Opt-in `LZ4_RA` is **~2× slower** (44.7/76.4/
  144.3 s vs 23.1/34.8/69.2 s at 4 cores) for only **~1.3× smaller** files —
  methylation intensity/beta matrices barely compress. `compress` threads
  `analyze()` → `runPreprocess()` → `processMethArrayExp()` to the node writes and
  compaction, is value-identical and bigmelon-readable either way, and is
  deliberately **not** a build-key field.
- **The structural win is NOT compression — don't re-investigate that.** Both
  fastMethyl and the upstream minfi benchmark write the same matrices, so they pay
  ~equal *absolute* compression cost — but that is ~7% of minfi's serial runtime
  and ~50% of fastMethyl's. Removing it from both *widens* the gap (4 cores,
  N=200: **4.6× compressed → 8.9× uncompressed**; minfi 612.8 s vs fastMethyl
  69.2 s). The lead is the parallel read + fused streaming, not a faster codec.
- **Core scaling depends on the codec.** Uncompressed (the default) is
  *read-bound* at large N, so more cores help: N=200 goes 69 → 65 → **50 s** over
  4 → 8 → 12 workers (~27%); at small N per-worker overhead dominates and extra
  cores regress (12 is slowest at N=50). Under opt-in `LZ4_RA` the serial
  master-side compression caps it — 4 → 8 → 12 buys almost nothing (144 → 135 →
  136 s). Either way peak RAM climbs with worker count (**COW inflation**, below).
- **Implemented and kept — parallel probe-QC compaction**
  (`.streamingCompactParallel`): the four matrix nodes are independent, so each is
  compacted to a temp GDS by a worker and the master splices it back with
  `copyto.gdsn` (a cheap **raw block copy, not a re-encode**). ~10% end-to-end at
  `LZ4_RA` (near-free once uncompressed). Trade-off: the four concurrent node
  buffers raise the N=200 peak, so peak is not strictly flat in N.
  `.streamingCompactRows` (serial) is the `nworkers == 1` fallback.
- **Tried and REVERTED — do not re-attempt without new evidence:**
  - *Betas-to-workers* — byte-identical but neutral-to-negative under
    `MulticoreParam`: shipping a third column back costs more IPC than the
    master's cheap `pmax(M,0)/(pmax(M,0)+pmax(U,0))` saves. **Master computes betas.**
  - *Compute/write overlap in the streaming loop* — achievable ceiling ~9 s, not
    worth the complexity.
  - *Sex-probe address cache* — a no-op (the lookup is already negligible).
  - *Unconditional `gc()` per block* — now `if (forking) gc(FALSE)`; for
    `MulticoreParam` effectively a no-op, kept only to bound the non-forking path.
  - *Parallel `dasen`* — serial quantile barrier, net-negative prototype.
  - **Streaming `dasen` — prototyped and validated, the right lever (`dev/benchmark_streamdasen.R`).**
    A two-pass between-array QN that reproduces `wateRmelon::dasen` **bit-exactly**
    (max diff ~1e-15 = float epsilon, at N = 50/100/200), without holding the full
    matrix and without `dasenrank`'s `temp.gds`. It works because the QN reference
    is `rowMeans` of the sorted columns (mean of order statistics) — accumulable in
    one streaming pass; a second pass dfsfits each column in RAM (cache the dfs2
    offset from pass 1, do not re-density) and maps ranks onto the reference.
    Results (450k): time **44.9/77.9/153.6 s** (0.99×/1.38×/1.06× vs bigmelon
    `dasen`'s 44.3/107.6/162.4 s — competitive-to-faster, no penalty); **peak
    bounded ~0.8–1.2 GiB and flat in N** vs bigmelon's 2.2→3.4 GiB *growing* (~4.3×
    less at N=200, gap widens). Because it is exact, the `dasenrank` subsample
    (`perc<1`) is an *optional* speed knob, never a necessity. **Crucially it must
    run AFTER outlier removal** (the reference must exclude outliers; it cannot be
    fused into the build, which runs before outliers are known — only within-sample
    correction like noob can fuse into the build). Pending promotion into the
    package as the `analyze(normalize=)` between-array path.
  - *FactoMineR / irlba PCA instead of bigmelon's GDS-native `prcomp`* — both
    **lose decisively** (`dev/benchmark_pca.R`, uncompressed, N = 50/100/200,
    `normbetas`, ncp = 10). bigmelon `prcomp(method = "quick")` runs in
    **0.2/0.3/0.6 s at ~0–2 MiB peak**; irlba (`prcomp_irlba` over a materialised
    matrix) is **8.8/13.7/18.4 s at ~1–2 GiB**, `stats::prcomp` worse
    (17/51/85 s), FactoMineR worst (**27.5/45.2/73.7 s at ~1.8–2.4 GiB**). The
    structural reason is the data shape: methylation is `p ≫ n`, so the GDS-native
    route only ever forms the small `n × n` cross-product and reads straight from
    the node, while every matrix route must first `gdsBetaMatrix()` the whole
    `samples × probes` matrix into RAM (3.5–7.7 s on its own) before decomposing.
    The `gdsBetaMatrix()` / `gdsSummarizedExperiment()` adapters exist for tools
    that genuinely need an in-RAM matrix or a `SummarizedExperiment` (FactoMineR
    diagnostics, limma, sva), **not** as a faster PCA — `pipeline.R` stays on
    `prcomp`. Do not swap the PCA route for speed without new evidence.

### Analysis-phase cost breakdown (`pipeline.R`'s `FUN`)

`dev/benchmark_pipeline.R` times each stage of the shipped pipeline separately on
a real 450k GDS (uncompressed, 4 cores, N = 50/100/200), in run order. Time (s) /
peak (MiB):

| phase | N=50 | N=100 | N=200 | scaling |
|---|---|---|---|---|
| build (`analyze(FUN=NULL)`) | 28.7 / 3313 | 47.6 / 5094 | 80.0 / 6584 | ~linear in N |
| `outlyx` (raw betas) | 1.0 | 1.4 | 2.1 | negligible, ~flat |
| **`dasen` → normbetas** | 47.7 / 2642 | 93.7 / 3852 | **176.2 / 4953** | **~linear in N, serial** |
| `prcomp` quick | 0.2 | 0.3 | 0.7 | free |
| `estimateCellCounts.gds` | 75.6 / 3976 | 79.7 / 3987 | 89.0 / 4451 | large but ~**fixed** in N |

The actionable conclusions, which **invert with cohort size**:

- **`dasen` is the scaling bottleneck and the single largest phase at scale**
  (176 s at N=200, > the build itself). It is serial — the between-sample
  quantile reference is a hard barrier (every sample maps onto the average of all
  samples' sorted intensities, so no sample's normbetas is final until the last is
  read; this is why *parallel `dasen`* was net-negative and stays reverted). The
  only real lever is bigmelon's `dasenrank` (subsampled reference, an
  approximation → opt-in with an equivalence check). This is the highest-ROI
  analysis-phase target for the large cohorts fastMethyl exists for.
- **`estimateCellCounts.gds` is a large but nearly *fixed* cost** (~75–89 s, +18%
  over a 4× cohort): it is dominated by loading the FlowSorted reference + probe
  picking, not the per-sample QP projection. So it dwarfs everything at small N
  but is not a scaling concern; `pipeline.R` already caches it (the
  `forceCellCounts` gate), which is the right mitigation. Its `perc` knob barely
  helps because the cost is the fixed part.
- **`outlyx` (1–2 s) and `prcomp` (<1 s) are noise.** A `detectOutliers` skip is a
  cleanliness nicety (the default `dropOutliers=FALSE` path computes a report it
  never uses), not a speed win; fusing outlyx's per-sample stats into the
  streaming write is **not** worth it for a 2 s phase.

### Memory / copy-on-write inflation (multi-worker forks)

`MulticoreParam` forks; a forked worker COW-breaks shared master pages two ways:
(1) the **GC mark phase** writes a mark bit into every reachable object, copying
~the whole live master heap per worker that GCs; (2) under R ≥ 4.0 **refcounting**,
even *reading* a shared object writes its header — an irreducible floor that GC
suppression cannot beat. The master heap is already minimised at every fork point
(`rm()`+`gc(FALSE)` before and within the block loop; the ~1.8 GB annotation frame
is freed before the reading pass, the fData frame built *after* it), so the live
heap at a fork is only ~40–55 MB (the index tables, which workers need). **Caveat:
the benchmark overstates production COW** — its pre-warmed ~2 GB annotation parent
is forked against, which a real `analyze()` run avoids. So **measure real
production COW (`Private_Dirty` per worker under `analyze()`, not the harness)
before building any mitigation.** Ranked options, none yet shipped: worker
GC-suppression (soft, 450k-only); a budget-capped worker count (OOM-safety, not
speed — best near-term ROI); `SnowParam` (no shared heap → a *hard* memory bound,
but never faster than fork — reserve for EPICv2 / Windows / memory-critical).
Dropped: mmap'd indices (~1–3% of peak), fork-once pools (inflation is
per-running-worker-GC, not per-fork-event).

## Project-specific rules

- **`inst/scripts/pipeline.R` is the shippable worked example**, built around
  `analyze()`. To add a pipeline, copy it and edit the CONFIG block.
- **Internal-only files** (not user-facing, excluded from the tarball via
  `.Rbuildignore`): `CLAUDE.md`, the entire `dev/` directory. Don't name them in
  shipped source (`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`) or the README.
  `dev/scratch/` is gitignored.
- **REUSE / SPDX.** Every file needs an `SPDX-FileCopyrightText` +
  `SPDX-License-Identifier: Artistic-2.0` pair — inline (`#` for R, `<!-- -->`
  for Markdown) or via a `REUSE.toml` glob for files that can't carry comments.
  Always use **Artistic-2.0** (the project license); never guess.
- **camelCase public API.** Exported names are camelCase by design; the
  `object_name_linter` is disabled to allow it. Keep new exports camelCase.
- **No emojis** in committed files unless explicitly requested.
- **Comments describe the current code** — intent and invariants, not the
  history of a change (see the global comment-style guidance).

## Container details

`dev/Containerfile` is based on `bioconductor/bioconductor_docker:latest` (the
current Bioconductor **release**, matching the environment users install into)
and pre-installs every Bioc dependency the tests need (minfi, minfiData,
minfiDataEPIC, bigmelon, wateRmelon, illuminaio, BiocParallel, the 450k/EPIC
manifest + annotation packages, FlowSorted.Blood.450k, covr, lintr, ...). The
cold build takes ~30 min; afterwards `ilo` reuses the cached image.
