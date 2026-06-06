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

**`analyze()` (`R/analyze.R`) is the sole exported function** (`NAMESPACE` is
hand-maintained, not roxygen-generated). It is a loan-pattern driver: one call
runs the preprocessing pipeline, optionally hands the open GDS handle to a
user-supplied `FUN`, and **guarantees the GDS is closed** on normal return and
on error (via `on.exit`, not the closure). `FUN` is optional —
`analyze(..., FUN = NULL)` just builds the QC'd GDS and returns its path +
metadata, subsuming the old `runPreprocess()` public role. It is deliberately
unopinionated about the analysis: **normalisation (dasen) and outlier detection
(outlyx) both live in `FUN`**, where the caller controls their order (correct
order: outlyx on raw betas first, then dasen). `analyze()` itself does no
normalisation. The heavy operations are dependency-injected
(`.preprocess`/`.open`/`.close`) so the control flow is unit-testable without
IDATs/minfi/gdsfmt/bigmelon.

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
  QC'd-GDS-out step `analyze()`'s default `.preprocess` calls. `.validateArgs()`
  checks every input up front (types, ranges, file existence, annotation-package
  availability, samplesheet structure, IDAT presence, cross-reactive CSV
  structure) before any IDAT byte is read, and a build-key sidecar caches
  results so re-runs that change thresholds/inputs rebuild instead of silently
  reusing a stale GDS.

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
  compaction rewrite rather than an in-place edit. The `layout = "methylset"`
  branch has no GDS to stream into and still assembles matrices — streaming is the
  `gds` path only. Do not reintroduce whole-matrix assembly.
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
  mistaken for a valid cache.

## Performance: what moves the needle, and the dead ends

Recorded so the same experiments are not repeated; this is the operative record
(distilled from `dev/` scratch design docs that have since been folded in here).
The harness is `dev/benchmark_analysis.R`, run **only** via `dev/run-benchmark.sh`
(a `systemd-run --user --scope` envelope — no resource caps belong in `.ilo.rc`).
Figures are 450k, the fixed ~14/17 GiB envelope, N = 50/100/200.

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
  - *Parallel `dasen`* — serial quantile barrier, net-negative prototype. The one
    un-pulled lever on the analysis phase is bigmelon's `dasenrank` (subsamples the
    reference — an approximation, so opt-in only, with an equivalence check).

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
