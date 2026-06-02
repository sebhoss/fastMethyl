<!--
SPDX-FileCopyrightText: 2026 The fastMethyl authors
SPDX-License-Identifier: Artistic-2.0
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## What this repo is

`fastMethyl` is a standalone R/Bioconductor package — a **performance layer over
minfi and bigmelon** for large Illumina methylation cohorts. It does not
reimplement or vendor either dependency: it imports them and adds faster paths
for the steps that dominate runtime and memory. The public API is camelCase
(Bioconductor convention), deliberately shadowing the minfi reader names.

Exported functions (`NAMESPACE` is hand-maintained, not roxygen-generated):

- `readMethArray()` / `readMethArrayExp()` / `readMethArraySheet()`
  (`R/read-meth.R`) — minfi-compatible IDAT readers with a `BPPARAM=` argument
  for parallel reading and an integer-indexed worker assembly hot path
  (~11× faster than the serial upstream reader on large cohorts). Return minfi
  `RGChannelSet` objects.
- `processMethArray()` / `processMethArrayExp()` (`R/read-meth-streaming.R`) —
  the fused, column-streaming pipeline: read + detectionP + preprocessRaw +
  es2gds + sample/probe QC in one pass that writes a bigmelon-compatible GDS
  **without ever assembling a full cohort-sized matrix in RAM**.
  `processMethArrayExp()` is the experiment-level wrapper (samplesheet + base
  dir + annotation package).
- `runPreprocess()` (`R/pipeline-preprocess.R`) — the one validated entry
  point: config in, QC'd GDS out. `.validateArgs()` checks every input up front
  (types, ranges, file existence, annotation-package availability, samplesheet
  structure, IDAT presence, cross-reactive CSV structure) before any IDAT byte
  is read, and a build-key sidecar caches results so re-runs that change
  thresholds/inputs rebuild instead of silently reusing a stale GDS.
- `analyze()` (`R/analyze.R`) — the recommended high-level driver: a loan
  pattern that runs `runPreprocess`, opens the GDS, optionally normalises it
  (dasen), hands the open handle to a user-supplied `FUN`, and **guarantees the
  GDS is closed** on normal return and on error (via `on.exit`, not the
  closure). The heavy operations are dependency-injected so the control flow is
  unit-testable without IDATs/minfi/gdsfmt/bigmelon.

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
| **Unit** | `test-unit-*.R` | yes | mocked mixed-cohort dispatch (needs real minfi classes) |
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
- **test-integration** — `R CMD INSTALL .` then
  `testthat::test_dir("tests/testthat", filter = "integration")` with
  `FASTMETHYL_RUN_INTEGRATION=1`; runs only after the fast checks pass.
- **bioc-check** — `BiocCheck::BiocCheck(".")`, advisory.
- **reuse** (`reuse.yml`) — REUSE/SPDX compliance over every file.

## Architecture notes worth knowing up front

- **Column-streaming GDS write (the core of `processMethArray`).** Samples are
  processed in column blocks: workers compute each sample's
  methylated/unmethylated/detection-p columns in parallel, and the master
  **appends** the surviving columns straight into the extendable,
  LZ4_RA-compressed `betas`/`methylated`/`unmethylated`/`pvals` nodes. Full
  matrices are never held in RAM, so **peak memory is flat in cohort size** — a
  1000-sample EPIC run peaks at roughly the same memory as a 100-sample one.
  Do not reintroduce whole-matrix assembly. `dev/streaming-writes-design.md`
  documents the design and its trade-offs.
- **QC interacts with streaming asymmetrically.** Sample QC
  (`sample_detP_threshold`) is per-column and decided as each sample streams
  past (dropped samples are simply not appended). Probe QC
  (`probe_detP_threshold`) is a reduction over all surviving samples, so its
  mask is known only after the full pass; it is accumulated into a per-probe
  boolean and, when it drops any probe, the matrix nodes are compacted to their
  final probe set by a **chunked row-subset rewrite** (`.streamingCompactRows`).
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

`dev/Containerfile` is based on `bioconductor/bioconductor_docker:devel` and
pre-installs every Bioc dependency the tests need (minfi, minfiData,
minfiDataEPIC, bigmelon, wateRmelon, illuminaio, BiocParallel, the 450k/EPIC
manifest + annotation packages, FlowSorted.Blood.450k, covr, lintr, ...). The
cold build takes ~30 min; afterwards `ilo` reuses the cached image.
