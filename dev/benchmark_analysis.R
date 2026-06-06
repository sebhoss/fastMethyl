# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Benchmarks fastMethyl's analyze() preprocessing against the upstream minfi
# pipeline that produces the same QC'd data, across cohort size N.
#
# Scope: this measures IDAT -> QC'd bigmelon GDS, the artifact both ecosystems
# produce before any analysis. Upstream does it as read.metharray.exp +
# detectionP + preprocessRaw + QC + a GDS write; fastMethyl does it as analyze()'s
# single fused streaming pass. Both write the same four matrix nodes with the SAME
# codec (BENCH_COMPRESS, default LZ4_RA; "" = uncompressed), so both pay the
# identical write/compression cost (bigmelon's es2gds writer is not used on the
# upstream side because it would impose its own node layout). BENCH_COMPRESS=""
# removes compression symmetrically, isolating the structural read+streaming win.
# What IS excluded from both columns is the downstream bigmelon *analysis*
# (dasen, outlyx, prcomp, estimateCellCounts.gds): that code is identical
# regardless of which reader produced the GDS, so timing it would only add the
# same work to both.
#
# Both sides apply the same QC: drop samples whose mean detection-p is at or
# above sampleThreshold, drop probes that fail detection in any surviving
# sample, and drop sex-chromosome + cross-reactive probes. Both then end with a
# QC'd bigmelon GDS on disk, so the two columns are the same task by two routes.
#
# Cohorts are scaled by symlinking the six real minfiData 450k IDAT pairs into
# an N-sample on-disk layout, so both pipelines do real on-disk IDAT reads (the
# read cost fastMethyl parallelises is genuinely exercised, not faked from an
# in-memory matrix).
#
# Run it through the wrapper -- the only supported entry point. It sizes the
# cgroup limits to the cohort and launches this script inside a systemd-run
# scope:
#   dev/run-benchmark.sh                # cohorts 50,100,200
#   BENCH_N=50 dev/run-benchmark.sh     # one cohort
# Env knobs (read by the wrapper):
#   BENCH_N        comma-separated cohort sizes
#   BENCH_WORKERS  fastMethyl reader workers (= CPU cores; minfi reads serially).
#                  Also sets the streaming block size, so fewer workers means
#                  smaller write bursts.
#   BENCH_SCRATCH  where cohorts + GDS files are written (default tempdir())
#
# Why a systemd-run scope: a large cohort can trigger a GDS-write storm that
# wedges the host, so the scope bounds CPU/memory/IO and makes runs reproducible.
# The rootless container nests under the scope via --cgroups=split (.ilo.rc), so
# every cap reaches the reader workers too. Hard-won gotchas, all handled by the
# wrapper:
#
#   * Delegate=yes is required: with --cgroups=split the container leaf otherwise
#     exposes only the `pids` controller, so cpu/io/memory limits still enforce
#     but are invisible from inside -- and the harness cannot read memory.current
#     to report peak RAM. Delegate=yes pushes the controllers down into the leaf.
#   * IOWriteBandwidthMax needs the real DEVICE NODE of the whole disk, never a
#     btrfs path (which stat()s to a synthetic dev_t systemd cannot map) and
#     never a partition (io.max keys on the disk request queue). The wrapper
#     resolves it from the scratch dir.
#   * The write cap only throttles buffered writeback (the GDS write) when the
#     cgroup is ALSO under memory pressure: MemoryHigh forces dirty pages to be
#     written back in the cgroup's own reclaim context, where the IO is
#     attributed to it. Alone, IOWriteBandwidthMax leaks past the cap.
#   * nproc inside the container still reports every logical CPU, ignoring
#     CPUQuota, so detectCores() over-subscribes -- BENCH_WORKERS is set to the
#     CPU core count explicitly. cpu/io/memory/pids are delegated to --user
#     scopes but cpuset is NOT, so CPU is capped with CPUQuota, not AllowedCPUs.
#     No nice/ionice: the cgroup caps are absolute, so priority hints add nothing
#     (and ionice is a no-op on NVMe's non-BFQ scheduler).

suppressPackageStartupMessages({
  library(fastMethyl); library(minfi); library(minfiData); library(gdsfmt)
})

annoPkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "12,24,48"), ",")[[1]])
workers <- as.integer(Sys.getenv("BENCH_WORKERS", "4"))
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
# Which sides to run: "both" (default), or "fastmethyl" to time only fastMethyl.
# The reader's worker count is the only thing CPU_CORES changes, and it affects
# *only* fastMethyl -- minfi reads serially, so its numbers are core-independent.
# A core-scaling sweep therefore re-runs only fastMethyl and reuses the minfi
# column from the matching both-sides run.
sides <- Sys.getenv("BENCH_SIDES", "both")
# fastMethyl reader backend: "multicore" (fork, default) or "snow" (socket
# cluster -- independent worker processes, no COW inflation; ships the index
# tables over a socket instead of inheriting them).
bpKind <- Sys.getenv("BENCH_BPPARAM", "multicore")
# GDS codec for the fastMethyl side: "LZ4_RA" (default) or "" (uncompressed).
bpCompress <- Sys.getenv("BENCH_COMPRESS", "LZ4_RA")
naMeasure <- list(time = NA_real_, peak = NA_real_, avg = NA_real_, value = NULL)
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)

# Remove every benchmark artifact this run writes into `scratch` -- cohort
# symlink dirs and the (multi-GB, uncompressed) GDS files. Each pipeline unlinks
# its own GDS on success, but a cohort that is OOM-killed mid-write leaves a
# partial file behind; this sweep reclaims those too, so a sweep never lingers on
# disk regardless of how a run ends. Called after each cohort and once at exit.
.sweepScratch <- function() {
  unlink(list.files(scratch, pattern = "^(bench_fm_|minfi_gds_|bench_cohort_)",
                    full.names = TRUE),
         recursive = TRUE, force = TRUE)
}

# Peak/average memory is read from the cgroup, not from a single process's RSS:
# fastMethyl forks MulticoreParam workers, so a master-only RSS read would miss
# them. memory.current accounts for the whole cgroup (master + workers + buffers)
# and is readable even though the cgroup fs is mounted read-only here; a tiny
# background poller samples it while a pipeline runs, and the figures are taken
# net of the pre-call baseline so each pipeline's own growth is isolated.
cgCurrent <- "/sys/fs/cgroup/memory.current"
hasCg <- file.exists(cgCurrent)
readCur <- function() as.numeric(readLines(cgCurrent, n = 1L))

# Run fn(), returning its value plus elapsed seconds and the peak + mean cgroup
# memory growth (bytes) observed during the call. Memory figures are NA when the
# cgroup memory controller is absent (an uncontained run).
measure <- function(fn) {
  base <- if (hasCg) readCur() else NA_real_
  sf <- tempfile("mem_"); pf <- tempfile("pid_")
  if (hasCg) {
    # $! is expanded by the container shell R spawns, not by ilo, so it resolves
    # to the poller's PID for a clean kill afterwards.
    system(paste0(
      "sh -c 'while true; do cat ", cgCurrent, "; sleep 0.05; done' > ",
      sf, " 2>/dev/null & echo $! > ", pf), wait = FALSE)
    Sys.sleep(0.1)
  }
  # Run the pipeline in a forked child so the parent heap stays clean between
  # cohorts. Otherwise R retains heap across the loop, the per-cohort baseline
  # drifts upward, and the "peak over baseline" delta shrinks for later cohorts
  # (and the first cohort's annotation-frame transient inflates its peak) -- which
  # made fastMethyl's peak non-monotonic. The child's allocations and its own
  # forked workers are reaped on exit, so every cohort measures from the same
  # baseline. The cgroup sampler sees the whole tree (idle parent + child +
  # workers), so the peak still accounts for the workers.
  out <- NULL
  t <- unname(system.time({
    job <- parallel::mcparallel(fn(), detached = FALSE)
    out <- parallel::mccollect(job)
  })[3])
  # The forked child's return value (fastMethyl returns its GDS size in bytes).
  val <- if (length(out)) out[[1L]] else NULL
  peak <- NA_real_
  avg <- NA_real_
  if (hasCg) {
    Sys.sleep(0.05)
    pid <- tryCatch(readLines(pf, n = 1L), error = function(e) NA_character_)
    if (length(pid) && !is.na(pid)) {
      system(paste("kill", pid), ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    s <- suppressWarnings(as.numeric(readLines(sf, warn = FALSE)))
    s <- s[is.finite(s)]
    if (length(s)) {
      peak <- max(s) - base
      avg  <- mean(s) - base
    }
  }
  unlink(c(sf, pf))
  list(time = t, peak = peak, avg = avg, value = val)
}
mib <- function(x) if (is.na(x)) NA_integer_ else as.integer(round(x / 1048576))

sampleThreshold <- 0.01
probeThreshold  <- 0.01
# A handful of real 450k CpG names to exercise the cross-reactive drop filter on
# both sides; the exact list is immaterial to timing, only that both pipelines
# remove the same probes.
xreactIds <- c("cg00000029", "cg00000108", "cg00000109",
               "cg00000165", "cg00000236")

# The six real minfiData 450k IDAT pairs, as <dir>/<barcode> basenames.
mdDir <- system.file("extdata", package = "minfiData")
srcBase <- sub("_Grn\\.idat$", "",
               list.files(mdDir, pattern = "_Grn\\.idat$",
                          recursive = TRUE, full.names = TRUE))
stopifnot(length(srcBase) >= 1L)

# Fabricate an N-sample on-disk cohort by symlinking the real IDAT pairs in a
# cycle, plus the samplesheet + cross-reactive CSV analyze() expects.
makeCohort <- function(n) {
  dir <- tempfile("bench_cohort_", tmpdir = scratch)
  dir.create(dir)
  bns <- sprintf("S%05d", seq_len(n))
  for (i in seq_len(n)) {
    src <- srcBase[(i - 1L) %% length(srcBase) + 1L]
    for (ch in c("_Grn.idat", "_Red.idat")) {
      file.symlink(paste0(src, ch), file.path(dir, paste0(bns[i], ch)))
    }
  }
  write.csv(
    data.frame(Sample_Name = bns, Basename = bns, stringsAsFactors = FALSE),
    file.path(dir, "samplesheet_bench.csv"), row.names = FALSE)
  xreact <- file.path(dir, "xreact.csv")
  write.csv(data.frame(TargetID = xreactIds, stringsAsFactors = FALSE),
            xreact, row.names = FALSE)
  list(dir = dir, bns = bns, xreact = xreact)
}

# Upstream minfi: read + detection-p + raw preprocess + the same QC, then write
# the QC'd matrices to a GDS that matches fastMethyl's output -- same four nodes,
# same codec. bigmelon's es2gds() is deliberately NOT used: it is the standard
# MethylSet->GDS writer but stores the matrices UNCOMPRESSED, so it would let the
# upstream column skip the compression work fastMethyl performs. Writing the
# nodes through gdsfmt with the SAME `compress` codec as the fastMethyl side
# (BENCH_COMPRESS, default LZ4_RA) makes both columns produce the same artifact
# and pay the identical write/compression cost -- so BENCH_COMPRESS="" removes it
# symmetrically from both sides, isolating the structural (read + streaming) win.
# The downstream bigmelon *analysis* (dasen/outlyx/...) is identical for both and
# stays out of both columns.
minfiPipeline <- function(coh) {
  targets <- data.frame(Basename = file.path(coh$dir, coh$bns),
                        stringsAsFactors = FALSE)
  rg   <- read.metharray.exp(targets = targets)
  detP <- detectionP(rg)
  mset <- preprocessRaw(rg)
  keepS <- colMeans(detP) < sampleThreshold
  detP  <- detP[, keepS, drop = FALSE]
  mset  <- mset[, keepS]
  keepP <- rowSums(detP >= probeThreshold) == 0L
  ann   <- getAnnotation(mset)
  sexP  <- rownames(mset) %in% ann$Name[ann$chr %in% c("chrX", "chrY")]
  xrP   <- rownames(mset) %in% xreactIds
  probeMask <- keepP & !sexP & !xrP
  msetQC <- mset[probeMask, ]
  detPQC <- detP[probeMask, , drop = FALSE]
  gp <- tempfile("minfi_gds_", tmpdir = scratch, fileext = ".gds")
  gf <- createfn.gds(gp)
  add.gdsn(gf, "betas",        getBeta(msetQC),   compress = bpCompress,
           closezip = TRUE)
  add.gdsn(gf, "methylated",   getMeth(msetQC),   compress = bpCompress,
           closezip = TRUE)
  add.gdsn(gf, "unmethylated", getUnmeth(msetQC), compress = bpCompress,
           closezip = TRUE)
  add.gdsn(gf, "pvals",        detPQC,            compress = bpCompress,
           closezip = TRUE)
  closefn.gds(gf)
  unlink(gp)
  invisible(NULL)   # NULL, not the MethylSet: avoids serialising it back from
                    # the forked measurement child
}

# fastMethyl: analyze() with FUN = NULL -- the "just build the QC'd GDS" path,
# which is exactly the preprocessing under test. forceRebuild = TRUE so the
# build-key cache never short-circuits a timed run.
fmPipeline <- function(coh) {
  cls <- tempfile("bench_fm_", tmpdir = scratch)
  args <- list(
    dataDirectory         = coh$dir,
    nonSpecificProbesPath = coh$xreact,
    targetPattern         = "bench",
    datasetClass          = cls,
    annotationPackage     = annoPkg,
    sampleDetPThreshold   = sampleThreshold,
    probeDetPThreshold    = probeThreshold,
    forceRebuild          = TRUE,
    compress              = bpCompress,
    verbose               = 0L,
    FUN                   = NULL)
  if (bpKind == "snow") {
    # One chunk per worker (tasks = workers) so the index tables are serialised
    # once per worker per dispatch, not once per sample. readerWorkers is ignored
    # when BPPARAM is supplied; 1L avoids runPreprocess's "ignored" warning.
    args$BPPARAM <- BiocParallel::SnowParam(workers = workers, tasks = workers)
    args$readerWorkers <- 1L
  } else {
    args$readerWorkers <- workers
  }
  do.call(analyze, args)
  gds <- paste0(cls, ".gds")
  sz <- if (file.exists(gds)) file.size(gds) else NA_real_
  unlink(paste0(cls, c(".gds", ".gds.buildkey.rds")))
  invisible(sz)   # GDS size in bytes, surfaced through measure()'s `value`
}

cat(sprintf("probes: 450k   workers (fastMethyl reader): %d   sides: %s   (minfi reads serially)\n",
            workers, sides))
cat(sprintf("memory: cgroup peak/mean over baseline%s\n",
            if (hasCg) "" else " UNAVAILABLE (no cgroup memory controller)"))
cat(sprintf("cohort sizes: %s\n\n", paste(Ns, collapse = ", ")))

# Warm the annotation package once, OUTSIDE any measured window. Both pipelines
# load it, but it is a one-time, cohort-independent session cost: whichever
# pipeline/cohort runs first would otherwise absorb the cold library load into
# its time and peak (e.g. a fastMethyl-only sweep would charge it to the first
# cohort, inflating that row and breaking comparability with a both-sides run
# where minfi warms it first). Loading it up front makes every measured run warm.
suppressPackageStartupMessages(library(annoPkg, character.only = TRUE))
invisible(suppressMessages(getAnnotation(get(annoPkg))))

rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  coh <- makeCohort(n)
  up <- if (sides == "fastmethyl") naMeasure else
    measure(function() suppressMessages(minfiPipeline(coh)))
  gc(FALSE)
  fm <- if (sides == "minfi") naMeasure else
    measure(function() suppressMessages(fmPipeline(coh)))
  gc(FALSE)
  unlink(coh$dir, recursive = TRUE)
  .sweepScratch()   # reclaim any GDS a crashed cohort left behind
  speedup <- if (is.na(up$time) || is.na(fm$time)) NA_real_ else
    round(up$time / fm$time, 1)
  fm_gds <- if (is.null(fm$value) || is.na(fm$value)) NA_integer_ else
    as.integer(round(fm$value / 1048576))
  rows[[length(rows) + 1L]] <- data.frame(
    N = n,
    upstream_s = round(up$time, 1),
    fastMethyl_s = round(fm$time, 1),
    speedup = speedup,
    upstream_peak_MiB = mib(up$peak),
    fastMethyl_peak_MiB = mib(fm$peak),
    fastMethyl_gds_MiB = fm_gds,
    upstream_avg_MiB = mib(up$avg),
    fastMethyl_avg_MiB = mib(fm$avg))
}
.sweepScratch()   # final safety net, in case the last cohort left anything
cat(sprintf("\n=== IDAT -> QC'd %s GDS: time + peak memory (downstream analysis excluded) ===\n",
            if (nzchar(bpCompress)) bpCompress else "uncompressed"))
print(do.call(rbind, rows), row.names = FALSE)
