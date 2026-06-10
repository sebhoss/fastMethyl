#!/usr/bin/env bash
# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Runs a fastMethyl benchmark inside a `systemd-run --user --scope` cgroup v2
# scope whose CPU / memory / IO limits are sized to the cohort. This is the ONLY
# supported way to run a benchmark: the limits live here, set per call, rather
# than in .ilo.rc -- so they can be tuned without editing a committed file, and
# every run is bounded and reproducible.
#
# BENCH_SCRIPT selects which benchmark to run (default the preprocessing one):
#   dev/benchmark_analysis.R  IDAT -> QC'd GDS, fastMethyl vs upstream minfi
#   dev/benchmark_pca.R       analysis-phase PCA routes on a normalised GDS
# The same envelope bounds both; the PCA benchmark is far lighter (it reuses one
# built GDS), so the envelope sized for the preprocessing run is ample for it.
#
# The rootless podman container nests under the scope (via --cgroups=split in
# .ilo.rc), so the caps reach the R master and its forked reader workers.
#
# Limits are derived from BENCH_MAX_N -- a FIXED reference cohort, not the sizes
# in this run -- so every sample size is bounded by the identical envelope and
# cross-size numbers stay comparable (different limits would confound any
# scaling claim). The pieces:
#   * CPUQuota   -- fixed core count (default 4); BENCH_WORKERS matches it.
#   * MemoryHigh -- throttle/reclaim just above the expected peak; this is also
#                   what makes the IO write cap actually bite (it forces dirty
#                   pages to be written back in the cgroup's own reclaim context,
#                   where the IO is attributed to the cgroup).
#   * MemoryMax  -- hard ceiling with margin, so a run up to N=200 cannot OOM.
#   * IOWriteBandwidthMax -- caps the GDS write storm; the device is the whole
#                   disk backing the scratch dir (io.max keys on the disk's
#                   request queue, not a partition), given as the real device
#                   node (a btrfs path stat()s to a synthetic dev_t systemd
#                   cannot map).
#
# Memory model (measured on 450k, net of the container baseline), evaluated at
# BENCH_MAX_N and CPU_CORES:
#   upstream peak   ~= 678 + 51 * N                MiB  (single-threaded; N=6 ->
#                                                        984, N=50 -> 3228)
#   fastMethyl peak ~= 1200 + 6 * N + 120 * cores  MiB  (block + per-worker
#                                                        buffers scale with cores;
#                                                        4c/N50 -> 1606, 8c -> 2087)
#   The envelope is sized to the LARGER of the two (so it holds at any core
#   count), plus ~2 GiB container/R/annotation resident set; MemoryHigh = 1.05x,
#   MemoryMax = 1.30x.
#
# Override anything via env: BENCH_N, BENCH_MAX_N, BENCH_WORKERS, CPU_CORES,
# IO_RATE, BENCH_SCRATCH, IO_DEVICE. Keep BENCH_MAX_N identical across every run
# you intend to compare.
#
#   dev/run-benchmark.sh                 # default cohorts 50,100,200
#   BENCH_N=50 dev/run-benchmark.sh      # one cohort
set -euo pipefail

cd "$(dirname "$0")/.."

BENCH_N="${BENCH_N:-50,100,200}"
# The resource envelope is sized from a FIXED reference cohort, NOT from the
# cohorts in this particular run, so every sample size -- whether run alone or
# together -- is bounded by the identical CPU/memory/IO limits. Comparing numbers
# across sample sizes (e.g. "does fastMethyl scale linearly?") is only valid if
# the limits do not change between them. Keep BENCH_MAX_N the same for every run
# whose results you intend to compare; it must be >= the largest cohort.
BENCH_MAX_N="${BENCH_MAX_N:-200}"
CPU_CORES="${CPU_CORES:-4}"
BENCH_WORKERS="${BENCH_WORKERS:-$CPU_CORES}"   # == cores; CPU is capped externally
BENCH_SIDES="${BENCH_SIDES:-both}"             # "both" or "fastmethyl" (minfi is
                                               # core-independent, so a CPU sweep
                                               # need only re-run fastMethyl)
BENCH_BPPARAM="${BENCH_BPPARAM:-multicore}"    # fastMethyl reader backend:
                                               # "multicore" (fork) or "snow"
# `-` (not `:-`): an explicitly empty BENCH_COMPRESS means "no compression" and
# must survive; only a genuinely unset var falls back to the LZ4_RA default.
BENCH_COMPRESS="${BENCH_COMPRESS-LZ4_RA}"       # GDS codec: "LZ4_RA" or "" (none)
BENCH_SCRIPT="${BENCH_SCRIPT:-dev/benchmark_analysis.R}"  # which benchmark to run
BENCH_PCA_NCP="${BENCH_PCA_NCP:-10}"            # components for the PCA benchmark
IO_RATE="${IO_RATE:-120M}"
scratch="${BENCH_SCRATCH:-$PWD/dev/scratch}"

# Validate the requested cohorts and ensure none exceeds the reference envelope.
maxN=0
IFS=',' read -ra _ns <<<"$BENCH_N"
for n in "${_ns[@]}"; do
  n="${n//[[:space:]]/}"
  [ -n "$n" ] && [ "$n" -gt "$maxN" ] && maxN="$n"
done
[ "$maxN" -gt 0 ] || { echo "no valid cohort sizes in BENCH_N='$BENCH_N'" >&2; exit 1; }
if [ "$maxN" -gt "$BENCH_MAX_N" ]; then
  echo "BENCH_N has N=$maxN > BENCH_MAX_N=$BENCH_MAX_N; raise BENCH_MAX_N and" \
       "re-run EVERY cohort so the envelope stays identical across sizes." >&2
  exit 1
fi

ceil_gib() { echo $(( ( $1 + 1023 ) / 1024 )); }
# Two peaks bound the run, both evaluated at BENCH_MAX_N: the single-threaded
# upstream pipeline (grows ~linearly in N, independent of cores) and fastMethyl
# (its streaming block + per-worker buffers grow with CPU_CORES, plus ~linear GDS
# write-cache in N). Size to whichever is larger so the envelope holds at any
# core count. At BENCH_MAX_N=200 upstream dominates; the fastMethyl term binds
# only for a small reference cohort run with many cores.
up_mib=$(( 678 + 51 * BENCH_MAX_N ))
fm_mib=$(( 1200 + 6 * BENCH_MAX_N + 120 * CPU_CORES ))
peak_mib=$(( up_mib > fm_mib ? up_mib : fm_mib ))
abs_mib=$(( 2048 + peak_mib ))
mem_high=$(ceil_gib $(( abs_mib * 105 / 100 )))
mem_max=$(ceil_gib $(( abs_mib * 130 / 100 )))
[ "$mem_high" -lt 2 ] && mem_high=2
[ "$mem_max" -le "$mem_high" ] && mem_max=$(( mem_high + 1 ))
# Explicit envelope override (whole GiB). The auto-size assumes the minfi side
# runs (its peak is ~51*N MiB); a fastMethyl-only large-N sweep needs only a few
# GiB, so the computed cap can exceed host RAM. Set MEM_MAX/MEM_HIGH to pin a
# host-safe ceiling -- the cgroup then OOM-kills an over-budget run instead of
# wedging the host.
[ -n "${MEM_HIGH:-}" ] && mem_high="$MEM_HIGH"
[ -n "${MEM_MAX:-}" ] && mem_max="$MEM_MAX"

# Whole disk backing the scratch dir (io.max keys on the disk request queue).
if [ -n "${IO_DEVICE:-}" ]; then
  iodev="$IO_DEVICE"
else
  src=$(findmnt -no SOURCE -T "$scratch" | sed 's/\[.*//')
  disk=$(lsblk -no PKNAME "$src" 2>/dev/null | head -1)
  iodev="/dev/${disk:-$(basename "$src")}"
fi

echo "cohorts $BENCH_N | envelope sized for N=$BENCH_MAX_N (fixed) | workers $BENCH_WORKERS | CPUQuota ${CPU_CORES}00%"
echo "MemoryHigh=${mem_high}G MemoryMax=${mem_max}G | IOWriteBandwidthMax='$iodev $IO_RATE'"

# Delegate=yes pushes the cpu/io/memory controllers down into the container's own
# cgroup (with --cgroups=split the leaf otherwise sees only `pids`), so the
# harness can read memory.current from inside to report peak RAM. The scope still
# enforces every cap hierarchically.
exec systemd-run --user --scope --unit="fm-bench-$$" -p Delegate=yes \
  -p "CPUQuota=${CPU_CORES}00%" \
  -p "MemoryHigh=${mem_high}G" \
  -p "MemoryMax=${mem_max}G" \
  -p "IOWriteBandwidthMax=$iodev $IO_RATE" \
  ilo bash -c "env BENCH_N='$BENCH_N' BENCH_WORKERS='$BENCH_WORKERS' BENCH_SIDES='$BENCH_SIDES' BENCH_BPPARAM='$BENCH_BPPARAM' BENCH_COMPRESS='$BENCH_COMPRESS' BENCH_PCA_NCP='$BENCH_PCA_NCP' R_LIBS=dev/scratch/rlib:/usr/local/lib/R/site-library Rscript $BENCH_SCRIPT"
