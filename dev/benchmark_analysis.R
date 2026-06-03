# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Profiler for the GDS analysis phase (bigmelon::dasen and bigmelon::outlyx)
# across cohort size N, on a synthetic 450k cohort fabricated from minfiData.
#
# History: this harness was written to benchmark prototype fastMethyl
# accelerations (a BPPARAM-parallel `dasen` and a Gram-matrix `outlyx`). The
# measurements showed both were net-negative versus the stock bigmelon GDS
# paths, so the accelerations were dropped. What
# remains is a profiler of the bigmelon baseline -- useful for understanding
# the real cost drivers before attempting any future acceleration.
#
# Measured (15 workers, 485512 probes):
#   outlyx: bigmelon 0.66s at N=96 (deterministic); the cost is gdsfmt's
#     C-level apply.gdsn over probes, not the SVD.
#   dasen:  bigmelon 66.5s (N=96), 127.9s (N=192); dominated by quantile
#     normalisation + I/O, not the per-sample density fits.
#
# Run (always nice/ionice it; the synthetic matrices are large):
#   nice -n 19 ionice -c 3 Rscript dev/benchmark_analysis.R
# Env knobs: BENCH_N (default "48,192"), BENCH_SCRATCH (default tempdir()).

suppressPackageStartupMessages({
  library(minfi); library(minfiData); library(bigmelon); library(gdsfmt)
})

Ns <- as.integer(strsplit(Sys.getenv("BENCH_N", "48,192"), ",")[[1]])
scratch <- Sys.getenv("BENCH_SCRATCH", unset = tempdir())
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
tt <- function(expr) unname(system.time(expr)[3])

message("Reading minfiData base (6 samples) ...")
base_mset <- preprocessRaw(
  minfi::read.metharray.exp(system.file("extdata", package = "minfiData"),
                            recursive = TRUE)
)
baseM <- getMeth(base_mset); baseU <- getUnmeth(base_mset)
ann <- annotation(base_mset); loci <- rownames(base_mset); nb <- ncol(base_mset)

# Fabricate an N-sample bigmelon GDS; es2gds returns it open.
synth_gds <- function(n) {
  idx <- sample.int(nb, n, replace = TRUE)
  jit <- function(x) round(x * matrix(runif(length(x), 0.85, 1.15), nrow(x)))
  M <- jit(baseM[, idx, drop = FALSE]); U <- jit(baseU[, idx, drop = FALSE])
  storage.mode(M) <- "integer"; storage.mode(U) <- "integer"
  bc <- sprintf("%010d_R%02dC%02d", 2e9 + seq_len(n),
                (seq_len(n) - 1L) %% 8L + 1L, (seq_len(n) - 1L) %% 6L + 1L)
  ms <- MethylSet(Meth = M, Unmeth = U, annotation = ann)
  rownames(ms) <- loci; colnames(ms) <- bc
  f <- tempfile(tmpdir = scratch, fileext = ".gds")
  g <- NULL
  invisible(utils::capture.output(g <- suppressMessages(es2gds(ms, f))))
  list(g = g, f = f)
}

cat(sprintf("\nprobes: %d   sample counts: %s\n\n", length(loci),
            paste(Ns, collapse = ", ")))
rows <- list()
for (n in Ns) {
  message(sprintf("=== N = %d ===", n))
  sg <- synth_gds(n); g <- sg$g
  t_dasen  <- tt(suppressMessages(dasen(g, node = "nb")))
  t_outlyx <- tt(suppressMessages(outlyx(g, plot = FALSE)))
  closefn.gds(g); unlink(sg$f)
  rows[[length(rows) + 1L]] <- data.frame(
    N = n, dasen_bigmelon = round(t_dasen, 2),
    outlyx_bigmelon = round(t_outlyx, 2)
  )
}
cat("\n=== bigmelon analysis-phase timings (seconds elapsed) ===\n")
print(do.call(rbind, rows), row.names = FALSE)
