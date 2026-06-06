# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# End-to-end installation test. Builds the source tarball, installs it into a
# *fresh, empty* library, and loads the package in a clean R process -- the
# closest local reproduction of what a user gets from
# `BiocManager::install("sebhoss/fastMethyl")`. Run inside the dev container:
#
#     ilo bash dev/install_test.sh
#
# It catches the install-breaking failures the testthat suite cannot, because
# the suite runs against an already-installed package: a NAMESPACE import of a
# symbol the resolved minfi does not export, a dependency used but missing from
# DESCRIPTION, or a file referenced by the package but excluded by .Rbuildignore.
# Dependencies (minfi, gdsfmt, BiocParallel, ...) are taken from the existing
# library; only fastMethyl itself is installed fresh, isolating its own
# installability from the ~30-package Bioconductor stack.
set -euo pipefail

cd "$(dirname "$0")/.."

LIB="$(mktemp -d "${BENCH_SCRATCH:-/tmp}/fm_install_lib_XXXXXX")"
BUILD_DIR="$(mktemp -d "${BENCH_SCRATCH:-/tmp}/fm_install_build_XXXXXX")"
trap 'rm -rf "$LIB" "$BUILD_DIR"' EXIT

echo ">> building source tarball"
R CMD build --no-build-vignettes --no-manual .
TARBALL="$(ls -1t fastMethyl_*.tar.gz | head -n 1)"
echo ">> built ${TARBALL}"

# Internal-only paths must not reach the shipped tarball (.Rbuildignore). A
# tarball that carries dev/ or CLAUDE.md would still install, so assert it here.
echo ">> checking tarball excludes internal-only files"
LEAKED="$(tar tzf "$TARBALL" | grep -E 'fastMethyl/(dev/|CLAUDE\.md|\.lintr|REUSE\.toml)' || true)"
if [ -n "$LEAKED" ]; then
    echo "FAIL: internal-only files leaked into the tarball:" >&2
    echo "$LEAKED" >&2
    exit 1
fi

echo ">> installing into a fresh empty library: ${LIB}"
R CMD INSTALL --library="$LIB" "$TARBALL"

echo ">> loading from the fresh library in a clean R process"
Rscript -e "
  .libPaths(c('${LIB}', .libPaths()))
  suppressPackageStartupMessages(ok <- requireNamespace('fastMethyl', quietly = TRUE))
  stopifnot('fastMethyl namespace loads' = ok)
  library(fastMethyl)
  stopifnot(
    'analyze() is exported'      = 'analyze' %in% getNamespaceExports('fastMethyl'),
    'analyze() is a function'    = is.function(fastMethyl::analyze),
    'minfi loaded transitively'  = 'minfi' %in% loadedNamespaces()
  )
  cat('>> loaded fastMethyl', as.character(packageVersion('fastMethyl')),
      'against minfi', as.character(packageVersion('minfi')), '\n')
"

echo ">> PASS: fastMethyl builds, installs into a clean library, and loads"
mv -f "$TARBALL" "$BUILD_DIR/" 2>/dev/null || true
