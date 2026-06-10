# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("fastMethyl")
  msg <- c(
    sprintf("fastMethyl %s", v),
    "  analyze(): one call -- build a QC'd, bigmelon-compatible GDS, then run",
    "    outlierRemoval -> normalization -> your analysis, GDS closed for you.",
    "    Hooks take a method name (\"outlyx\"/\"dasen\"), a function, or \"none\".",
    "    See system.file('scripts', 'pipeline.R', package = 'fastMethyl')",
    "    for a worked example.",
    "  Built on minfi + bigmelon."
  )
  packageStartupMessage(paste(msg, collapse = "\n"))
}
