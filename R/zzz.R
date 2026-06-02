# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("fastMethyl")
  msg <- c(
    sprintf("fastMethyl %s", v),
    "  runPreprocess(): config -> QC'd, LZ4_RA-compressed GDS.",
    "    See system.file('scripts', 'pipeline.R', package = 'fastMethyl')",
    "    for a worked example.",
    "  processMethArray() / processMethArrayExp(): fused, parallel",
    "    read + detectionP + preprocessRaw streamed to a GDS in one pass.",
    "  readMethArray() / readMethArrayExp(): minfi-compatible readers",
    "    with a BPPARAM= argument for parallel IDAT reading.",
    "  Built on minfi + bigmelon."
  )
  packageStartupMessage(paste(msg, collapse = "\n"))
}
