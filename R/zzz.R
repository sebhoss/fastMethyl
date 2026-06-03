# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion("fastMethyl")
  msg <- c(
    sprintf("fastMethyl %s", v),
    "  analyze(): one call -- preprocess (config -> QC'd, LZ4_RA-compressed",
    "    GDS) -> your analysis function on the open GDS -> closed for you.",
    "    Omit FUN to just build the GDS; normalise/outlier-detect inside FUN.",
    "    See system.file('scripts', 'pipeline.R', package = 'fastMethyl')",
    "    for a worked example.",
    "  Built on minfi + bigmelon."
  )
  packageStartupMessage(paste(msg, collapse = "\n"))
}
