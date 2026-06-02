# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Small minfi internal helpers that fastMethyl depends on but minfi does not
# export. They are pure (no S4 classes, no heavy dependencies) and reproduced
# here verbatim so the package needs only minfi's *public* API. `.guessArrayTypes`
# itself lives in read-meth.R; it consumes the default-annotation constants
# defined below.
#
# Keep these in sync with minfi if its array/annotation defaults ever change.

# Default annotation package suffix per array type.
.default.27k.annotation <- "ilmn12.hg19"
.default.450k.annotation <- "ilmn12.hg19"
.default.epic.annotation <- "ilm10b4.hg19"
.default.allergy.annotation <- "ilm12.hg19"
.default.epicv2.annotation <- "20a1.hg38"

# annotation tuple -> annotation package name (e.g. "...EPICanno.ilm10b4.hg19").
.getAnnotationString <- function(annotation) {
  if (length(annotation) == 1) {
    if (annotation == "Unknown") {
      stop("Cannot get Annotation object for an 'Unknown' array")
    }
    return(sprintf("%sanno", annotation))
  }
  if (all(c("array", "annotation") %in% names(annotation))) {
    if (annotation["array"] == "Unknown") {
      stop("Cannot get Annotation object for an 'Unknown' array")
    }
    return(
      sprintf("%sanno.%s", annotation["array"], annotation["annotation"])
    )
  }
  stop("unable to get the annotation string for this object")
}
