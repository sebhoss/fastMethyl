# SPDX-FileCopyrightText: 2026 The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Expose fastMethyl's non-exported helpers by bare name so the test files can
# call them directly (the package is loaded by tests/testthat.R, so its
# namespace is always available here).

.fm_internals <- c(
  ".normalize_verbose",
  ".guessArrayTypes",
  ".alignProbeAddresses",
  ".streamingResolveManifest",
  ".streamingBuildIndices",
  ".streamingChunkKeepSamples",
  ".streamingChunkProbeFails",
  ".streamingAnnoKeep",
  ".streamingQCfinalize",
  ".bytes_human",
  ".bytes_human_count",
  ".mem_report",
  ".path_size",
  ".bpmapply_try",
  ".bpmapply_with_progress",
  ".default.27k.annotation",
  ".default.450k.annotation",
  ".default.epic.annotation",
  ".default.allergy.annotation",
  ".default.epicv2.annotation",
  ".getAnnotationString",
  ".validateArgs",
  ".inspectExistingGDS",
  ".buildKeyPath",
  ".currentBuildKey",
  ".compareBuildKey",
  ".withGds",
  ".dasenNormalize",
  "processMethArrayExp"
)

.fm_ns <- asNamespace("fastMethyl")
for (.fm_n in .fm_internals) {
  if (exists(.fm_n, envir = .fm_ns, inherits = FALSE)) {
    assign(.fm_n, get(.fm_n, envir = .fm_ns), envir = environment())
  }
}
rm(.fm_ns, .fm_internals)
if (exists(".fm_n", inherits = FALSE)) rm(.fm_n)
