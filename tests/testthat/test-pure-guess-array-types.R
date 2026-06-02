# .guessArrayTypes() maps a raw IDAT probe count to an (array, annotation)
# tuple. Every supported array type has a probe-count band; outside any
# band the function returns "Unknown". This test pins one probe count per
# band so a future re-shuffle of the band boundaries can't silently
# orphan an array type that the pipeline used to handle.

test_that(".guessArrayTypes resolves the 450k band", {
    out <- .guessArrayTypes(622399L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylation450k")
})

test_that(".guessArrayTypes resolves the current-EPIC band", {
    out <- .guessArrayTypes(1052641L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylationEPIC")
})

test_that(".guessArrayTypes resolves the early-access EPIC band", {
    out <- .guessArrayTypes(1032641L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylationEPIC")
})

test_that(".guessArrayTypes resolves the EPIC v2 band", {
    out <- .guessArrayTypes(1105209L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylationEPICv2")
})

test_that(".guessArrayTypes resolves the 27k band", {
    out <- .guessArrayTypes(55300L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylation27k")
})

test_that(".guessArrayTypes resolves the 54.7k Allergy band", {
    out <- .guessArrayTypes(54750L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylationAllergy")
})

test_that(".guessArrayTypes resolves the HorvathMammalMethylChip40 band", {
    out <- .guessArrayTypes(41050L)
    expect_equal(unname(out[["array"]]), "HorvathMammalMethylChip40")
})

test_that(".guessArrayTypes resolves the 43.6k Allergy band", {
    out <- .guessArrayTypes(43665L)
    expect_equal(unname(out[["array"]]), "IlluminaHumanMethylationAllergy")
})

test_that(".guessArrayTypes falls back to Unknown for out-of-band counts", {
    out <- .guessArrayTypes(9999L)
    expect_equal(unname(out[["array"]]),      "Unknown")
    expect_equal(unname(out[["annotation"]]), "Unknown")
})

test_that(".streamingResolveManifest refuses to guess for an unknown probe count", {
    expect_error(.streamingResolveManifest(9999L),
                 regexp = "could not infer array type")
})
