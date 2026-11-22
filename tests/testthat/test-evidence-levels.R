# The duplicated-level bug that broke 2.1.0 came from this list being written
# out inline in the print method and hand-edited each time a layer was added.
# It now lives in one constant; these guard that it stays one.

test_that("the evidence level constant has no duplicates", {
    expect_false(as.logical(anyDuplicated(cpgdirection:::CPGD_EVIDENCE)))
    expect_silent(factor("measured", levels = cpgdirection:::CPGD_EVIDENCE))
})

test_that("every layer the ladder can emit is a known level", {
    expect_true(all(c(
        "measured", "smr_high", "catalogue_consensus", "smr_moderate",
        "catalogue_single", "smr_weak", "distance_only",
        "distance_targeted", "catalogue_targeted",
        "targeted_last_resort", "distance_tissue_conflict",
        "distance_uninformative", "tissue_conflict",
        "no_evidence"
    ) %in% cpgdirection:::CPGD_EVIDENCE))
})

test_that("the levels are ordered strongest first", {
    lv <- cpgdirection:::CPGD_EVIDENCE
    expect_equal(lv[1], "measured")
    expect_lt(match("smr_high", lv), match("catalogue_consensus", lv))
    expect_lt(match("catalogue_consensus", lv), match("smr_moderate", lv))
    expect_lt(match("smr_moderate", lv), match("catalogue_single", lv))
    expect_equal(lv[length(lv)], "no_evidence")
})
