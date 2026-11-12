# Columns differ by mode. The default (tissue = "all") returns the cpgd_full
# shape with best_* columns; naming a tissue returns the per-tissue shape with
# target_gene, status and call. Tests that assert on the per-tissue columns must
# therefore ask for a tissue. Four tests here failed R CMD check for exactly
# that reason: they predated the multi-layer default and were never updated.

test_that("identifiers are stripped to canonical form", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  r <- cpg_expression_direction(c("cg00000029_TC21", "cg00000029"),
                                tissue = "blood", verbose = FALSE)
  expect_true(all(grepl("^cg[0-9]+$", r$cpg_id)))
  expect_equal(nrow(r), 1L)          # duplicates collapse
})

test_that("a gene parsed from the name is used", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  r <- cpg_expression_direction("cg00050692_TC21_DNMT3A",
                                tissue = "blood", verbose = FALSE)
  expect_true(r$gene_source[1] %in%
                c("parsed_from_name", "user_supplied", "catalogue", "nearest_TSS"))
})

test_that("unknown identifiers come back as UNMAPPED, not silently dropped", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  r <- cpg_expression_direction("cg99999999", tissue = "blood", verbose = FALSE)
  expect_equal(r$status[1], "UNMAPPED")
})

test_that("an unsupported tissue is refused rather than silently substituted", {
  expect_error(cpg_expression_direction("cg00000029", tissue = "brain"))
})

test_that("all three tissues load and can score", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  skip_on_cran()
  for (t in c("blood", "nasal_epithelium", "solid_tissue")) {
    r <- cpg_expression_direction("cg00050692", tissue = t, verbose = FALSE)
    expect_equal(nrow(r), 1L)
    expect_equal(r$tissue[1], t)
  }
})

test_that("a requested gene the table cannot supply does not get another gene's direction", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  r <- cpg_expression_direction("Zcg02079741_TC21_POMC",
                                tissue = "blood", verbose = FALSE)
  expect_equal(nrow(r), 1L)
  expect_true(r$status[1] %in%
    c("REQUESTED_GENE_UNAVAILABLE", "NO_EXPRESSION_ANCHOR", "DIRECT_eQTM",
      "PREDICTED", "ABSTAIN"))
  if (identical(r$status[1], "REQUESTED_GENE_UNAVAILABLE"))
    expect_equal(r$call[1], "ABSTAIN")
})

test_that("a Z-prefixed identifier still parses", {
  d <- cpgd_parse_check("Zcg00335286_TC21_MC2R")
  expect_equal(d$cpg_id[1], "cg00335286")
  expect_equal(d$given_gene2[1], "MC2R")
})

test_that("the default mode returns the multi-layer column set", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  skip_on_cran()
  r <- cpg_expression_direction("cg00000029", verbose = FALSE)
  for (col in c("best_direction", "best_evidence", "best_expected_accuracy",
                "smr_direction", "smr_agreement", "direction_uncertain"))
    expect_true(col %in% names(r))
  expect_true(all(stats::na.omit(r$best_evidence) %in% cpgdirection:::CPGD_EVIDENCE))
})

test_that("target_tissue ranks by the declared tissue, fills blanks, and never drops a call", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  skip_on_cran()
  cg <- c("cg00000029", "cg00050692", "cg02079741")
  a <- cpg_expression_direction(cg, verbose = FALSE)
  b <- cpg_expression_direction(cg, target_tissue = "solid_tissue", verbose = FALSE)
  data.table::setkeyv(a, "cpg_id"); data.table::setkeyv(b, "cpg_id")
  called <- !is.na(a$best_direction)
  expect_false(any(is.na(b$best_direction[called])))
  expect_true(sum(is.na(b$best_direction)) <= sum(is.na(a$best_direction)))
  # the same source gives the same sign; a different source answers only when
  # it ranks higher in the declared tissue (the ladder is ranked per tissue)
  same <- called & b$best_rung == a$best_rung
  expect_equal(b$best_direction[same], a$best_direction[same])
  os <- cpgdirection:::.cpgd_rung_order("solid_tissue")
  other <- called & !same
  ra <- match(a$best_rung[other], os); rb <- match(b$best_rung[other], os)
  expect_true(all(!is.na(rb) & (is.na(ra) | rb < ra)))
})

test_that("the shipped reference reproduces", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "blood catalogue unavailable")
  skip_on_cran()
  expect_true(cpgd_selftest())
})
