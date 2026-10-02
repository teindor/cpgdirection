# Distal links (include = "distal") and the tumour eQTM consensus
# (include_onco = TRUE): both are opt-in, both propose targets, neither moves
# best_direction. Fixture-driven and hermetic.

library(data.table)

fx_manifest <- function(...) {
  d <- data.table(...)
  d[, array := "EPICv2"]
  d[, annotation_source := "Illumina_EPICv2"]
  if (!"refgene_group" %in% names(d)) d[, refgene_group := "TSS200"]
  d
}
fx_measured_row <- function(cpg, gene, direction = 1, tissue = "blood") {
  data.table(cpg_id = cpg, target_gene = gene, tissue = tissue,
             direction = direction, tss_dist = 500)
}
fx_distal <- function(cpg, gene, source = "rE2G", biosample = "whole_blood",
                      tissue_class = "blood", score = 0.5) {
  data.table(cpg_id = cpg, target_gene = gene, source = source,
             biosample = biosample, tissue_class = tissue_class, score = score,
             element_chr = "1", element_start = 1000L, element_end = 1500L,
             element_class = "distal")
}
fx_onco <- function(cpg, gene, direction = -1, tier = "O1", n = 7L, agree = 1) {
  data.table(cpg_id = cpg, target_gene = gene, direction = direction,
             onco_tier = tier, n_cancers = n, n_positive = 0L, n_negative = n,
             sign_agreement = agree, median_r = -0.4, max_abs_r = 0.6,
             region = "DISTAL", tss_dist_approx = 40000L, cancers = "BRCA;LUAD",
             source = "fixture")
}
CG1 <- "cg10000001"
CG2 <- "cg10000002"
base_src <- list(
  manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
  lookup   = NULL, measured = NULL, smr = NULL, probe_qc = NULL)
# c() would leave two entries named "measured" and [[ ]] would find the NULL
src_with <- function(...) utils::modifyList(base_src, list(...))

test_that("distal links are off by default and never counted as annotation", {
  src <- src_with(distal = fx_distal(CG1, "FARGENE"))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = src)
  expect_false("FARGENE" %in% p$target_gene)
  expect_false("distal_sources" %in% names(p))
})

test_that("include = 'distal' proposes the far gene with provenance, no direction", {
  src <- src_with(distal = rbind(
    fx_distal(CG1, "FARGENE", biosample = "whole_blood", score = 0.4),
    fx_distal(CG1, "FARGENE", source = "scE2G", biosample = "CD14_mono", score = 0.7),
    fx_distal(CG1, "ISLETGENE", biosample = "islet_beta", tissue_class = "other")))
  p <- cpg_gene_pairs(CG1, include = c("manifest", "distal"), universal = FALSE,
                      verbose = FALSE, sources = src)
  far <- p[p$target_gene == "FARGENE", ]
  expect_equal(nrow(far), 1L)
  expect_true(far$has_distal)
  expect_equal(far$mapping_primary, "distal_link")
  expect_equal(far$distal_sources, "rE2G;scE2G")
  expect_equal(far$distal_n_biosamples, 2L)
  expect_equal(far$distal_score_max, 0.7)
  # a distal link is a target proposal: with no direction evidence the pair
  # abstains rather than inheriting a sign from the link
  expect_true(is.na(far$best_direction))
  # blood question -> islet-only link is not consulted
  expect_false("ISLETGENE" %in% p$target_gene)
  # the manifest pair is untouched and carries no distal provenance
  a <- p[p$target_gene == "GENEA", ]
  expect_false(a$has_distal)
  expect_true(is.na(a$distal_sources))
  expect_equal(attr(p, "source_counts")[["distal"]], 1L)
})

test_that("non-blood tissue consults every biosample class", {
  src <- src_with(
    distal = fx_distal(CG1, "ISLETGENE", biosample = "islet_beta", tissue_class = "other"))
  p <- cpg_gene_pairs(CG1, tissue = "solid_tissue", include = c("manifest", "distal"),
                      universal = FALSE, verbose = FALSE, sources = src)
  expect_true("ISLETGENE" %in% p$target_gene)
})

test_that("a distal target that also has measured evidence resolves normally", {
  src <- src_with(
    measured = fx_measured_row(CG1, "FARGENE", direction = -1),
    distal   = fx_distal(CG1, "FARGENE"))
  p <- cpg_gene_pairs(CG1, include = c("manifest", "measured", "distal"),
                      universal = FALSE, verbose = FALSE, sources = src)
  far <- p[p$target_gene == "FARGENE", ]
  expect_equal(far$best_evidence, "measured")
  expect_equal(far$best_direction, -1)
  expect_equal(far$mapping_primary, "measured_eQTM")   # measured outranks the link
  expect_true(far$has_distal)
})

test_that("include_onco reports tumour evidence beside, never inside, best_direction", {
  src <- src_with(
    measured = fx_measured_row(CG1, "GENEA", direction = 1),
    onco     = fx_onco(CG1, "GENEA", direction = -1))
  p <- cpg_gene_pairs(CG1, include_onco = TRUE, universal = FALSE,
                      verbose = FALSE, sources = src)
  a <- p[p$target_gene == "GENEA", ]
  expect_equal(a$best_direction, 1)          # peripheral measured sign wins
  expect_equal(a$onco_direction, -1)
  expect_equal(a$onco_tier, "O1")
  expect_equal(a$onco_n_cancers, 7L)
  expect_false(a$onco_agreement)             # tumours disagree with blood: visible
  expect_true(a$has_onco)
  # provenance: the tumour layer never outranks a peripheral source, and every
  # mapping source the function can emit has a rank (an unranked source sorts
  # first and silently becomes mapping_primary)
  expect_equal(a$mapping_primary, "measured_eQTM")
  expect_true(all(unlist(strsplit(p$mapping_sources, ";")) %in%
                    names(.CPGD_MAPPING_PRIORITY)))
  srcm <- src_with(onco = fx_onco(CG1, "GENEA"))   # manifest + tumour only
  pm <- cpg_gene_pairs(CG1, include_onco = TRUE, universal = FALSE,
                       verbose = FALSE, sources = srcm)
  expect_equal(pm[pm$target_gene == "GENEA", ]$mapping_primary, "EPICv2_manifest")

  # a gene known only to the tumour layer is discovered, flagged, and abstains
  src2 <- src_with(onco = fx_onco(CG1, "TUMGENE", direction = -1))
  p2 <- cpg_gene_pairs(CG1, include_onco = TRUE, universal = FALSE,
                       verbose = FALSE, sources = src2)
  t2 <- p2[p2$target_gene == "TUMGENE", ]
  expect_equal(nrow(t2), 1L)
  expect_equal(t2$mapping_primary, "tumour_eQTM")
  expect_true(is.na(t2$best_direction))
  expect_true(is.na(t2$onco_agreement))
})

test_that("include_onco = FALSE leaves no tumour columns", {
  src <- src_with(onco = fx_onco(CG1, "GENEA"))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = src)
  expect_false("onco_direction" %in% names(p))
  expect_false(any(p$has_onco))
})

test_that(".cpgd_distal_summary collapses per pair and handles empties", {
  D <- rbind(fx_distal(CG1, "G1", score = 0.2),
             fx_distal(CG1, "G1", biosample = "CD4_T", score = 0.9),
             fx_distal(CG2, "G2", score = NA_real_))
  s <- cpgdirection:::.cpgd_distal_summary(D, c(CG1, CG2))
  expect_equal(nrow(s), 2L)
  expect_equal(s[s$cpg_id == CG1, ]$distal_score_max, 0.9)
  expect_true(is.na(s[s$cpg_id == CG2, ]$distal_score_max))
  e <- cpgdirection:::.cpgd_distal_summary(D, "cg99999999")
  expect_equal(nrow(e), 0L)
  expect_true(all(c("distal_sources", "distal_score_max") %in% names(e)))
})
