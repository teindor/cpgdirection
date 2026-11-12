# The measured layer after the eQTM Atlas rebuild (2.99.8): cohort provenance
# columns, the replication-based accuracy string, and the pair resolver carrying
# both through.

test_that("the measured table carries cohort provenance and no tumour class", {
  skip_if(!cpgd_has_data("measured_eqtms"), "measured table unavailable")
  M <- cpgd_measured_eqtms()
  expect_true(all(c("cpg_id", "target_gene", "tissue", "direction", "tss_dist") %in% names(M)))
  expect_true(all(M$direction %in% c(-1, 1)))
  expect_true(all(abs(M$tss_dist) <= 250000, na.rm = TRUE))
  expect_false("tumour" %in% M$tissue)
  skip_if(!"n_cohorts" %in% names(M), "measured table predates the 2.99.8 rebuild")
  expect_true(all(c("n_cohorts", "cohort_agreement", "cohorts", "source") %in% names(M)))
  expect_true(all(M$n_cohorts >= 1L))
  expect_true(all(M$cohort_agreement > 0.5 & M$cohort_agreement <= 1))
  expect_true(any(M$n_cohorts >= 2L))
  expect_true(all(M$source %in% c("eqtm_atlas", "packaged_v1", "helix_dryad_2022", "wang2021_S2",
                                  "kreitmaier2026_msk", "liu2022_kidney", "delahaye2018_placenta",
                                  "voisin2016_adipose", "kennedy2018_gtp")))
})

test_that("kidney, placenta, adipose and GTP blood are single-cohort donor groups in their tissue", {
  co <- c("Liu2022_kidney", "Delahaye2018_placenta", "Voisin2016_adipose", "Kennedy2018_GTP_blood")
  expect_equal(cpgdirection:::.cpgd_donor_group(co), co)
  skip_if(!cpgd_has_data("measured_eqtms"), "measured table unavailable")
  M <- cpgd_measured_eqtms()
  skip_if(!"liu2022_kidney" %in% M$source, "measured table predates the kidney source")
  has <- function(x) M[grepl(x, M$cohorts, fixed = TRUE), ]
  expect_true(all(has("Liu2022_kidney")$tissue == "solid_tissue"))
  expect_true(all(has("Delahaye2018_placenta")$tissue == "solid_tissue"))
  expect_true(all(has("Voisin2016_adipose")$tissue == "solid_tissue"))
  expect_true(all(has("Kennedy2018_GTP_blood")$tissue == "blood"))
  expect_gt(nrow(has("Liu2022_kidney")), 29000L)
  expect_gt(nrow(has("Delahaye2018_placenta")), 2000L)
  expect_gt(nrow(has("Voisin2016_adipose")), 350L)
  expect_gt(nrow(has("Kennedy2018_GTP_blood")), 3500L)
  # no pair of the table is left split between cohorts
  expect_true(all(M$n_cohorts == lengths(strsplit(M$cohorts, "|", fixed = TRUE)) | M$cohorts == "packaged_v1"))
  # a kidney pair also measured, with the same sign, in another donor group is replicated; a kidney-only pair is not
  K <- has("Liu2022_kidney"); K <- K[!grepl(";", K$target_gene), ]
  only <- K[K$cohorts == "Liu2022_kidney", ][1, ]
  both <- K[grepl("FUSION|GTEx", K$cohorts) & K$cohort_agreement >= 0.99, ][1, ]
  p <- cpg_gene_pairs(c(only$cpg_id, both$cpg_id), genes = c(only$target_gene, both$target_gene),
                      gene_mode = "pairwise", probe_qc = "flag", tissue = "solid_tissue", verbose = FALSE)
  d <- as.data.frame(p)
  r1 <- d[d$cpg_id == only$cpg_id & d$target_gene == only$target_gene, ]
  r2 <- d[d$cpg_id == both$cpg_id & d$target_gene == both$target_gene, ]
  expect_equal(r1$measured_n_cohorts, 1L)
  expect_equal(r1$best_rung, "measured_solid_single")
  expect_equal(r1$best_direction, only$direction)
  expect_equal(r2$best_rung, "measured_solid_replicated")
})

test_that("the osteoarthritis joint tissues are one donor group in solid tissue", {
  g <- cpgdirection:::.cpgd_donor_group(c("Kreitmaier2026_cartilage", "Kreitmaier2026_synovium", "FUSION",
                                           "GTEx_Lung", "MESA_Tcells"))
  expect_equal(g, c("Kreitmaier2026", "Kreitmaier2026", "FUSION", "GTEx", "MESA"))
  skip_if(!cpgd_has_data("measured_eqtms"), "measured table unavailable")
  M <- cpgd_measured_eqtms()
  skip_if(!"kreitmaier2026_msk" %in% M$source, "measured table predates the osteoarthritis source")
  K <- M[grepl("Kreitmaier2026", M$cohorts, fixed = TRUE), ]
  expect_true(all(K$tissue == "solid_tissue"))
  expect_gt(nrow(K), 15000L)
  # a pair seen only in cartilage and synovium is measured in one study, not replicated
  only <- K[K$cohorts == "Kreitmaier2026_cartilage|Kreitmaier2026_synovium" & K$cohort_agreement >= 0.99 &
              !grepl(";", K$target_gene), ][1, ]
  both <- K[grepl("FUSION|GTEx", K$cohorts) & K$cohort_agreement >= 0.99 & !grepl(";", K$target_gene), ][1, ]
  p <- cpg_gene_pairs(c(only$cpg_id, both$cpg_id), genes = c(only$target_gene, both$target_gene),
                      gene_mode = "pairwise", probe_qc = "flag", tissue = "solid_tissue", verbose = FALSE)
  d <- as.data.frame(p)
  r1 <- d[d$cpg_id == only$cpg_id & d$target_gene == only$target_gene, ]
  r2 <- d[d$cpg_id == both$cpg_id & d$target_gene == both$target_gene, ]
  expect_equal(r1$measured_n_cohorts, 2L)
  expect_false(grepl("replicated", r1$best_rung))
  expect_equal(r2$best_rung, "measured_solid_replicated")
})

test_that("a measured pair reports its cohort count and a replication-based accuracy", {
  skip_if(!cpgd_has_data("measured_eqtms"), "measured table unavailable")
  M <- cpgd_measured_eqtms()
  skip_if(!"n_cohorts" %in% names(M), "measured table predates the 2.99.8 rebuild")
  # blood pairs measured in exactly one table row and without an SMR sign, so
  # the measured rung answers in the blood order
  key <- paste(M$cpg_id, M$target_gene)
  S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
  smr_key <- if (is.null(S)) character(0) else paste(S$cpg_id, S$target_gene)
  cand <- M[key %in% names(which(table(key) == 1L)) & !grepl(";", M$target_gene) &
              M$tissue == "blood" & !(key %in% smr_key), ]
  ngrp <- vapply(strsplit(cand$cohorts, "|", fixed = TRUE), function(x)
    length(unique(cpgdirection:::.cpgd_donor_group(x))), integer(1))
  rep <- cand[ngrp >= 2L & cand$cohort_agreement >= 0.99, ][1, ]
  one <- cand[cand$n_cohorts == 1L, ][1, ]
  # two cohorts of one study (MESA's two cell types) are one donor group
  same <- cand[cand$cohorts == "MESA_monocytes|MESA_Tcells" & cand$cohort_agreement >= 0.99, ][1, ]
  p <- cpg_gene_pairs(c(rep$cpg_id, one$cpg_id, same$cpg_id),
                      genes = c(rep$target_gene, one$target_gene, same$target_gene),
                      gene_mode = "pairwise", probe_qc = "flag", verbose = FALSE)
  d <- as.data.frame(p)
  r1 <- d[d$cpg_id == rep$cpg_id & d$target_gene == rep$target_gene, ]
  r2 <- d[d$cpg_id == one$cpg_id & d$target_gene == one$target_gene, ]
  r3 <- d[d$cpg_id == same$cpg_id & d$target_gene == same$target_gene, ]
  expect_equal(nrow(r1), 1L); expect_equal(nrow(r2), 1L); expect_equal(nrow(r3), 1L)
  expect_equal(r1$best_rung, "measured_blood_replicated")
  expect_gte(r1$measured_n_cohorts, 2L)
  expect_match(r1$best_expected_accuracy, "measured in blood, 2\\+ independent cohorts agree")
  expect_equal(r2$best_rung, "measured_blood_single")
  expect_equal(r2$measured_n_cohorts, 1L)
  expect_match(r2$best_expected_accuracy, "measured in blood, one study")
  expect_equal(r3$measured_n_cohorts, 2L)
  expect_equal(r3$best_rung, "measured_blood_single")
})

test_that("an older measured table (no cohort columns) still resolves as one cohort", {
  M <- data.table::data.table(cpg_id = "cg00000001", target_gene = "GENEX", tissue = "blood",
                              direction = -1, tss_dist = 100)
  p <- cpg_gene_pairs("cg00000001", genes = "GENEX", gene_mode = "pairwise",
                      sources = list(measured = M, smr = NULL, lookup = NULL, onco = NULL),
                      verbose = FALSE)
  d <- as.data.frame(p)
  expect_equal(d$best_evidence, "measured")
  expect_equal(d$measured_n_cohorts, 1L)
  expect_equal(d$measured_cohort_agreement, 1)
  expect_match(d$best_expected_accuracy, "measured in blood, one study")
})
