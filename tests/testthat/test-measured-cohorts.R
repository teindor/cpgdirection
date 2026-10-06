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
  expect_true(all(M$source %in% c("eqtm_atlas", "packaged_v1")))
})

test_that("a measured pair reports its cohort count and a replication-based accuracy", {
  skip_if(!cpgd_has_data("measured_eqtms"), "measured table unavailable")
  M <- cpgd_measured_eqtms()
  skip_if(!"n_cohorts" %in% names(M), "measured table predates the 2.99.8 rebuild")
  # pairs measured in exactly one tissue row, so no cross-tissue sign conflict
  key <- paste(M$cpg_id, M$target_gene)
  cand <- M[key %in% names(which(table(key) == 1L)) & !grepl(";", M$target_gene), ]
  rep <- cand[cand$n_cohorts >= 2L & cand$cohort_agreement >= 0.99, ][1, ]
  one <- cand[cand$n_cohorts == 1L, ][1, ]
  p <- cpg_gene_pairs(c(rep$cpg_id, one$cpg_id), genes = c(rep$target_gene, one$target_gene),
                      gene_mode = "pairwise", probe_qc = "flag", verbose = FALSE)
  d <- as.data.frame(p)
  r1 <- d[d$cpg_id == rep$cpg_id & d$target_gene == rep$target_gene, ]
  r2 <- d[d$cpg_id == one$cpg_id & d$target_gene == one$target_gene, ]
  expect_equal(nrow(r1), 1L); expect_equal(nrow(r2), 1L)
  expect_equal(r1$best_evidence, "measured")
  expect_gte(r1$measured_n_cohorts, 2L)
  expect_match(r1$best_expected_accuracy, "2\\+ independent cohorts")
  expect_equal(r2$best_evidence, "measured")
  expect_equal(r2$measured_n_cohorts, 1L)
  expect_match(r2$best_expected_accuracy, "^0\\.96-0\\.99")
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
  expect_match(d$best_expected_accuracy, "^0\\.96-0\\.99")
})
