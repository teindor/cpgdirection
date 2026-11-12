# Reported tissue contexts and region flags (1.1.0): information only. They
# never vote and never supply or change best_direction.

test_that("the reported contexts table is well formed and disjoint from the voting table's tissues", {
  X <- cpgd_reported_contexts()
  expect_true(all(c("cpg_id", "target_gene", "context", "direction", "tss_dist", "n_cohorts",
                    "cohort_agreement", "cohorts", "source", "evidence", "flags") %in% names(X)))
  expect_setequal(unique(X$context), c("fibroblasts", "liver", "pancreatic_islets", "retina", "saliva"))
  expect_true(all(X$direction %in% c(-1, 1)))
  expect_true(all(abs(X$tss_dist) <= 250000))
  expect_false(anyDuplicated(X[, c("context", "cpg_id", "target_gene")]) > 0)
  expect_true(all(X$evidence[X$context == "saliva"] == "smr"))
  expect_true(all(X$evidence[X$context != "saliva"] == "measured_eqtm"))
  expect_true(all(stats::na.omit(X$flags) %in% c("mhc", "inv17", "mhc|inv17")))
  # liver: two independent cohorts; a pair both hold agrees
  lv <- X[X$context == "liver"]
  expect_true(any(lv$n_cohorts == 2L))
  expect_true(all(lv$cohort_agreement[lv$n_cohorts == 2L] == 1))
  expect_identical(nrow(cpgd_reported_contexts("liver")), nrow(lv))
  expect_error(cpgd_reported_contexts("brain"), "should be one of")
  # the description table: none clears the bar for a voting source with its interval
  D <- attr(X, "contexts")
  expect_setequal(D$context, unique(X$context))
  expect_true(all(D$agreement_lo < 0.85))
  # none of these cohorts is in the measured table
  M <- cpgd_measured_eqtms()
  co <- unique(unlist(strsplit(X$cohorts, "|", fixed = TRUE)))
  expect_false(any(vapply(co, function(z) any(grepl(z, M$cohorts, fixed = TRUE)), logical(1))))
  expect_false(any(M$tissue %in% unique(X$context)))
  # the accessor returns a copy
  X[, direction := 0L]
  expect_true(all(cpgd_reported_contexts()$direction %in% c(-1, 1)))
})

test_that("region flags: one row per flagged CpG, logical columns", {
  f <- cpgd_region_flags()
  expect_identical(names(f), c("cpg_id", "in_corsiv", "in_imprinted_or_parental_asm", "in_asm"))
  expect_false(anyDuplicated(f$cpg_id) > 0)
  expect_true(all(vapply(f[, -1], is.logical, logical(1))))
  expect_true(all(f$in_corsiv | f$in_imprinted_or_parental_asm | f$in_asm))
  expect_gt(sum(f$in_corsiv), 4000L); expect_gt(sum(f$in_asm), 50000L)
})

test_that("cpgd_annotate() adds columns and changes nothing else", {
  X <- cpgd_reported_contexts(); f <- cpgd_region_flags()
  lv <- X[X$context == "liver"][1]
  x <- data.frame(cpg_id = c(lv$cpg_id, f$cpg_id[f$in_corsiv][1], "cg_not_there"),
                  target_gene = c(lv$target_gene, "NOGENE", "NOGENE"), best_direction = c(1, -1, NA))
  a <- cpgd_annotate(x)
  expect_identical(a$cpg_id, x$cpg_id); expect_identical(a$best_direction, x$best_direction)
  expect_true(grepl(paste0("liver:", if (lv$direction > 0) "+1" else "-1"), a$context_directions[1], fixed = TRUE))
  expect_identical(a$context_n[3], 0L); expect_true(is.na(a$context_directions[3]))
  expect_true(a$in_corsiv[2]); expect_false(a$in_corsiv[3])
  # without a gene column: every context sign of the CpG, with its gene
  b <- cpgd_annotate(data.frame(cpg_id = lv$cpg_id), region_flags = FALSE)
  expect_true(grepl(paste0("liver:", lv$target_gene, ":"), b$context_directions, fixed = TRUE))
  expect_false("in_asm" %in% names(b))
  expect_false("context_n" %in% names(cpgd_annotate(x, contexts = FALSE)))
  expect_error(cpgd_annotate(data.frame(a = 1)), "cpg_id")
})

test_that("a reported context never reaches best_direction", {
  X <- cpgd_reported_contexts()
  M <- cpgd_measured_eqtms()
  only <- X[!(X$cpg_id %chin% M$cpg_id) & !grepl(";", X$target_gene, fixed = TRUE)][1:20]
  only <- only[!is.na(only$cpg_id)]
  p <- cpg_gene_pairs(only$cpg_id, genes = only$target_gene, gene_mode = "pairwise", universal = FALSE,
                      verbose = FALSE, probe_qc = "ignore")
  expect_false(any(grepl("liver|retina|islet|fibroblast|saliva",
                         paste(p$best_rung, p$best_evidence, p$measured_tissues))))
})
