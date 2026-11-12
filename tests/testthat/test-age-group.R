# The age_group argument (1.1.0, spec approved 4 Nov 2026). Age moves no
# direction and no rung; age_group changes which validation figures are
# quoted, which age-specific columns are returned, and (require_detected) an
# optional filter on the SMR layer. The default "adult" must reproduce the
# earlier output exactly.

library(data.table)

# per-pair accuracy by rung: the band is what the ladder ranks (replicated ->
# smr_acc_95, high -> smr_acc_90, moderate -> smr_acc_80, weak -> smr_acc_65)
ag_acc <- c(replicated = 0.97, high = 0.93, moderate = 0.85, weak = 0.70, reported = 0.60)
ag_smr <- function(cpg, gene, direction, rung, tier = "S1", child = "sig",
                   outcome = "eQTLGen", replication = "concordant") {
  data.table(cpg_id = cpg, target_gene = gene, direction = direction,
             smr_tier = tier, p_SMR = 1e-10, n_instruments = 3L,
             instrument_agreement = 1, cpg_gene_dist = 1000,
             top_instrument = "rs1", p_HEIDI = 0.5, nsnp_HEIDI = 3L,
             heidi_status = "pass", outcome = outcome,
             replication = replication, smr_rung = rung,
             smr_accuracy = ag_acc[[rung]], child_blood_helix = child)
}
ag_fx_smr <- function() rbindlist(list(
  ag_smr("cg20000001", "GENEA", -1, "replicated", "S1", "sig"),
  ag_smr("cg20000001", "GENEB",  1, "moderate",   "S2", "null", replication = "eqtlgen_only"),
  ag_smr("cg20000002", "GENEC", -1, "high",       "S1", "untested", replication = "eqtlgen_only"),
  ag_smr("cg20000003", "GENED",  1, "weak",       "S3", "nominal", replication = "eqtlgen_only"),
  ag_smr("cg20000003", "GENEE", -1, "reported",   "S2", "sig", replication = "discordant")))
ag_fx_measured <- function() data.table(cpg_id = "cg20000002", target_gene = "GENEC",
                                        tissue = "blood", direction = 1, tss_dist = 500)
ag_src <- function(smr = ag_fx_smr()) list(manifest = NULL, lookup = NULL,
                                           measured = ag_fx_measured(), smr = smr,
                                           probe_qc = NULL, onco = NULL)
ag_cpgs <- c("cg20000001", "cg20000002", "cg20000003")
ag_pairs <- function(...) cpg_gene_pairs(ag_cpgs, universal = FALSE, verbose = FALSE, ...)

# ---- constants -------------------------------------------------------------

test_that("child blood quotes HELIX where HELIX scores the source on 100+ pairs, adult otherwise", {
  R <- CPGD_RUNG_ACCURACY
  for (src in c("smr_acc_95", "smr_acc_80", "smr_acc_65", "catalogue_consensus_A")) {
    ch <- R[R$source == src & R$target == "blood" & R$age == "child", ]
    k <- cpgdirection:::.cpgd_rung_label(src, "blood", "child")
    a <- cpgdirection:::.cpgd_rung_label(src, "blood", "adult")
    if (nrow(ch) && ch$n >= 100) {
      expect_true(grepl("in child blood", k, fixed = TRUE), info = src)
      expect_identical(cpgdirection:::.cpgd_rung_label(src, "blood", "any"), paste0(a, " | ", k))
    } else {
      expect_identical(k, a, info = src)
      expect_identical(cpgdirection:::.cpgd_rung_label(src, "blood", "any"), a)
    }
  }
  expect_true(grepl("in adult blood", cpgdirection:::.cpgd_rung_label("smr_acc_90", "blood", "adult")))
  expect_setequal(names(CPGD_SMR_DETECTION_CHILD),
                  c("replicated", "high", "moderate", "weak", "reported"))
  expect_identical(CPGD_CHILD_BLOOD_LEVELS, c("sig", "nominal", "null", "untested"))
})

# ---- hermetic: the pair API ------------------------------------------------

test_that("the default equals age_group = 'adult' and carries no age-specific column", {
  p0 <- ag_pairs(sources = ag_src())
  p1 <- ag_pairs(sources = ag_src(), age_group = "adult")
  expect_identical(p0, p1)
  expect_false("smr_child_blood_helix" %in% names(p0))
  expect_null(attr(p0, "age_group"))
  # the adult strings are the adult-blood validation's
  r <- p0[p0$target_gene == "GENEA", ]
  expect_equal(r$best_evidence, "smr_replicated")
  expect_identical(r$best_expected_accuracy,
                   cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "adult"))
  expect_identical(r$best_rung, "smr_acc_95")
})

test_that("the default output: as at 8896789 except the ranked ladder and the two new columns", {
  # At 8896789 this fixture gave directions -1, 1, 1, 1, NA with evidence
  # smr_replicated, smr_moderate, measured, smr_weak, no_evidence. 1.1.0 ranks
  # the ladder by the lower bound of each source's figure in the call's tissue,
  # and ranks SMR pairs by the band of their per-pair accuracy: in blood a
  # single measured study outranks the 0.90-0.95 band, so GENEC keeps the
  # measured sign. The columns are those of 8896789 plus best_rung (after
  # best_evidence), best_accuracy (after best_expected_accuracy) and
  # smr_accuracy (after smr_outcome).
  p <- ag_pairs(sources = ag_src())
  expect_identical(setdiff(names(p), c("best_rung", "best_accuracy")), c(
    "input_id", "input", "cpg_id", "target_gene", "pair_id", "mapping_sources", "mapping_primary",
    "mapping_strength", "best_direction", "best_evidence", "best_confidence", "direction_tier",
    "probability_plus1", "best_expected_accuracy", "usable", "abstain_reason", "n_targets_for_cpg",
    "is_coeffect", "has_manifest", "has_lookup", "has_measured_eqtm", "has_smr", "has_brain",
    "has_distal", "has_onco", "has_input_annotation", "measured_direction", "measured_tissues",
    "measured_n_cohorts", "measured_cohort_agreement", "lookup_direction", "lookup_probability_plus1",
    "lookup_confidence", "lookup_evidence_tier", "n_tissues_calling", "tissue_agreement",
    "smr_direction", "smr_tier", "p_SMR", "n_instruments", "instrument_agreement", "p_HEIDI",
    "heidi_status", "smr_rung", "smr_replication", "smr_outcome", "smr_accuracy", "abs_dist", "p_universal",
    "dir_universal", "dist_unanimous", "onco_direction", "onco_tier", "onco_n_cancers",
    "onco_sign_agreement", "onco_median_r", "onco_agreement", "input_ids", "inputs"))
  expect_identical(names(p)[match("best_evidence", names(p)) + 1L], "best_rung")
  expect_identical(names(p)[match("best_expected_accuracy", names(p)) + 1L], "best_accuracy")
  expect_identical(p$pair_id, c("cg20000001|GENEA", "cg20000001|GENEB", "cg20000002|GENEC",
                                "cg20000003|GENED", "cg20000003|GENEE"))
  ob <- cpgdirection:::.cpgd_rung_order("blood")
  expect_true(match("measured_blood_single", ob) < match("smr_acc_90", ob))
  expect_identical(p$best_direction, c(-1, 1, 1, 1, NA))
  expect_identical(p$best_evidence, c("smr_replicated", "smr_moderate", "measured", "smr_weak",
                                      "no_evidence"))
  rk <- c("smr_acc_95", "smr_acc_80", "measured_blood_single", "smr_acc_65")
  expect_identical(p$best_rung, c(rk, "no_evidence"))
  expect_identical(p$best_expected_accuracy,
                   c(cpgdirection:::.cpgd_rung_label(rk, "blood", "adult"), NA))
  expect_true(all(grepl("in adult blood", p$best_expected_accuracy[1:4], fixed = TRUE)))
  # the number: the lower Wilson bound of the same validation, floored at 0.5
  R <- CPGD_RUNG_ACCURACY[CPGD_RUNG_ACCURACY$target == "blood" & CPGD_RUNG_ACCURACY$age == "adult", ]
  lo <- R$lo[match(rk, R$source)]
  expect_identical(p$best_accuracy, c(pmax(lo, 0.5), NA))
  # in solid tissue too the single measured study outranks that SMR band:
  # GENEC is +1
  q <- ag_pairs(sources = ag_src(), tissue = "solid_tissue")
  expect_identical(q$best_rung[3], "measured_blood_single")
  expect_identical(q$best_direction[3], 1)
})

test_that("an SMR resource without child_blood_helix: warned, NA, and require_detected ignored", {
  old <- ag_fx_smr()[, child_blood_helix := NULL]
  expect_warning(k <- ag_pairs(sources = ag_src(old), age_group = "child"), "predates")
  expect_true(all(is.na(k$smr_child_blood_helix)))
  expect_warning(d <- ag_pairs(sources = ag_src(old), age_group = "child",
                               require_detected = TRUE), "cannot be applied")
  # the SMR layer is kept, not silently emptied
  expect_identical(d$best_evidence, ag_pairs(sources = ag_src())$best_evidence)
})

test_that("the default output does not depend on whether the table carries child_blood_helix", {
  # a pre-1.1.0 table (no child column) and the new table give the identical
  # default output: the column never leaks into the adult answer
  old <- ag_fx_smr()[, child_blood_helix := NULL]
  expect_identical(ag_pairs(sources = ag_src(old)), ag_pairs(sources = ag_src()))
})

test_that("age_group = 'child' adds child_blood_helix and the child figures, nothing else", {
  a <- ag_pairs(sources = ag_src())
  k <- ag_pairs(sources = ag_src(), age_group = "child")
  expect_true("smr_child_blood_helix" %in% names(k))
  expect_true(all(stats::na.omit(k$smr_child_blood_helix) %in% CPGD_CHILD_BLOOD_LEVELS))
  expect_equal(k$smr_child_blood_helix[k$target_gene == "GENEB"], "null")
  expect_equal(attr(k, "age_group"), "child")
  # directions, evidence and the pair set are age-independent
  expect_identical(k$pair_id, a$pair_id)
  expect_identical(k$best_direction, a$best_direction)
  expect_identical(k$best_evidence, a$best_evidence)
  # each row quotes its source's child-blood figure
  expect_identical(k$best_expected_accuracy[1:4],
                   cpgdirection:::.cpgd_rung_label(k$best_rung[1:4], "blood", "child"))
  expect_true(grepl("in child blood", k$best_expected_accuracy[k$target_gene == "GENEA"], fixed = TRUE))
  out <- utils::capture.output(print(k))
  expect_true(any(grepl("child blood (HELIX, 6-11 y", out, fixed = TRUE)))
})

test_that("age_group = 'any' quotes adult and child figures side by side", {
  y <- ag_pairs(sources = ag_src(), age_group = "any")
  s <- y$best_expected_accuracy[y$target_gene == "GENEA"]
  expect_identical(s, paste0(
    cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "adult"), " | ",
    cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "child")))
  expect_true("smr_child_blood_helix" %in% names(y))
})

test_that("require_detected drops only the SMR rows HELIX did not detect", {
  pw <- function(...) cpg_gene_pairs(
    c("cg20000001", "cg20000001", "cg20000002", "cg20000003", "cg20000003"),
    genes = c("GENEA", "GENEB", "GENEC", "GENED", "GENEE"), gene_mode = "pairwise",
    universal = FALSE, verbose = FALSE, sources = ag_src(), ...)
  a <- pw(age_group = "child")
  expect_message(k <- pw(age_group = "child", require_detected = TRUE),
                 "require_detected: 2 SMR pair")
  g <- function(x, gene, col) x[[col]][x$target_gene == gene]
  # sig / nominal rows are kept, unchanged
  for (gn in c("GENEA", "GENED", "GENEE"))
    for (col in c("smr_direction", "smr_rung", "best_direction", "best_evidence"))
      expect_identical(g(k, gn, col), g(a, gn, col), info = paste(gn, col))
  # a null SMR row loses its SMR evidence; it does not fall back to anything else
  expect_true(is.na(g(k, "GENEB", "smr_direction")))
  expect_equal(g(k, "GENEB", "best_evidence"), "no_evidence")
  # measured evidence is never filtered: with its (outranking) SMR row dropped,
  # GENEC takes its measured sign
  expect_equal(g(k, "GENEC", "best_evidence"), "measured")
  expect_identical(g(k, "GENEC", "best_direction"), g(a, "GENEC", "measured_direction"))
  expect_identical(g(k, "GENEC", "measured_direction"), g(a, "GENEC", "measured_direction"))
  expect_true(is.na(g(k, "GENEC", "smr_direction")))
  # in discovery mode a dropped SMR row proposes no target
  d <- suppressMessages(ag_pairs(sources = ag_src(), age_group = "child",
                                 require_detected = TRUE))
  expect_false("GENEB" %in% d$target_gene)
  expect_true(all(c("GENEA", "GENEC", "GENED", "GENEE") %in% d$target_gene))
})

test_that("require_detected without age_group = 'child' warns and is ignored", {
  expect_warning(p <- ag_pairs(sources = ag_src(), require_detected = TRUE),
                 "applies only with age_group")
  expect_identical(p, ag_pairs(sources = ag_src()))
  expect_warning(ag_pairs(sources = ag_src(), age_group = "any", require_detected = TRUE),
                 "applies only with age_group")
  expect_error(ag_pairs(sources = ag_src(), require_detected = NA), "TRUE or FALSE")
})

test_that("match.arg rejects any other age_group", {
  expect_error(ag_pairs(sources = ag_src(), age_group = "teen"), "should be one of")
  expect_error(cpgd_smr_directions(age_group = "teen"), "should be one of")
  expect_error(cpg_expression_direction("cg00000029", age_group = "teen", verbose = FALSE),
               "should be one of")
  expect_error(cpg_direction_all_tissues("cg00000029", age_group = "teen", verbose = FALSE),
               "should be one of")
})

# ---- the shipped SMR resource ---------------------------------------------

.ag_child_table <- function() {
  S <- tryCatch(cpgdirection:::.cpgd_smr_table(), error = function(e) NULL)
  if (is.null(S) || !nrow(S) || !"child_blood_helix" %in% names(S)) return(NULL)
  S
}

test_that("cpgd_smr_directions(): adult is the table without child_blood_helix", {
  skip_if(!cpgd_has_data("smr_directions"), "SMR table unavailable")
  a <- cpgd_smr_directions()
  expect_identical(a, cpgd_smr_directions(age_group = "adult"))
  expect_false("child_blood_helix" %in% names(a))
  expect_identical(data.table::key(a), c("cpg_id", "target_gene"))
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  expect_identical(names(a), setdiff(names(S), "child_blood_helix"))
  expect_equal(nrow(a), nrow(S))
})

test_that("cpgd_smr_directions(): child adds child_blood_helix with the four values only", {
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  k <- cpgd_smr_directions(age_group = "child")
  expect_true("child_blood_helix" %in% names(k))
  expect_true(all(k$child_blood_helix %in% CPGD_CHILD_BLOOD_LEVELS))
  expect_true(all(CPGD_CHILD_BLOOD_LEVELS %in% k$child_blood_helix))
  # informational only: the rung rules do not look at it
  expect_true(any(k$smr_rung == "replicated" & k$child_blood_helix == "null"))
  expect_true(any(k$smr_rung == "reported" & k$child_blood_helix == "sig"))
  expect_identical(cpgd_smr_directions(age_group = "any")$child_blood_helix, k$child_blood_helix)
})

test_that("cpgd_smr_directions(): require_detected keeps exactly the sig/nominal pairs", {
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  expect_message(d <- cpgd_smr_directions(age_group = "child", require_detected = TRUE),
                 "require_detected")
  expect_true(all(d$child_blood_helix %in% CPGD_CHILD_DETECTED))
  expect_equal(nrow(d), sum(S$child_blood_helix %in% CPGD_CHILD_DETECTED))
  expect_warning(w <- cpgd_smr_directions(require_detected = TRUE), "applies only")
  expect_identical(w, cpgd_smr_directions())
})

test_that("the age views of the SMR table are copies: editing them never reaches the cache", {
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  d0 <- cpgdirection:::.cpgd_smr_table()$direction[1]
  for (ag in c("adult", "child", "any")) {
    k <- cpgd_smr_directions(age_group = ag)
    k[1L, direction := 99L]
    expect_identical(cpgdirection:::.cpgd_smr_table()$direction[1], d0, info = ag)
  }
})

test_that("require_detected keeps smr_in_table honest and says why SMR is missing", {
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  det <- S$child_blood_helix %in% CPGD_CHILD_DETECTED
  gone <- setdiff(unique(S$cpg_id[!det]), S$cpg_id[det])[1:3]
  r <- suppressMessages(cpg_expression_direction(gone, universal = FALSE, verbose = FALSE,
                                                  age_group = "child", require_detected = TRUE))
  expect_true(all(r$smr_in_table))
  expect_true(all(is.na(r$smr_direction)))
  out <- paste(utils::capture.output(print(r)), collapse = " ")
  expect_true(grepl("require_detected removed them", out, fixed = TRUE))
  expect_false(grepl("never assayed", out, fixed = TRUE))
})

test_that("cpg_expression_direction(): default equals adult; child adds the column, not a direction", {
  S <- .ag_child_table(); skip_if(is.null(S), "SMR table predates child_blood_helix")
  cp <- unique(S$cpg_id[S$smr_rung %in% c("replicated", "moderate")])[1:30]
  a0 <- cpg_expression_direction(cp, universal = FALSE, verbose = FALSE)
  a1 <- cpg_expression_direction(cp, universal = FALSE, verbose = FALSE, age_group = "adult")
  expect_identical(a0, a1)
  expect_false("smr_child_blood_helix" %in% names(a0))
  k <- cpg_expression_direction(cp, universal = FALSE, verbose = FALSE, age_group = "child")
  expect_true(all(stats::na.omit(k$smr_child_blood_helix) %in% CPGD_CHILD_BLOOD_LEVELS))
  expect_identical(k$best_direction, a0$best_direction)
  expect_identical(k$best_evidence, a0$best_evidence)
  smr <- k$best_evidence %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak")
  skip_if(!any(smr), "no SMR-won row in this sample")
  expect_identical(k$best_expected_accuracy[smr],
                   cpgdirection:::.cpgd_rung_label(k$best_rung[smr], "blood", "child"))
  expect_identical(a0$best_expected_accuracy[smr],
                   cpgdirection:::.cpgd_rung_label(a0$best_rung[smr], "blood", "adult"))
  expect_true(any(grepl("child blood (HELIX", utils::capture.output(print(k)), fixed = TRUE)))
  # require_detected removes undetected SMR evidence and nothing measured
  expect_message(kd <- cpg_expression_direction(cp, universal = FALSE, verbose = FALSE,
                                                age_group = "child", require_detected = TRUE),
                 "require_detected")
  expect_true(all(kd$smr_child_blood_helix %in% c(CPGD_CHILD_DETECTED, NA)))
  meas <- a0$best_evidence == "measured"
  expect_identical(kd$best_evidence == "measured", meas)
  expect_identical(kd$best_direction[meas], a0$best_direction[meas])
})

test_that("the brain bridge caveat differs between adult and child", {
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  a <- cpg_brain_bridge("cg06846259")
  k <- cpg_brain_bridge("cg06846259", age_group = "child")
  expect_false(identical(attr(a, "brain_caveat"), attr(k, "brain_caveat")))
  expect_identical(attr(a, "brain_caveat"), CPGD_BRAIN_CAVEAT[["adult"]])
  expect_match(attr(k, "brain_caveat"), "1-25 y")
  expect_match(attr(k, "brain_caveat"), "fetal/< 1 y untested", fixed = TRUE)
  expect_match(attr(a, "brain_caveat"), "fetal and infant (< 1 y) cortex is untested", fixed = TRUE)
  expect_identical(attr(cpg_brain_bridge("cg06846259", age_group = "any"), "brain_caveat"),
                   CPGD_BRAIN_CAVEAT[["any"]])
  pk <- paste(utils::capture.output(print(k)), collapse = " ")
  expect_true(grepl("Age group child", pk, fixed = TRUE))
  expect_error(cpg_brain_bridge("cg06846259", age_group = "teen"), "should be one of")
  expect_warning(cpg_brain_bridge("cg06846259", require_detected = TRUE), "applies only")
})
