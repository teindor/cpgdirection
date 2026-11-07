# Measured brain cortex pairs (1.1.0) and the brain-target ladder of
# cpg_brain_bridge(): measured brain > brain T1 > brain T2 > measured blood >
# blood SMR (ruling 3, tools/validate_ladder_transport.R, 4 Nov 2026). The
# counts pinned below are the ones the documentation states, so a drifting
# data file fails loudly.

library(data.table)

test_that("the cohort vote: one vote per cohort, majority sign, ties dropped", {
  D <- data.table(
    cpg_id      = c("cgA", "cgA", "cgB", "cgB", "cgC", "cgC", "cgC", "cgD", "cgD", "cgD"),
    target_gene = c("G1",  "g1",  "G2",  "G2",  "G3",  "G3",  "G3",  "G4",  "G4",  "G4"),
    direction   = c(1,     1,     1,     -1,    -1,    -1,    1,     1,     -1,    -1),
    cohort      = c("c1",  "c2",  "c1",  "c2",  "c1",  "c2",  "c3",  "c1",  "c1",  "c2"),
    source      = "s", tss_dist = 10)
  M <- cpgdirection:::.cpgd_cohort_vote(D)
  expect_identical(data.table::key(M), c("cpg_id", "target_gene"))
  expect_false(anyDuplicated(M[, c("cpg_id", "target_gene"), with = FALSE]) > 0)
  a <- M[M$cpg_id == "cgA"]
  expect_equal(a$target_gene, "G1"); expect_equal(a$direction, 1L)
  expect_equal(a$n_cohorts, 2L); expect_equal(a$cohort_agreement, 1)
  expect_equal(a$cohorts, "c1|c2")
  expect_false("cgB" %in% M$cpg_id)                       # 1:1 tie dropped
  cc <- M[M$cpg_id == "cgC"]
  expect_equal(cc$direction, -1L); expect_equal(cc$n_cohorts, 3L)
  expect_equal(cc$cohort_agreement, 0.667)
  # cohort c1 disagrees with itself on cgD and casts no vote: c2 alone decides
  d <- M[M$cpg_id == "cgD"]
  expect_equal(d$direction, -1L); expect_equal(d$n_cohorts, 1L); expect_equal(d$cohorts, "c2")
})

test_that("brain_measured_eqtms holds the two documented cohorts, one row per pair per cohort", {
  R <- cpgd_brain_measured(collapse = FALSE)
  expect_true(all(c("cpg_id", "target_gene", "direction", "cohort", "source") %in% names(R)))
  expect_true(all(R$direction %in% c(-1, 1)))
  expect_false(anyDuplicated(R[, c("cpg_id", "target_gene", "cohort"), with = FALSE]) > 0)
  n <- table(R$cohort)
  expect_equal(as.integer(n[["ROSMAP_DLPFC_aged"]]), 8815L)
  expect_equal(as.integer(n[["Jaffe_DLPFC_ageadj"]]), 383L)
  expect_setequal(unique(R$source), c("wang2021_S2", "lead18b_GSE74193_GSE30272"))
  # the cohorts agree on 174 of 180 shared pairs (0.967, above the 0.85 bar)
  w <- dcast(R, cpg_id + target_gene ~ cohort, value.var = "direction")
  sh <- w[!is.na(w$ROSMAP_DLPFC_aged) & !is.na(w$Jaffe_DLPFC_ageadj)]
  expect_equal(nrow(sh), 180L)
  expect_equal(sum(sh$ROSMAP_DLPFC_aged == sh$Jaffe_DLPFC_ageadj), 174L)
  # the per-cohort copy is a copy: editing it leaves the cache alone
  R[, direction := 0]
  expect_true(all(cpgd_brain_measured(collapse = FALSE)$direction %in% c(-1, 1)))
})

test_that("cpgd_brain_measured() collapses by the cohort vote", {
  R <- cpgd_brain_measured(collapse = FALSE)
  M <- cpgd_brain_measured()
  expect_identical(names(M), c("cpg_id", "target_gene", "direction", "n_cohorts",
                               "cohort_agreement", "cohorts", "source", "tss_dist"))
  expect_false(anyDuplicated(M[, c("cpg_id", "target_gene"), with = FALSE]) > 0)
  # 9,018 pairs less the 6 cross-cohort ties
  expect_equal(nrow(M), data.table::uniqueN(R[, c("cpg_id", "target_gene"), with = FALSE]) - 6L)
  expect_equal(sum(M$n_cohorts == 2L), 174L)
  expect_true(all(M$cohort_agreement[M$n_cohorts == 2L] == 1))
  expect_error(cpgd_brain_measured(collapse = NA), "TRUE or FALSE")
  p <- M[M$cpg_id == "cg06846259" & M$target_gene == "POMC"]
  expect_equal(p$direction, -1L); expect_equal(p$n_cohorts, 2L)
  # a copy: editing it leaves the table the bridge reads alone
  M[, direction := 0L]
  expect_true(all(cpgd_brain_measured()$direction %in% c(-1L, 1L)))
})

test_that("the bridge reports measured cortex pairs per CpG", {
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  r <- cpg_brain_bridge("cg06846259")
  expect_true(all(c("measured_brain_genes", "measured_brain_directions",
                    "measured_brain_n_cohorts", "measured_brain_cohorts") %in% names(r)))
  g <- strsplit(r$measured_brain_genes, ";")[[1]]
  d <- strsplit(r$measured_brain_directions, ";")[[1]]
  expect_true("POMC" %in% g)
  expect_equal(d[g == "POMC"], "-1")
  expect_match(r$measured_brain_cohorts, "ROSMAP_DLPFC_aged")
  out <- paste(utils::capture.output(print(r)), collapse = " ")
  expect_true(grepl("brain (measured", out, fixed = TRUE))
  # a CpG measured in no cortex cohort comes back NA, not dropped
  M <- cpgd_brain_measured()
  B <- cpgd_brain_directions()
  none <- setdiff(unique(B$cpg_id), M$cpg_id)[1]
  r2 <- cpg_brain_bridge(c("cg06846259", none))
  expect_equal(nrow(r2), 2L)
  expect_true(is.na(r2$measured_brain_genes[r2$cpg_id == none]))
})

test_that("the brain-target ladder: measured brain, then T1, then T2, then blood", {
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  M <- cpgd_brain_measured()
  B <- cpgd_brain_directions()
  B <- B[B$tier %in% c("T1", "T2") & !grepl("_", B$target_gene)]
  mkey <- paste(M$cpg_id, M$target_gene)
  bkey <- paste(B$cpg_id, toupper(B$target_gene))

  # measured brain outranks the brain SMR layer; the brain SMR columns keep
  # their own answer
  r <- cpg_brain_bridge("cg06846259_POMC")
  expect_equal(r$requested_gene_best_evidence, "measured_brain")
  expect_equal(as.integer(r$requested_gene_best_direction), -1L)
  expect_equal(r$requested_gene_tier, "T1")
  expect_identical(r$requested_gene_best_accuracy,
                   CPGD_BRAIN_TARGET_ACCURACY[["measured_brain_replicated"]])
  expect_identical(attr(r, "accuracy"), "lower")
  expect_equal(r$requested_gene_accuracy, 0.929)
  expect_equal(cpg_brain_bridge("cg06846259_POMC", accuracy = "median")$requested_gene_accuracy,
               0.983)
  # a column subset keeps the class and prints as a plain table
  sub <- r[, c("requested_gene", "requested_gene_best_evidence"), with = FALSE]
  expect_true(any(grepl("measured_brain", utils::capture.output(print(sub)), fixed = TRUE)))

  for (tr in c("T1", "T2")) {
    x <- B[B$tier == tr & !(bkey %in% mkey)]
    x <- x[!duplicated(x$cpg_id)][1:3]          # the bridge answers one gene per CpG
    x <- x[!is.na(x$cpg_id)]
    skip_if(!nrow(x), paste("no unmeasured", tr, "pair"))
    q <- cpg_brain_bridge(paste0(x$cpg_id, "_", x$target_gene))
    q <- q[match(x$cpg_id, q$cpg_id)]
    expect_true(all(q$requested_gene_best_evidence == paste0("brain_", tr)), info = tr)
    expect_equal(as.integer(q$requested_gene_best_direction), as.integer(x$direction), info = tr)
    expect_true(all(q$requested_gene_best_accuracy == CPGD_BRAIN_ACCURACY[[tr]]), info = tr)
    expect_true(all(grepl("in-sample", q$requested_gene_best_accuracy)), info = tr)
    # the overlap-adjusted figure, as one number
    expect_true(all(q$requested_gene_accuracy == c(T1 = 0.926, T2 = 0.899)[[tr]]), info = tr)
  }
})

test_that("blood evidence is a labelled fallback and never outranks brain", {
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  M <- cpgd_brain_measured()
  B <- cpgd_brain_directions()
  B12 <- B[B$tier %in% c("T1", "T2")]
  brain_any <- unique(c(paste(M$cpg_id, M$target_gene), paste(B12$cpg_id, toupper(B12$target_gene))))
  Mb <- cpgd_measured_eqtms()
  Mb <- Mb[Mb$tissue == "blood" & !grepl("_", Mb$target_gene)]
  bl <- paste(Mb$cpg_id, toupper(Mb$target_gene))

  # a pair with both brain SMR and measured blood: brain wins
  both <- Mb[bl %in% paste(B12$cpg_id, toupper(B12$target_gene)) &
               !(bl %in% paste(M$cpg_id, M$target_gene))]
  both <- both[!duplicated(both$cpg_id)][1:3]
  both <- both[!is.na(both$cpg_id)]
  if (nrow(both)) {
    q <- cpg_brain_bridge(paste0(both$cpg_id, "_", both$target_gene))
    expect_true(all(q$requested_gene_best_evidence %in% c("brain_T1", "brain_T2")))
  }

  # measured blood only
  mb <- Mb[!(bl %in% brain_any)]
  mb <- mb[!duplicated(mb$cpg_id)][1:3]
  mb <- mb[!is.na(mb$cpg_id)]
  skip_if(!nrow(mb), "no blood-only measured pair")
  q <- cpg_brain_bridge(paste0(mb$cpg_id, "_", mb$target_gene))
  q <- q[match(mb$cpg_id, q$cpg_id)]
  expect_true(all(q$requested_gene_best_evidence == "measured_blood"))
  expect_equal(as.integer(q$requested_gene_best_direction), as.integer(sign(mb$direction)))
  expect_true(all(grepl("fallback only", q$requested_gene_best_accuracy)))
  expect_true(all(q$requested_gene_accuracy == 0.822))
  expect_true(all(is.na(q$requested_gene_direction)))     # the brain SMR answer stays NA

  # blood SMR only, a calling rung: fallback; a reported rung: no answer
  S <- cpgdirection:::.cpgd_smr_table()
  sk <- paste(S$cpg_id, toupper(S$target_gene))
  free <- !(sk %in% brain_any) & !(sk %in% bl) & !grepl("_", S$target_gene)
  sm <- S[free & S$smr_rung == "replicated"]
  sm <- sm[!duplicated(sm$cpg_id)][1:3]
  sm <- sm[!is.na(sm$cpg_id)]
  skip_if(!nrow(sm), "no SMR-only pair")
  q <- cpg_brain_bridge(paste0(sm$cpg_id, "_", sm$target_gene))
  q <- q[match(sm$cpg_id, q$cpg_id)]
  expect_true(all(q$requested_gene_best_evidence == "blood_SMR"))
  expect_true(all(q$requested_gene_best_accuracy ==
                    CPGD_BRAIN_TARGET_ACCURACY[["blood_SMR_replicated"]]))
  expect_true(all(q$requested_gene_accuracy == 0.745))
  rp <- S[free & S$smr_rung == "reported" &
            !(S$cpg_id %chin% S$cpg_id[S$smr_rung != "reported"])][1]
  if (nrow(rp) && !is.na(rp$cpg_id)) {
    q <- cpg_brain_bridge(paste0(rp$cpg_id, "_", rp$target_gene))
    expect_true(is.na(q$requested_gene_best_evidence))
  }
})

test_that("require_detected filters the blood-SMR fallback only", {
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  S <- cpgdirection:::.cpgd_smr_table()
  skip_if(!"child_blood_helix" %in% names(S), "SMR table predates child_blood_helix")
  M <- cpgd_brain_measured()
  B <- cpgd_brain_directions()
  B12 <- B[B$tier %in% c("T1", "T2")]
  Mb <- cpgd_measured_eqtms()
  taken <- unique(c(paste(M$cpg_id, M$target_gene), paste(B12$cpg_id, toupper(B12$target_gene)),
                    paste(Mb$cpg_id[Mb$tissue == "blood"], toupper(Mb$target_gene[Mb$tissue == "blood"]))))
  sk <- paste(S$cpg_id, toupper(S$target_gene))
  ok <- !(sk %in% taken) & !grepl("_", S$target_gene) & S$smr_rung == "replicated"
  nul <- S[ok & S$child_blood_helix == "null"][1]
  sig <- S[ok & S$child_blood_helix == "sig" & S$cpg_id != nul$cpg_id][1]
  skip_if(is.na(nul$cpg_id) || is.na(sig$cpg_id), "no suitable pairs")
  inp <- c(paste0(nul$cpg_id, "_", nul$target_gene), paste0(sig$cpg_id, "_", sig$target_gene))
  a <- cpg_brain_bridge(inp, age_group = "child")
  expect_message(k <- cpg_brain_bridge(inp, age_group = "child", require_detected = TRUE),
                 "require_detected: [0-9,]+ blood-SMR pair")
  a <- a[match(c(nul$cpg_id, sig$cpg_id), a$cpg_id)]
  k <- k[match(c(nul$cpg_id, sig$cpg_id), k$cpg_id)]
  expect_equal(a$requested_gene_best_evidence, c("blood_SMR", "blood_SMR"))
  expect_equal(a$requested_gene_child_blood_helix, c("null", "sig"))
  expect_true(is.na(k$requested_gene_best_evidence[1]))
  expect_equal(k$requested_gene_best_evidence[2], "blood_SMR")
  expect_false("requested_gene_child_blood_helix" %in%
                 names(cpg_brain_bridge(inp)))
})
