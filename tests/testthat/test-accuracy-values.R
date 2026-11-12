# The ranked ladder and its numbers (1.1.0). Every evidence source is scored
# in each tissue against an independent yardstick (CPGD_RUNG_ACCURACY,
# tools/validate_rung_accuracy.R); the ladder order in a tissue and the number
# a call carries both come from that table.

library(data.table)

# per-pair accuracy by rung: the band is what the ladder ranks
acc_acc <- c(replicated = 0.97, high = 0.93, moderate = 0.85, weak = 0.70, reported = 0.60)
acc_smr <- function(cpg, gene, direction, rung, tier = "S1", child = "sig",
                    outcome = "eQTLGen", replication = "concordant") {
  data.table(cpg_id = cpg, target_gene = gene, direction = direction,
             smr_tier = tier, p_SMR = 1e-10, n_instruments = 3L,
             instrument_agreement = 1, cpg_gene_dist = 1000,
             top_instrument = "rs1", p_HEIDI = 0.5, nsnp_HEIDI = 3L,
             heidi_status = "pass", outcome = outcome,
             replication = replication, smr_rung = rung,
             smr_accuracy = acc_acc[[rung]], child_blood_helix = child)
}
acc_src <- function() list(
  manifest = NULL, lookup = NULL, probe_qc = NULL, onco = NULL,
  measured = data.table(cpg_id = "cg20000002", target_gene = "GENEC",
                        tissue = "blood", direction = 1, tss_dist = 500),
  smr = rbindlist(list(
    acc_smr("cg20000001", "GENEA", -1, "replicated", "S1", "sig"),
    acc_smr("cg20000001", "GENEB",  1, "moderate",   "S2", "null", replication = "eqtlgen_only"),
    acc_smr("cg20000002", "GENEC", -1, "high",       "S1", "untested", replication = "eqtlgen_only"),
    acc_smr("cg20000003", "GENED",  1, "weak",       "S3", "nominal", replication = "eqtlgen_only"),
    acc_smr("cg20000003", "GENEE", -1, "reported",   "S2", "sig", replication = "discordant"))))
acc_pairs <- function(...) cpg_gene_pairs(c("cg20000001", "cg20000002", "cg20000003"),
                                          universal = FALSE, verbose = FALSE,
                                          sources = acc_src(), ...)
num <- function(...) cpgdirection:::.cpgd_accuracy_number(...)
order_in <- function(t) cpgdirection:::.cpgd_rung_order(t)
R <- CPGD_RUNG_ACCURACY
# every quoted string carries a number, and nothing without a string has one
expect_mapped <- function(str, val, info = NULL) {
  str <- as.character(str)
  quoted <- !is.na(str)
  expect_false(any(quoted & is.na(val)), info = info)
  expect_true(all(is.na(val[!quoted])), info = info)
}

# ---- the validation table ----------------------------------------------------

test_that("the rank table is well formed and its intervals are the Wilson intervals of its counts", {
  expect_identical(names(R), c("source", "target", "age", "yardstick", "n", "n_cpg", "k",
                               "agree", "lo", "hi"))
  expect_false(anyDuplicated(paste(R$source, R$target, R$age)) > 0)
  expect_true(all(R$k <= R$n & R$n_cpg <= R$n & R$n > 0))
  expect_setequal(unique(paste(R$target, R$age)),
                  c("blood adult", "blood child", "nasal_epithelium any", "solid_tissue any"))
  w <- t(mapply(function(k, n) cpgdirection:::.cpgd_wilson(k, n), R$k, R$n))
  expect_equal(R$agree, round(w[, "agree"], 4), tolerance = 1e-4)
  expect_equal(R$lo, round(w[, "lo"], 4), tolerance = 1e-4)
  expect_equal(R$hi, round(w[, "hi"], 4), tolerance = 1e-4)
  # every source the ladders use has a description; all but two have a row:
  # the replicated nasal sign (nasal epithelium has one donor group) and the
  # solid-tissue model's tier-B single call (since the 0.99.7 lookups the solid
  # model's tier-B calls abstain, so no yardstick pair)
  expect_true(all(CPGD_LADDER_SOURCES %in% names(CPGD_RUNG_DESC)))
  expect_identical(setdiff(CPGD_LADDER_SOURCES, R$source),
                   c("measured_nasal_replicated", "catalogue_single_solid_B", "distance_only_b6"))
})

test_that("cpgd_accuracy_values(): one row per quoted string, probabilities, lower <= median", {
  V <- cpgd_accuracy_values()
  expect_identical(names(V), c("rung", "population", "median", "lower", "upper",
                               "n", "interval", "basis", "label"))
  expect_false(anyDuplicated(V$label) > 0)
  expect_false(any(grepl(" | ", V$label, fixed = TRUE)))
  q <- !is.na(V$median)
  expect_true(all(V$lower[q] >= 0 & V$lower[q] <= V$median[q] & V$median[q] <= 1))
  u <- q & !is.na(V$upper)
  expect_true(all(V$median[u] <= V$upper[u]))
  # every label has a number; a source no yardstick scores has the neutral 0.5
  expect_true(all(q))
  uq <- grepl("^unquantified", V$label)
  expect_true(all(V$median[uq] == 0.5 & V$lower[uq] == 0.5))
  expect_true(all(nzchar(V$basis) & nzchar(V$interval)))
  # every ladder label and every brain label is in it
  for (t in CPGD_TISSUES) for (a in c("adult", "child")) {
    lab <- cpgdirection:::.cpgd_rung_label(names(CPGD_RUNG_DESC), t, a)
    expect_true(all(stats::na.omit(lab) %in% V$label), info = paste(t, a))
  }
  expect_true(all(CPGD_BRAIN_TARGET_ACCURACY %in% V$label))
})

test_that("a figure is the source's own in its tissue from 30 pairs, cross-tissue below", {
  fig <- cpgdirection:::.cpgd_rung_figure1
  cross <- function(s) {
    all <- R[R$source == s, ]
    list(n = sum(all$n), agree = round(sum(all$k) / sum(all$n), 4))
  }
  for (i in seq_len(nrow(R))) {
    s <- R$source[i]; t <- R$target[i]; info <- paste(s, t, R$age[i])
    if (t == "blood" && R$age[i] == "child") {
      # child blood: the HELIX row where it has 100+ pairs, the adult figure otherwise
      f <- fig(s, "blood", "child")
      if (R$n[i] >= CPGD_CHILD_MIN_N) {
        expect_identical(f$basis, "own", info = info)
        expect_identical(f$population, "child blood", info = info)
        expect_equal(f$agree, R$agree[i], info = info)
        expect_equal(f$lo, R$lo[i], info = info)
      } else {
        expect_identical(f, fig(s, "blood", "adult"), info = info)
      }
      next
    }
    f <- fig(s, t, if (t == "blood") "adult" else "any")
    if (R$n[i] >= CPGD_FIGURE_MIN_N) {
      expect_identical(f$basis, "own", info = info)
      expect_equal(f$agree, R$agree[i], info = info)
      expect_equal(f$lo, R$lo[i], info = info)
    } else if (cross(s)$n > 0) {
      expect_identical(f$basis, "cross-tissue", info = info)
      expect_equal(f$n, cross(s)$n, info = info)
      expect_equal(f$agree, cross(s)$agree, info = info)
    } else {
      expect_null(f, info = info)
    }
  }
})

# ---- the order -------------------------------------------------------------------

test_that("the order is the lower bound alone; a point below 0.5 never answers", {
  fig <- cpgdirection:::.cpgd_rung_figure1
  for (t in CPGD_TISSUES) {
    age <- if (t == "blood") "adult" else "any"
    o <- order_in(t)
    F <- lapply(o, function(s) fig(s, t, age))
    lo <- vapply(F, `[[`, numeric(1), "lo")
    pt <- vapply(F, `[[`, numeric(1), "agree")
    # non-increasing lower bound; ties broken by the point estimate
    expect_false(is.unsorted(rev(lo)), info = t)
    tie <- which(diff(lo) == 0)
    expect_true(all(pt[tie] >= pt[tie + 1L]), info = t)
    # in: every source with a figure whose point is >= 0.5 (a distance band:
    # >= 0.55, the distance floor); out: the rest
    for (s in CPGD_LADDER_SOURCES) {
      f <- fig(s, t, age)
      floor_s <- if (startsWith(s, "distance_")) cpgdirection:::CPGD_DISTANCE_FLOOR else 0.5
      expect_identical(s %in% o, !is.null(f) && f$agree >= floor_s, info = paste(t, s))
    }
  }
  ob <- order_in("blood")
  # a small sample does not outrank a large one on a noisy point: the blood
  # model's tier-B single calls (40 pairs) sit below the catalogue consensus
  # (2,145 pairs) although their point estimates are level
  expect_gt(match("catalogue_single_blood_B", ob), match("catalogue_consensus_A", ob))
  # the SMR bands are ranked like any source: the top band above a single
  # measured blood study, the bands in the order of their accuracy
  expect_lt(match("smr_acc_95", ob), match("measured_blood_single", ob))
  expect_false(is.unsorted(match(paste0("smr_acc_", c(95, 90, 80, 65, 55)), ob)))
  # distance is ranked by band: the nearest band above the solid model's and
  # the nasal model's single calls, and the bands in the order of their distance
  expect_gt(match("catalogue_single_solid_A", ob), match("distance_only_b1", ob))
  expect_lt(match("distance_only_b1", ob), match("catalogue_single_nasal_A", ob))
  expect_lt(match("distance_only_b1", ob), match("distance_only_b2", ob))
  expect_lt(match("distance_only_b2", ob), match("distance_only_b3", ob))
  # the lowest SMR band is below 0.5 in solid tissue: it never answers there
  expect_false("smr_acc_55" %in% order_in("solid_tissue"))
  # nasal epithelium has one donor group (EVA-PR): no nasal sign is replicated
  for (t in CPGD_TISSUES) expect_false("measured_nasal_replicated" %in% order_in(t), info = t)
  # the order takes no age: age_group changes numbers, never rungs (tested on
  # the pairs below and on real data)
  expect_identical(names(formals(cpgdirection:::.cpgd_rung_order)), c("target", "pool"))
})

test_that("a point below 0.5 is left out of the order, whatever its lower bound", {
  # pinned on made-up figures, so that it does not depend on the data (with
  # the 0.99.7 lookups, catalogue_single_blood_A in nasal epithelium is one:
  # 0.468 [0.333])
  local_mocked_bindings(.cpgd_rung_figure1 = function(source, target = "blood", age = "adult")
    switch(source,
           a = list(agree = 0.49, lo = 0.48),
           b = list(agree = 0.60, lo = 0.40),
           c = list(agree = 0.90, lo = 0.70),
           NULL))
  expect_identical(cpgdirection:::.cpgd_rung_order("blood", pool = c("a", "b", "c", "d")),
                   c("c", "b"))
})

test_that("the ladder answers with the first source in the order; conflicts stop at distance", {
  n <- 4L
  cand <- list(a = c(1, NA, NA, NA), b = c(-1, -1, NA, NA),
               distance_only = c(1, 1, 1, NA), c = c(NA, NA, -1, -1))
  r <- cpgdirection:::.cpgd_apply_ladder(cand, c("a", "b", "distance_only", "c"),
                                         conflict = c(FALSE, FALSE, TRUE, TRUE), n)
  expect_identical(r$best, c(1, -1, NA, NA))
  expect_identical(r$key, c("a", "b", "tissue_conflict", "tissue_conflict"))
  r <- cpgdirection:::.cpgd_apply_ladder(cand, c("c", "a", "b", "distance_only"),
                                         conflict = rep(FALSE, n), n)
  expect_identical(r$key, c("a", "b", "c", "c"))
  # with distance in bands, conflicts stop at the highest-ranked band
  cand <- list(a = c(1, NA, NA, NA), distance_only_b1 = c(NA, 1, NA, NA), b = c(NA, NA, -1, -1),
               distance_only_b3 = c(NA, NA, NA, 1))
  r <- cpgdirection:::.cpgd_apply_ladder(cand, c("a", "distance_only_b1", "b", "distance_only_b3"),
                                         conflict = c(FALSE, FALSE, TRUE, FALSE), n)
  expect_identical(r$key, c("a", "distance_only_b1", "tissue_conflict", "b"))
  # the coarse label: an SMR band takes the rung of the row that supplied it
  expect_identical(cpgdirection:::.cpgd_rung_evidence(
    c("measured_solid_single", "smr_acc_95", "smr_acc_80", "smr_acc_65",
      "catalogue_single_nasal_B", "catalogue_consensus_A", "onco_consensus", "no_evidence",
      "distance_only_b3"),
    c(NA, "smr_replicated", "smr_moderate", "smr_weak", NA, NA, NA, NA, NA)),
    c("measured", "smr_replicated", "smr_moderate", "smr_weak", "catalogue_single",
      "catalogue_consensus", "onco_consensus", "no_evidence", "distance_only"))
})

test_that("in every tissue's order, the first source with a sign answers", {
  set.seed(20261104)
  n <- 300L
  for (t in CPGD_TISSUES) {
    o <- order_in(t)
    cand <- lapply(o, function(s) ifelse(stats::runif(n) < 0.12, sample(c(-1, 1), n, TRUE), NA_real_))
    names(cand) <- o
    r <- cpgdirection:::.cpgd_apply_ladder(cand, o, conflict = rep(FALSE, n), n)
    M <- do.call(cbind, cand)
    first <- apply(M, 1L, function(v) if (any(!is.na(v))) o[which(!is.na(v))[1]] else NA_character_)
    expect_identical(r$key, unname(first), info = t)
    expect_identical(r$best, vapply(seq_len(n), function(i)
      if (is.na(first[i])) NA_real_ else M[i, first[i]], numeric(1)), info = t)
  }
})

test_that("replicated means two donor groups: one study's strata or cell types are one cohort", {
  dg <- cpgdirection:::.cpgd_donor_group
  expect_identical(dg(c("GTEx: Lung", "EVA-PR: Controls", "Framingham Heart Study",
                        "HELIX_2022_Dryad_adj", "MESA_Tcells", "FUSION", "packaged_v1")),
                   c("GTEx", "EVAPR", "FHS", "HELIX", "MESA", "FUSION", "packaged_v1"))
  M <- data.table(cpg_id = c("cg1", "cg2", "cg3", "cg3", "cg4", "cg5"),
                  target_gene = "G", tissue = c("blood", "blood", "blood", "blood",
                                                "nasal_epithelium", "solid_tissue"),
                  direction = c(1, -1, 1, 1, 1, -1),
                  n_cohorts = c(2L, 2L, 1L, 1L, 3L, 2L), cohort_agreement = 1,
                  cohorts = c("MESA_monocytes|MESA_Tcells", "HELIX_2022_Dryad_adj|MESA_monocytes",
                              "packaged_v1", "MESA_Tcells",
                              "EVA-PR: All Subjects|EVA-PR: Controls|EVA-PR: Atopic Asthma Only",
                              "GTEx: Lung|GTEx: Ovary"))
  q <- c("cg1", "cg2", "cg3", "cg4", "cg5")
  r <- cpgdirection:::.cpgd_measured_by_tissue(M, q, rep("G", 5))
  expect_identical(r$blood$d, c(1, -1, 1, NA, NA))
  expect_identical(r$blood$rep, c(FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(r$nasal$rep[4], FALSE)
  expect_identical(r$solid$rep[5], FALSE)
  # a table without cohort names falls back to n_cohorts
  r2 <- cpgdirection:::.cpgd_measured_by_tissue(M[, !"cohorts"], q, rep("G", 5))
  expect_identical(r2$blood$rep[1:2], c(TRUE, TRUE))
})

test_that("measured signs: only the queried CpGs are read, and a split sign abstains", {
  f <- cpgdirection:::.cpgd_measured_by_tissue
  M <- data.table(cpg_id = c("cg1", "cg1", "cg2", "cg2", "cg9"),
                  target_gene = "G", tissue = "blood",
                  direction = c(1, -1, -1, -1, 1),
                  n_cohorts = 1L, cohort_agreement = 1,
                  cohorts = c("FUSION", "MESA_Tcells", "FUSION",
                              "HELIX_2022_Dryad_adj", "GTEx: Lung"))
  q <- c("cg1", "cg2", "cg3")
  r <- f(M, q, rep("G", 3))
  expect_identical(r$blood$d, c(NA, -1, NA))
  expect_identical(r$blood$rep, c(FALSE, TRUE, FALSE))
  # records of CpGs nobody asked about change nothing
  junk <- data.table(cpg_id = sprintf("cgX%05d", 1:5000), target_gene = "G",
                     tissue = "blood", direction = 1, n_cohorts = 1L,
                     cohort_agreement = 1, cohorts = "FUSION")
  expect_identical(f(rbind(M, junk), q, rep("G", 3)), r)
  # a query with no measured record returns the empty answer
  r0 <- f(M, c("cg7", "cg8"), rep("G", 2))
  expect_true(all(is.na(unlist(lapply(r0, `[[`, "d")))))
})

test_that("catalogue candidates: DIRECT rows are measured signs; consensus and conflict need one gene", {
  cc <- cpgdirection:::.cpgd_catalogue_cand
  m <- function(...) matrix(c(...), nrow = 1L)
  # two predicted calls agreeing on one gene, a DIRECT call aside
  r <- cc(m(1, 1, -1), m("A", "B", "A"), m("PREDICTED", "PREDICTED", "DIRECT_eQTM"),
          m("G", "G", "G"))
  expect_identical(r$cons, 1); expect_identical(r$worst, "B"); expect_false(r$conflict)
  expect_identical(as.vector(r$direct), c(NA, NA, -1))
  expect_true(all(is.na(r$single)))
  # two predicted calls on different genes: two single calls, no consensus, no conflict
  r <- cc(m(1, -1, NA), m("A", "A", NA), m("PREDICTED", "PREDICTED", NA), m("G", "H", NA))
  expect_true(is.na(r$cons)); expect_false(r$conflict)
  expect_identical(as.vector(r$single), c(1, -1, NA))
  # the same gene, opposite signs: a conflict
  r <- cc(m(1, -1, NA), m("A", "A", NA), m("PREDICTED", "PREDICTED", NA), m("G", "G", NA))
  expect_true(r$conflict); expect_true(is.na(r$cons))
  # one predicted call: a single call, also with min_tissues = 1
  r <- cc(m(NA, -1, NA), m(NA, "A", NA), m(NA, "PREDICTED", NA), NULL, min_tissues = 1L)
  expect_true(is.na(r$cons)); expect_identical(as.vector(r$single), c(NA, -1, NA))
})

test_that("without genes, a CpG-level consensus or conflict never spans two genes", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19"),
              "lookup tables unavailable")
  set.seed(20261106)
  ids <- sample(cpgd_cpg_positions()$cpg_id, 600)
  f <- as.data.frame(suppressWarnings(suppressMessages(cpg_expression_direction(ids, verbose = FALSE))))
  g <- f[, paste0("gene_", c("blood", "nasal", "solid"))]
  d <- f[, paste0("dir_", c("blood", "nasal", "solid"))]
  st <- f[, paste0("status_", c("blood", "nasal", "solid"))]
  pred <- !is.na(d) & st == "PREDICTED"
  ngene <- vapply(seq_len(nrow(f)), function(i) length(unique(toupper(unlist(g[i, pred[i, ]])))), 1L)
  hit <- grepl("^catalogue_consensus", f$best_rung) | f$best_rung %in% "tissue_conflict"
  expect_true(all(ngene[hit] == 1L))
  # a single call names the tissue whose model made it, and that model called
  s1 <- regmatches(f$best_rung, regexpr("(?<=^catalogue_single_)[a-z]+", f$best_rung, perl = TRUE))
  i <- grep("^catalogue_single_", f$best_rung)
  expect_true(all(pred[cbind(i, match(s1, c("blood", "nasal", "solid")))]))
})

test_that("every direction has a weight, in every output; only an abstention has none", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19"),
              "lookup tables unavailable")
  set.seed(20261107)
  ids <- sample(cpgd_cpg_positions()$cpg_id, 300)
  q <- function(x) as.data.frame(suppressWarnings(suppressMessages(x)))
  both <- function(d, dir, w, info) {
    expect_identical(is.na(d[[dir]]), is.na(d[[w]]), info = info)
  }
  for (a in c("lower", "median")) {
    both(q(cpg_expression_direction(ids, verbose = FALSE, accuracy = a)),
         "best_direction", "best_accuracy", paste("ced", a))
    for (t in CPGD_TISSUES) {
      # a declared target fills every row it can; the last resort carries its curve's figure
      f <- q(cpg_expression_direction(ids, target_tissue = t, verbose = FALSE, accuracy = a))
      both(f, "best_direction", "best_accuracy", paste("target", t, a))
      lr <- f$best_evidence %in% "targeted_last_resort"
      expect_true(all(startsWith(f$best_rung[lr], "distance_curve_")), info = t)
      # the single-tissue result: every call has a number (solid tier B: 0.5)
      g <- q(cpg_expression_direction(ids, tissue = t, verbose = FALSE, accuracy = a))
      expect_identical(!g$call %in% "ABSTAIN", !is.na(g$accuracy), info = paste("single", t, a))
    }
  }
})

test_that("min_tissues = 1 never turns a lone model's call into a consensus", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19"),
              "lookup tables unavailable")
  set.seed(20261105)
  ids <- sample(cpgd_cpg_positions()$cpg_id, 400)
  f <- suppressWarnings(suppressMessages(cpg_expression_direction(ids, min_tissues = 1L,
                                                                  verbose = FALSE)))
  one <- f$n_tissues_calling %in% 1L
  expect_gt(sum(one), 0)
  expect_false(any(grepl("^catalogue_consensus", f$best_rung[one])))
})

test_that("SMR keys: the band of the per-pair accuracy; below the floor and off the ladder give none", {
  key <- cpgdirection:::.cpgd_smr_key
  band <- cpgdirection:::.cpgd_smr_band
  expect_identical(band(c(0.999, 0.95, 0.949, 0.90, 0.80, 0.799, 0.65, 0.649, 0.55, 0.549, 0.2, NA)),
                   c("acc_95", "acc_95", "acc_90", "acc_90", "acc_80", "acc_65", "acc_65", "acc_55",
                     "acc_55", NA, NA, NA))
  expect_identical(cpgdirection:::CPGD_SMR_FLOOR, 0.55)
  # the rung does not change the key; a row that is not on the ladder has none
  expect_identical(key(c("smr_replicated", "smr_high", "smr_moderate", "smr_weak", NA, "x"),
                       c(0.97, 0.97, 0.97, 0.60, 0.97, 0.97)),
                   c("acc_95", "acc_95", "acc_95", "acc_55", NA, NA))
  expect_identical(key("smr_moderate", c(0.91, 0.30)), c("acc_90", NA))
  # the bands are ladder sources, and in blood their figures fall with the band
  expect_true(all(paste0("smr_acc_", c(95, 90, 80, 65, 55)) %in% CPGD_LADDER_SOURCES))
  f <- vapply(paste0("smr_acc_", c(95, 90, 80, 65, 55)),
              function(s) cpgdirection:::.cpgd_rung_figure1(s, "blood", "adult")$agree, numeric(1))
  expect_false(is.unsorted(rev(f)))
  # a table without the accuracy column: said once, and nothing is supplied
  env <- cpgdirection:::.cpgd_env
  old <- env$smr_accuracy_warned
  env$smr_accuracy_warned <- FALSE
  on.exit(env$smr_accuracy_warned <- old, add = TRUE)
  S0 <- acc_src(); S0$smr <- data.table::copy(S0$smr)[, smr_accuracy := NULL]
  expect_warning(p0 <- cpg_gene_pairs(c("cg20000001", "cg20000002", "cg20000003"), universal = FALSE,
                                      verbose = FALSE, sources = S0), "predates the per-pair accuracy")
  expect_false(any(startsWith(p0$best_evidence, "smr_")))
  expect_true(all(!is.na(p0$smr_direction)))
  expect_true(all(is.na(p0$smr_accuracy)))
})

test_that("distance keys: the band of |CpG - TSS|; a band below the floor supplies nothing", {
  band <- cpgdirection:::.cpgd_distance_band
  expect_identical(band(c(0, -499, 500, 1499, 1500, 4999, 5000, 9999, 10000, 49999, 50000, 250000, 250001, NA)),
                   c("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4", "b5", "b5", "b6", "b6", "b7", NA))
  expect_identical(cpgdirection:::.cpgd_distance_key("distance_curve_solid", c(100, 60000, NA)),
                   c("distance_curve_solid_b1", "distance_curve_solid_b6", NA))
  expect_identical(cpgdirection:::CPGD_DISTANCE_FLOOR, 0.55)
  sup <- cpgdirection:::.cpgd_distance_supplies
  fig <- cpgdirection:::.cpgd_rung_figure1
  for (t in CPGD_TISSUES) for (k in c("distance_curve_blood_b1", "distance_curve_solid_b6", "distance_only_b7")) {
    f <- fig(k, t, if (t == "blood") "adult" else "any")
    expect_identical(sup(k, t), !is.null(f) && f$agree >= 0.55, info = paste(t, k))
  }
  expect_identical(sup(c(NA, "distance_curve_blood_b1"), "blood"), c(FALSE, TRUE))
  # every band is described, so every emitted key has a quoted string
  keys <- c(paste0("distance_only_b", 1:7),
            paste0("distance_curve_", rep(c("blood", "nasal", "solid"), each = 7), "_b", 1:7))
  expect_true(all(keys %in% names(cpgdirection:::CPGD_RUNG_DESC)))
  expect_false(any(is.na(cpgdirection:::.cpgd_rung_label(keys, "blood", "adult"))))
  # near the TSS the distance source is worth more than far from it
  expect_gt(fig("distance_only_b1", "blood", "adult")$agree, fig("distance_only_b5", "blood", "adult")$agree)
})

# ---- the numbers ---------------------------------------------------------------------

test_that("strings map to numbers; lower floors at 0.5; 'any' takes the smaller", {
  la <- cpgdirection:::.cpgd_rung_label("smr_acc_65", "blood", "adult")
  fa <- cpgdirection:::.cpgd_rung_figure1("smr_acc_65", "blood", "adult")
  expect_equal(num(la, "median"), fa$agree)
  expect_equal(num(la), max(fa$lo, 0.5))
  lk <- cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "child")
  ly <- cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "any")
  expect_equal(num(ly), min(num(cpgdirection:::.cpgd_rung_label("smr_acc_95", "blood", "adult")),
                            num(lk)))
  expect_identical(num(c(NA, "not a quoted string")), rep(NA_real_, 2))
  uq <- cpgdirection:::.cpgd_rung_label("catalogue_model_solid_B", "solid_tissue", "adult")
  expect_match(uq, "^unquantified .*neutral weight 0.5")
  expect_identical(c(num(uq), num(uq, "median")), c(0.5, 0.5))
  expect_identical(num(character(0)), numeric(0))
})

test_that("the accuracy argument and the session option", {
  arg <- cpgdirection:::.cpgd_accuracy_arg
  expect_identical(arg("lower"), "lower")
  expect_identical(arg("median"), "median")
  expect_identical(arg(NULL), "lower")
  expect_error(arg("mean"), "lower")
  expect_error(arg(c("lower", "median")), "lower")
  expect_error(acc_pairs(accuracy = "upper"), "lower")
})

test_that("pairs: best_accuracy is the number of best_expected_accuracy, lower by default", {
  p <- acc_pairs()
  expect_identical(attr(p, "accuracy"), "lower")
  expect_identical(p$best_accuracy, num(p$best_expected_accuracy))
  expect_mapped(p$best_expected_accuracy, p$best_accuracy)
  m <- acc_pairs(accuracy = "median")
  expect_identical(attr(m, "accuracy"), "median")
  expect_identical(m$best_accuracy, num(m$best_expected_accuracy, "median"))
  expect_true(all(m$best_accuracy >= p$best_accuracy, na.rm = TRUE))
  cols <- setdiff(names(p), "best_accuracy")
  expect_identical(names(m), names(p))
  expect_equal(as.data.frame(m)[cols], as.data.frame(p)[cols], ignore_attr = TRUE)
  # the tissue sets the order and the figures
  s <- acc_pairs(tissue = "solid_tissue")
  expect_true(all(grepl("in solid tissue", stats::na.omit(s$best_expected_accuracy), fixed = TRUE)))
})

test_that("options(cpgdirection.accuracy =) sets the default for every call", {
  old <- options(cpgdirection.accuracy = "median")
  on.exit(options(old))
  expect_identical(acc_pairs()$best_accuracy, acc_pairs(accuracy = "median")$best_accuracy)
  expect_identical(attr(acc_pairs(), "accuracy"), "median")
  expect_identical(attr(acc_pairs(accuracy = "lower"), "accuracy"), "lower")
  options(cpgdirection.accuracy = "nonsense")
  expect_error(acc_pairs(), "lower")
})

test_that("age_group picks the child numbers; 'any' the smaller of the two; never the order", {
  a <- acc_pairs()
  k <- acc_pairs(age_group = "child")
  y <- acc_pairs(age_group = "any")
  expect_identical(k$best_rung, a$best_rung)
  expect_identical(k$best_direction, a$best_direction)
  expect_identical(k$best_expected_accuracy,
                   cpgdirection:::.cpgd_rung_label(ifelse(is.na(k$best_direction), NA, k$best_rung),
                                                   "blood", "child"))
  expect_identical(y$best_accuracy, pmin(a$best_accuracy, k$best_accuracy))
})

test_that("the weighting an engine derives stays within [-1, 1] and keeps the sign", {
  p <- acc_pairs()
  ok <- !is.na(p$best_direction)
  p_plus <- ifelse(p$best_direction[ok] > 0, p$best_accuracy[ok], 1 - p$best_accuracy[ok])
  mult <- 2 * p_plus - 1
  expect_equal(mult, p$best_direction[ok] * (2 * p$best_accuracy[ok] - 1))
  expect_true(all(abs(mult) <= 1 & (mult == 0 | sign(mult) == p$best_direction[ok])))
})

# ---- with the packaged data --------------------------------------------------------

test_that("on real data every emitted string has a number, and the chosen source is the best present", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19"),
              "lookup tables unavailable")
  P <- cpgd_cpg_positions()
  set.seed(20261104)
  ids <- c("cg00050692", "cg00000029", "cg06846259", "cg16867657", "cg09809672",
           sample(P$cpg_id, 400))
  q <- function(...) suppressWarnings(suppressMessages(...))
  for (ag in c("adult", "child", "any")) {
    for (tt in list(NULL, "blood", "nasal_epithelium", "solid_tissue")) {
      f <- q(cpg_expression_direction(ids, verbose = FALSE, age_group = ag, target_tissue = tt))
      expect_mapped(f$best_expected_accuracy, f$best_accuracy, info = paste(ag, tt))
      expect_identical(attr(f, "rank_tissue"), if (is.null(tt)) "blood" else tt)
      fm <- q(cpg_expression_direction(ids, verbose = FALSE, age_group = ag, target_tissue = tt,
                                       accuracy = "median"))
      expect_identical(fm$best_accuracy, num(fm$best_expected_accuracy, "median"))
      expect_identical(fm$best_rung, f$best_rung)
    }
  }
  # every answer is a ladder source of the call's tissue (or a declared
  # target's own model or curve); child and adult give the same answers
  o <- order_in("blood")
  f <- q(cpg_expression_direction(ids, verbose = FALSE))
  ans <- f$best_rung[!is.na(f$best_direction)]
  expect_gt(length(ans), 0)
  expect_true(all(ans %in% o))
  fk <- q(cpg_expression_direction(ids, verbose = FALSE, age_group = "child"))
  expect_identical(fk$best_rung, f$best_rung)
  expect_identical(fk$best_direction, f$best_direction)
  # measured pairs named with their gene: the measured rungs answer, by tissue
  if (cpgd_has_data("measured_eqtms")) {
    Ms <- as.data.frame(cpgd_measured_eqtms())
    Ms <- Ms[!duplicated(Ms$cpg_id), ][sample(sum(!duplicated(Ms$cpg_id)), 80), ]
    g <- q(cpg_expression_direction(paste0(Ms$cpg_id, "_", Ms$target_gene), verbose = FALSE))
    ga <- g$best_rung[!is.na(g$best_direction)]
    expect_true(any(startsWith(ga, "measured_")))
    expect_true(all(ga %in% o))
    expect_mapped(g$best_expected_accuracy, g$best_accuracy)
  }
  expect_identical(names(f)[match("best_expected_accuracy", names(f)) + 1L], "best_accuracy")
  expect_identical(names(f)[match("best_evidence", names(f)) + 1L], "best_rung")

  for (t in CPGD_TISSUES) {
    s1 <- q(cpg_expression_direction(ids, tissue = t, verbose = FALSE))
    expect_mapped(s1$expected_accuracy, s1$accuracy, info = t)
    s2 <- q(cpg_expression_direction(ids, tissue = t, verbose = FALSE, accuracy = "median"))
    expect_identical(s2$accuracy, num(s2$expected_accuracy, "median"))
  }
  a <- q(cpg_direction_all_tissues(ids, verbose = FALSE, accuracy = "median"))
  expect_mapped(a$best_expected_accuracy, a$best_accuracy)
  pp <- q(cpg_gene_pairs(ids[1:300], include_onco = TRUE, verbose = FALSE, age_group = "any"))
  expect_mapped(pp$best_expected_accuracy, pp$best_accuracy)
  old <- options(cpgdirection.accuracy = "median"); on.exit(options(old), add = TRUE)
  s3 <- q(cpg_expression_direction(ids, tissue = "blood", verbose = FALSE))
  expect_identical(s3$accuracy, num(s3$expected_accuracy, "median"))
})

test_that("cpg_direction_universal() carries the number", {
  skip_if_not(cpgd_has_data("cpg_positions_hg19"), "positions unavailable")
  P <- cpgd_cpg_positions()
  p <- P[match(c("cg03405789", "cg08215831"), P$cpg_id), ]
  tss <- data.frame(gene = "CRH", chr = p$chr[1], tss = p$pos[1] + 2000L)
  u <- suppressMessages(cpg_direction_universal(c("cg03405789", "cg08215831"),
                                                genes = c("CRH", "CRH"), gene_tss = tss,
                                                verbose = FALSE))
  expect_gt(nrow(u), 0)
  expect_identical(names(u)[match("expected_accuracy", names(u)) + 1L], "accuracy")
  expect_identical(u$accuracy, num(u$expected_accuracy))
  expect_true(all(grepl("distance curves agree", u$expected_accuracy, fixed = TRUE)))
  um <- suppressMessages(cpg_direction_universal(c("cg03405789", "cg08215831"),
                                                 genes = c("CRH", "CRH"), gene_tss = tss,
                                                 tissue = "solid_tissue", verbose = FALSE,
                                                 accuracy = "median"))
  expect_true(all(grepl("solid-tissue distance curve", um$expected_accuracy, fixed = TRUE)))
  expect_identical(um$accuracy, num(um$expected_accuracy, "median"))
  expect_identical(attr(um, "accuracy"), "median")
})

test_that("the bridge: option, both values, age, and the cortex order of the blood SMR fallback", {
  expect_identical(CPGD_BRAIN_SMR_ORDER, c("acc_95", "acc_90", "acc_80", "acc_55", "acc_65"))
  V <- cpgd_accuracy_values()
  cx <- V$median[match(paste0("brain_fallback_blood_SMR_", CPGD_BRAIN_SMR_ORDER), V$rung)]
  expect_false(is.unsorted(rev(cx)))
  skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
  ids <- c("cg06846259_POMC", "cg06846259_CREBBP", "cg00000029_RBL2", "cg16867657_ELOVL2")
  q <- function(...) suppressWarnings(suppressMessages(...))
  b <- q(cpg_brain_bridge(ids))
  expect_mapped(b$requested_gene_best_accuracy, b$requested_gene_accuracy)
  bm <- q(cpg_brain_bridge(ids, accuracy = "median"))
  expect_identical(bm$requested_gene_accuracy, num(bm$requested_gene_best_accuracy, "median"))
  old <- options(cpgdirection.accuracy = "median"); on.exit(options(old), add = TRUE)
  expect_identical(q(cpg_brain_bridge(ids))$requested_gene_accuracy, bm$requested_gene_accuracy)
  options(old)
  for (ag in c("child", "any"))
    expect_identical(q(cpg_brain_bridge(ids, age_group = ag))$requested_gene_accuracy,
                     b$requested_gene_accuracy, info = ag)
})
