# The ranked ladder and its numbers (1.1.0). Every evidence source is scored
# in each tissue against an independent yardstick (CPGD_RUNG_ACCURACY,
# tools/validate_rung_accuracy.R); the ladder order in a tissue and the number
# a call carries both come from that table.

library(data.table)

acc_smr <- function(cpg, gene, direction, rung, tier = "S1", child = "sig",
                    outcome = "eQTLGen", replication = "concordant") {
  data.table(cpg_id = cpg, target_gene = gene, direction = direction,
             smr_tier = tier, p_SMR = 1e-10, n_instruments = 3L,
             instrument_agreement = 1, cpg_gene_dist = 1000,
             top_instrument = "rs1", p_HEIDI = 0.5, nsnp_HEIDI = 3L,
             heidi_status = "pass", outcome = outcome,
             replication = replication, smr_rung = rung,
             child_blood_helix = child)
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
# every quoted string carries a number (only a declared target's last-resort
# sign has none, by design), and nothing without a string has one
expect_mapped <- function(str, val, info = NULL) {
  str <- as.character(str)
  quoted <- !is.na(str) & !grepl("^unquantified", str)
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
  # solid-tissue model's tier-B single call (no yardstick pair)
  expect_true(all(CPGD_LADDER_SOURCES %in% names(CPGD_RUNG_DESC)))
  expect_identical(setdiff(CPGD_LADDER_SOURCES, R$source),
                   c("measured_nasal_replicated", "catalogue_single_solid_B"))
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
  expect_true(all(grepl("^unquantified", V$label[!q])))
  expect_true("targeted_last_resort" %in% V$rung[!q])
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
    # in: every source with a figure whose point is >= 0.5; out: the rest
    for (s in CPGD_LADDER_SOURCES) {
      f <- fig(s, t, age)
      expect_identical(s %in% o, !is.null(f) && f$agree >= 0.5, info = paste(t, s))
    }
  }
  ob <- order_in("blood")
  # a small sample does not outrank a large one on a noisy point: the blood
  # model's tier-B single calls (12 pairs, 0.83 [0.55, 0.95]) sit below the
  # catalogue consensus (534 pairs, 0.80 [0.76, 0.83])
  expect_gt(match("catalogue_single_blood_B", ob), match("catalogue_consensus_A", ob))
  expect_lt(match("smr_replicated_S1", ob), match("measured_blood_single", ob))   # 0.89 vs 0.87
  expect_false("catalogue_single_solid_A" %in% ob)                               # 0.39 in blood
  expect_lt(match("distance_only", ob), match("catalogue_single_nasal_A", ob))
  # nasal epithelium has one donor group (EVA-PR): no nasal sign is replicated
  for (t in CPGD_TISSUES) expect_false("measured_nasal_replicated" %in% order_in(t), info = t)
  # the order takes no age: age_group changes numbers, never rungs (tested on
  # the pairs below and on real data)
  expect_identical(names(formals(cpgdirection:::.cpgd_rung_order)), c("target", "pool"))
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
  expect_identical(cpgdirection:::.cpgd_rung_evidence(
    c("measured_solid_single", "smr_replicated_S2", "smr_moderate_finngen", "smr_weak_discordant",
      "catalogue_single_nasal_B", "catalogue_consensus_A", "onco_consensus", "no_evidence")),
    c("measured", "smr_replicated", "smr_moderate", "smr_weak", "catalogue_single",
      "catalogue_consensus", "onco_consensus", "no_evidence"))
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

test_that("SMR keys: the discordant-cell-type stratum is its own source", {
  key <- cpgdirection:::.cpgd_smr_key
  expect_identical(key(rep("smr_moderate", 4), c("S2", "S2", "S3", "S2"),
                       c("eQTLGen", "eQTLGen", "eQTLGen", "FinnGen_PBMC"),
                       c("eqtlgen_only", "discordant_celltype", "concordant", "finngen_only")),
                   c("moderate", "moderate_discordant_celltype", "moderate_s3_replicated",
                     "moderate_finngen"))
  expect_identical(key(c("smr_replicated", "smr_replicated", "smr_high", "smr_weak", "smr_weak", "x"),
                       c("S1", "S2", "S1", "S1", "S3", "S1"), "eQTLGen",
                       c("concordant", "concordant", "eqtlgen_only", "discordant", "eqtlgen_only", "x")),
                   c("replicated_S1", "replicated_S2", "high", "weak_discordant", "weak", NA))
  f <- cpgdirection:::.cpgd_rung_figure1("smr_moderate_discordant_celltype", "blood", "adult")
  m <- cpgdirection:::.cpgd_rung_figure1("smr_moderate", "blood", "adult")
  expect_lt(f$agree, m$agree)
})

# ---- the numbers ---------------------------------------------------------------------

test_that("strings map to numbers; lower floors at 0.5; 'any' takes the smaller", {
  la <- cpgdirection:::.cpgd_rung_label("smr_weak", "blood", "adult")
  fa <- cpgdirection:::.cpgd_rung_figure1("smr_weak", "blood", "adult")
  expect_equal(num(la, "median"), fa$agree)
  expect_equal(num(la), max(fa$lo, 0.5))
  lk <- cpgdirection:::.cpgd_rung_label("smr_replicated_S1", "blood", "child")
  ly <- cpgdirection:::.cpgd_rung_label("smr_replicated_S1", "blood", "any")
  expect_equal(num(ly), min(num(cpgdirection:::.cpgd_rung_label("smr_replicated_S1", "blood", "adult")),
                            num(lk)))
  expect_identical(num(c(NA, "not a quoted string", CPGD_LAST_RESORT_ACCURACY)), rep(NA_real_, 3))
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
  expect_identical(CPGD_BRAIN_SMR_ORDER, c("replicated", "high", "weak", "moderate"))
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
