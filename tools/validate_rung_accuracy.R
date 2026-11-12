#!/usr/bin/env Rscript
# Every evidence source of the direction ladder, scored in each target tissue
# against ONE yardstick per tissue that is independent of the source it scores
# (4 Nov 2026). These figures order the ladder per tissue (ladder rank =
# point estimate in the call's tissue) and are the numbers best_accuracy
# carries.
#
# Yardsticks (measured cis eQTMs, one sign per pair; a pair's yardstick sign is
# the majority of the yardstick group's cohorts, ties dropped):
#   blood, adult  MESA monocytes + T cells (Wang et al. 2021 S2; n = 1,201 / 214).
#                 Not in the catalogue training cohorts (HELIX, Framingham,
#                 GTEx whole blood, EVA-PR white blood cells), not an SMR arm.
#                 With --blueprint, the BLUEPRINT monocyte, neutrophil and
#                 T-cell eQTMs (Chen et al. 2016; n = 197) are a second
#                 adult-blood fold: a yardstick only, never a measured source.
#   blood, child  HELIX Dryad cell-adjusted (Ruiz-Arenas et al. 2022; 6-11 y,
#                 n = 823), restricted to pairs the blood catalogue does NOT
#                 hold as DIRECT_eQTM, so no catalogue training label is reused.
#   nasal         EVA-PR nasal epithelium (all Atlas subsets), restricted to
#                 pairs the nasal catalogue does not hold as DIRECT_eQTM. The
#                 nasal catalogue was trained on EVA-PR: its figure is
#                 same-study, flagged in the report.
#   solid         four folds pooled, each donor group held out in turn and
#                 scored by the measured solid pairs of the others: FUSION
#                 skeletal muscle, GTEx solid tissues, the osteoarthritis joint
#                 tissues (Kreitmaier 2026) and kidney tubule (Liu 2022); a
#                 pair measured in several counts once (the first fold in that
#                 order). Independent of the TCGA-trained catalogue and tumour
#                 layer.
# Fold rule (PI ruling, 12 Nov 2026): every admitted donor group with at least
# 5,000 pairs in its tissue is a held-out fold there, provided it is
# independent of the sources it scores. Framingham and EVA-PR white blood
# cells pass the size rule in blood but trained the blood catalogue, so they
# are not folds; placenta, adipose and GTP whole blood are below 5,000 pairs:
# they vote as measured sources and are not folds.
# Measured data are what the package ships (cpgd_measured_eqtms()): the Wang
# 2021 liver rows are not shipped and are not used. Cohorts are counted by
# donor group (.cpgd_donor_group(): EVA-PR's strata, MESA's two cell types and
# GTEx's tissues are one group each). A measured source never includes the
# yardstick's own donors: GTEx whole blood is excluded where GTEx is the
# yardstick, EVA-PR white blood cells where EVA-PR nasal is, packaged_v1
# (Framingham + HELIX Atlas) where HELIX is.
#
# Sources (each scored standalone, on every yardstick pair where it calls):
#   measured_<tissue>_<replicated|single>  sign from non-yardstick cohorts of
#     that tissue (replicated: 2+ donor groups, all cohorts agree)
#   smr_acc_<band>     the accuracy bands of .cpgd_smr_key(): pairs whose rung
#     is not `reported`, grouped by their held-out per-pair accuracy
#     (smr_accuracy; 0.55-0.65, 0.65-0.80, 0.80-0.90, 0.90-0.95, 0.95+)
#   catalogue_consensus_<A|B>, catalogue_single_<tissue>_<A|B>  the package's
#     call rule (tier A conf >= 0.30, tier B >= 0.50, every row alike),
#     consensus = 2+ tissue models calling and agreeing, worst tier, formed
#     from PREDICTED calls only, as the ladders form it: a DIRECT_eQTM row is a
#     measured sign from a training cohort, not a prediction, and the ladders
#     treat it as one.
#   catalogue_model_<tissue>_<A|B>  one tissue model's PREDICTED call
#   onco_consensus     Onco-eQTM O1
#   distance_only_<band>, distance_curve_<tissue>_<band>   the three distance
#     curves unanimous (pooled |p - 0.5| * 2 >= 0.02), and each tissue's own
#     curve, at the yardstick's |CpG - TSS|, scored separately in six distance
#     bands (0-0.5, 0.5-1.5, 1.5-5, 5-10, 10-50, 50-250 kb). The blood and
#     nasal curves were fitted on the Atlas Framingham, HELIX and EVA-PR pairs,
#     so distance is scored against MESA and the solid folds only.
# Agreement with Wilson 95% intervals. The yardstick carries its own error, so
# every figure is an agreement with an independent measurement, slightly below
# the probability of being right; all sources in a tissue face the same
# yardstick, so their order is unaffected.
#
# Usage (package root):
#   Rscript tools/validate_rung_accuracy.R --atlas-long <measured_eqtms_atlas_long.csv.gz> [--data-dir ~/cpgd_local]
# Inputs: tools/out/measured_eqtms_add_helix.csv.gz, tools/out/measured_eqtms_add_wang2021.csv.gz,
#   tools/out/measured_eqtms_add_kreitmaier2026.csv.gz (osteoarthritis cartilage and synovium, one donor group),
#   tools/out/measured_eqtms_add_nov2026_sources.csv.gz (kidney, placenta, adipose, GTP whole blood),
#   optionally --blueprint tools/out/blueprint_eqtm_pooled.csv.gz (tools/build_blueprint_yardstick.py: cell,
#   cpg_id, target_gene, direction, tss_dist; one row per cell type and pair),
#   the Atlas long table written by inst/scripts/build_measured_eqtms.R, package data.
# Writes tools/reports/validate_rung_accuracy.md and tools/out/rung_accuracy.csv.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = path.expand(data_dir))
atlas_long <- path.expand(get_arg("--atlas-long", "tools/out/measured_eqtms_atlas_long.csv.gz"))
add_helix  <- get_arg("--helix", "tools/out/measured_eqtms_add_helix.csv.gz")
add_wang   <- get_arg("--wang", "tools/out/measured_eqtms_add_wang2021.csv.gz")
add_oa     <- get_arg("--kreitmaier", "tools/out/measured_eqtms_add_kreitmaier2026.csv.gz")
add_nov    <- get_arg("--nov-sources", "tools/out/measured_eqtms_add_nov2026_sources.csv.gz")
add_bp     <- get_arg("--blueprint")   # BLUEPRINT eQTM pairs: an adult-blood yardstick fold only
out_md     <- get_arg("--out-md", "tools/reports/validate_rung_accuracy.md")
out_csv    <- get_arg("--out-csv", "tools/out/rung_accuracy.csv")
out_r      <- get_arg("--write-r")   # e.g. R/rung_accuracy_table.R: the table the package ranks by
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)),
          file.exists(atlas_long), file.exists(add_helix), file.exists(add_wang), file.exists(add_oa),
          file.exists(add_nov), is.null(add_bp) || file.exists(add_bp))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
rd <- function(f) cpgdirection:::.cpgd_fread(f)
gk <- function(x) gsub("_", "-", toupper(trimws(x)))
wilson <- function(k, n, z = stats::qnorm(0.975)) {
  p <- k / n; d <- 1 + z^2 / n
  c0 <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(p = p, lo = pmax(0, c0 - h), hi = pmin(1, c0 + h))
}

# ---- per-cohort measured signs, with donor groups ---------------------------
A <- rd(atlas_long)
A <- A[!is.na(direction) & direction != 0,
       .(cohort, tissue, cpg_id, g = gk(target_gene), d = sign(as.numeric(direction)),
         tss_dist = as.numeric(tss_dist))]
A[, grp := cpgdirection:::.cpgd_donor_group(cohort)]
H <- rd(add_helix)[, .(cohort = "HELIX_2022_Dryad_adj", tissue = "blood", cpg_id,
                      g = gk(target_gene), d = sign(as.numeric(direction)),
                      tss_dist = as.numeric(tss_dist), grp = "HELIX")]
W <- rd(add_wang)
# the package ships the MESA blood rows of this file, not the liver rows
W <- W[tissue == "blood"]
W <- W[, .(cohort = cohorts, tissue, cpg_id, g = gk(target_gene), d = sign(as.numeric(direction)),
           tss_dist = as.numeric(tss_dist),
           grp = cpgdirection:::.cpgd_donor_group(sub("[|].*$", "", cohorts)))]
# osteoarthritis cartilage and synovium (Kreitmaier 2026): one row per cohort; where both cohorts measured a
# pair the add file keeps it only when they agree, so each cohort carries the pair's sign
K <- rd(add_oa)
K <- K[, .(cohort = unlist(strsplit(cohorts, "|", fixed = TRUE))),
       by = .(tissue, cpg_id, g = gk(target_gene), d = sign(as.numeric(direction)), tss_dist = as.numeric(tss_dist))]
K[, grp := cpgdirection:::.cpgd_donor_group(cohort)]
stopifnot(all(K$grp == "Kreitmaier2026"), all(K$tissue == "solid_tissue"))
Mp <- as.data.table(cpgd_measured_eqtms())
Mp <- Mp[source == "packaged_v1", .(cohort = "packaged_v1", tissue, cpg_id, g = gk(target_gene),
                                   d = sign(as.numeric(direction)),
                                   tss_dist = as.numeric(tss_dist), grp = "packaged_v1")]
# kidney, placenta, adipose, GTP whole blood (12 Nov 2026): one row per cohort and pair, each its own donor group
N <- rd(add_nov)
N <- N[, .(cohort = cohorts, tissue, cpg_id, g = gk(target_gene), d = sign(as.numeric(direction)),
           tss_dist = as.numeric(tss_dist), grp = cpgdirection:::.cpgd_donor_group(cohorts))]
stopifnot(all(N$grp == N$cohort))
C <- unique(rbindlist(list(A, H, W, K, N, Mp), use.names = TRUE)[!is.na(d) & d != 0],
            by = c("cohort", "cpg_id", "g"))
# BLUEPRINT: never a measured source (grp "BLUEPRINT" is excluded from every fold's sources)
BP <- NULL
if (!is.null(add_bp)) {
  # one row per cell type and pair (tools/build_blueprint_yardstick.py); the yardstick sign of a pair is the
  # majority over the cell types that hold it, ties dropped (collapse())
  BP <- rd(add_bp)
  BP <- unique(BP[, .(cohort = paste0("BLUEPRINT_", cell), tissue = "blood", cpg_id, g = gk(target_gene),
                      d = sign(as.numeric(direction)), tss_dist = as.numeric(tss_dist), grp = "BLUEPRINT")][
    !is.na(d) & d != 0], by = c("cohort", "cpg_id", "g"))
  BP[, pk := paste(cpg_id, g)]
}
C[, pk := paste(cpg_id, g)]
cat(sprintf("measured cohort x pair rows: %s\n", format(nrow(C), big.mark = ",")))

# majority sign of a set of cohort rows per pair (ties dropped), with the
# number of cohorts and whether they all agree
collapse <- function(X) {
  X[, .(d = { s <- sum(d); if (s > 0) 1 else if (s < 0) -1 else NA_real_ },
        nco = .N, unanimous = length(unique(d)) == 1L,
        tss_dist = stats::median(tss_dist, na.rm = TRUE)), by = .(pk, cpg_id, g)][!is.na(d)]
}

# ---- yardstick folds ----------------------------------------------------------
L <- lapply(setNames(CPGD_TISSUES, CPGD_TISSUES), function(t)
  as.data.table(cpgdirection:::.cpgd_lookup(t)))
direct <- lapply(L, function(x) unique(paste(x$cpg_id[x$status == "DIRECT_eQTM"],
                                             gk(x$target_gene[x$status == "DIRECT_eQTM"]))))
folds <- list(
  list(target = "blood", age = "adult", name = "MESA", y = C[grp == "MESA" & tissue == "blood"],
       exclude = c("MESA"), drop_direct = NULL, distance = TRUE),
  list(target = "blood", age = "child", name = "HELIX (non-catalogue pairs)",
       y = C[grp == "HELIX"], exclude = c("HELIX", "packaged_v1"), drop_direct = "blood",
       distance = FALSE),
  list(target = "nasal_epithelium", age = "any", name = "EVA-PR nasal (non-catalogue pairs)",
       y = C[grp == "EVAPR" & tissue == "nasal_epithelium"], exclude = c("EVAPR"),
       drop_direct = "nasal_epithelium", distance = FALSE),
  list(target = "solid_tissue", age = "any", name = "FUSION", y = C[grp == "FUSION"],
       exclude = c("FUSION"), drop_direct = NULL, distance = TRUE),
  list(target = "solid_tissue", age = "any", name = "GTEx solid",
       y = C[grp == "GTEx" & tissue == "solid_tissue"], exclude = c("GTEx"), drop_direct = NULL,
       distance = TRUE),
  list(target = "solid_tissue", age = "any", name = "osteoarthritis joint tissues", y = C[grp == "Kreitmaier2026"],
       exclude = c("Kreitmaier2026"), drop_direct = NULL, distance = TRUE),
  list(target = "solid_tissue", age = "any", name = "kidney", y = C[grp == "Liu2022_kidney"],
       exclude = c("Liu2022_kidney"), drop_direct = NULL, distance = TRUE))
if (!is.null(BP)) folds <- append(folds, list(
  list(target = "blood", age = "adult", name = "BLUEPRINT", y = BP,
       exclude = c("BLUEPRINT"), drop_direct = NULL, distance = TRUE)), after = 1L)
for (f in folds) stopifnot(nrow(f$y) >= 5000L)

# ---- the non-measured sources, keyed by pair ---------------------------------
S <- cpgdirection:::.cpgd_smr_table()
SS <- data.table(pk = paste(S$cpg_id, gk(S$target_gene)), d = as.numeric(S$direction),
                 tier = as.character(S$smr_tier), p = as.numeric(S$p_SMR),
                 rung = cpgdirection:::.cpgd_smr_rung(S),
                 acc = cpgdirection:::.cpgd_smr_accuracy(S))
setorderv(SS, c("pk", "tier", "p")); SS <- SS[, .SD[1L], by = pk]
SS[, src := fcase(rung == "replicated", "smr_replicated", rung == "high", "smr_high",
                  rung == "moderate", "smr_moderate", rung == "weak", "smr_weak",
                  default = NA_character_)]
SS[, skey := cpgdirection:::.cpgd_smr_key(src, acc)]
SS <- SS[!is.na(skey), .(pk, d, source = paste0("smr_", skey))]

O <- cpgdirection:::.cpgd_onco_keyed(NULL)
OO <- if (is.null(O)) NULL else
  O[tier %in% "O1" & !is.na(d), .(pk = paste(cpg_id, .gkey), d, source = "onco_consensus")]

short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")
catcall <- function(t) {
  x <- L[[t]]
  x <- data.table(pk = paste(x$cpg_id, gk(x$target_gene)), d = as.numeric(x$direction),
                  cf = as.numeric(x$confidence), ti = as.character(x$evidence_tier),
                  st = as.character(x$status))
  setorderv(x, c("pk", "cf"), order = c(1L, -1L)); x <- x[, .SD[1L], by = pk]
  # the package's call rule (.cpgd_one_tissue() and the pair resolver)
  x[, call := fcase(is.na(d), NA_real_,
                    ti == "A" & !is.na(cf) & cf >= 0.30, d,
                    ti == "B" & !is.na(cf) & cf >= 0.50, d, default = NA_real_)]
  x[!is.na(call), .(pk, call, ti, st, tissue = t)]
}
CC <- rbindlist(lapply(CPGD_TISSUES, catcall))
# the ladders form consensus and single calls from predicted calls only (a
# DIRECT_eQTM call is a measured sign: .cpgd_catalogue_cand())
CAT <- CC[st == "PREDICTED", .(n_call = .N, n_pos = sum(call == 1), n_neg = sum(call == -1),
              worst = max(ti), single_tissue = if (.N == 1L) tissue[1] else NA_character_,
              single_tier = if (.N == 1L) ti[1] else NA_character_,
              single_d = if (.N == 1L) call[1] else NA_real_), by = pk]
CAT[, source := fcase(n_call >= 2L & (n_pos == n_call | n_neg == n_call),
                      paste0("catalogue_consensus_", worst),
                      n_call == 1L, paste0("catalogue_single_", short[single_tissue], "_", single_tier),
                      default = NA_character_)]
CAT[, d := fifelse(n_call == 1L, single_d, fifelse(n_pos > n_neg, 1, -1))]
CAT <- CAT[!is.na(source) & source %like% "_[AB]$", .(pk, d, source)]
# each tissue model's own PREDICTED call, whatever the others say (the
# single-tissue result labels these; its DIRECT rows are measured signs)
MOD <- CC[st == "PREDICTED",
          .(pk, d = call, source = paste0("catalogue_model_", short[tissue], "_", ti))][
  source %like% "_[AB]$"]

K <- fread(system.file("extdata", "distance_curves.csv", package = "cpgdirection"))
dist_call <- function(dist) {
  ad <- log10(pmax(abs(dist), 1))
  P3 <- sapply(c("blood", "nasal_epithelium", "solid_tissue"), function(t) {
    k <- K[tissue == t]
    stats::approx(log10(pmax(k$dist_mid, 1)), k$p_positive, xout = ad, rule = 2)$y
  })
  P3 <- matrix(P3, ncol = 3)
  npos <- rowSums(P3 >= 0.5); pm <- rowMeans(P3)
  list(unanimous = ifelse((npos == 3 | npos == 0) & abs(pm - 0.5) * 2 >= 0.02,
                          ifelse(npos == 3, 1, -1), NA_real_),
       curve = ifelse(P3 >= 0.5, 1, -1))
}

fixed_sources <- rbindlist(list(SS, OO, CAT, MOD), use.names = TRUE)

# ---- score every source on every fold ----------------------------------------
score <- function(f) {
  Y <- collapse(f$y)
  if (!is.null(f$drop_direct)) Y <- Y[!(pk %in% direct[[f$drop_direct]])]
  setnames(Y, "d", "y")
  # measured sources: non-yardstick, non-overlapping cohorts, by tissue
  R <- C[!(grp %in% f$exclude) & pk %chin% Y$pk]
  if (f$target == "solid_tissue" && "GTEx" %in% f$exclude) R <- R[grp != "GTEx"]
  if (f$target == "nasal_epithelium") R <- R[grp != "EVAPR"]
  MS <- R[, .(d = { s <- sum(d); if (s > 0) 1 else if (s < 0) -1 else NA_real_ },
              rep = uniqueN(grp) >= 2L & length(unique(d)) == 1L), by = .(pk, tissue)][!is.na(d)]
  MS[, source := paste0("measured_", short[tissue], "_", fifelse(rep, "replicated", "single"))]
  DS <- NULL
  if (isTRUE(f$distance)) {
    # by distance band (.cpgd_distance_band(): 0-0.5, 0.5-1.5, 1.5-5, 5-10, 10-50, 50-250 kb)
    dc <- dist_call(Y$tss_dist)
    bd <- cpgdirection:::.cpgd_distance_band(Y$tss_dist)
    DS <- rbindlist(c(
      list(data.table(pk = Y$pk, d = dc$unanimous, source = paste0("distance_only_", bd))),
      lapply(1:3, function(j) data.table(pk = Y$pk, d = dc$curve[, j],
                                        source = paste0("distance_curve_", short[[j]], "_", bd)))))[!is.na(d)]
  }
  allsrc <- rbindlist(list(fixed_sources[pk %chin% Y$pk], MS[, .(pk, d, source)], DS),
                      use.names = TRUE)
  m <- merge(allsrc, Y[, .(pk, cpg_id, y)], by = "pk")
  m[, `:=`(target = f$target, age = f$age, yardstick = f$name)]
}
pairs <- rbindlist(lapply(folds, score))
res <- pairs[, .(n = .N, k = sum(d == y), n_cpg = uniqueN(cpg_id)),
             by = .(target, age, yardstick, source)]
w <- wilson(res$k, res$n); res[, `:=`(agree = w$p, lo = w$lo, hi = w$hi)]

# pooled per target tissue and age (solid: both folds, each pair once, the
# first fold's yardstick sign)
pooled_pairs <- unique(pairs, by = c("target", "age", "source", "pk"))
pool <- pooled_pairs[, .(n = .N, k = sum(d == y), n_cpg = uniqueN(cpg_id)),
                     by = .(target, age, source)]
ynames <- rbindlist(lapply(folds, function(f) data.table(target = f$target, age = f$age,
                                                         name = f$name)))[
  , .(yardstick = paste(unique(name), collapse = " + ")), by = .(target, age)]
pool <- merge(pool, ynames, by = c("target", "age"), sort = FALSE)
w <- wilson(pool$k, pool$n); pool[, `:=`(agree = w$p, lo = w$lo, hi = w$hi)]
setorderv(pool, c("target", "age", "agree"), order = c(1L, 1L, -1L))
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(rbind(res[, level := "fold"], pool[, level := "pooled"], fill = TRUE), out_csv)

# ---- report -------------------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
md <- c("# Every ladder source, scored in each target tissue against an independent yardstick",
        "",
        sprintf("Run %s. Agreement with the yardstick's measured sign, Wilson 95%% CI. Sources are scored standalone, on every yardstick pair where they call (catalogue sources: formed from predicted calls only, as the ladders form them). `n_cpg`: distinct CpGs.",
                format(Sys.time(), "%Y-%m-%d %H:%M")),
        "",
        "Yardsticks: adult blood = MESA monocytes + T cells; child blood = HELIX Dryad pairs outside the blood catalogue's DIRECT set; nasal = EVA-PR nasal pairs outside the nasal catalogue's DIRECT set (the nasal catalogue is trained on EVA-PR: its rows are same-study); solid = FUSION, GTEx solid, the osteoarthritis joint tissues and kidney, each held out in turn, pooled, each pair once (with --blueprint, adult blood pools MESA and BLUEPRINT the same way). Measured cohorts are counted by donor group. SMR pairs are scored by accuracy band (rung not `reported`). Distance is not scored against HELIX or EVA-PR (its blood and nasal curves were fitted on them).",
        "")
for (tg in unique(pool$target)) for (ag in unique(pool$age[pool$target == tg])) {
  x <- pool[target == tg & age == ag]
  md <- c(md, sprintf("## %s (%s): %s", tg, ag, x$yardstick[1]), "",
          "| rank | source | n | n_cpg | agree | 95% CI |", "|---|---|---|---|---|---|",
          sprintf("| %d | %s | %s | %s | %s | [%s, %s] |", seq_len(nrow(x)), x$source,
                  format(x$n, big.mark = ","), format(x$n_cpg, big.mark = ","),
                  f3(x$agree), f3(x$lo), f3(x$hi)), "")
}
md <- c(md, "## By fold", "", "| target | age | yardstick | source | n | agree | 95% CI |",
        "|---|---|---|---|---|---|---|",
        res[order(target, age, yardstick, -agree),
            sprintf("| %s | %s | %s | %s | %s | %s | [%s, %s] |", target, age, yardstick, source,
                    format(n, big.mark = ","), f3(agree), f3(lo), f3(hi))])
dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
writeLines(md, out_md)
cat(sprintf("wrote %s and %s\n", out_md, out_csv))

# ---- the package table ---------------------------------------------------------
if (!is.null(out_r)) {
  x <- pool[order(target, age, source)]
  # one element per line for text (a wrap must never fall inside a string);
  # numbers packed, broken only between elements
  q <- function(v) paste0('"', v, '"')
  num <- function(v, d = 4) sprintf(paste0("%.", d, "f"), v)
  int <- function(v) sprintf("%dL", as.integer(v))
  wrap <- function(el, per = 1L) {
    g <- split(el, ceiling(seq_along(el) / per))
    paste0("    ", vapply(g, paste, character(1), collapse = ", "), collapse = ",\n")
  }
  body <- c(
    "# Generated by tools/validate_rung_accuracy.R --write-r; do not edit by hand.",
    sprintf("# %s. Every ladder source scored in each target tissue against an", format(Sys.Date())),
    "# independent measured yardstick (see tools/reports/validate_rung_accuracy.md).",
    "# Pooled rows only: target x age x source, Wilson 95% interval.",
    "CPGD_RUNG_ACCURACY <- data.frame(",
    paste0("  source = c(\n", wrap(q(x$source)), "),"),
    paste0("  target = c(\n", wrap(q(x$target), 4L), "),"),
    paste0("  age = c(\n", wrap(q(x$age), 6L), "),"),
    paste0("  yardstick = c(\n", wrap(q(x$yardstick)), "),"),
    paste0("  n = c(\n", wrap(int(x$n), 8L), "),"),
    paste0("  n_cpg = c(\n", wrap(int(x$n_cpg), 8L), "),"),
    paste0("  k = c(\n", wrap(int(x$k), 8L), "),"),
    paste0("  agree = c(\n", wrap(num(x$agree), 8L), "),"),
    paste0("  lo = c(\n", wrap(num(x$lo), 8L), "),"),
    paste0("  hi = c(\n", wrap(num(x$hi), 8L), "),"),
    "  stringsAsFactors = FALSE)")
  writeLines(body, out_r)
  cat(sprintf("wrote %s (%d rows)\n", out_r, nrow(x)))
}
