#!/usr/bin/env Rscript
# Same-person overlap between the brain SMR layer and its cortex references (4 Nov 2026).
#
# Brain-mMeta (exposure) pools ROSMAP (n=468), Hannon et al. (n=166) and Jaffe et al. (n=526) cortex mQTLs;
# BrainMeta v2 (outcome) holds ROSMAP (n=832) and LIBD (n=128) among its 2,865 samples. Every cortex eQTM
# reference the package has (ROSMAP aged DLPFC; Jaffe DLPFC) therefore shares donors with an SMR arm, so a
# direct agreement figure can be inflated by shared sampling noise. This script estimates the brain layer's
# accuracy WITHOUT using any comparison between donor-sharing datasets:
#
#   1. ratio_T = rho(T, X) / rho(ROSMAP, X) on the SAME pairs, pooled over comparators X that share no donors
#      with the brain arms, with ROSMAP or with each other (measured blood without EVA-PR/GTEx; EVA-PR nasal;
#      solid tissue without GTEx). Under conditional independence this is rho_T / rho_ROSMAP, whatever X's own
#      accuracy or its tissue's concordance with cortex. rho = 2 * agreement - 1.
#   2. ROSMAP's own sign accuracy from ROSMAP vs Jaffe (no shared donors) and the same calibration for Jaffe:
#      rho_R^2 = rho(R, J) / (rho_J / rho_R).
#   3. accuracy vs cortex truth = (1 + ratio * rho_R) / 2; agreement expected against an INDEPENDENT cortex cohort
#      of ROSMAP's quality = (1 + ratio * rho_R^2) / 2. Compare with the in-sample agreement with ROSMAP.
#   4. log-odds weights for the brain-target ladder and a check that a weighted vote would not beat it.
# Cluster bootstrap by CpG for every interval. Read-only; changes nothing in the package.
#
# Usage (package root):
#   Rscript tools/validate_brain_overlap.R [--data-dir ~/cpgd_local] [--boot 1000] [--jaffe-p 0.01]
# Inputs: tools/out/brain_measured_eqtms_wang2021.csv.gz, tools/out/brain_dev_eqtm_ageadj.csv.gz, package data.
# Writes tools/reports/validate_brain_overlap.md and tools/out/brain_overlap_estimates.csv.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
data_dir <- get_arg("--data-dir"); if (!is.null(data_dir)) options(cpgdirection.data_dir = path.expand(data_dir))
nboot <- as.integer(get_arg("--boot", "1000")); jp <- as.numeric(get_arg("--jaffe-p", "0.01"))
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
rd <- function(f) cpgdirection:::.cpgd_fread(f)
gk <- function(x) gsub("_", "-", toupper(trimws(x)))
key <- function(c, g) paste(c, gk(g))
set.seed(20261104)

# ---- per-pair raters ------------------------------------------------------------------------------
B <- as.data.table(cpgd_brain_directions())[tier %in% c("T1", "T2")]
B <- unique(B[, .(k = key(cpg_id, target_gene), A = as.integer(sign(direction)), tier, p_SMR)], by = "k")
W <- rd("tools/out/brain_measured_eqtms_wang2021.csv.gz")
R <- W[group == "brain_aged", .(k = key(cpg_id, target_gene), R = as.integer(sign(direction)), FDR)][order(FDR)]
R <- unique(R, by = "k")[, FDR := NULL]
J <- rd("tools/out/brain_dev_eqtm_ageadj.csv.gz")[set == "fetal_to_25"]
J <- J[, .(k = key(cpg_id, target_gene), J = as.integer(sign(beta_adj)), Jp = p_adj)][order(Jp)]
J <- unique(J, by = "k")[Jp < jp][, Jp := NULL]
M <- as.data.table(cpgd_measured_eqtms())
mm <- function(t, keep, nm) { x <- M[tissue == t & keep(cohorts), .(k = key(cpg_id, target_gene), d = as.integer(sign(direction)))]
  x <- unique(x, by = "k"); setnames(x, "d", nm); x }
everything <- function(z) rep(TRUE, length(z))
Bl  <- mm("blood", everything, "Bl")
BlI <- mm("blood", function(z) !grepl("EVA-PR|GTEx", z), "BlI")
Na  <- mm("nasal_epithelium", everything, "Na")
SoI <- mm("solid_tissue", function(z) !grepl("GTEx", z), "SoI")
S <- as.data.table(tryCatch(cpgdirection:::.cpgd_smr_table(), error = function(e) cpgd_smr_directions()))
S <- S[smr_rung %in% c("replicated", "high", "moderate", "weak")]
S <- unique(S[, .(k = key(cpg_id, target_gene), Sb = as.integer(sign(direction)), hi = smr_rung %in% c("replicated", "high"))], by = "k")
U <- Reduce(function(a, b) merge(a, b, by = "k", all = TRUE), list(B, R, J, Bl, BlI, Na, SoI, S))
U[, `:=`(A1 = fifelse(tier == "T1", A, NA_integer_), A2 = fifelse(tier == "T2", A, NA_integer_),
         SbH = fifelse(hi %in% TRUE, Sb, NA_integer_), SbL = fifelse(hi %in% FALSE, Sb, NA_integer_), cpg = sub(" .*", "", k))]
rat <- c("A1", "A2", "R", "J", "Bl", "BlI", "Na", "SoI", "SbH", "SbL")
U <- U[rowSums(!is.na(as.matrix(U[, ..rat]))) >= 2]
message(sprintf("pairs with >= 2 raters: %s", format(nrow(U), big.mark = ",")))

# ---- estimators ---------------------------------------------------------------------------------
Xs <- c("BlI", "Na", "SoI")
through <- function(D, a, b) { num <- 0; den <- 0
  for (X in Xs) { i <- !is.na(D[[a]]) & !is.na(D[[b]]) & !is.na(D[[X]]); if (!any(i)) next
    num <- num + sum(2 * (D[[a]][i] == D[[X]][i]) - 1); den <- den + sum(2 * (D[[b]][i] == D[[X]][i]) - 1) }
  num / den }
rr <- function(D, a, b) { i <- !is.na(D[[a]]) & !is.na(D[[b]]); 2 * mean(D[[a]][i] == D[[b]][i]) - 1 }
est <- function(D) {
  rJR <- through(D, "J", "R"); rR2 <- min(1, rr(D, "R", "J") / rJR); rR <- sqrt(rR2)
  t1 <- through(D, "A1", "R"); t2 <- through(D, "A2", "R")
  c(ROSMAP_sign_accuracy = (1 + rR) / 2, Jaffe_sign_accuracy = (1 + min(1, rR * rJR)) / 2,
    T1_ratio_to_ROSMAP = t1, T2_ratio_to_ROSMAP = t2,
    T1_accuracy = (1 + min(1, t1 * rR)) / 2, T2_accuracy = (1 + min(1, t2 * rR)) / 2,
    T1_vs_independent_cortex = (1 + min(1, t1 * rR2)) / 2, T2_vs_independent_cortex = (1 + min(1, t2 * rR2)) / 2,
    T1_vs_ROSMAP_insample = (1 + rr(D, "A1", "R")) / 2, T2_vs_ROSMAP_insample = (1 + rr(D, "A2", "R")) / 2,
    measured_blood_accuracy = (1 + min(1, rr(D, "Bl", "R") / rR)) / 2, measured_blood_vs_ROSMAP = (1 + rr(D, "Bl", "R")) / 2,
    blood_SMR_hi_accuracy = (1 + rr(D, "SbH", "R") / rR) / 2, blood_SMR_lo_accuracy = (1 + rr(D, "SbL", "R") / rR) / 2) }
pt <- est(U)
cp <- unique(U$cpg); idx <- split(seq_len(nrow(U)), U$cpg)
bs <- t(replicate(nboot, est(U[unlist(idx[sample(cp, length(cp), TRUE)], use.names = FALSE)])))
ci <- apply(bs, 2, quantile, c(.025, .975), na.rm = TRUE)
E <- data.table(quantity = names(pt), estimate = round(pt, 3), lo = round(ci[1, ], 3), hi = round(ci[2, ], 3))
print(E)

# ---- dose-response: overlap bias would concentrate where the SMR instrument is weak -------------
dr <- rbindlist(lapply(c("A1", "A2"), function(a) { d <- U[!is.na(get(a)) & !is.na(R) & !is.na(p_SMR)]
  br <- unique(quantile(d$p_SMR, c(0, 1/3, 2/3, 1))); d[, ps := cut(p_SMR, br, include.lowest = TRUE, labels = FALSE)]
  rbindlist(lapply(sort(unique(d$ps)), function(q) { e <- d[d$ps == q]
    data.table(tier = sub("A", "T", a), p_SMR_tertile = q, n = nrow(e), insample_vs_ROSMAP = round(mean(e[[a]] == e$R), 3),
               ratio_via_independent = round(through(e, a, "R"), 3)) })) }))
print(dr)

# ---- weights and the ladder ---------------------------------------------------------------------
lg <- function(p) log(p / (1 - p))
w <- c(measured_brain = lg(pt[["ROSMAP_sign_accuracy"]]), brain_T1 = lg(pt[["T1_accuracy"]]), brain_T2 = lg(pt[["T2_accuracy"]]),
       measured_blood = lg(pt[["measured_blood_accuracy"]]), blood_SMR_hi = lg(pt[["blood_SMR_hi_accuracy"]]),
       blood_SMR_lo = lg(pt[["blood_SMR_lo_accuracy"]]))
Wt <- data.table(source = names(w), log_odds_weight = round(w, 2))
cmp <- rbindlist(lapply(c("R", "J"), function(judge) {
  V <- copy(U); V[, Mb := get(setdiff(c("R", "J"), judge))]
  V[, lad := fcoalesce(Mb, A1, A2, Bl, SbH, SbL)]
  s <- with(V, fcoalesce(Mb * w[[1]], 0) + fcoalesce(A1 * w[[2]], 0) + fcoalesce(A2 * w[[3]], 0) +
              fcoalesce(Bl * w[[4]], 0) + fcoalesce(SbH * w[[5]], 0) + fcoalesce(SbL * w[[6]], 0))
  V[, vote := fifelse(s > 0, 1L, fifelse(s < 0, -1L, lad))]
  d <- V[!is.na(lad) & !is.na(get(judge))]
  data.table(judged_by = c(R = "ROSMAP", J = "Jaffe")[[judge]], pairs = nrow(d), ladder_correct = round(mean(d$lad == d[[judge]]), 4),
             weighted_vote_correct = round(mean(d$vote == d[[judge]]), 4), decisions_differ = sum(d$lad != d$vote)) }))
e <- U[!is.na(A2) & !is.na(Bl) & !is.na(SbH) & Bl == SbH & Bl != A2 & !is.na(R)]
f <- U[!is.na(Bl) & !is.na(SbH) & !is.na(R)]
edge <- sprintf(paste0("T2 opposed by both measured blood and blood SMR (replicated/high), ROSMAP-judged: %d pairs, ROSMAP sides with T2 on %d. ",
                       "The two blood sources err together on %d pairs against %.1f expected if independent, so their weights do not add."),
                nrow(e), sum(e$A2 == e$R), sum(f$Bl != f$R & f$SbH != f$R), nrow(f) * mean(f$Bl != f$R) * mean(f$SbH != f$R))
message(edge)

md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
dir.create("tools/reports", showWarnings = FALSE); dir.create("tools/out", showWarnings = FALSE)
fwrite(E, "tools/out/brain_overlap_estimates.csv")
writeLines(c("# Same-person overlap: brain SMR layer vs its cortex references", "",
  sprintf("Run %s. Jaffe signs at p < %s (all-ages age-adjusted model). Cluster bootstrap by CpG, %d resamples.", format(Sys.time(), "%Y-%m-%d %H:%M"), jp, nboot), "",
  "## Overlap-adjusted estimates (no comparison between donor-sharing datasets is used)", "", md(E), "",
  "## Dose-response by SMR strength", "", md(dr), "",
  "## Log-odds weights for the brain-target ladder", "", md(Wt), "", md(cmp), "", edge, ""),
  "tools/reports/validate_brain_overlap.md")
message("wrote tools/reports/validate_brain_overlap.md and tools/out/brain_overlap_estimates.csv")
