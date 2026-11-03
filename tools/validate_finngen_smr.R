#!/usr/bin/env Rscript
# tools/validate_finngen_smr.R -- validate the cell-type SMR contexts.
# (1) concordance with the measured blood eQTMs shipped in the package, by
#     context and tier, with Wilson 95% CIs and the majority-class baseline;
# (2) agreement with the bulk GoDMC x eQTLGen layer (cpgd_smr_directions())
#     on shared CpG-gene pairs, by tier;
# (3) cross-context consistency: per pair, number of contexts calling and the
#     fraction agreeing; PBMC pseudobulk vs each cell type.
# Usage: Rscript tools/validate_finngen_smr.R --pattern 'smr/smr_blood_*.directions.csv.gz' \
#          --out tools/reports/validate_finngen_smr.md
suppressPackageStartupMessages({ library(data.table); library(cpgdirection) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
pat <- get_arg("--pattern", "smr/smr_blood_*.directions.csv.gz")
out <- get_arg("--out", "validate_finngen_smr.md")
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n; c <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  sprintf("%.1f%% (%.1f-%.1f)", 100 * p, 100 * (c - h) / d, 100 * (c + h) / d)
}
files <- Sys.glob(pat); stopifnot(length(files) > 0)
X <- rbindlist(lapply(files, fread), fill = TRUE)
X[, direction := as.integer(direction)]
L <- c(sprintf("# FinnGen cell-type SMR validation (%s)", format(Sys.Date())), "",
       sprintf("Files: %d; context-pairs: %d; contexts: %s", length(files), nrow(X), paste(sort(unique(X$context)), collapse = ", ")), "")

# (1) measured blood eQTMs
M <- cpgd_measured_eqtms()[tissue == "blood", .(cpg_id, target_gene, m_dir = as.integer(direction))]
J <- merge(X, M, by = c("cpg_id", "target_gene"))
base <- max(mean(J$m_dir == 1), mean(J$m_dir == -1))
L <- c(L, "## 1. Concordance with measured blood eQTMs", "",
       sprintf("Overlap %d context-pairs (%d unique pairs); majority-class baseline %.1f%%.", nrow(J), uniqueN(J[, .(cpg_id, target_gene)]), 100 * base), "",
       "| context | tier | n | concordance (95% CI) |", "|---|---|---:|---|")
T1 <- J[, .(n = .N, k = sum(direction == m_dir)), by = .(context, smr_tier)][order(context, smr_tier)]
L <- c(L, T1[, sprintf("| %s | %s | %d | %s |", context, smr_tier, n, wilson(k, n))])
T1a <- J[, .(n = .N, k = sum(direction == m_dir)), by = smr_tier][order(smr_tier)]
L <- c(L, "", "All contexts pooled:", "", "| tier | n | concordance (95% CI) |", "|---|---:|---|",
       T1a[, sprintf("| %s | %d | %s |", smr_tier, n, wilson(k, n))], "")

# (2) bulk SMR layer
B <- cpgd_smr_directions()[, .(cpg_id, target_gene, b_dir = as.integer(direction), b_tier = smr_tier)]
K <- merge(X, B, by = c("cpg_id", "target_gene"))
L <- c(L, "## 2. Agreement with the bulk GoDMC x eQTLGen layer", "",
       sprintf("Shared pairs: %d of %d context-pairs (%.1f%%).", nrow(K), nrow(X), 100 * nrow(K) / nrow(X)), "",
       "| context | n shared | sign agreement (95% CI) | S1 both | S1 agreement |", "|---|---:|---|---:|---|")
T2 <- K[, .(n = .N, k = sum(direction == b_dir), n1 = sum(smr_tier == "S1" & b_tier == "S1"),
            k1 = sum(smr_tier == "S1" & b_tier == "S1" & direction == b_dir)), by = context][order(context)]
L <- c(L, T2[, sprintf("| %s | %d | %s | %d | %s |", context, n, wilson(k, n), n1, ifelse(n1 > 0, wilson(k1, n1), "-"))], "")

# (3) cross-context consistency
C <- X[context != "blood_PBMC", .(n_ctx = .N, agree = max(mean(direction == 1), mean(direction == -1))), by = .(cpg_id, target_gene)]
L <- c(L, "## 3. Cross-context consistency (cell types, PBMC excluded)", "",
       sprintf("Unique pairs: %d; called in >= 2 cell types: %d; of those, unanimous: %.1f%%, any disagreement: %.1f%%.",
               nrow(C), sum(C$n_ctx >= 2), 100 * mean(C[n_ctx >= 2, agree == 1]), 100 * mean(C[n_ctx >= 2, agree < 1])), "",
       "| n cell types calling | pairs | unanimous |", "|---:|---:|---|")
T3 <- C[, .(pairs = .N, k = sum(agree == 1)), by = n_ctx][order(n_ctx)]
L <- c(L, T3[, sprintf("| %d | %d | %s |", n_ctx, pairs, wilson(k, pairs))], "")
P <- X[context == "blood_PBMC", .(cpg_id, target_gene, p_dir = direction)]
Q <- merge(X[context != "blood_PBMC"], P, by = c("cpg_id", "target_gene"))
T4 <- Q[, .(n = .N, k = sum(direction == p_dir)), by = context][order(context)]
L <- c(L, "PBMC pseudobulk vs each cell type on shared pairs:", "", "| cell type | n shared | agreement (95% CI) |", "|---|---:|---|",
       T4[, sprintf("| %s | %d | %s |", context, n, wilson(k, n))], "")
if ("heidi_status" %in% names(X)) {
  T5 <- X[, .(pass = sum(heidi_status == "pass"), fail = sum(heidi_status == "fail"), nt = sum(heidi_status == "not_tested")), by = context][order(context)]
  L <- c(L, "## 4. HEIDI", "", "| context | pass | fail | not tested |", "|---|---:|---:|---:|",
         T5[, sprintf("| %s | %d | %d | %d |", context, pass, fail, nt)], "")
}
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(L, out); cat(L, sep = "\n")
