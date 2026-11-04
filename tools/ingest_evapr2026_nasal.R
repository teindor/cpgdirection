#!/usr/bin/env Rscript
# Lead 2b (3 Nov 2026): add the EVA-PR 2026 cis-eQTM pairs (Kim et al. 2026, Dryad
# 10.5061/dryad.rxwdbrvpq, FDR < 0.01, n = 455, 450K, hg19, TPM outcome, M-value
# predictor) to the measured nasal layer. Same cohort as the shipped EVA-PR nasal
# pairs (Atlas exports of Kim 2020/2023): this is COVERAGE from a fuller export of
# the same analysis, not replication, and the rows are tagged accordingly.
#
# Facts from tools/kim2026_nasal_dryad.R (3 Nov 2026): 23,953 cis pairs in the file;
# 10,501 of the 11,045 shipped nasal pairs are in it with sign agreement 1.000;
# 13,452 pairs are not in the shipped table; new-pair signs agree with measured
# blood at 0.926 and with measured solid tissue at 0.904.
#
# Usage (package root):
#   Rscript tools/ingest_evapr2026_nasal.R [--dryad ~/kim2026_dryad] [--window 250000]
#          [--existing inst/extdata/measured_eqtms.csv.gz] [--out tools/out]
# Writes <out>/measured_eqtms_v3.csv.gz (candidate replacement for inst/extdata) and
# <out>/ingest_evapr2026_nasal.md. Nothing in the package is changed by running it.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
dryad    <- path.expand(get_arg("--dryad", "~/kim2026_dryad"))
window   <- as.numeric(get_arg("--window", "250000"))
existing <- get_arg("--existing", "inst/extdata/measured_eqtms.csv.gz")
out_dir  <- get_arg("--out", "tools/out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(existing)) stop("packaged measured table not found at ", existing, " (run from the package root)")

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
big <- function(x) format(x, big.mark = ",")

# ---- 1. the Dryad cis-eQTM ---------------------------------------------------------
Q <- fread(file.path(dryad, "cis-eQTM.nocelltypes.fdr0.01.Original.TPM.csv"), showProgress = FALSE)
setnames(Q, c("chromosome cpg", "position cpg", "chromosome gene"), c("chr_cpg", "pos_cpg", "chr_gene"))
n_file <- nrow(Q)
Q <- Q[!is.na(beta) & beta != 0 & !is.na(fdr) & fdr < 0.01 & analysis == "cis"]
Q[, target_gene := toupper(trimws(gene))]
Q <- Q[nzchar(target_gene) & !grepl("^ENSG", target_gene)]
# package convention (inst/scripts/build_measured_eqtms.R): tss_dist = TSS - CpG position, same chromosome
Q[, tss_dist := as.numeric(TSS) - as.numeric(pos_cpg)]
Q[, same_chr := sub("^chr", "", as.character(chr_cpg)) == sub("^chr", "", as.character(chr_gene))]
n_sig <- nrow(Q)
Q <- Q[same_chr & abs(tss_dist) <= window]
n_win <- nrow(Q)
setorder(Q, cpg, target_gene, pvalue)
Q <- unique(Q, by = c("cpg", "target_gene"))
N <- Q[, .(cpg_id = cpg, target_gene, tissue = "nasal_epithelium", direction = sign(beta), tss_dist,
           pvalue, fdr, dist_file = dist)]
message(sprintf("Dryad cis-eQTM: %s rows; %s signed at FDR<0.01; %s within %s bp on the same chromosome; %s unique pairs",
                big(n_file), big(n_sig), big(n_win), big(window), big(nrow(N))))
# sanity on the distance convention: the file's own dist should equal +/- tss_dist
conv <- N[, .(same_sign = mean(sign(dist_file) == sign(tss_dist), na.rm = TRUE), same_abs = mean(abs(dist_file) == abs(tss_dist), na.rm = TRUE))]
message(sprintf("file 'dist' vs TSS - pos: |equal| %.3f, same sign %.3f (the package stores TSS - CpG)", conv$same_abs, conv$same_sign))

# ---- 2. the packaged measured table ------------------------------------------------
E <- fread(existing)
E[, target_gene := toupper(trimws(target_gene))]
En <- E[tissue == "nasal_epithelium"]
Eb <- E[tissue == "blood"]
Es <- E[tissue == "solid_tissue"]
j <- merge(N, En[, .(cpg_id, target_gene, ship_dir = direction, cohorts, n_cohorts, cohort_agreement)], by = c("cpg_id", "target_gene"))
new <- N[!En, on = c("cpg_id", "target_gene")]
jb <- merge(new, Eb[, .(cpg_id, target_gene, b_dir = direction)], by = c("cpg_id", "target_gene"))
js <- merge(new, Es[, .(cpg_id, target_gene, s_dir = direction)], by = c("cpg_id", "target_gene"))
message(sprintf("shipped nasal %s; shared %s, agreement %s; NEW pairs %s (%s CpGs); new vs blood %s; new vs solid %s",
                big(nrow(En)), big(nrow(j)), fmt(sum(j$direction == j$ship_dir), nrow(j)), big(nrow(new)), big(uniqueN(new$cpg_id)),
                fmt(sum(jb$direction == jb$b_dir), nrow(jb)), fmt(sum(js$direction == js$s_dir), nrow(js))))
dist_tab <- rbind(new[, .(set = "new", n = .N, neg_share = round(mean(direction == -1), 3)), by = .(bin = cut(abs(tss_dist), c(-Inf, 1500, 10000, 100000, Inf), labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb")))],
                  En[, .(set = "shipped", n = .N, neg_share = round(mean(direction == -1), 3)), by = .(bin = cut(abs(tss_dist), c(-Inf, 1500, 10000, 100000, Inf), labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb")))])[order(set, bin)]
print(dist_tab)

# ---- 3. merged table in the packaged schema -------------------------------------------
# shared pairs: keep the packaged row, append the Dryad cohort tag so provenance shows both exports
E2 <- copy(E)
key_shared <- paste(j$cpg_id, j$target_gene)
E2[tissue == "nasal_epithelium" & paste(cpg_id, target_gene) %in% key_shared & !grepl("EVA-PR_2026", cohorts),
   cohorts := paste(cohorts, "EVA-PR_2026_Dryad", sep = "|")]
add <- new[, .(cpg_id, target_gene, tissue, direction, tss_dist, n_cohorts = 1L, cohort_agreement = 1,
               cohorts = "EVA-PR_2026_Dryad", source = "evapr_2026_dryad")]
miss <- setdiff(names(E2), names(add)); for (m in miss) add[, (m) := NA]
add <- add[, names(E2), with = FALSE]
M <- rbind(E2, add)
setorder(M, tissue, cpg_id, target_gene)
stopifnot(!anyDuplicated(M[, .(cpg_id, target_gene, tissue)]))
fwrite(M, file.path(out_dir, "measured_eqtms_v3.csv.gz"))
fwrite(new, file.path(out_dir, "evapr2026_new_nasal_pairs.csv.gz"))

rep <- c("# EVA-PR 2026 (Dryad) cis-eQTM pairs added to the measured nasal layer", "",
         sprintf("Run %s. Source: Kim et al. 2026 Sci Adv, Dryad 10.5061/dryad.rxwdbrvpq, cis-eQTM FDR < 0.01 (n = 455 youth, 450K, hg19).", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
         sprintf("- file rows %s; signed at FDR < 0.01 %s; same chromosome within %s bp %s; unique pairs %s.", big(n_file), big(n_sig), big(window), big(n_win), big(nrow(N))),
         sprintf("- file `dist` vs TSS - CpG: |equal| %.3f, same sign %.3f.", conv$same_abs, conv$same_sign),
         sprintf("- shipped nasal pairs %s; shared %s with sign agreement %s (same cohort, same analysis family).", big(nrow(En)), big(nrow(j)), fmt(sum(j$direction == j$ship_dir), nrow(j))),
         sprintf("- NEW nasal pairs %s over %s CpGs; vs measured blood %s; vs measured solid %s.", big(nrow(new)), big(uniqueN(new$cpg_id)), fmt(sum(jb$direction == jb$b_dir), nrow(jb)), fmt(sum(js$direction == js$s_dir), nrow(js))), "",
         "## Distance and sign mix, new vs shipped nasal pairs", "", md(dist_tab), "",
         sprintf("## Merged table: %s rows (%s nasal, was %s); written to %s", big(nrow(M)), big(sum(M$tissue == "nasal_epithelium")), big(nrow(En)), file.path(out_dir, "measured_eqtms_v3.csv.gz")), "",
         "Provenance: shared pairs keep their Atlas cohort string with `|EVA-PR_2026_Dryad` appended; new pairs carry",
         "`cohorts = EVA-PR_2026_Dryad`, `source = evapr_2026_dryad`, `n_cohorts = 1`. Same cohort as the Atlas EVA-PR",
         "exports, so no pair gains replication from this; it is coverage. Ship by copying measured_eqtms_v3.csv.gz to",
         "inst/extdata/measured_eqtms.csv.gz (1.1.0), re-running the distance curves is NOT required (they are fitted,",
         "not recomputed at load), and updating inst/scripts/make-extdata.md.", "")
writeLines(rep, file.path(out_dir, "ingest_evapr2026_nasal.md"))
message("wrote ", file.path(out_dir, "measured_eqtms_v3.csv.gz"), ", evapr2026_new_nasal_pairs.csv.gz and ingest_evapr2026_nasal.md")
