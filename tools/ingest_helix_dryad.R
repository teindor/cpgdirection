#!/usr/bin/env Rscript
# Lead 14 (3 Nov 2026): HELIX children's blood cis-eQTM catalogue, full and signed
# (Ruiz-Arenas et al. 2022, eLife; Dryad 10.5061/dryad.fxpnvx0t0, v2021-11-30).
# 823 children aged 6-11, 450K + Affymetrix HTA 2.0, hg19, cis = TSS +/- 500 kb; two models
# (cell-composition adjusted = the paper's main model; unadjusted); sigPair = the paper's
# significance call (CpG-level multiple testing, then nominal p per pair).
#
# What this script does
#   1. reads the significant cell-adjusted pairs (eQTM_autosome_adj.cells_SIG.txt.gz), keeps
#      pairs with a gene symbol within the package window (|TSS - CpG| <= 250 kb), collapses
#      transcript clusters to genes by majority sign (ties dropped), sign = sign(log2FC);
#   2. checks: adjusted vs unadjusted model (sign agreement on shared significant pairs);
#      HELIX vs the shipped blood pairs (all cohorts; the packaged_v1 HELIX rows specifically);
#      HELIX child blood vs the blood SMR layer by rung (a validation of the ladder on children);
#   3. with --full, reads the 13.6 M-row cell-adjusted file and keeps every TESTED pair that
#      the SMR layer covers (within the window), significant or not -- the tested-but-null set;
#   4. writes tools/out/measured_eqtms_add_helix.csv.gz (packaged schema), optional
#      tools/out/helix_tested_vs_smr.csv.gz, and tools/out/ingest_helix_dryad.md.
#   Nothing in the package is changed.
#
# Usage (package root):
#   Rscript tools/ingest_helix_dryad.R [--dir ~/helix_dryad] [--window 250000] [--full] [--data-dir <dir>]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
dir <- path.expand(get_arg("--dir", "~/helix_dryad"))
window <- as.numeric(get_arg("--window", "250000"))
do_full <- "--full" %in% args
out_dir <- get_arg("--out", "tools/out")
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else suppressPackageStartupMessages(library(cpgdirection))

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
big <- function(x) format(x, big.mark = ",")
rd <- function(f, ...) fread(cmd = sprintf("gzip -dc %s", shQuote(f)), showProgress = FALSE, ...)

# ---- 1. significant pairs, both models -------------------------------------------------------
prep <- function(f) {
  d <- rd(f)
  d <- d[sigPair == TRUE & !is.na(log2FC) & log2FC != 0 & !is.na(TC_gene) & nzchar(TC_gene)]
  d[, target_gene := toupper(trimws(TC_gene))]
  d[, tss_dist := as.numeric(TC_gene_TSS) - as.numeric(CpG_pos)]
  d <- d[abs(tss_dist) <= window]
  d[, sgn := sign(log2FC)]
  # transcript clusters -> gene: majority sign, ties dropped; keep the most significant TC's distance
  setorder(d, CpG, target_gene, p.value)
  d[, .(direction = { s <- sum(sgn); if (s == 0) NA_real_ else sign(s) }, n_tc = .N, tc_agreement = round(max(mean(sgn == 1), mean(sgn == -1)), 3),
        tss_dist = tss_dist[1], p_min = p.value[1], log2FC_top = log2FC[1]), by = .(cpg_id = CpG, target_gene)][!is.na(direction)]
}
A <- prep(file.path(dir, "eQTM_autosome_adj.cells_SIG.txt.gz"))
U <- prep(file.path(dir, "eQTM_autosome_unadj.cells_SIG.txt.gz"))
message(sprintf("cell-adjusted significant pairs within %s bp with a symbol: %s (%s CpGs, %s genes), negative %.3f; unadjusted: %s",
                big(window), big(nrow(A)), big(uniqueN(A$cpg_id)), big(uniqueN(A$target_gene)), mean(A$direction == -1), big(nrow(U))))
au <- merge(A[, .(cpg_id, target_gene, a = direction)], U[, .(cpg_id, target_gene, u = direction)], by = c("cpg_id", "target_gene"))
rows <- list(data.table(check = "adjusted vs unadjusted model (shared significant pairs)", agree = fmt(sum(au$a == au$u), nrow(au))))

# ---- 2. against the shipped blood layer and the SMR rungs -------------------------------------
M <- fread("inst/extdata/measured_eqtms.csv.gz"); M[, target_gene := toupper(trimws(target_gene))]
Mb <- M[tissue == "blood"]
x <- merge(A[, .(cpg_id, target_gene, a = direction)], Mb[, .(cpg_id, target_gene, b = direction, cohorts, source)], by = c("cpg_id", "target_gene"))
rows[[length(rows) + 1]] <- data.table(check = "HELIX adj vs shipped blood (all cohorts)", agree = fmt(sum(x$a == x$b), nrow(x)))
h1 <- x[grepl("HELIX", cohorts) | source == "packaged_v1"]
rows[[length(rows) + 1]] <- data.table(check = "HELIX adj vs shipped rows tagged HELIX / packaged_v1", agree = fmt(sum(h1$a == h1$b), nrow(h1)))
h2 <- x[!grepl("HELIX", cohorts) & source != "packaged_v1"]
rows[[length(rows) + 1]] <- data.table(check = "HELIX adj vs shipped blood from OTHER cohorts (independent)", agree = fmt(sum(h2$a == h2$b), nrow(h2)))
new <- A[!Mb, on = c("cpg_id", "target_gene")]
S <- as.data.table(cpgd_smr_directions()); S[, target_gene := toupper(trimws(target_gene))]
if (!"smr_rung" %in% names(S)) S[, smr_rung := fcase(smr_tier == "S1", "high", smr_tier == "S2", "moderate", default = "weak")]
if (!"outcome" %in% names(S)) S[, outcome := "eQTLGen"]
sj <- merge(A[, .(cpg_id, target_gene, a = direction)], S[, .(cpg_id, target_gene, s = direction, smr_rung, smr_tier, replication = if ("replication" %in% names(S)) replication else NA_character_)], by = c("cpg_id", "target_gene"))
by_rung <- sj[, .(n = .N, agree = fmt(sum(a == s), .N)), by = smr_rung][order(smr_rung)]
by_rr <- sj[, .(n = .N, agree = fmt(sum(a == s), .N)), by = .(smr_rung, replication, smr_tier)][order(smr_rung, replication, smr_tier)]
val <- rbindlist(rows); message("\n== sign agreement"); print(val)
message("\n== HELIX child blood (adj) vs blood SMR, by rung"); print(by_rung); print(by_rr)
# distance / sign profile of what is new
prof <- rbind(new[, .(set = "new", n = .N, neg = round(mean(direction == -1), 3)), by = .(bin = cut(abs(tss_dist), c(-Inf, 1500, 10000, 100000, Inf), labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb")))],
              Mb[, .(set = "shipped_blood", n = .N, neg = round(mean(direction == -1), 3)), by = .(bin = cut(abs(tss_dist), c(-Inf, 1500, 10000, 100000, Inf), labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb")))])[order(set, bin)]
print(prof)

# ---- 3. optional: tested-but-null set for the SMR-covered pairs --------------------------------
null_note <- "not built (run with --full)"
if (do_full) {
  message("\nreading the full cell-adjusted file (13.6 M rows) ...")
  Fz <- rd(file.path(dir, "eQTM_autosome_adj.cells.txt.gz"), select = c("CpG", "TC", "log2FC", "SE", "p.value", "sigPair", "CpG_pos", "TC_gene_TSS", "TC_gene"))
  Fz <- Fz[!is.na(TC_gene) & nzchar(TC_gene)]
  Fz[, target_gene := toupper(trimws(TC_gene))]
  Fz[, tss_dist := as.numeric(TC_gene_TSS) - as.numeric(CpG_pos)]
  Fz <- Fz[abs(tss_dist) <= window]
  setkey(Fz, CpG, target_gene)
  keep <- unique(S[, .(CpG = cpg_id, target_gene)])
  T <- Fz[keep, nomatch = 0L]
  setorder(T, CpG, target_gene, p.value)
  T <- T[, .(helix_log2FC = log2FC[1], helix_se = SE[1], helix_p = p.value[1], helix_sig = any(sigPair == TRUE), n_tc = .N, tss_dist = tss_dist[1]), by = .(cpg_id = CpG, target_gene)]
  T <- merge(T, S[, .(cpg_id, target_gene, smr_direction = direction, smr_rung, smr_tier)], by = c("cpg_id", "target_gene"))
  fwrite(T, file.path(out_dir, "helix_tested_vs_smr.csv.gz"))
  # how often does a significant SMR pair have ANY detectable HELIX association, by rung; and sign agreement at nominal p < 0.05
  nt <- T[, .(pairs = .N, helix_sig = round(mean(helix_sig), 3), nominal_p05 = round(mean(helix_p < 0.05), 3),
              agree_at_p05 = fmt(sum(sign(helix_log2FC)[helix_p < 0.05] == smr_direction[helix_p < 0.05]), sum(helix_p < 0.05))), by = smr_rung][order(smr_rung)]
  message("\n== SMR pairs tested in HELIX: detection and sign agreement by rung"); print(nt)
  null_note <- c(md(nt), "", sprintf("%s SMR pairs were tested in HELIX within the window; table in helix_tested_vs_smr.csv.gz.", big(nrow(T))))
  rm(Fz); invisible(gc())
}

# ---- 4. write ----------------------------------------------------------------------------------
add <- A[, .(cpg_id, target_gene, tissue = "blood", direction, tss_dist, n_cohorts = 1L, cohort_agreement = 1,
             cohorts = "HELIX_2022_Dryad_adj", source = "helix_dryad_2022")]
fwrite(add, file.path(out_dir, "measured_eqtms_add_helix.csv.gz"))
rep <- c("# HELIX (Dryad fxpnvx0t0) child-blood cis-eQTMs -- adapter report", "",
         sprintf("Run %s. Cell-adjusted significant pairs within %s bp with a gene symbol: %s (%s CpGs, %s genes; negative %.3f). New to the shipped blood layer: %s pairs (%s CpGs).",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), big(window), big(nrow(A)), big(uniqueN(A$cpg_id)), big(uniqueN(A$target_gene)), mean(A$direction == -1), big(nrow(new)), big(uniqueN(new$cpg_id))), "",
         "## Sign agreement", "", md(val), "",
         "## HELIX child blood vs the blood SMR layer, by rung", "", md(by_rung), "", "by rung x replication x tier:", "", md(by_rr), "",
         "## Distance and sign profile, new HELIX pairs vs shipped blood", "", md(prof), "",
         "## Tested-but-null set", "", null_note, "",
         "## Reading", "",
         "Admission to the measured blood layer follows the cohort rule used for every Atlas cohort: agreement with the OTHER blood cohorts on",
         "shared pairs in the 0.85+ range. The SMR-by-rung table is the ladder validated on 6-11-year-olds; if the rungs keep their order and",
         "their blood figures within a few points, the adult-blood validation transfers to children. The tested-but-null set is for",
         "calibration work (which SMR rungs have the highest share of pairs with no detectable association at all) and is not shipped.", "")
writeLines(rep, file.path(out_dir, "ingest_helix_dryad.md"))
message("\nwrote ", file.path(out_dir, "measured_eqtms_add_helix.csv.gz"), if (do_full) ", helix_tested_vs_smr.csv.gz" else "", " and ingest_helix_dryad.md")
