#!/usr/bin/env Rscript
# November 2026 data build (3 Nov 2026): fold the validated new sources into the bundled layers NOW.
#
#   inst/extdata/measured_eqtms.csv.gz   + HELIX child blood (Dryad fxpnvx0t0; signed, cell-adjusted; 0.988 vs other cohorts)
#                                         + MESA sorted monocytes / T cells (Wang 2021 S2; 0.93 vs shipped blood)
#                                         - the 10 duplicated blood rows (DDTL, HLA-DQA2 double entries)
#                                         packaged_v1 blood rows that HELIX now covers are REPLACED by the signed HELIX row
#                                         liver (0.841, below the 0.85 bar) only with --include-liver
#   inst/extdata/epicv2_probe_qc.csv.gz  + T2T-CHM13 / pangenome uniqueness flags (Dong 2025), informational only
#   inst/extdata/brain_measured_eqtms.csv.gz (NEW) ROSMAP aged DLPFC measured pairs (Wang 2021 S2; covariates incl. age
#                                         regressed; T1 0.995, T2 0.939 against them). Jaffe groups excluded (age not regressed).
#
# Cohort rule for merging into the measured table (as in inst/scripts/build_measured_eqtms.R): one vote per cohort
# per pair, sign by majority, ties dropped; n_cohorts / cohort_agreement / cohorts recomputed.
#
# Usage (package root):
#   Rscript tools/build_nov2026_extdata.R [--uniq ~/probe_uniqueness/epicv2_probe_uniqueness.csv.gz] [--include-liver]
# Inputs from earlier runs: tools/out/measured_eqtms_add_helix.csv.gz, tools/out/measured_eqtms_add_wang2021.csv.gz,
# tools/out/brain_measured_eqtms_wang2021.csv.gz. Writes the three extdata files and tools/reports/build_nov2026_extdata.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
uniq_f <- path.expand(get_arg("--uniq", "~/probe_uniqueness/epicv2_probe_uniqueness.csv.gz"))
incl_liver <- "--include-liver" %in% args
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)))
need <- c("inst/extdata/measured_eqtms.csv.gz", "inst/extdata/epicv2_probe_qc.csv.gz", "tools/out/measured_eqtms_add_helix.csv.gz",
          "tools/out/measured_eqtms_add_wang2021.csv.gz", "tools/out/brain_measured_eqtms_wang2021.csv.gz", uniq_f)
miss <- need[!file.exists(need)]; if (length(miss)) stop("missing inputs: ", paste(miss, collapse = ", "))
big <- function(x) format(x, big.mark = ",")
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
rep <- c("# November 2026 data build", "", sprintf("Run %s.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

# ---- 1. measured table ----------------------------------------------------------------------
E0 <- fread("inst/extdata/measured_eqtms.csv.gz")
cols <- names(E0)
E <- copy(E0); E[, target_gene := toupper(trimws(target_gene))]
n_dup <- sum(duplicated(E[, .(cpg_id, target_gene, tissue)]))
E <- unique(E, by = c("cpg_id", "target_gene", "tissue"))
H <- fread("tools/out/measured_eqtms_add_helix.csv.gz"); H[, target_gene := toupper(trimws(target_gene))]
W <- fread("tools/out/measured_eqtms_add_wang2021.csv.gz"); W[, target_gene := toupper(trimws(target_gene))]
if (!incl_liver) W <- W[tissue == "blood"]
W <- W[!is.na(tss_dist) & abs(tss_dist) <= 250000]
H <- H[!is.na(tss_dist) & abs(tss_dist) <= 250000]
hkey <- paste(H$cpg_id, H$target_gene)
repl <- E$source == "packaged_v1" & E$tissue == "blood" & paste(E$cpg_id, E$target_gene) %in% hkey
n_repl <- sum(repl); E <- E[!repl]
votes <- function(D, existing) {
  D <- copy(D)
  D[, n := as.integer(n_cohorts)]
  D[, k := as.integer(round(n * cohort_agreement))]
  D[, `:=`(plus = fifelse(direction > 0, k, n - k), minus = fifelse(direction > 0, n - k, k))]
  D[, .(cpg_id, target_gene, tissue, plus, minus, tss_dist, cohorts, source, existing = existing)]
}
V <- rbindlist(list(votes(E, TRUE), votes(H, FALSE), votes(W, FALSE)), use.names = TRUE)
setorder(V, cpg_id, target_gene, tissue, -existing)
old_dir <- E[, .(cpg_id, target_gene, tissue, old = direction)]
M <- V[, .(plus = sum(plus), minus = sum(minus), tss_dist = tss_dist[!is.na(tss_dist)][1],
           cohorts = paste(unique(unlist(strsplit(cohorts, "|", fixed = TRUE))), collapse = "|"), source = source[1]),
       by = .(cpg_id, target_gene, tissue)]
M[, n_cohorts := plus + minus]
n_tie <- sum(M$plus == M$minus)
M <- M[plus != minus]
M[, direction := as.integer(sign(plus - minus))]
M[, cohort_agreement := round(pmax(plus, minus) / n_cohorts, 3)]
chk <- merge(M[, .(cpg_id, target_gene, tissue, direction)], old_dir, by = c("cpg_id", "target_gene", "tissue"))
n_flip <- sum(chk$direction != chk$old)
M <- M[, intersect(cols, names(M)), with = FALSE]
setcolorder(M, intersect(cols, names(M)))
setorder(M, tissue, cpg_id, target_gene)
stopifnot(!anyDuplicated(M[, .(cpg_id, target_gene, tissue)]))
fwrite(M, "inst/extdata/measured_eqtms.csv.gz", compress = "gzip")
tab <- merge(E0[, .(before = .N), by = tissue], M[, .(after = .N, multi_cohort = sum(n_cohorts > 1)), by = tissue], by = "tissue", all = TRUE)
message("measured table:"); print(tab)
rep <- c(rep, "## measured_eqtms", "", md(tab), "",
         sprintf("- duplicated rows removed: %s; packaged_v1 blood rows replaced by signed HELIX rows: %s", big(n_dup), big(n_repl)),
         sprintf("- pairs dropped as cohort ties: %s; existing pairs whose sign changed: %s", big(n_tie), big(n_flip)),
         sprintf("- liver included: %s", incl_liver), "")

# ---- 2. probe QC -----------------------------------------------------------------------------
Q <- fread("inst/extdata/epicv2_probe_qc.csv.gz")
newc <- c("t2t_ambiguous", "pan_t2t_ambiguous", "pangenome_xpop_ambiguous", "pan_t2t_eur_ambiguous", "t2t_snp_maf05", "chm13_chr", "chm13_pos")
Q[, intersect(newc, names(Q)) := NULL]
U <- fread(uniq_f)[, c("cpg_id", newc), with = FALSE]
Q <- merge(Q, U, by = "cpg_id", all.x = TRUE, sort = FALSE)
if ("qc_sources" %in% names(Q)) Q[!is.na(t2t_ambiguous) & !grepl("Dong2025", qc_sources), qc_sources := paste(qc_sources, "Dong2025_T2T_pangenome", sep = ";")]
setorder(Q, cpg_id)
fwrite(Q, "inst/extdata/epicv2_probe_qc.csv.gz", compress = "gzip")
qs <- Q[, .(cpgs = .N, t2t_ambiguous = sum(t2t_ambiguous %in% TRUE), pangenome_xpop_ambiguous = sum(pangenome_xpop_ambiguous %in% TRUE),
            of_which_not_qc_exclude = sum(pangenome_xpop_ambiguous %in% TRUE & !(qc_exclude %in% TRUE)))]
message("probe QC:"); print(qs)
rep <- c(rep, "## epicv2_probe_qc (flags informational; qc_exclude unchanged)", "", md(qs), "")

# ---- 3. brain measured -----------------------------------------------------------------------
B <- fread("tools/out/brain_measured_eqtms_wang2021.csv.gz")[group == "brain_aged"]
B <- B[, .(cpg_id, target_gene = toupper(trimws(target_gene)), direction = as.integer(direction), corr = round(corr, 4),
           FDR = signif(FDR, 3), tss_dist, cohort = "ROSMAP_DLPFC_aged", source = "wang2021_S2")]
B <- unique(B, by = c("cpg_id", "target_gene"))
setorder(B, cpg_id, target_gene)
fwrite(B, "inst/extdata/brain_measured_eqtms.csv.gz", compress = "gzip")
message(sprintf("brain measured: %s pairs over %s CpGs", big(nrow(B)), big(uniqueN(B$cpg_id))))
rep <- c(rep, "## brain_measured_eqtms (new)", "",
         sprintf("%s ROSMAP aged-DLPFC pairs over %s CpGs and %s genes; negative %.3f.", big(nrow(B)), big(uniqueN(B$cpg_id)), big(uniqueN(B$target_gene)), mean(B$direction == -1)), "")

dir.create("tools/reports", showWarnings = FALSE)
writeLines(rep, "tools/reports/build_nov2026_extdata.md")
message("wrote inst/extdata/{measured_eqtms,epicv2_probe_qc,brain_measured_eqtms}.csv.gz and tools/reports/build_nov2026_extdata.md")
