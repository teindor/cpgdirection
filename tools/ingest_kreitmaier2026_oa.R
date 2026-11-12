#!/usr/bin/env Rscript
# Lead 8 (11 Nov 2026; PI ruling "admit as a measured solid-tissue source"): Kreitmaier ... Zeggini 2026,
# Ann Rheum Dis (10.1016/j.ard.2026.01.020), eQTM maps of osteoarthritis joint tissues, full summary statistics
# from the Musculoskeletal Knowledge Portal (msk.hugeamp.org/downloads.html, "eQTMs ALL" = eqtm_full.tar.gz, 3.16 GB).
#   eqtm_lg_full.txt.gz   low-grade cartilage  (n = 227)   eqtm_hg_full.txt.gz  high-grade cartilage (n = 127)
#   eqtm_syn_full.txt.gz  synovium             (n = 226)   eqtm_fp_full.txt.gz  fat pad (n = 45; no significant pair, not used)
# EPIC array, RNA-seq, MatrixEQTL, CpG within 1 Mb of the TSS; every tested pair is in the file. A pair is
# significant (eQTM_cpg = 1) when its gene passes 5% FDR (Storey-Tibshirani) and the pair passes the gene's
# Bonferroni threshold. All tissues come from the same patients: ONE donor group.
#
# What this does (nothing in the package is changed):
#   1. reads the significant pairs of the three tissues straight from the tar;
#   2. gene = the Ensembl ID of the file, named as the catalogue names it (cpgdirectionData 0.99.7 gene table),
#      else by its Ensembl GRCh37 symbol; TSS = the same table's strand-aware hg19 TSS, else Ensembl GRCh37;
#   3. keeps pairs on one chromosome within 250 kb; direction = sign(t); tss_dist = TSS - CpG;
#   4. two cohorts, Kreitmaier2026_cartilage (low- and high-grade, one vote, opposite signs dropped) and
#      Kreitmaier2026_synovium; one row per pair, majority sign, ties dropped;
#   5. checks the signs against the shipped measured table and previews the merge;
#   6. writes <out>/measured_eqtms_add_kreitmaier2026.csv.gz and <out>/ingest_kreitmaier2026_oa.md.
#
# Usage (package root):
#   Rscript tools/ingest_kreitmaier2026_oa.R --tar ~/Downloads/eqtm_full.tar.gz --gene-table <gene_table.csv.gz> \
#     --ens <ens_grch37.csv.gz> [--out tools/out] [--data-dir <dir>] [--window 250000]
# gene_table.csv.gz: pools/gene_table.csv.gz of catalogue_build_2.99.10.tar.gz (Zenodo record 24091940).
# ens_grch37.csv.gz: written by cpgdirectionData inst/scripts/make-lookups/export_inputs.R (annotables grch37).

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
tar_f <- path.expand(get_arg("--tar", "~/Downloads/eqtm_full.tar.gz"))
gt_f <- get_arg("--gene-table"); ens_f <- get_arg("--ens")
out_dir <- get_arg("--out", "tools/out"); data_dir <- get_arg("--data-dir")
window <- as.numeric(get_arg("--window", "250000"))
stopifnot(file.exists(tar_f), !is.null(gt_f), file.exists(gt_f), !is.null(ens_f), file.exists(ens_f))
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) suppressPackageStartupMessages(library(cpgdirection))
} else suppressPackageStartupMessages(library(cpgdirection))

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(trimws(x), collapse = " | "), " |"))) }
big <- function(x) format(x, big.mark = ",")
SRC <- "kreitmaier2026_msk"

# ---- 1. significant pairs from the tar ---------------------------------------------------
read_sig <- function(t) {
  cmd <- sprintf("tar -xzOf %s eqtm_%s_full.txt.gz | gunzip -c | awk -F'\\t' 'NR==1 || $8==1'", shQuote(tar_f), t)
  d <- fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE)
  stopifnot(all(c("msite", "gene", "statistic", "pvalue", "eQTM_cpg") %in% names(d)))
  d[, .(group = t, cpg_id = as.character(msite), gene_id = sub("\\..*$", "", as.character(gene)), t = as.numeric(statistic), p = as.numeric(pvalue),
        gene_name = if ("gene_name" %in% names(d)) as.character(gene_name) else NA_character_)]
}
message("reading the significant pairs from ", tar_f, " ...")
W <- rbindlist(lapply(c("lg", "hg", "syn"), read_sig))
W <- W[grepl("^cg\\d{8}$", cpg_id) & grepl("^ENSG\\d+$", gene_id) & is.finite(t) & t != 0]
raw <- W[, .(significant_pairs = .N, cpgs = uniqueN(cpg_id), genes = uniqueN(gene_id)), by = group]
print(raw)

# ---- 2. gene name and TSS ------------------------------------------------------------------
GT <- fread(gt_f); stopifnot(all(c("gene_id", "chr", "strand", "display", "tss") %in% names(GT)))
EN <- unique(fread(ens_f), by = "ensgene"); stopifnot(all(c("ensgene", "symbol", "chr", "start", "end", "strand") %in% names(EN)))
EN[, `:=`(chr = paste0("chr", sub("^chr", "", as.character(chr))), sym = toupper(trimws(as.character(symbol))))]
EN[, tss := fifelse(strand == 1, as.numeric(start), as.numeric(end))]
EN <- EN[nzchar(sym) & sym != "NA" & chr %in% paste0("chr", c(1:22, "X", "Y"))]
taken <- toupper(GT$display)
EN[, n_sym := .N, by = sym]
E2 <- EN[!ensgene %in% GT$gene_id & n_sym == 1L & !sym %in% taken, .(gene_id = ensgene, target_gene = sym, gchr = chr, tss, gene_source = "ensembl_grch37")]
GM <- rbind(GT[, .(gene_id, target_gene = toupper(display), gchr = chr, tss = as.numeric(tss), gene_source = "catalogue_0.99.7")], E2)
stopifnot(!anyDuplicated(GM$gene_id), !anyDuplicated(GM$target_gene))
W <- merge(W, GM, by = "gene_id", all.x = TRUE)
n_nogene <- W[is.na(target_gene), .N]
name_diff <- W[!is.na(target_gene) & !is.na(gene_name) & nzchar(gene_name), mean(toupper(gene_name) != target_gene)]
W <- W[!is.na(target_gene)]

# ---- 3. position, window, sign -----------------------------------------------------------
P <- as.data.table(cpgd_cpg_positions())
W <- merge(W, P[, .(cpg_id, chr, pos = as.numeric(pos))], by = "cpg_id", all.x = TRUE)
n_nopos <- W[is.na(pos), .N]; W <- W[!is.na(pos)]
n_chr <- W[chr != gchr, .N]; W <- W[chr == gchr]
W[, tss_dist := as.integer(round(tss - pos))]
n_far <- W[abs(tss_dist) > window, .N]; W <- W[abs(tss_dist) <= window]
W[, `:=`(direction = sign(t), cohort = fifelse(group == "syn", "Kreitmaier2026_synovium", "Kreitmaier2026_cartilage"))]
kept <- W[, .(pairs_kept = .N, neg_share = sprintf("%.3f", mean(direction == -1)), median_abs_dist = big(as.integer(median(abs(tss_dist))))), by = group]
steps <- data.table(step = c("no gene in the catalogue table or Ensembl GRCh37 (or an ambiguous symbol)", "CpG without an hg19 position in the package",
                             "CpG and gene on different chromosomes", sprintf("farther than %s bp from the TSS", big(window))),
                    rows_dropped = big(c(n_nogene, n_nopos, n_chr, n_far)))
x <- merge(W[group == "lg", .(cpg_id, target_gene, a = direction)], W[group == "hg", .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
lg_hg <- fmt(sum(x$a == x$b), nrow(x))
# one vote per cohort
C <- W[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, tss_dist = tss_dist[1], gene_id = gene_id[1], gene_source = gene_source[1],
           max_abs_t = max(abs(t)), min_p = min(p)), by = .(cohort, cpg_id, target_gene)]
n_tie_cart <- C[is.na(direction), .N]; C <- C[!is.na(direction)]
x <- merge(C[cohort == "Kreitmaier2026_cartilage", .(cpg_id, target_gene, a = direction)], C[cohort == "Kreitmaier2026_synovium", .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
cart_syn <- fmt(sum(x$a == x$b), nrow(x))
add <- C[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, tss_dist = tss_dist[1], n_cohorts = .N,
             cohort_agreement = round(max(mean(direction == 1), mean(direction == -1)), 3), cohorts = paste(sort(unique(cohort)), collapse = "|"),
             source = SRC, gene_id = gene_id[1], gene_source = gene_source[1], max_abs_t = round(max(max_abs_t), 3), min_p = signif(min(min_p), 3)),
         by = .(cpg_id, target_gene)]
n_tie <- add[is.na(direction), .N]; add <- add[!is.na(direction)]
add[, `:=`(tissue = "solid_tissue", direction = as.integer(direction))]
setcolorder(add, c("cpg_id", "target_gene", "tissue", "direction", "tss_dist", "n_cohorts", "cohort_agreement", "cohorts", "source"))
setorder(add, cpg_id, target_gene)
message(sprintf("pairs to add: %s (%s CpGs, %s genes)", big(nrow(add)), big(uniqueN(add$cpg_id)), big(uniqueN(add$target_gene))))

# ---- 4. validation against the shipped measured table ---------------------------------------
M <- as.data.table(cpgd_measured_eqtms()); M[, target_gene := toupper(trimws(target_gene))]
stopifnot(!SRC %in% M$source)
rows <- list()
chk <- function(lab, A, B) { x <- merge(A[, .(cpg_id, target_gene, d = direction)], B[, .(cpg_id, target_gene, s = direction)], by = c("cpg_id", "target_gene"))
  data.table(check = lab, agreement = fmt(sum(x$d == x$s), nrow(x))) }
for (ti in c("solid_tissue", "blood", "nasal_epithelium")) {
  rows[[length(rows) + 1L]] <- chk(paste("all added pairs vs shipped", ti), add, M[tissue == ti])
  for (co in c("Kreitmaier2026_cartilage", "Kreitmaier2026_synovium"))
    rows[[length(rows) + 1L]] <- chk(paste(sub("Kreitmaier2026_", "", co), "vs shipped", ti), C[cohort == co], M[tissue == ti])
}
Ms <- M[tissue == "solid_tissue"]
one <- Ms[n_cohorts == 1L]
for (co in sort(unique(one$cohorts))) rows[[length(rows) + 1L]] <- chk(paste("all added pairs vs solid pairs measured only in", co), add, one[cohorts == co])
rows[[length(rows) + 1L]] <- chk("all added pairs vs solid pairs measured in 2+ cohorts", add, Ms[n_cohorts >= 2L])
rows[[length(rows) + 1L]] <- data.table(check = "low-grade vs high-grade cartilage (both significant)", agreement = lg_hg)
rows[[length(rows) + 1L]] <- data.table(check = "cartilage vs synovium (both significant)", agreement = cart_syn)
val <- rbindlist(rows); print(val)

# ---- 5. merge preview: one vote per cohort, majority sign, ties dropped ------------------------
dg <- function(co) fcase(startsWith(co, "GTEx"), "GTEx", startsWith(co, "Kreitmaier2026"), "Kreitmaier2026", default = co)
votes <- function(D) { n <- as.integer(D$n_cohorts); k <- as.integer(round(n * D$cohort_agreement))
  data.table(cpg_id = D$cpg_id, target_gene = D$target_gene, plus = fifelse(D$direction > 0, k, n - k), minus = fifelse(D$direction > 0, n - k, k), cohorts = D$cohorts) }
V <- rbind(votes(Ms), votes(add))
R <- V[, .(plus = sum(plus), minus = sum(minus), cohorts = paste(unique(unlist(strsplit(cohorts, "|", fixed = TRUE))), collapse = "|")), by = .(cpg_id, target_gene)]
R[, `:=`(n = plus + minus, direction = sign(plus - minus))]
R[, agr := pmax(plus, minus) / n]
ngrp <- function(co) vapply(strsplit(co, "|", fixed = TRUE), function(z) uniqueN(dg(z)), integer(1))
R[, groups := ngrp(cohorts)]; Ms[, groups := ngrp(cohorts)]
old <- Ms[, .(cpg_id, target_gene, old = direction)]
R <- merge(R, old, by = c("cpg_id", "target_gene"), all.x = TRUE)
inboth <- merge(add[, .(cpg_id, target_gene, d = direction)], Ms[, .(cpg_id, target_gene, s = direction)], by = c("cpg_id", "target_gene"))
prev <- data.table(item = c("solid pairs in the shipped measured table", "pairs added by this source", "of them new to the solid table", "already there, same sign", "already there, opposite sign",
                            "pairs dropped as ties after the vote", "shipped pairs whose sign changes", "solid pairs after the merge",
                            "solid pairs with 2+ independent donor groups agreeing (>= 0.99), before", "the same, after"),
                   n = big(c(nrow(Ms), nrow(add), nrow(add) - nrow(inboth), sum(inboth$d == inboth$s), sum(inboth$d != inboth$s),
                             R[direction == 0, .N], R[!is.na(old) & direction != 0 & direction != old, .N], R[direction != 0, .N],
                             Ms[groups >= 2L & cohort_agreement >= 0.99, .N], R[direction != 0 & groups >= 2L & agr >= 0.99, .N])))
print(prev)

# ---- 6. write -----------------------------------------------------------------------------------
fwrite(add, file.path(out_dir, "measured_eqtms_add_kreitmaier2026.csv.gz"))
rep <- c("# Kreitmaier et al. 2026 osteoarthritis eQTMs as measured solid-tissue pairs (lead 8)", "",
  sprintf("Run %s. Source: %s (MSK portal, full summary statistics). Nothing in the package is changed.", format(Sys.time(), "%Y-%m-%d %H:%M"), basename(tar_f)), "",
  "Low- and high-grade cartilage and synovium come from the same patients: one donor group (`Kreitmaier2026`).", "",
  "## Significant pairs in the files (gene at 5% FDR, pair under the gene's Bonferroni threshold)", "", md(raw), "",
  "## Rows dropped while placing the pairs", "", md(steps), "",
  sprintf("Kept, by tissue (the file's gene name differs from the name used here for %.1f%% of rows; the Ensembl ID decides):", 100 * name_diff), "", md(kept), "",
  sprintf("Cartilage pairs dropped because low- and high-grade disagree: %s. Pairs dropped because cartilage and synovium disagree: %s.", big(n_tie_cart), big(n_tie)), "",
  "## Sign agreement", "", md(val), "",
  "## Merge preview (one vote per cohort, majority sign, ties dropped; GTEx tissues are one donor group)", "", md(prev), "",
  "## Output", "",
  sprintf("- `measured_eqtms_add_kreitmaier2026.csv.gz`: %s pairs over %s CpGs and %s genes; columns of the measured table plus `gene_id`, `gene_source`, `max_abs_t`, `min_p`.",
          big(nrow(add)), big(uniqueN(add$cpg_id)), big(uniqueN(add$target_gene))),
  sprintf("- cohorts: `Kreitmaier2026_cartilage` %s pairs, `Kreitmaier2026_synovium` %s pairs, both %s; source `%s`.",
          big(add[cohorts == "Kreitmaier2026_cartilage", .N]), big(add[cohorts == "Kreitmaier2026_synovium", .N]), big(add[n_cohorts == 2L, .N]), SRC), "",
  "## For the integration pass", "",
  "- `tools/build_*_extdata.R`: fold the file into `inst/extdata/measured_eqtms.csv.gz` by the cohort rule.",
  "- `.cpgd_donor_group()` in `R/rung_accuracy.R`: `startsWith(cohort, \"Kreitmaier2026\") -> \"Kreitmaier2026\"`.",
  sprintf("- `tests/testthat/test-measured-cohorts.R`: add `\"%s\"` to the allowed sources.", SRC),
  "- Roxygen of `cpgd_measured_eqtms()` and `inst/scripts/make-extdata.md`: the new cohorts.",
  "- `tools/validate_rung_accuracy.R`: FUSION and GTEx solid gain an independent solid comparator.", "")
writeLines(rep, file.path(out_dir, "ingest_kreitmaier2026_oa.md"))
message("wrote ", file.path(out_dir, "measured_eqtms_add_kreitmaier2026.csv.gz"), " and ingest_kreitmaier2026_oa.md")
