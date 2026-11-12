#!/usr/bin/env Rscript
# Four measured eQTM sources admitted on 12 Nov 2026 (PI rulings of that day; evidence bar: sign agreement
# >= 0.85 with what the package already measures). Each is one cohort and one donor group.
#
#   kidney    Liu et al. 2022, Nat Genet (10.1038/s41588-022-01097-w), Supplementary Table 12
#             (41588_2022_1097_MOESM4_ESM.xlsx, sheet "Table S12"); kidney tubule, n = 414, EPIC, RNA-seq,
#             CpG within 1 Mb of the TSS, global FDR < 0.05.
#             THE PUBLISHED SHEET IS MISALIGNED: its columns "CpG" and "Gene" (Ensembl ID) form one block and
#             the columns Chr, Gene Symbol, CpG Position, TSS Position, Distance, Beta, P value form another;
#             each block is internally consistent, but the two are sorted differently, so a row's CpG ID does
#             not belong to its statistics. The CpG of a statistics row is recovered from Chr + CpG Position
#             (the hg19 position of the package); the Ensembl gene is the one paired with that CpG in the
#             other block whose TSS lies nearest the row's TSS. The recovery is checked against Table S11 of
#             the same workbook, whose CpG IDs match their positions.
#   placenta  Delahaye et al. 2018, PLoS Genet (10.1371/journal.pgen.1007785), S20 Table
#             (pgen.1007785.s032.xlsx); n = 74, 450K, RNA-seq, 100 kb window, permutation-selected pairs.
#   blood     Kennedy et al. 2018, BMC Genomics (10.1186/s12864-018-4842-3), Additional file 1
#             (12864_2018_4842_MOESM1_ESM.txt); Grady Trauma Project whole blood, all African American,
#             n = 333, 450K, HT-12. The file lists pairs at p < 1e-5; same-chromosome pairs within the window
#             are kept (the paper's own threshold, p <= 1e-11, is recorded in `strict`).
#             Additional file 2 of that paper is MESA monocytes, already in the package: not read.
#   adipose   Voisin et al. 2016, PLoS One (10.1371/journal.pone.0157776), S20 Table
#             (pone.0157776.s020.xlsx); adipose tissue, n = 118, 450K, FDR < 5%.
#
# What this does (nothing in the package is changed):
#   1. reads the published tables; 2. names each gene as the catalogue names it (cpgdirectionData gene table),
#   else by its unique Ensembl GRCh37 symbol; 3. places each CpG at the package's hg19 position, keeps pairs on
#   one chromosome within the window, direction = sign of the reported effect, tss_dist = TSS - CpG;
#      A CpG the package holds no position for (the position table follows EPIC v2, which dropped about
#      91,000 probes of the 450K array) is placed by the Illumina 450K manifest, build 37 (--hm450; the
#      manifest's MAPINFO equals the package position on 99.9% of the probes both hold). The shipped
#      measured table already carries such CpGs; `pos_source` records which table placed each pair;
#   4. one vote per cohort per pair (opposite signs inside a cohort dropped); 5. checks signs against the
#   shipped measured table; 6. writes <out>/measured_eqtms_add_nov2026_sources.csv.gz (one row per cohort and
#   pair) and <out>/ingest_nov2026_measured_sources.md.
#
# Usage (package root):
#   Rscript tools/ingest_nov2026_measured_sources.R --kidney <MOESM4.xlsx> --placenta <s032.xlsx> \
#     --gtp <MOESM1_ESM.txt> --adipose <s020.xlsx> --gene-table <gene_table.csv.gz> --ens <ens_grch37.csv.gz> \
#     [--hm450 <GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz>] [--out tools/out] [--data-dir <dir>] [--window 250000]

suppressPackageStartupMessages({ library(data.table); library(readxl) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
f_kid <- get_arg("--kidney"); f_pla <- get_arg("--placenta"); f_gtp <- get_arg("--gtp"); f_adi <- get_arg("--adipose")
gt_f <- get_arg("--gene-table"); ens_f <- get_arg("--ens"); hm_f <- get_arg("--hm450")
out_dir <- get_arg("--out", "tools/out"); data_dir <- get_arg("--data-dir")
window <- as.numeric(get_arg("--window", "250000"))
for (f in list(f_kid, f_pla, f_gtp, f_adi, gt_f, ens_f)) stopifnot(!is.null(f), file.exists(path.expand(f)))
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
xl <- function(f, sheet, skip) suppressMessages(as.data.table(read_excel(path.expand(f), sheet = sheet, skip = skip, guess_max = 100000, .name_repair = "minimal")))

# ---- genes and positions --------------------------------------------------------------------
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
P <- as.data.table(cpgd_cpg_positions())[, .(cpg_id, chr, pos = as.numeric(pos), pos_source = "package")]
hm_note <- "No 450K manifest given (--hm450): CpGs the package holds no position for are dropped."
if (!is.null(hm_f)) {
  stopifnot(file.exists(path.expand(hm_f)))
  message("reading the 450K manifest ...")
  H <- suppressWarnings(fread(path.expand(hm_f), skip = "IlmnID,", fill = TRUE, select = c("Name", "Genome_Build", "CHR", "MAPINFO"), colClasses = "character"))
  H <- H[grepl("^c[gh]\\d", Name) & Genome_Build == "37" & CHR %in% c(1:22, "X", "Y") & grepl("^\\d+$", MAPINFO),
         .(cpg_id = Name, chr = paste0("chr", CHR), pos = as.numeric(MAPINFO), pos_source = "hm450_manifest")]
  H <- unique(H, by = "cpg_id")
  cal <- merge(H, P, by = "cpg_id")
  hm_same <- cal[, mean(chr.x == chr.y & pos.x == pos.y)]
  stopifnot(nrow(H) > 400000, nrow(cal) > 300000, hm_same > 0.99)
  hm_note <- sprintf("450K manifest (build 37): %s probes; %s also in the package table, identical position in %.2f%%; %s probes the package does not hold are placed by the manifest.",
                     big(nrow(H)), big(nrow(cal)), 100 * hm_same, big(sum(!H$cpg_id %in% P$cpg_id)))
  message(hm_note)
  P <- rbind(P, H[!cpg_id %in% P$cpg_id])
}
Pk <- unique(P[grepl("^cg", cpg_id)], by = c("chr", "pos"))
by_symbol <- function(D) {                       # D: cohort, cpg_id, sym, stat, p (a symbol field may hold several names)
  D <- D[!is.na(sym) & nzchar(sym)]
  D <- D[, .(sym = toupper(trimws(unlist(strsplit(sym, "\\s*///\\s*|\\s*;\\s*")))), cpg_id, stat, p, strict), by = .(cohort, row = seq_len(nrow(D)))]
  merge(D[, row := NULL], GM, by.x = "sym", by.y = "target_gene")[, target_gene := sym][, sym := NULL]
}

# ---- kidney (Liu 2022, Table S12; recovered) --------------------------------------------------
message("reading the kidney workbook (large) ...")
k12 <- xl(f_kid, "Table S12", 1)
stopifnot(all(c("Chr", "CpG", "Gene", "Gene Symbol", "CpG Position", "TSS Position", "Beta", "P value") %in% names(k12)))
A <- unique(k12[, .(cpg_id = as.character(CpG), ens = sub("\\..*$", "", as.character(Gene)))])
B <- k12[, .(chr = as.character(Chr), pos = as.numeric(`CpG Position`), sym_file = toupper(trimws(as.character(`Gene Symbol`))),
             tss_file = as.numeric(`TSS Position`), stat = as.numeric(Beta), p = as.numeric(`P value`))]
listed <- merge(k12[, .(cpg_id = as.character(CpG), chr_row = as.character(Chr))], P[, .(cpg_id, chr)], by = "cpg_id")
kid_listed_chr <- listed[, mean(chr == chr_row)]
B <- merge(B, Pk[, .(chr, pos, cpg_id)], by = c("chr", "pos"), all.x = TRUE)
kid_nopos <- B[is.na(cpg_id), .N]; B <- B[!is.na(cpg_id)]
# gene: the Ensembl ID paired with this CpG in the other block whose GRCh37 TSS is nearest the row's TSS
EA <- merge(A, EN[, .(ens = ensgene, ens_tss = tss, ens_chr = chr)], by = "ens")
B[, row := .I]
cand <- merge(B[, .(row, cpg_id, chr, tss_file)], EA, by = "cpg_id", allow.cartesian = TRUE)[chr == ens_chr]
cand[, gap := abs(ens_tss - tss_file)]
setorder(cand, row, gap); cand <- cand[, .SD[1L], by = row][gap <= 5000]
B <- merge(B, cand[, .(row, gene_id = ens)], by = "row", all.x = TRUE)
kid_by_block <- B[!is.na(gene_id), .N]
B1 <- merge(B[!is.na(gene_id)], GM, by = "gene_id")
B2 <- merge(B[is.na(gene_id)][, gene_id := NULL], GM, by.x = "sym_file", by.y = "target_gene")[, target_gene := sym_file]
KID <- rbind(B1[, .(cohort = "Liu2022_kidney", cpg_id, gene_id, target_gene, gchr, tss, gene_source, stat, p, strict = TRUE)],
             B2[, .(cohort = "Liu2022_kidney", cpg_id, gene_id, target_gene, gchr, tss, gene_source, stat, p, strict = TRUE)])
kid_nogene <- nrow(B) - nrow(KID)
# check on Table S11 (top target gene per mQTL CpG, with the eQTM beta)
k11 <- xl(f_kid, "Table S11", 1)
s11 <- k11[, .(cpg_id = as.character(`meQTL CpG`), sym_file = toupper(trimws(as.character(.SD[[13]]))), beta11 = suppressWarnings(as.numeric(`eQTM beta`)),
               fdr = suppressWarnings(as.numeric(`eQTM global FDR`)))][!is.na(fdr) & fdr < 0.05]
chk11 <- merge(s11, B[, .(cpg_id, sym_file, stat)], by = c("cpg_id", "sym_file"))
kid_chk <- sprintf("%s of %s Table S11 pairs at global FDR < 0.05 found; identical beta in %s; same sign in %s",
                   big(nrow(chk11)), big(nrow(s11)), big(sum(abs(chk11$beta11 - chk11$stat) < 5e-5)), big(sum(sign(chk11$beta11) == sign(chk11$stat))))
message("kidney: ", kid_chk)
stopifnot(nrow(chk11) > 10000, mean(sign(chk11$beta11) == sign(chk11$stat)) > 0.999)

# ---- placenta, blood (GTP), adipose -------------------------------------------------------------
pl <- xl(f_pla, 1, 1); stopifnot(all(c("Gene", "CpG", "P.value", "linear_coefficient") %in% names(pl)))
PLA <- by_symbol(pl[, .(cohort = "Delahaye2018_placenta", cpg_id = as.character(CpG), sym = as.character(Gene), stat = as.numeric(linear_coefficient), p = as.numeric(P.value), strict = TRUE)])
g <- fread(path.expand(f_gtp)); stopifnot(all(c("CpG.probe", "annot.gene", "status", "p.val", "beta") %in% names(g)))
g <- g[status != "TRANS"]
GTP <- by_symbol(g[, .(cohort = "Kennedy2018_GTP_blood", cpg_id = as.character(CpG.probe), sym = as.character(annot.gene), stat = as.numeric(beta), p = as.numeric(p.val), strict = as.numeric(p.val) <= 1e-11)])
ad <- xl(f_adi, 1, 0); stopifnot(all(c("TargetID", "Gene Symbol", "Beta", "P-value") %in% names(ad)))
ADI <- by_symbol(ad[, .(cohort = "Voisin2016_adipose", cpg_id = as.character(TargetID), sym = as.character(`Gene Symbol`), stat = as.numeric(Beta), p = as.numeric(`P-value`), strict = TRUE)])
raw <- data.table(source = c("kidney (Liu 2022, Table S12)", "placenta (Delahaye 2018, S20)", "blood (Kennedy 2018, Additional file 1, not trans)", "adipose (Voisin 2016, S20)"),
                  rows_in_file = big(c(nrow(k12), nrow(pl), nrow(g), nrow(ad))), rows_with_a_gene = big(c(nrow(KID), nrow(PLA), nrow(GTP), nrow(ADI))))

# ---- position, window, sign, one vote per cohort -------------------------------------------------
W <- rbind(KID, PLA[, names(KID), with = FALSE], GTP[, names(KID), with = FALSE], ADI[, names(KID), with = FALSE])
W <- W[grepl("^c[gh]\\d", cpg_id) & is.finite(stat) & stat != 0]
W <- merge(W, P, by = "cpg_id", all.x = TRUE)
drops <- W[, .(no_position = sum(is.na(pos)), other_chromosome = sum(!is.na(pos) & chr != gchr),
               beyond_window = sum(!is.na(pos) & chr == gchr & abs(tss - pos) > window)), by = cohort]
W <- W[!is.na(pos) & chr == gchr]; W[, tss_dist := as.integer(round(tss - pos))]; W <- W[abs(tss_dist) <= window]
W[, direction := sign(stat)]
C <- W[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, tss_dist = tss_dist[1], gene_id = gene_id[1], gene_source = gene_source[1], pos_source = pos_source[1],
           max_abs_stat = max(abs(stat)), min_p = min(p), strict = any(strict)), by = .(cohort, cpg_id, target_gene)]
ties <- C[is.na(direction), .N, by = cohort]; C <- C[!is.na(direction)]
src_of <- c(Liu2022_kidney = "liu2022_kidney", Delahaye2018_placenta = "delahaye2018_placenta", Kennedy2018_GTP_blood = "kennedy2018_gtp", Voisin2016_adipose = "voisin2016_adipose")
tis_of <- c(Liu2022_kidney = "solid_tissue", Delahaye2018_placenta = "solid_tissue", Kennedy2018_GTP_blood = "blood", Voisin2016_adipose = "solid_tissue")
add <- C[, .(cpg_id, target_gene, tissue = unname(tis_of[cohort]), direction = as.integer(direction), tss_dist, n_cohorts = 1L, cohort_agreement = 1,
             cohorts = cohort, source = unname(src_of[cohort]), gene_id, gene_source, max_abs_stat = signif(max_abs_stat, 4), min_p = signif(min_p, 3), strict, pos_source)]
setorder(add, source, cpg_id, target_gene)
stopifnot(!anyDuplicated(add[, .(source, cpg_id, target_gene)]))
kept <- add[, .(pairs = big(.N), cpgs = big(uniqueN(cpg_id)), genes = big(uniqueN(target_gene)), negative = sprintf("%.3f", mean(direction == -1)),
                within_10kb = big(sum(abs(tss_dist) <= 1e4)), median_abs_dist = big(as.integer(median(abs(tss_dist))))), by = .(cohort = cohorts, tissue)]
print(kept)

# ---- validation against the shipped measured table -------------------------------------------------
M <- as.data.table(cpgd_measured_eqtms()); M[, target_gene := toupper(trimws(target_gene))]
stopifnot(!any(unname(src_of) %in% M$source))
anyM <- M[, .(s = sum(direction)), by = .(cpg_id, target_gene)][s != 0][, direction := sign(s)]
chk <- function(A, Bt) { x <- merge(A[, .(cpg_id, target_gene, d = direction)], unique(Bt[, .(cpg_id, target_gene, s = direction)], by = c("cpg_id", "target_gene")), by = c("cpg_id", "target_gene"))
  fmt(sum(x$d == x$s), nrow(x)) }
val <- rbindlist(lapply(names(src_of), function(co) { a <- add[cohorts == co]
  data.table(cohort = co, all_shared_measured = chk(a, anyM), blood = chk(a, M[tissue == "blood"]), solid = chk(a, M[tissue == "solid_tissue"]), nasal = chk(a, M[tissue == "nasal_epithelium"]),
             solid_2plus_cohorts = chk(a, M[tissue == "solid_tissue" & n_cohorts >= 2L]),
             new_pairs_in_its_tissue = big(nrow(a) - nrow(merge(a[, .(cpg_id, target_gene)], M[tissue == unname(tis_of[co]), .(cpg_id, target_gene)], by = c("cpg_id", "target_gene"))))) }))
print(val)
by_pos <- add[, { a <- .SD; data.table(pairs = big(.N), all_shared_measured = chk(a, anyM)) }, by = .(cohort = cohorts, pos_source), .SDcols = c("cpg_id", "target_gene", "direction")]
print(by_pos)
gt <- add[cohorts == "Kennedy2018_GTP_blood"]
gtp_strata <- sprintf("pairs at the paper's p <= 1e-11: %s, against shipped blood %s; pairs at 1e-11 < p < 1e-5: %s, against shipped blood %s",
                      big(sum(gt$strict)), chk(gt[strict == TRUE], M[tissue == "blood"]), big(sum(!gt$strict)), chk(gt[strict == FALSE], M[tissue == "blood"]))
cross <- rbindlist(lapply(combn(names(src_of), 2, simplify = FALSE), function(z) {
  x <- merge(add[cohorts == z[1], .(cpg_id, target_gene, a = direction)], add[cohorts == z[2], .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
  data.table(pair = paste(z, collapse = " x "), agreement = fmt(sum(x$a == x$b), nrow(x))) }))

# ---- write -----------------------------------------------------------------------------------------
fwrite(add, file.path(out_dir, "measured_eqtms_add_nov2026_sources.csv.gz"))
rep <- c("# Four measured sources admitted on 12 Nov 2026: kidney, placenta, GTP whole blood, adipose", "",
  sprintf("Run %s. Nothing in the package is changed. Window %s bp; tss_dist = TSS - CpG; one row per cohort and pair.", format(Sys.time(), "%Y-%m-%d %H:%M"), big(window)), "",
  "## Rows read", "", md(raw), "",
  "## Kidney table recovery (Liu 2022, Table S12)", "",
  sprintf("- The chromosome of the listed CpG ID matches its row in %.1f%% of rows: the CpG ID / Ensembl ID block is sorted differently from the statistics block.", 100 * kid_listed_chr),
  sprintf("- CpG recovered from Chr + CpG Position; %s rows sit on positions the package does not hold and are dropped.", big(kid_nopos)),
  sprintf("- Gene taken from the Ensembl ID paired with the recovered CpG in the other block (TSS within 5 kb of the row's TSS): %s rows; by the row's gene symbol otherwise; %s rows without a usable gene.", big(kid_by_block), big(kid_nogene)),
  sprintf("- Check on Table S11: %s.", kid_chk), "",
  "## CpG positions", "", hm_note, "", "## Rows dropped while placing the pairs", "", md(drops), "",
  if (nrow(ties)) c("Pairs dropped for opposite signs inside a cohort:", "", md(ties), "") else "No pair was dropped for opposite signs inside a cohort.", "",
  "## Kept", "", md(kept), "",
  "## Sign agreement with the shipped measured table (bar: 0.85 on all shared measured pairs)", "", md(val), "",
  paste0("GTP whole blood, ", gtp_strata, "."), "",
  "By the table that placed the CpG:", "", md(by_pos), "",
  "## Between the new sources", "", md(cross), "",
  "## Output", "",
  sprintf("- `measured_eqtms_add_nov2026_sources.csv.gz`: %s rows; columns of the measured table plus `gene_id`, `gene_source`, `max_abs_stat`, `min_p`, `strict`, `pos_source`.", big(nrow(add))),
  "- sources `liu2022_kidney`, `delahaye2018_placenta`, `kennedy2018_gtp`, `voisin2016_adipose`; each cohort is its own donor group.", "")
writeLines(rep, file.path(out_dir, "ingest_nov2026_measured_sources.md"))
message("wrote ", file.path(out_dir, "measured_eqtms_add_nov2026_sources.csv.gz"), " and ingest_nov2026_measured_sources.md")
