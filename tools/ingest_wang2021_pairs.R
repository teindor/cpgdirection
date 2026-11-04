#!/usr/bin/env Rscript
# Lead 18 (3 Nov 2026, v2): Wang ... Liu 2021, Epigenetics 16:1336 (10.1080/15592294.2021.1993607),
# Supplemental Table S2 -- correlated gene-CpG pairs (cGCPs; Spearman, BH FDR < 0.05, CpG within the
# gene +/- 10 kb, 450K, hg19) in six groups. S2's expression identifier differs by dataset:
#   brain (aged group)          ROSMAP, n=502      gene symbols
#   brain developmental group   Jaffe DLPFC, n=129 HEEBO oligo IDs (GEO GSE30272 platform)
#   brain aging group           Jaffe DLPFC, n=118 HEEBO oligo IDs
#   monocyte / Tcell            MESA, n=1,201/214  Illumina HT-12 v4 probe IDs (ILMN_)
#   liver                       Horvath, n=75      Affymetrix HuGene 1.1 ST transcript clusters
# v2 maps each to a symbol (illuminaHumanv4.db, hugene11sttranscriptcluster.db, GEOquery platform
# table cached in ~/wang2021/heebo_map.csv), then does what v1 did: brain -> separate brain resource;
# MESA -> measured blood rows; liver -> measured solid rows; sign = sign(corr); tss_dist = TSS - CpG
# from the package hg19 TSS table and S2's MAPINFO; validation against the shipped layers and the
# brain SMR layer by tier. Nothing in the package is changed.
#
# Usage (package root):
#   Rscript tools/ingest_wang2021_pairs.R [--xlsx ~/wang2021/Supplemental_Table_S2.xlsx] [--out tools/out] [--data-dir <dir>]
# Needs: readxl, AnnotationDbi, illuminaHumanv4.db, hugene11sttranscriptcluster.db, GEOquery (installed on demand).

suppressPackageStartupMessages({ library(data.table); library(readxl) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
xlsx <- path.expand(get_arg("--xlsx", "~/wang2021/Supplemental_Table_S2.xlsx"))
out_dir <- get_arg("--out", "tools/out")
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else suppressPackageStartupMessages(library(cpgdirection))
need_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) { message("installing ", p); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(p, ask = FALSE, update = FALSE) }
  suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE))
}

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
big <- function(x) format(x, big.mark = ",")

# ---- 1. read, group, map identifiers to symbols -------------------------------------------
W <- as.data.table(read_excel(xlsx, sheet = 1, skip = 1))
stopifnot(all(c("Tissue", "Ename", "Mname", "corr", "pvalue", "FDR", "MAPINFO") %in% names(W)))
W[, `:=`(cpg_id = Mname, eid = trimws(as.character(Ename)), corr = as.numeric(corr), FDR = as.numeric(FDR), pvalue = as.numeric(pvalue))]
W <- W[grepl("^cg\\d{8}$", cpg_id) & nzchar(eid) & !is.na(corr) & corr != 0 & !is.na(FDR) & FDR < 0.05]
tl <- tolower(W$Tissue)
W[, group := fcase(grepl("brain", tl) & grepl("develop", tl), "brain_developmental",
                   grepl("brain", tl) & grepl("aged", tl), "brain_aged",
                   grepl("brain", tl) & grepl("ag(e)?ing", tl), "brain_ageing",
                   grepl("t ?cell", tl), "mesa_tcell", grepl("mono", tl), "mesa_monocyte", grepl("liver", tl), "liver", default = "unmapped")]
stopifnot(!any(W$group == "unmapped"))
W[, target_gene := NA_character_]
W[group == "brain_aged", target_gene := toupper(eid)]
# MESA: Illumina HT-12 v4 probes
ilmn <- W[grepl("^ILMN_", eid), unique(eid)]
if (length(ilmn) && need_pkg("illuminaHumanv4.db")) {
  m <- suppressMessages(AnnotationDbi::mapIds(illuminaHumanv4.db::illuminaHumanv4.db, keys = ilmn, keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
  W[grepl("^ILMN_", eid), target_gene := toupper(unname(m[eid]))]
}
# liver: HuGene 1.1 ST transcript clusters
affy <- W[group == "liver" & grepl("^\\d+$", eid), unique(eid)]
if (length(affy) && need_pkg("hugene11sttranscriptcluster.db")) {
  m <- suppressMessages(AnnotationDbi::mapIds(hugene11sttranscriptcluster.db::hugene11sttranscriptcluster.db, keys = affy, keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
  W[group == "liver" & grepl("^\\d+$", eid), target_gene := toupper(unname(m[eid]))]
}
# Jaffe DLPFC: HEEBO oligos -> symbols from the GSE30272 platform table (cached)
heebo_cache <- file.path(dirname(xlsx), "heebo_map.csv")
hb <- W[grepl("^HEEBO", eid), unique(eid)]
if (length(hb)) {
  if (!file.exists(heebo_cache) && need_pkg("GEOquery")) {
    message("fetching the GSE30272 platform annotation (HEEBO) from GEO ...")
    g <- tryCatch(GEOquery::getGEO("GSE30272", GSEMatrix = TRUE, getGPL = TRUE), error = function(e) { message("GEO fetch failed: ", conditionMessage(e)); NULL })
    if (!is.null(g)) {
      fd <- as.data.table(Biobase::fData(g[[1]]))
      idc <- names(fd)[vapply(fd, function(v) mean(grepl("^HEEBO", as.character(v))) > 0.5, logical(1))][1]
      syc <- names(fd)[grepl("symbol", names(fd), ignore.case = TRUE)][1]
      if (is.na(syc)) syc <- names(fd)[grepl("^gene$|gene.?name|gene_assignment", names(fd), ignore.case = TRUE)][1]
      message("HEEBO platform columns used: id = ", idc, ", symbol = ", syc)
      hm <- unique(fd[, .(eid = as.character(get(idc)), sym = toupper(trimws(sub("[ /;,|].*$", "", as.character(get(syc))))))])[nzchar(sym) & sym != "NA"]
      fwrite(hm, heebo_cache)
    }
  }
  if (file.exists(heebo_cache)) { hm <- fread(heebo_cache); W <- merge(W, hm[, .(eid, hsym = sym)], by = "eid", all.x = TRUE); W[grepl("^HEEBO", eid), target_gene := hsym]; W[, hsym := NULL] }
}
W <- W[!is.na(target_gene) & nzchar(target_gene) & target_gene != "NA"]
message("rows with a gene symbol after mapping:"); print(W[, .N, by = group])
W[, tissue := fcase(group %like% "^brain", "brain", group %like% "^mesa", "blood", group == "liver", "solid_tissue")]
W[, direction := sign(corr)]
W[, cohort := fcase(group == "brain_developmental", "Jaffe_DLPFC_developmental", group == "brain_ageing", "Jaffe_DLPFC_ageing", group == "brain_aged", "ROSMAP_DLPFC_aged",
                    group == "mesa_monocyte", "MESA_monocytes", group == "mesa_tcell", "MESA_Tcells", group == "liver", "Horvath_liver")]
G <- as.data.table(cpgd_gene_tss()); G[, .gkey := gsub("_", "-", toupper(gene))]
W[, cpg_pos := as.numeric(MAPINFO)]
W[, .gkey := gsub("_", "-", target_gene)]
W <- merge(W, unique(G[, .(.gkey, tss = as.numeric(tss))], by = ".gkey"), by = ".gkey", all.x = TRUE)
W[, tss_dist := tss - cpg_pos]
message(sprintf("rows %s; tss_dist available for %.1f%%", big(nrow(W)), 100 * mean(!is.na(W$tss_dist))))
# several probes/oligos per gene: collapse to one row per (group, CpG, gene) by majority sign, ties dropped
W <- W[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, n_ids = .N, corr = corr[which.min(FDR)], pvalue = min(pvalue), FDR = min(FDR),
           tss_dist = tss_dist[1], tissue = tissue[1], cohort = cohort[1], gene_anno = GENEanno[1], cgi = CGIanno[1], same_tad = TAD[1]), by = .(group, cpg_id, target_gene)][!is.na(direction)]
grp <- W[, .(pairs = .N, cpgs = uniqueN(cpg_id), genes = uniqueN(target_gene), neg_share = round(mean(direction == -1), 3), median_abs_corr = round(median(abs(corr)), 3)), by = .(group, tissue)][order(tissue, group)]
print(grp)

# ---- 2. validation ---------------------------------------------------------------------------
M <- fread("inst/extdata/measured_eqtms.csv.gz"); M[, target_gene := toupper(trimws(target_gene))]
Mb <- M[tissue == "blood"]; Ms <- M[tissue == "solid_tissue"]
rows <- list()
for (g in c("mesa_monocyte", "mesa_tcell")) {
  x <- merge(W[group == g, .(cpg_id, target_gene, d = direction)], Mb[, .(cpg_id, target_gene, b = direction, cohorts)], by = c("cpg_id", "target_gene"))
  rows[[length(rows) + 1]] <- data.table(check = paste(g, "vs shipped blood"), agree = fmt(sum(x$d == x$b), nrow(x)))
  e <- x[grepl("EvoImmunoPop", cohorts)]
  rows[[length(rows) + 1]] <- data.table(check = paste(g, "vs EvoImmunoPop monocytes"), agree = fmt(sum(e$d == e$b), nrow(e)))
}
x <- merge(W[group == "mesa_monocyte", .(cpg_id, target_gene, a = direction)], W[group == "mesa_tcell", .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
rows[[length(rows) + 1]] <- data.table(check = "mesa_monocyte vs mesa_tcell", agree = fmt(sum(x$a == x$b), nrow(x)))
x <- merge(W[group == "liver", .(cpg_id, target_gene, d = direction)], Ms[, .(cpg_id, target_gene, s = direction)], by = c("cpg_id", "target_gene"))
rows[[length(rows) + 1]] <- data.table(check = "liver vs shipped solid_tissue", agree = fmt(sum(x$d == x$s), nrow(x)))
B <- W[tissue == "brain"]; bg <- unique(B$group)
for (i in seq_along(bg)) for (j in seq_along(bg)) if (i < j) {
  x <- merge(B[group == bg[i], .(cpg_id, target_gene, a = direction)], B[group == bg[j], .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
  rows[[length(rows) + 1]] <- data.table(check = paste(bg[i], "vs", bg[j]), agree = fmt(sum(x$a == x$b), nrow(x)))
}
x <- merge(B[, .(cpg_id, target_gene, group, d = direction)], Mb[, .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
for (g in bg) { y <- x[group == g]; rows[[length(rows) + 1]] <- data.table(check = paste(g, "vs measured blood"), agree = fmt(sum(y$d == y$b), nrow(y))) }
val <- rbindlist(rows); message("\n== sign agreement"); print(val)
BD <- tryCatch(as.data.table(cpgd_brain_directions()), error = function(e) { message("brain_directions: ", conditionMessage(e)); NULL })
bt <- bt_all <- near <- NULL
if (!is.null(BD)) {
  BD[, target_gene := toupper(trimws(target_gene))]
  x <- merge(B[, .(cpg_id, target_gene, group, d = direction)], BD[, .(cpg_id, target_gene, s = direction, tier, cpg_gene_dist)], by = c("cpg_id", "target_gene"))
  bt_all <- x[, .(n = .N, agree = fmt(sum(d == s), .N), baseline = sprintf("%.3f", max(mean(d == 1), mean(d == -1)))), by = tier][order(tier)]
  bt <- x[, .(n = .N, agree = fmt(sum(d == s), .N)), by = .(group, tier)][order(group, tier)]
  near <- x[abs(cpg_gene_dist) <= 100000, .(n = .N, agree = fmt(sum(d == s), .N)), by = tier][order(tier)]
  message("\n== brain measured (Wang S2) vs brain SMR layer, by tier"); print(bt_all); print(bt); message("within 100 kb:"); print(near)
}

# ---- 3. write ----------------------------------------------------------------------------------
add <- W[tissue %in% c("blood", "solid_tissue"), .(cpg_id, target_gene, tissue, direction, tss_dist, cohorts = cohort)]
add <- add[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, tss_dist = tss_dist[1], n_cohorts = .N,
               cohort_agreement = round(max(mean(direction == 1), mean(direction == -1)), 3), cohorts = paste(sort(unique(cohorts)), collapse = "|"), source = "wang2021_S2"),
           by = .(cpg_id, target_gene, tissue)][!is.na(direction)]
newb <- add[tissue == "blood"][!Mb, on = c("cpg_id", "target_gene")]
news <- add[tissue == "solid_tissue"][!Ms, on = c("cpg_id", "target_gene")]
brain <- B[, .(cpg_id, target_gene, group, cohort, direction, corr, pvalue, FDR, tss_dist, gene_anno, cgi, same_tad)]
fwrite(add, file.path(out_dir, "measured_eqtms_add_wang2021.csv.gz"))
fwrite(brain, file.path(out_dir, "brain_measured_eqtms_wang2021.csv.gz"))
rep <- c("# Wang et al. 2021 Supplemental Table S2 -- measured CpG-gene pairs for brain, sorted blood cells and liver (v2, identifiers mapped)", "",
         sprintf("Run %s. Source: %s (450K, hg19, Spearman, BH FDR < 0.05, CpG within gene +/- 10 kb).", format(Sys.time(), "%Y-%m-%d %H:%M"), xlsx), "",
         "## Groups (after mapping to symbols and collapsing probes/oligos per gene)", "", md(grp), "",
         "## Sign agreement checks", "", md(val), "",
         if (!is.null(bt_all)) c("## Brain measured vs brain SMR layer (cpgd_brain_directions), by tier", "", md(bt_all), "", "by group x tier:", "", md(bt), "", "within 100 kb of the gene:", "", md(near), "") else "## Brain vs SMR: brain_directions not available in this run",
         "", "## Outputs", "",
         sprintf("- `measured_eqtms_add_wang2021.csv.gz`: %s blood rows (MESA; %s new to the shipped blood table) and %s solid rows (liver; %s new).",
                 big(sum(add$tissue == "blood")), big(nrow(newb)), big(sum(add$tissue == "solid_tissue")), big(nrow(news))),
         sprintf("- `brain_measured_eqtms_wang2021.csv.gz`: %s brain rows over %s CpGs (developmental %s, ageing %s, aged %s) -- separate resource for cpg_brain_bridge().",
                 big(nrow(brain)), big(uniqueN(brain$cpg_id)), big(sum(brain$group == "brain_developmental")), big(sum(brain$group == "brain_ageing")), big(sum(brain$group == "brain_aged"))), "")
writeLines(rep, file.path(out_dir, "ingest_wang2021_pairs.md"))
message("\nwrote ", file.path(out_dir, "measured_eqtms_add_wang2021.csv.gz"), ", brain_measured_eqtms_wang2021.csv.gz and ingest_wang2021_pairs.md")
