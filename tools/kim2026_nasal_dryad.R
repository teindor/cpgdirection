#!/usr/bin/env Rscript
# Lead 2 (Nov 2026 scan, a1): Kim et al. 2026 Sci Adv 12:eady3555 -- EVA-PR nasal
# epithelium eQTM / eQTL / meQTL summary statistics on Dryad (CC0),
# https://doi.org/10.5061/dryad.rxwdbrvpq (version 6, 25 Mar 2026).
#
# What this script settles, on the Mac (the cloud sandbox cannot reach Dryad):
#   1. downloads the small files through the Dryad API (README, cis-eQTM 3.4 MB,
#      cis-eQTL 6.5 MB, cis-meQTL 54 MB, TableS12rev, TableS13rev; the 1.1 GB
#      trans-eQTM only with --trans);
#   2. prints the README and every file's header, so the column meanings are on
#      record (the sandbox never saw them);
#   3. looks for the per-pair causal classification (DNAm->expression /
#      expression->DNAm / bidirectional, 23,948 / 9,638 / 1,305). The paper cites
#      NO table or data file for it (tables S1-S13 are eQTM/eCpG/eGene lists,
#      adjusted QTLs and the Soliai replication), so expect to find nothing; if a
#      column of that kind exists anywhere, it is reported with its value counts;
#   4. compares the 2026 cis-eQTM signs (FDR < 0.01, n = 455, M-value predictor,
#      TPM outcome) with the nasal pairs already shipped in cpgdirection
#      (cpgd_measured_eqtms(), tissue == "nasal_epithelium", built from the eQTM
#      Atlas exports of the SAME cohort, Kim 2020/2023). Same cohort => this is a
#      re-analysis, not replication; the comparison measures how much the stricter
#      threshold and the TPM model move the shipped signs, nothing more;
#   5. compares the 2026 nasal signs with the blood SMR layer (eQTLGen outcome,
#      and FinnGen rows if the rebuilt table is in place) as a cross-tissue check.
#
# Usage (package root):
#   Rscript tools/kim2026_nasal_dryad.R [--dir ~/kim2026_dryad] [--trans] [--data-dir <dir>]
# Writes <dir>/kim2026_nasal_dryad.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
dir      <- path.expand(get_arg("--dir", "~/kim2026_dryad"))
want_trans <- "--trans" %in% args
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)
dir.create(dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else suppressPackageStartupMessages(library(cpgdirection))

# ---- 1. download ---------------------------------------------------------------------
files <- data.table(
  name = c("README.md", "cis-eQTM.nocelltypes.fdr0.01.Original.TPM.csv",
           "cis-eQTL.nocelltypes.fdr0.01.Original.TPM.csv",
           "cis-meQTL.nocelltypes.fdr0.01.Original.TPM.csv",
           "TableS12rev.csv", "TableS13rev.csv",
           "trans-eQTL.nocelltypes.fdr0.01.Original.TPM.csv",
           "trans-meQTL.nocelltypes.fdr0.01.Original.TPM.csv",
           "trans-eQTM.nocelltypes.fdr0.01.Original.TPM.csv"),
  id   = c(4672010, 4671961, 4671957, 4671960, 4671963, 4671964, 4671958, 4671959, 4671962),
  size = c(4730, 3394542, 6494968, 53993290, 362811, 546810, 235473, 2429426, 1138819072))
if (!want_trans) files <- files[!grepl("^trans-eQTM", name)]
for (k in seq_len(nrow(files))) {
  dest <- file.path(dir, files$name[k])
  if (file.exists(dest) && file.size(dest) > 1000) next
  url <- sprintf("https://datadryad.org/downloads/file_stream/%d", files$id[k])
  message("downloading ", files$name[k], " (", format(files$size[k], big.mark = ","), " bytes)")
  st <- system2("curl", c("-sS", "-L", "--fail", "--retry", "3", "-A", shQuote("Mozilla/5.0 (Macintosh; Intel Mac OS X 14_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130 Safari/537.36"), "-o", shQuote(dest), shQuote(url)))
  if (st != 0 || !file.exists(dest)) stop("download failed: ", url)
  if (file.size(dest) != files$size[k])
    warning(sprintf("%s: got %s bytes, Dryad lists %s", files$name[k], file.size(dest), files$size[k]))
}
rep <- c("# Kim et al. 2026 (EVA-PR nasal) Dryad files -- inspection and comparison", "",
         sprintf("Run %s; files in `%s`.", format(Sys.time(), "%Y-%m-%d %H:%M"), dir), "")

# ---- 2. README and headers -------------------------------------------------------------
readme <- readLines(file.path(dir, "README.md"), warn = FALSE)
cat(readme, sep = "\n"); rep <- c(rep, "## README.md (verbatim)", "", "```", readme, "```", "")
peek <- function(f, n = 3) {
  d <- fread(f, nrows = n, showProgress = FALSE)
  wc <- trimws(system2("wc", c("-l", shQuote(f)), stdout = TRUE))
  nr <- as.integer(sub(" .*", "", wc)) - 1L
  message(sprintf("\n== %s: %s rows, %d columns", basename(f), format(nr, big.mark = ","), ncol(d)))
  print(d)
  c(sprintf("### %s -- %s rows", basename(f), format(nr, big.mark = ",")), "",
    paste0("columns: `", paste(names(d), collapse = "`, `"), "`"), "")
}
for (f in files$name[-1]) rep <- c(rep, peek(file.path(dir, f)))

# ---- 3. any causal-classification column anywhere? -----------------------------------
class_pat <- "mediat|direction|causal|bidirect|DNAm.*(->|to).*(expr|gene)|JS|AB_|class"
found <- character()
for (f in files$name[-1]) {
  d <- fread(file.path(dir, f), nrows = 50000, showProgress = FALSE)
  hits <- names(d)[grepl(class_pat, names(d), ignore.case = TRUE)]
  for (h in hits) {
    v <- d[[h]]
    if (is.character(v) || is.factor(v) || uniqueN(v) <= 6) {
      tab <- sort(table(v), decreasing = TRUE)[1:min(6, uniqueN(v))]
      line <- sprintf("%s / %s: %s", basename(f), h, paste(names(tab), tab, sep = "=", collapse = "; "))
      message("candidate classification column: ", line); found <- c(found, line)
    }
  }
}
rep <- c(rep, "## Causal-classification columns", "",
         if (length(found)) paste0("- ", found) else
           paste("None. The 23,948 / 9,638 / 1,305 JS-AB classification is not in the Dryad deposit",
                 "(the paper cites no table for it either): request it from S. Kim / J. Celedon."), "")

# ---- 4. cis-eQTM signs vs the shipped nasal pairs (same cohort) ------------------------
Q <- fread(file.path(dir, "cis-eQTM.nocelltypes.fdr0.01.Original.TPM.csv"), showProgress = FALSE)
pick <- function(d, pat, test = NULL) {
  nm <- names(d)[grepl(pat, names(d), ignore.case = TRUE)]
  if (!is.null(test)) nm <- nm[vapply(nm, function(x) test(d[[x]]), logical(1))]
  if (length(nm)) nm[1] else NA_character_
}
cpg_col  <- pick(Q, ".", function(v) is.character(v) && mean(grepl("^cg\\d{8}$", head(v, 1000))) > 0.9)
gene_col <- pick(Q, "gene|symbol|ensembl|feature|target", function(v) is.character(v))
beta_col <- pick(Q, "^(beta|estimate|coef|slope|b|effect)", function(v) is.numeric(v))
if (is.na(beta_col)) beta_col <- pick(Q, "beta|estimate|coef|slope|effect|statistic|^t$", function(v) is.numeric(v))
message(sprintf("\ncis-eQTM columns used: cpg=%s gene=%s effect=%s", cpg_col, gene_col, beta_col))
stopifnot(!is.na(cpg_col), !is.na(gene_col), !is.na(beta_col))
K <- Q[, .(cpg_id = get(cpg_col), gene = as.character(get(gene_col)), kim_dir = sign(get(beta_col)))]
K <- K[!is.na(kim_dir) & kim_dir != 0]
# gene identifiers: the package keys on symbols; map Ensembl ids if that is what we got
if (mean(grepl("^ENSG", head(K$gene, 1000))) > 0.5) {
  K[, ensembl := sub("\\..*$", "", gene)]
  sym <- tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(K$ensembl), keytype = "ENSEMBL",
                          column = "SYMBOL", multiVals = "first")
  }, error = function(e) NULL)
  if (is.null(sym)) stop("gene column is Ensembl; install org.Hs.eg.db to map to symbols")
  K[, gene := unname(sym[ensembl])]
  K <- K[!is.na(gene)]
} else if (any(grepl("_", K$gene))) {
  # compound labels: the informative symbol is the last token (package convention)
  K[, gene := sub(".*_", "", gene)]
}
K[, gene := toupper(trimws(gene))]
K <- K[, .(kim_dir = if (uniqueN(kim_dir) == 1L) kim_dir[1] else NA_real_), by = .(cpg_id, gene)][!is.na(kim_dir)]
rep <- c(rep, "## cis-eQTM (FDR < 0.01) as a sign table", "",
         sprintf("%s CpG-gene pairs with a sign (columns `%s`, `%s`, `%s`); %s CpGs, %s genes; negative %.1f%%.",
                 format(nrow(K), big.mark = ","), cpg_col, gene_col, beta_col,
                 format(uniqueN(K$cpg_id), big.mark = ","), format(uniqueN(K$gene), big.mark = ","),
                 100 * mean(K$kim_dir < 0)), "")

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }

M <- as.data.table(cpgd_measured_eqtms())
M[, gene := toupper(trimws(target_gene))]
Mn <- M[tissue == "nasal_epithelium", .(cpg_id, gene, ship_dir = direction, cohorts, n_cohorts)]
Mb <- M[tissue == "blood", .(cpg_id, gene, blood_dir = direction)]
Ms <- M[tissue == "solid_tissue", .(cpg_id, gene, solid_dir = direction)]
j <- merge(K, Mn, by = c("cpg_id", "gene"))
message(sprintf("\n2026 nasal vs shipped nasal (same cohort): overlap %s of %s shipped nasal pairs; sign agreement %s",
                format(nrow(j), big.mark = ","), format(nrow(Mn), big.mark = ","), fmt(sum(j$kim_dir == j$ship_dir), nrow(j))))
new_pairs <- K[!Mn, on = c("cpg_id", "gene")]
jb <- merge(K, Mb, by = c("cpg_id", "gene")); js <- merge(K, Ms, by = c("cpg_id", "gene"))
jbs <- merge(Mn, Mb, by = c("cpg_id", "gene"))
rep <- c(rep, "## Against the shipped catalogues", "",
         sprintf("- shipped nasal pairs: %s; 2026 cis-eQTM covers %s of them; sign agreement %s (same cohort, re-analysis).",
                 format(nrow(Mn), big.mark = ","), format(nrow(j), big.mark = ","), fmt(sum(j$kim_dir == j$ship_dir), nrow(j))),
         sprintf("- 2026 pairs absent from the shipped nasal table: %s (%s CpGs).",
                 format(nrow(new_pairs), big.mark = ","), format(uniqueN(new_pairs$cpg_id), big.mark = ",")),
         sprintf("- 2026 nasal vs shipped blood: %s; shipped nasal vs shipped blood, for reference: %s.",
                 fmt(sum(jb$kim_dir == jb$blood_dir), nrow(jb)), fmt(sum(jbs$ship_dir == jbs$blood_dir), nrow(jbs))),
         sprintf("- 2026 nasal vs shipped solid tissue: %s.", fmt(sum(js$kim_dir == js$solid_dir), nrow(js))), "")

# ---- 5. cross-tissue: 2026 nasal signs vs the blood SMR layer ------------------------
S <- as.data.table(cpgd_smr_directions())
S[, gene := toupper(trimws(target_gene))]
if (!"smr_rung" %in% names(S)) S[, smr_rung := fcase(smr_tier == "S1", "high", smr_tier == "S2", "moderate", default = "weak")]
if (!"outcome" %in% names(S)) S[, outcome := "eQTLGen"]
sj <- merge(K, S[, .(cpg_id, gene, smr_dir = direction, smr_tier, smr_rung, outcome)], by = c("cpg_id", "gene"))
by_rung <- sj[, .(n = .N, agree = fmt(sum(kim_dir == smr_dir), .N)), by = .(smr_rung, outcome = sub("_.*", "", outcome))][order(smr_rung, outcome)]
message("\n== 2026 nasal signs vs blood SMR, by rung"); print(by_rung)
# the same for the shipped nasal pairs, so the two references can be compared like for like
sn <- merge(Mn, S[, .(cpg_id, gene, smr_dir = direction, smr_rung)], by = c("cpg_id", "gene"))
by_rung_ship <- sn[, .(n = .N, agree = fmt(sum(ship_dir == smr_dir), .N)), by = smr_rung][order(smr_rung)]
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
rep <- c(rep, "## 2026 nasal signs vs the blood SMR layer (cross-tissue)", "", md(by_rung), "",
         "Shipped nasal pairs vs the same SMR rows, for comparison:", "", md(by_rung_ship), "")

fwrite(K, file.path(dir, "kim2026_cis_eqtm_signs.csv.gz"), compress = "gzip")
writeLines(rep, file.path(dir, "kim2026_nasal_dryad.md"))
message("\nwrote ", file.path(dir, "kim2026_nasal_dryad.md"), " and kim2026_cis_eqtm_signs.csv.gz")
