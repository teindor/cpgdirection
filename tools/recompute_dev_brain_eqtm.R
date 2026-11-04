#!/usr/bin/env Rscript
# Lead 18b (4 Nov 2026): recompute the Jaffe developmental-cortex eQTMs with age removed.
# Wang 2021 S2 developmental pairs were Spearman correlations across fetal-to-25-year DLPFC with age NOT regressed,
# so their signs can reflect two developmental trends rather than cis regulation. This script re-tests exactly
# those CpG-gene pairs (plus the Jaffe ageing pairs) in the source data with age modelled:
#   expr ~ meth + ns(age, df) + fetal + sex + race          (all developmental samples)
#   expr ~ meth + ns(age, df) + sex + race                  (postnatal 1-25 y only: the Project Alpha age range)
# and compares the age-adjusted signs to the unadjusted Wang signs, to ROSMAP aged cortex and to measured blood.
#
# Data (GEO): methylation GSE74193 (Jaffe 2016, 450K DLPFC); expression GSE30272 (Colantuoni 2011, HEEBO, BrainCloud).
# Samples are paired by BrainCloud donor number (BrNNNN) found in titles/characteristics; if that fails, pass a
# two-column CSV (meth_gsm, expr_gsm) with --pairing. Large download: GSE74193 processed betas are several GB;
# if the series matrix carries no values, download GSE74193_GEO_procData.csv.gz and pass --meth.
#
# Usage (package root):
#   Rscript tools/recompute_dev_brain_eqtm.R [--wang tools/out/brain_measured_eqtms_wang2021.csv.gz]
#          [--heebo ~/wang2021/heebo_map.csv] [--meth path/GSE74193_GEO_procData.csv.gz] [--pairing pairs.csv]
#          [--age-df 4] [--cache ~/lead18b]
# Writes tools/out/brain_dev_eqtm_ageadj.csv.gz and tools/reports/recompute_dev_brain_eqtm.md. Changes nothing shipped.

suppressPackageStartupMessages({ library(data.table); library(splines) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
wang_f  <- get_arg("--wang", "tools/out/brain_measured_eqtms_wang2021.csv.gz")
heebo_f <- path.expand(get_arg("--heebo", "~/wang2021/heebo_map.csv"))
meth_f  <- get_arg("--meth", NULL); pair_f <- get_arg("--pairing", NULL)
age_df  <- as.integer(get_arg("--age-df", "4"))
cache   <- path.expand(get_arg("--cache", "~/lead18b")); dir.create(cache, showWarnings = FALSE, recursive = TRUE)
for (p in c("GEOquery", "Biobase")) if (!requireNamespace(p, quietly = TRUE)) stop("install ", p, " (BiocManager::install('", p, "'))")
wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }

# ---- 1. pairs to test --------------------------------------------------------------------------
Wg <- fread(wang_f)[group %in% c("brain_developmental", "brain_ageing")]
Wg[, target_gene := toupper(trimws(target_gene))]
hm <- fread(heebo_f)[, .(eid = as.character(eid), sym = toupper(sym))]
cpgs <- unique(Wg$cpg_id); genes <- unique(Wg$target_gene)
message(sprintf("pairs to re-test: %d (%d CpGs, %d genes)", nrow(unique(Wg[, .(cpg_id, target_gene)])), length(cpgs), length(genes)))

# ---- 2. GEO: phenotypes, expression, methylation -----------------------------------------------
pheno <- function(es) { pd <- as.data.table(Biobase::pData(es), keep.rownames = "gsm")
  txt <- do.call(paste, c(pd[, lapply(.SD, as.character)], sep = " | "))
  getn <- function(re) suppressWarnings(as.numeric(sub(re, "\\1", ignore.case = TRUE, perl = TRUE, regmatches(txt, regexpr(re, txt, ignore.case = TRUE, perl = TRUE)))[match(seq_along(txt), which(grepl(re, txt, ignore.case = TRUE, perl = TRUE)))]))
  data.table(gsm = pd$gsm, title = if ("title" %in% names(pd)) as.character(pd$title) else NA_character_, donor = toupper(regmatches(txt, regexpr("Br[0-9]{3,5}", txt, ignore.case = TRUE))[match(seq_along(txt), which(grepl("Br[0-9]{3,5}", txt, ignore.case = TRUE)))]),
             age = getn("age[^|:]*:\\s*(-?[0-9.]+)"),
             sex = ifelse(grepl("(sex|gender)[^|:]*:", txt, ignore.case = TRUE), toupper(substr(sub(".*(sex|gender)[^|:]*:\\s*([A-Za-z]+).*", "\\2", txt, ignore.case = TRUE), 1, 1)), NA_character_),
             race = ifelse(grepl("race[^|:]*:", txt, ignore.case = TRUE), trimws(sub(".*race[^|:]*:\\s*([^|]+).*", "\\1", txt, ignore.case = TRUE)), NA_character_), txt = txt) }
ge <- GEOquery::getGEO("GSE30272", GSEMatrix = TRUE, getGPL = FALSE, destdir = cache)[[1]]
PE <- pheno(ge); X <- Biobase::exprs(ge)
ids <- intersect(rownames(X), hm[sym %in% genes, eid]); X <- X[ids, , drop = FALSE]
gm <- hm[eid %in% ids]
gm_g <- split(gm$eid, gm$sym)
message(sprintf("expression: %d samples, %d HEEBO probes for %d target genes", ncol(X), nrow(X), length(gm_g)))

gm_m <- GEOquery::getGEO("GSE74193", GSEMatrix = TRUE, getGPL = FALSE, destdir = cache)[[1]]
PM <- pheno(gm_m)
B <- Biobase::exprs(gm_m)
if (!nrow(B) && !is.null(meth_f)) {
  B <- fread(meth_f, showProgress = TRUE)
  idc <- names(B)[1]; B <- B[get(idc) %in% cpgs]; rn <- B[[idc]]
  bc <- grep("_Beta$", names(B), value = TRUE)
  B <- as.matrix(B[, bc, with = FALSE]); rownames(B) <- rn; colnames(B) <- sub("_Beta$", "", bc)
  message("procData columns are sample names, not GSMs: if pairing fails, map them with --pairing (meth column name, expr GSM)")
}
if (!nrow(B)) stop("GSE74193 series matrix has no beta values: download GSE74193_GEO_procData.csv.gz and pass --meth")
if (!any(colnames(B) %in% PM$gsm)) { k <- match(colnames(B), sub("_.*$", "", PM$title)); colnames(B)[!is.na(k)] <- PM$gsm[k[!is.na(k)]]
  message(sprintf("procData columns matched to GSMs by sample title: %d of %d", sum(!is.na(k)), length(k))) }
B <- B[intersect(rownames(B), cpgs), , drop = FALSE]
message(sprintf("methylation: %d samples, %d target CpGs", ncol(B), nrow(B)))

# ---- 3. pair samples ------------------------------------------------------------------------
# GSE74193 includes schizophrenia cases and QC-failed arrays: keep controls that pass QC, one array per donor (bestqc first).
PM[, `:=`(control = grepl("group: Control", txt, fixed = TRUE), drop = grepl("dropsample[^|]*: TRUE", txt), best = grepl("bestqc[^|]*: TRUE", txt))]
PMc <- PM[control & !drop & gsm %in% colnames(B)][order(-best)]
PMc <- unique(PMc, by = "donor")
message(sprintf("GSE74193 controls passing QC, one per donor: %d of %d arrays", nrow(PMc), nrow(PM)))
fp <- function(D) D[, fp := paste(sprintf("%.3f", age), sex, toupper(trimws(race)))][]
if (!is.null(pair_f)) { P <- fread(pair_f); setnames(P, 1:2, c("meth_gsm", "expr_gsm")) } else if (sum(!is.na(PE$donor)) >= 50) {
  P <- merge(PMc[!is.na(donor), .(meth_gsm = gsm, donor)], PE[!is.na(donor), .(expr_gsm = gsm, donor)], by = "donor")
} else {
  # GSE30272 carries no BrNum: pair by demographic fingerprint (age to 3 decimals, sex, race); keep only 1:1 matches
  a <- fp(copy(PMc[!is.na(age)])); b <- fp(copy(PE[!is.na(age)]))
  a1 <- a[, if (.N == 1) .SD, by = fp]; b1 <- b[, if (.N == 1) .SD, by = fp]
  P <- merge(a1[, .(fp, meth_gsm = gsm)], b1[, .(fp, expr_gsm = gsm)], by = "fp")
  agechk <- merge(P, a[, .(meth_gsm = gsm, am = age)], by = "meth_gsm")[, max(abs(am - as.numeric(sub(" .*", "", fp))))]
  message(sprintf("fingerprint pairing: %d unique-fingerprint donors in methylation, %d in expression, %d matched 1:1 (max age diff %.4f)",
                  nrow(a1), nrow(b1), nrow(P), agechk))
  fwrite(P, file.path(cache, "pairing_fingerprint.csv"))
}
P <- P[meth_gsm %in% colnames(B) & expr_gsm %in% colnames(X)]
P <- unique(P, by = "meth_gsm"); P <- unique(P, by = "expr_gsm")
P <- merge(P, PE[, .(expr_gsm = gsm, age, sex, race)], by = "expr_gsm")
P[is.na(age), age := PM$age[match(meth_gsm[is.na(age)], PM$gsm)]]
P <- P[!is.na(age)]
P[, fetal := as.integer(age < 0)]
message(sprintf("paired donors: %d (fetal %d, 0-25 y %d, >25 y %d)", nrow(P), sum(P$fetal), sum(P$age >= 0 & P$age <= 25), sum(P$age > 25)))
if (nrow(P) < 50) {
  fwrite(PM[, .(gsm, donor, age, txt)], file.path(cache, "GSE74193_pheno.csv")); fwrite(PE[, .(gsm, donor, age, txt)], file.path(cache, "GSE30272_pheno.csv"))
  stop("fewer than 50 paired samples: inspect ", cache, "/GSE*_pheno.csv, build a meth_gsm,expr_gsm CSV and pass --pairing")
}

# ---- 4. re-test every pair ---------------------------------------------------------------------
fitset <- function(S, label, use_fetal) {
  if (nrow(S) < 30) return(NULL)
  bm <- B[, S$meth_gsm, drop = FALSE]; xm <- X[, S$expr_gsm, drop = FALSE]
  cov <- data.frame(age = S$age, sex = factor(S$sex), race = factor(S$race), fetal = S$fetal)
  keep <- vapply(cov, function(v) length(unique(v[!is.na(v)])) > 1, logical(1)); cov <- cov[, keep, drop = FALSE]
  f_adj <- paste("~ m + ns(age, df =", age_df, ")", if (use_fetal && "fetal" %in% names(cov)) "+ fetal" else "",
                 if ("sex" %in% names(cov)) "+ sex" else "", if ("race" %in% names(cov)) "+ race" else "")
  out <- vector("list", nrow(Wp)); i <- 0L
  for (j in seq_len(nrow(Wp))) {
    cg <- Wp$cpg_id[j]; g <- Wp$target_gene[j]
    if (!cg %in% rownames(bm) || is.null(gm_g[[g]])) next
    y <- colMeans(xm[gm_g[[g]], , drop = FALSE], na.rm = TRUE); m <- bm[cg, ]
    d <- cbind(cov, y = y, m = m); d <- d[complete.cases(d), , drop = FALSE]; if (nrow(d) < 30 || sd(d$m) == 0) next
    a <- summary(lm(as.formula(paste("y", f_adj)), d))$coefficients["m", ]
    u <- summary(lm(y ~ m, d))$coefficients["m", ]
    i <- i + 1L; out[[i]] <- data.table(set = label, cpg_id = cg, target_gene = g, n = nrow(d),
      beta_unadj = u[1], p_unadj = u[4], beta_adj = a[1], p_adj = a[4])
  }
  R <- rbindlist(out[seq_len(i)]); R[, fdr_adj := p.adjust(p_adj, "BH")]; R
}
Wp <- unique(Wg[, .(cpg_id, target_gene)])
R <- rbindlist(list(fitset(P, "fetal_to_25", TRUE),
                    fitset(P[age >= 1 & age <= 25], "postnatal_1_25", FALSE),
                    fitset(P[age > 25], "adult_over_25", FALSE)))
R[, `:=`(dir_unadj = as.integer(sign(beta_unadj)), dir_adj = as.integer(sign(beta_adj)))]
dir.create("tools/out", showWarnings = FALSE, recursive = TRUE); dir.create("tools/reports", showWarnings = FALSE)
fwrite(R, "tools/out/brain_dev_eqtm_ageadj.csv.gz")

# ---- 5. compare --------------------------------------------------------------------------------
ref <- list(wang_dev_unadj = Wg[group == "brain_developmental", .(cpg_id, target_gene, r = direction)],
            wang_ageing    = Wg[group == "brain_ageing", .(cpg_id, target_gene, r = direction)])
bf <- "inst/extdata/brain_measured_eqtms.csv.gz"; if (file.exists(bf)) ref$rosmap_aged <- fread(bf)[, .(cpg_id, target_gene = toupper(target_gene), r = direction)]
mf <- "inst/extdata/measured_eqtms.csv.gz"; if (file.exists(mf)) ref$measured_blood <- fread(mf)[tissue == "blood", .(cpg_id, target_gene = toupper(target_gene), r = direction)]
rows <- list()
for (s in unique(R$set)) for (lev in c("all", "nominal_p05", "fdr05")) {
  Z <- R[set == s]; Z <- switch(lev, all = Z, nominal_p05 = Z[p_adj < 0.05], fdr05 = Z[fdr_adj < 0.05])
  rows[[length(rows) + 1]] <- data.table(set = s, filter = lev, pairs = nrow(Z),
    adj_vs_unadj_same_data = fmt(sum(Z$dir_adj == Z$dir_unadj), nrow(Z)),
    attenuation_median = sprintf("%.2f", median(abs(Z$beta_adj) / pmax(abs(Z$beta_unadj), 1e-12))))
  for (rn in names(ref)) { x <- merge(Z[, .(cpg_id, target_gene, d = dir_adj)], ref[[rn]], by = c("cpg_id", "target_gene"))
    rows[[length(rows)]][, (rn) := fmt(sum(x$d == x$r), nrow(x))] }
}
T <- rbindlist(rows, fill = TRUE); print(T)

rep <- c("# Lead 18b: developmental cortex eQTMs with age removed", "",
  sprintf("Run %s. Paired donors %d (fetal %d; 1-25 y %d; >25 y %d). Model: expr ~ meth + ns(age, df=%d) [+ fetal] + sex + race; HEEBO probes averaged per gene.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(P), sum(P$fetal), sum(P$age >= 1 & P$age <= 25), sum(P$age > 25), age_df), "",
  md(T), "",
  "## Decision rule", "",
  "- If postnatal 1-25 y age-adjusted signs (FDR < 0.05) agree with ROSMAP aged cortex at >= 0.85, the 0.73 figure was age confounding:",
  "  the adult brain directions transport to youth and the developmental caveat can be softened to fetal/infant only.",
  "- If agreement stays near 0.73 after adjustment, the developmental difference is real biology: keep the caveat and consider a",
  "  separate developmental measured rung built from these age-adjusted pairs (cohort Jaffe_DLPFC_dev_ageadj).",
  "- If few pairs survive adjustment, report the developmental layer as unvalidated and do not ship it.", "")
writeLines(rep, "tools/reports/recompute_dev_brain_eqtm.md")
message("wrote tools/out/brain_dev_eqtm_ageadj.csv.gz and tools/reports/recompute_dev_brain_eqtm.md")
