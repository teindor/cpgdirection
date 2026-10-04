#!/usr/bin/env Rscript
# Validate the Onco-eQTM consensus layer against the package's peripheral
# references, to decide whether any tier earns a rung on the evidence ladder.
#
# The question is narrow: when a tumour consensus sign exists for a CpG-gene
# pair that ALSO has an independently measured peripheral sign, how often do
# they agree, and does that beat the majority-class baseline and the rungs the
# ladder already has (catalogue_single 0.62-0.87, catalogue_consensus 0.77-0.87,
# smr_moderate 0.84-0.86)?
#
# Reference sets, each evaluated on its own (they are independent data sets; the
# measured table carries no finer cohort key, so "leave one out" here means
# leave one REFERENCE SET out):
#   blood    measured eQTMs, whole blood
#   nasal    measured eQTMs, nasal epithelium
#   solid    measured eQTMs, solid tissue
#   smr_S1   two-sample SMR, tier S1 (validated 0.959 against measured)
#
# Usage
#   Rscript tools/validate_onco.R --data-dir <hub_upload dir> \
#       --out tools/reports/validate_onco.md
#
# Prints the tables and writes them as markdown. Nothing in the package is
# changed by this script.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
data_dir <- get_arg("--data-dir")
out_md   <- get_arg("--out", "tools/reports/validate_onco.md")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)

# package code: an installed cpgdirection, else load_all() from the package root,
# else (no Bioconductor deps on this machine) source R/ directly -- the data
# accessors used here need only data.table.
# Prefer the package root (the code under validation), not whatever older
# cpgdirection happens to be installed in the library.
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else if (requireNamespace("cpgdirection", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cpgdirection))
} else stop("run from the cpgdirection package root or install cpgdirection")
if (!exists("cpgd_distal_links")) stop("cpgdirection < 2.99.6 loaded; run from the package root")

wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(c - h, c + h)
}
fmt_ci <- function(k, n) {
  if (n == 0) return("-")
  ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2])
}
md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  h <- paste0("| ", paste(names(d), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  r <- apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(h, s, r), collapse = "\n")
}

# ---- data --------------------------------------------------------------------------
O <- as.data.table(cpgd_onco_eqtm())
if (is.null(O) || !nrow(O)) stop("onco_eqtm_consensus layer not available; set --data-dir")
O[, target_gene := toupper(trimws(target_gene))]
# a usable tumour sign: the consensus direction where one exists, otherwise the
# sign of the median r (single-cancer O4 and the disagreeing O3 majority).
O[, onco_sign := ifelse(!is.na(direction), direction, sign(median_r))]
O <- O[!is.na(onco_sign) & onco_sign != 0]
O[, tss_bin := cut(abs(tss_dist_approx), c(-Inf, 1500, 10000, 100000, Inf),
                   labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb"))]

M <- as.data.table(cpgd_measured_eqtms())
# a measured gene field may hold "A;B": one reference row per gene
M <- M[, list(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
       by = c("cpg_id", "tissue", "direction", "tss_dist")]
M <- unique(M[nzchar(target_gene)])

S <- tryCatch(as.data.table(cpgd_smr_directions()), error = function(e) NULL)
refs <- list(
  blood  = M[tissue == "blood",            list(cpg_id, target_gene, ref_dir = direction)],
  nasal  = M[tissue == "nasal_epithelium", list(cpg_id, target_gene, ref_dir = direction)],
  solid  = M[tissue == "solid_tissue",     list(cpg_id, target_gene, ref_dir = direction)])
if (!is.null(S) && "smr_tier" %in% names(S)) {
  refs$smr_S1 <- S[smr_tier == "S1", list(cpg_id, target_gene = toupper(target_gene), ref_dir = direction)]
}
refs <- lapply(refs, function(r) unique(r[!is.na(ref_dir) & ref_dir != 0]))

# ---- evaluation ----------------------------------------------------------------------
eval_block <- function(J, by) {
  J[, list(n = .N,
           agree = fmt_ci(sum(onco_sign == ref_dir), .N),
           baseline = sprintf("%.3f", max(mean(ref_dir == 1), mean(ref_dir == -1))),
           onco_pos = sprintf("%.2f", mean(onco_sign == 1)),
           ref_pos  = sprintf("%.2f", mean(ref_dir == 1))),
    by = by][order(get(by[1]))]
}
ppv_block <- function(J, by) {
  # the value of each tumour sign separately: P(reference sign == s | tumour sign == s)
  J[, list(n_onco_pos = sum(onco_sign == 1),
           ppv_pos = fmt_ci(sum(onco_sign == 1 & ref_dir == 1), sum(onco_sign == 1)),
           n_onco_neg = sum(onco_sign == -1),
           ppv_neg = fmt_ci(sum(onco_sign == -1 & ref_dir == -1), sum(onco_sign == -1))),
    by = by][order(get(by[1]))]
}

report <- c("# Onco-eQTM consensus vs peripheral references", "",
            sprintf("Layer: %s pairs with a usable tumour sign (consensus direction, or sign of median r for O3/O4).",
                    format(nrow(O), big.mark = ",")),
            "Agreement = share of overlapping pairs where the tumour sign equals the reference sign; Wilson 95% CI.",
            "Baseline = majority-class share of the reference sign within the same overlap.",
            "Ladder rungs for comparison: catalogue_single 0.62-0.87, catalogue_consensus 0.77-0.87, smr_moderate 0.84-0.86, smr_high 0.95-0.97.",
            "")
summary_rows <- list()
for (nm in names(refs)) {
  R <- refs[[nm]]
  J <- merge(R, O, by = c("cpg_id", "target_gene"))
  cov_cpg <- mean(unique(R$cpg_id) %in% O$cpg_id)
  report <- c(report, sprintf("## Reference: %s", nm), "",
              sprintf("Reference pairs: %s; CpGs with any tumour record: %.1f%%; overlapping pairs: %s.",
                      format(nrow(R), big.mark = ","), 100 * cov_cpg, format(nrow(J), big.mark = ",")), "")
  if (!nrow(J)) { report <- c(report, "_no overlap_", ""); next }
  message(sprintf("[%s] overlap %d pairs", nm, nrow(J)))
  t1 <- eval_block(J, "onco_tier"); report <- c(report, "### By tier", "", md_table(t1), "")
  print(t1)
  t2 <- eval_block(J, c("onco_tier", "region")); report <- c(report, "### By tier and region", "", md_table(t2), "")
  t3 <- eval_block(J, c("onco_tier", "tss_bin")); report <- c(report, "### By tier and |distance to TSS|", "", md_table(t3), "")
  t4 <- ppv_block(J, "onco_tier"); report <- c(report, "### Value of each tumour sign, by tier", "", md_table(t4), "")
  # agreement-strength gradient inside O1/O2: more cancers, higher median |r|
  J[, r_bin := cut(abs(median_r), c(0, 0.3, 0.5, 1), labels = c("|r|<0.3", "0.3-0.5", ">=0.5"), include.lowest = TRUE)]
  t5 <- eval_block(J[onco_tier %in% c("O1", "O2")], c("onco_tier", "r_bin"))
  report <- c(report, "### O1/O2 by median |r|", "", md_table(t5), "")
  for (tier in c("O1", "O2")) {
    k <- J[onco_tier == tier]
    if (nrow(k)) summary_rows[[length(summary_rows) + 1L]] <-
      data.table(reference = nm, tier = tier, n = nrow(k),
                 agreement = sum(k$onco_sign == k$ref_dir) / nrow(k),
                 baseline = max(mean(k$ref_dir == 1), mean(k$ref_dir == -1)))
  }
}
if (length(summary_rows)) {
  Sm <- rbindlist(summary_rows)
  Sm[, verdict := fcase(
    agreement >= 0.84 & n >= 200, "at smr_moderate level",
    agreement >= 0.77 & n >= 200, "at catalogue_consensus level",
    agreement >= 0.62 & n >= 200 & agreement > baseline + 0.05, "at catalogue_single level",
    n < 200, "too few pairs",
    default = "below ladder / not above baseline")]
  Sm[, `:=`(agreement = sprintf("%.3f", agreement), baseline = sprintf("%.3f", baseline))]
  report <- c(report, "## Summary", "", md_table(Sm), "",
              "A rung is earned only if the agreement beats the majority baseline AND reaches an existing rung's",
              "validated range on at least two independent references with n >= 200. Promoter-proximal O1 pairs are",
              "the natural candidate; gene-body and distal pairs are where tumour and blood biology diverge.", "")
  print(Sm)
}
dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
writeLines(report, out_md)
message("wrote ", out_md)
