#!/usr/bin/env Rscript
# Validate the FinnGen immune-multiome cell-type SMR layer (GoDMC blood mQTL x
# FinnGen L1 cell-type cis-eQTL, HEIDI run against the dense mQTL BESD) against
# the package's independent references, to decide whether any context/tier
# earns a rung on the evidence ladder.
#
# Mirrors tools/validate_onco.R. The questions, in order:
#   1. Agreement with measured eQTMs (blood, nasal, solid; cohort-resolved) by
#      context x tier -- the only INDEPENDENT validation (different data).
#   2. The same by HEIDI stratum and by nsnp_HEIDI, so the record states for this
#      layer whether HEIDI discriminates (2.2.5 found it did not on GoDMC x eQTLGen).
#   3. Consistency with the existing SMR layer (GoDMC x eQTLGen): SAME exposure,
#      independent outcome study -- a pipeline check, not validation.
#   4. PBMC vs each cell type on shared pairs: how much the cell types add beyond
#      the pooled PBMC call (sign disagreement = cell-type-specific regulation or
#      noise; the measured check says which).
#
# Usage (from the cpgdirection package root)
#   Rscript tools/validate_celltype_smr.R \
#       --layer-glob "~/finngen_multiome/smr/smrd_blood_*.directions.csv.gz" \
#       --data-dir <hub_upload dir>   (optional; local cpgdirectionData files) \
#       --out tools/reports/validate_celltype_smr.md
#
# Prints the tables and writes them as markdown. Nothing in the package is changed.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
layer_glob <- get_arg("--layer-glob", "~/finngen_multiome/smr/smrd_blood_*.directions.csv.gz")
data_dir   <- get_arg("--data-dir")
out_md     <- get_arg("--out", "tools/reports/validate_celltype_smr.md")
min_n      <- as.integer(get_arg("--min-n", "200"))
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)

# package code: the package root under validation first, then an installed copy
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else if (requireNamespace("cpgdirection", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cpgdirection))
} else stop("run from the cpgdirection package root or install cpgdirection")

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
read_gz <- function(f) fread(cmd = sprintf("gzip -dc %s", shQuote(f)))

# ---- the candidate layer -------------------------------------------------------------
files <- Sys.glob(path.expand(layer_glob))
if (!length(files)) stop("no layer files match ", layer_glob)
L <- rbindlist(lapply(files, read_gz), fill = TRUE)
L[, target_gene := toupper(trimws(target_gene))]
L <- L[!is.na(direction) & direction != 0 & nzchar(target_gene)]
L[, heidi_status := fifelse(is.na(heidi_status) | heidi_status == "", "not_tested", heidi_status)]
L[, nsnp_bin := cut(nsnp_HEIDI, c(-Inf, 2, 5, 10, 20, Inf),
                    labels = c("<3", "3-5", "6-10", "11-20", ">20"))]
L[is.na(nsnp_bin), nsnp_bin := "<3"]
L[, inst_bin := cut(n_instruments, c(0, 1, 2, 4, Inf), labels = c("1", "2", "3-4", ">=5"))]
L[, tss_bin := cut(abs(cpg_gene_dist), c(-Inf, 1500, 10000, 100000, Inf),
                   labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb"))]
ctx_order <- c("blood_PBMC", "blood_CD4_T", "blood_CD8_T", "blood_Mono", "blood_NK",
               "blood_B", "blood_DC", "blood_other_T", "blood_other")
L[, context := factor(context, levels = c(ctx_order, setdiff(unique(context), ctx_order)))]
message(sprintf("layer: %s pairs, %s CpGs, %s genes, %d contexts",
                format(nrow(L), big.mark = ","), format(uniqueN(L$cpg_id), big.mark = ","),
                format(uniqueN(L$target_gene), big.mark = ","), uniqueN(L$context)))

# ---- references ----------------------------------------------------------------------
M <- as.data.table(cpgd_measured_eqtms())
if (!"cohorts" %in% names(M)) M[, cohorts := "packaged_v1"]
M <- M[, list(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
       by = c("cpg_id", "tissue", "direction", "tss_dist", "cohorts")]
M <- unique(M[nzchar(target_gene)])
refs <- list(
  blood = M[tissue == "blood",            list(cpg_id, target_gene, ref_dir = direction, cohorts)],
  nasal = M[tissue == "nasal_epithelium", list(cpg_id, target_gene, ref_dir = direction, cohorts)],
  solid = M[tissue == "solid_tissue",     list(cpg_id, target_gene, ref_dir = direction, cohorts)])
refs <- lapply(refs, function(r) unique(r[!is.na(ref_dir) & ref_dir != 0]))

S <- tryCatch(as.data.table(cpgd_smr_directions()), error = function(e) NULL)
if (!is.null(S)) S[, target_gene := toupper(trimws(target_gene))]

# ---- evaluation helpers --------------------------------------------------------------
eval_block <- function(J, by) {
  J[, list(n = .N,
           agree = fmt_ci(sum(direction == ref_dir), .N),
           baseline = sprintf("%.3f", max(mean(ref_dir == 1), mean(ref_dir == -1))),
           layer_pos = sprintf("%.2f", mean(direction == 1)),
           ref_pos = sprintf("%.2f", mean(ref_dir == 1))),
    by = by][do.call(order, lapply(by, function(b) get(b)))]
}
agree_rate <- function(J) if (nrow(J)) sum(J$direction == J$ref_dir) / nrow(J) else NA_real_

report <- c("# Cell-type SMR layer (GoDMC x FinnGen immune multiome L1) vs references", "",
            sprintf("Layer: %s pairs over %d contexts, files `%s`.",
                    format(nrow(L), big.mark = ","), uniqueN(L$context), layer_glob),
            "Agreement = share of overlapping pairs where the layer sign equals the reference sign; Wilson 95% CI.",
            "Baseline = majority-class share of the reference sign within the same overlap.",
            "Ladder rungs for comparison: catalogue_single 0.62-0.87, catalogue_consensus 0.77-0.87,",
            "smr_moderate 0.84-0.86 (GoDMC x eQTLGen S2), smr_high 0.95-0.97 (S1).",
            "HEIDI is REPORTED, never applied (2.2.5 precedent); the by-stratum tables below are the record for this layer.", "")

# layer composition
comp <- L[, list(pairs = .N, cpgs = uniqueN(cpg_id), genes = uniqueN(target_gene),
                 S1 = sum(smr_tier == "S1"), S2 = sum(smr_tier == "S2"), S3 = sum(smr_tier == "S3"),
                 heidi_pass = sum(heidi_status == "pass"), heidi_fail = sum(heidi_status == "fail"),
                 not_tested = sum(heidi_status == "not_tested"),
                 pos_share = sprintf("%.2f", mean(direction == 1))), by = context][order(context)]
report <- c(report, "## Layer composition", "", md_table(comp), "")
print(comp)

# ---- 1-2. independent references -----------------------------------------------------
summary_rows <- list()
for (nm in names(refs)) {
  R <- refs[[nm]]
  J <- merge(L, R, by = c("cpg_id", "target_gene"), allow.cartesian = TRUE)
  report <- c(report, sprintf("## Reference: measured %s eQTMs", nm), "",
              sprintf("Reference pairs: %s; overlapping layer rows: %s (%s distinct pairs).",
                      format(nrow(R), big.mark = ","), format(nrow(J), big.mark = ","),
                      format(uniqueN(J[, paste(cpg_id, target_gene)]), big.mark = ",")), "")
  if (!nrow(J)) { report <- c(report, "_no overlap_", ""); next }
  message(sprintf("[%s] overlap %d rows", nm, nrow(J)))
  t1 <- eval_block(J, c("context", "smr_tier"))
  report <- c(report, "### By context and tier", "", md_table(t1), ""); print(t1)
  t2 <- eval_block(J, c("context", "heidi_status"))
  report <- c(report, "### By context and HEIDI status", "", md_table(t2), "")
  t3 <- eval_block(J, c("smr_tier", "heidi_status"))
  report <- c(report, "### All contexts pooled: tier x HEIDI status", "", md_table(t3), ""); print(t3)
  t4 <- eval_block(J[heidi_status != "not_tested"], c("nsnp_bin"))
  report <- c(report, "### All contexts pooled: by nsnp_HEIDI (tested pairs)", "", md_table(t4), "")
  t5 <- eval_block(J, c("smr_tier", "inst_bin"))
  report <- c(report, "### All contexts pooled: tier x instruments", "", md_table(t5), "")
  t6 <- eval_block(J, c("context", "tss_bin"))
  report <- c(report, "### By context and |CpG - gene| distance", "", md_table(t6), "")
  if (uniqueN(J$cohorts) > 1L) {
    t7 <- eval_block(J[context == "blood_PBMC"], c("cohorts"))[n >= 50]
    report <- c(report, "### PBMC by measured cohort (n >= 50)", "", md_table(t7), "")
  }
  # HEIDI discrimination test, the 2.2.5 statistic: pass vs fail concordance, pooled
  pf <- J[heidi_status %in% c("pass", "fail")]
  if (nrow(pf) && uniqueN(pf$heidi_status) == 2L) {
    tab <- table(pf$heidi_status, pf$direction == pf$ref_dir)
    pv <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    report <- c(report, sprintf("HEIDI pass vs fail concordance (pooled): pass %s, fail %s, Fisher p = %.3g.",
                                fmt_ci(sum(pf$heidi_status == "pass" & pf$direction == pf$ref_dir), sum(pf$heidi_status == "pass")),
                                fmt_ci(sum(pf$heidi_status == "fail" & pf$direction == pf$ref_dir), sum(pf$heidi_status == "fail")),
                                pv), "")
  }
  for (cx in levels(droplevels(J$context))) for (tier in c("S1", "S2", "S3")) {
    k <- J[context == cx & smr_tier == tier]
    if (nrow(k)) summary_rows[[length(summary_rows) + 1L]] <-
      data.table(reference = nm, context = cx, tier = tier, n = nrow(k),
                 agreement = agree_rate(k),
                 baseline = max(mean(k$ref_dir == 1), mean(k$ref_dir == -1)))
  }
}

# ---- 3. consistency with the existing SMR layer (same exposure) ---------------------
if (!is.null(S) && "smr_tier" %in% names(S)) {
  R <- unique(S[!is.na(direction) & direction != 0,
                list(cpg_id, target_gene, ref_dir = direction, ref_tier = smr_tier, ref_heidi = heidi_status)])
  J <- merge(L, R, by = c("cpg_id", "target_gene"))
  report <- c(report, "## Consistency: existing SMR layer (GoDMC x eQTLGen) -- same exposure, independent outcome", "",
              sprintf("Overlapping rows: %s. This is a pipeline check, not independent validation: the instruments are shared.",
                      format(nrow(J), big.mark = ",")), "")
  if (nrow(J)) {
    t8 <- eval_block(J, c("context", "smr_tier"))
    report <- c(report, "### By context and layer tier", "", md_table(t8), ""); print(t8)
    t9 <- eval_block(J, c("smr_tier", "ref_tier"))
    report <- c(report, "### Layer tier x eQTLGen tier (pooled)", "", md_table(t9), "")
    cov_pairs <- L[, list(pairs = .N, in_eqtlgen_layer = sum(paste(cpg_id, target_gene) %in% paste(R$cpg_id, R$target_gene)),
                          new_pairs = sum(!paste(cpg_id, target_gene) %in% paste(R$cpg_id, R$target_gene))), by = context][order(context)]
    report <- c(report, "### Coverage: pairs already in the eQTLGen layer vs new", "", md_table(cov_pairs), ""); print(cov_pairs)
  }
}

# ---- 4. cell types vs pooled PBMC ---------------------------------------------------
P <- L[context == "blood_PBMC", list(cpg_id, target_gene, pbmc_dir = direction, pbmc_tier = smr_tier)]
if (nrow(P)) {
  J <- merge(L[context != "blood_PBMC"], P, by = c("cpg_id", "target_gene"))
  t10 <- J[, list(shared = .N, agree = fmt_ci(sum(direction == pbmc_dir), .N),
                  only_in_celltype = L[context == .BY$context, .N] - .N), by = context][order(context)]
  report <- c(report, "## Cell types vs pooled PBMC on shared pairs", "",
              "Disagreement here is either cell-type-specific regulation or noise; the measured-reference tables say which.", "",
              md_table(t10), "")
  print(t10)
  # where a cell type disagrees with PBMC and a measured blood sign exists, who is right?
  D <- merge(J[direction != pbmc_dir], refs$blood, by = c("cpg_id", "target_gene"))
  if (nrow(D)) {
    t11 <- D[, list(n = .N, celltype_right = fmt_ci(sum(direction == ref_dir), .N),
                    pbmc_right = fmt_ci(sum(pbmc_dir == ref_dir), .N)), by = context][order(context)]
    report <- c(report, "### Disagreements with PBMC that have a measured blood sign", "", md_table(t11), "")
    print(t11)
  }
}

# ---- summary and verdict -------------------------------------------------------------
if (length(summary_rows)) {
  Sm <- rbindlist(summary_rows)
  Sm[, verdict := fcase(
    agreement >= 0.95 & n >= min_n, "at smr_high level",
    agreement >= 0.84 & n >= min_n, "at smr_moderate level",
    agreement >= 0.77 & n >= min_n, "at catalogue_consensus level",
    agreement >= 0.62 & n >= min_n & agreement > baseline + 0.05, "at catalogue_single level",
    n < min_n, "too few pairs",
    default = "below ladder / not above baseline")]
  Sm[, `:=`(agreement = sprintf("%.3f", agreement), baseline = sprintf("%.3f", baseline))]
  setorder(Sm, reference, context, tier)
  report <- c(report, "## Summary", "", md_table(Sm), "",
              sprintf("A rung is earned only if the agreement beats the majority baseline AND reaches an existing rung's validated range on at least two independent references with n >= %d.", min_n),
              "Expected if the layer behaves like its parent: S1 ~0.96, S2 ~0.85, S3 ~0.70 against measured blood.",
              "The nasal and solid references test transfer out of blood; the brain bridge is not addressed by this layer.", "")
  print(Sm)
}
dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
writeLines(report, out_md)
message("wrote ", out_md)
