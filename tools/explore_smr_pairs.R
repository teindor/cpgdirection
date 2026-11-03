#!/usr/bin/env Rscript
# Lead 1 -- SMR-tested pairs as a second pair source (exploration, nothing registered).
#
# The FinnGen runs wrote, per context and chromosome, SMR's own result table
# (smr/bychr/smrd_<ctx>_chr<c>.smr): every exposure probe (CpG) with a cis-mQTL at
# p < 5e-8 tested against every outcome probe (gene) in cis, with the signed SMR
# effect b_SMR = b_outcome / b_exposure at the top mQTL SNP, p_SMR, p_HEIDI and
# nsnp_HEIDI. The directions tables used so far come from the clumped-instrument
# Wald ratios and only ATTACH these results; pairs with no clumped instrument were
# never considered. This script asks what those pairs are worth.
#
# Questions
#   1. How many SMR-tested pairs are there per context, how many pass p_SMR < 5e-8
#      (and 1e-5, 1e-3), and how many of those are NEW: absent from the instrument
#      table of the same context and absent from the eQTLGen layer?
#   2. On pairs covered by both, does sign(b_SMR) agree with the instrument-based
#      direction (pipeline sanity)?
#   3. Against the measured blood eQTMs, how accurate is sign(b_SMR) by p_SMR
#      stratum x HEIDI status x novelty, and does any new stratum reach a rung
#      (smr_moderate 0.84, smr_high 0.95) on n >= 200?
#
# Usage (from the cpgdirection package root)
#   Rscript tools/explore_smr_pairs.R --smr-dir ~/finngen_multiome/smr/bychr \
#       [--contexts PBMC,CD4_T,...] [--out tools/reports/explore_smr_pairs.md]
#       [--data-dir <hub_upload dir>]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
smr_dir  <- path.expand(get_arg("--smr-dir", "~/finngen_multiome/smr/bychr"))
ctx_arg  <- get_arg("--contexts", "PBMC,CD4_T,CD8_T,Mono,NK,B,DC,other_T,other")
out_md   <- get_arg("--out", "tools/reports/explore_smr_pairs.md")
data_dir <- get_arg("--data-dir")
p_thr    <- 5e-8
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)
contexts <- strsplit(ctx_arg, ",", fixed = TRUE)[[1]]

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
fmt_ci <- function(k, n) { if (n == 0) return("-"); ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2]) }
md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  paste(c(paste0("| ", paste(names(d), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
          apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))), collapse = "\n")
}
read_gz <- function(f) fread(cmd = sprintf("gzip -dc %s", shQuote(f)))
gkey <- function(g) gsub("_", "-", toupper(trimws(g)))
big <- function(x) format(x, big.mark = ",")

# ---- references -----------------------------------------------------------------------
E <- as.data.table(cpgd_smr_directions())
E <- E[if ("outcome" %in% names(E)) outcome == "eQTLGen" else TRUE]
E[, .gkey := gkey(target_gene)]
ekey <- unique(paste(E$cpg_id, E$.gkey))
M <- as.data.table(cpgd_measured_eqtms())[tissue == "blood"]
M <- unique(M[, list(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
              by = .(cpg_id, ref_dir = direction)])
M[, .gkey := gkey(target_gene)]
M <- unique(M[, .(cpg_id, .gkey, ref_dir)])
message(sprintf("eQTLGen layer %s pairs; measured blood %s pairs", big(length(ekey)), big(nrow(M))))

report <- c("# SMR-tested pairs as a second pair source -- exploration", "",
            sprintf("Source: `%s/smrd_<ctx>_chr<c>.smr` (SMR v1.03 two-BESD mode; exposure = GoDMC dense cis mQTL, outcome = FinnGen L1 cis-eQTL; HEIDI vs 1000G EUR).", smr_dir),
            "Direction = sign(b_SMR) at the top cis-mQTL SNP. HEIDI pass = p >= 0.05 (table contract). 'instrument' = pair present in the context's clumped-instrument directions table; 'eQTLGen' = pair present in the shipped GoDMC x eQTLGen layer.",
            "Rungs for comparison: smr_moderate 0.84-0.86, smr_high 0.95-0.97, replicated S2 0.958.", "")
summary_rows <- list(); strata_rows <- list()

for (ct in contexts) {
  files <- Sys.glob(file.path(smr_dir, sprintf("smrd_%s_chr*.smr", ct)))
  if (!length(files)) { message("no .smr files for ", ct); next }
  S <- rbindlist(lapply(files, function(f) {
    d <- fread(f, showProgress = FALSE)
    need <- c("Expo_ID", "Expo_bp", "Outco_ID", "Outco_Gene", "Outco_bp", "topSNP", "p_Expo",
              "b_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI")
    miss <- setdiff(need, names(d)); if (length(miss)) stop(f, ": missing ", paste(miss, collapse = ","))
    d[, need, with = FALSE]
  }), use.names = TRUE)
  S <- S[!is.na(b_SMR) & b_SMR != 0 & !is.na(p_SMR)]
  S[, `:=`(cpg_id = Expo_ID, target_gene = toupper(trimws(Outco_Gene)), ensembl_id = Outco_ID,
           direction = sign(b_SMR), cpg_gene_dist = Outco_bp - Expo_bp)]
  S <- S[nzchar(target_gene)]
  S[, .gkey := gkey(target_gene)]
  S[, heidi_status := fifelse(is.na(p_HEIDI), "not_tested", fifelse(p_HEIDI >= 0.05, "pass", "fail"))]
  # one row per pair (a gene can appear under several outcome probes): smallest p_SMR
  setorder(S, cpg_id, .gkey, p_SMR)
  S <- S[!duplicated(S[, .(cpg_id, .gkey)])]
  S[, p_bin := cut(p_SMR, c(0, 5e-8, 1e-5, 1e-3, 1), labels = c("<5e-8", "5e-8..1e-5", "1e-5..1e-3", ">1e-3"), right = FALSE)]

  # the instrument-based table of the same context
  inst_file <- file.path(dirname(smr_dir), sprintf("smrd_blood_%s.directions.csv.gz", ct))
  I <- if (file.exists(inst_file)) read_gz(inst_file) else NULL
  if (!is.null(I)) {
    I[, .gkey := gkey(target_gene)]
    I <- unique(I[, .(cpg_id, .gkey, inst_dir = as.numeric(direction), inst_tier = smr_tier)])
    S <- merge(S, I, by = c("cpg_id", ".gkey"), all.x = TRUE)
  } else S[, `:=`(inst_dir = NA_real_, inst_tier = NA_character_)]
  S[, in_instrument := !is.na(inst_dir)]
  S[, in_eqtlgen := paste(cpg_id, .gkey) %in% ekey]
  S[, novelty := fcase(in_instrument & in_eqtlgen, "both",
                       in_instrument, "instrument_only",
                       in_eqtlgen, "eqtlgen_only",
                       default = "NEW")]

  # 1. counts
  cnt <- S[, .(tested = .N,
               p5e8 = sum(p_SMR < 5e-8), p1e5 = sum(p_SMR < 1e-5),
               new_p5e8 = sum(p_SMR < 5e-8 & novelty == "NEW"),
               new_p5e8_heidi_pass = sum(p_SMR < 5e-8 & novelty == "NEW" & heidi_status == "pass"),
               new_p1e5 = sum(p_SMR < 1e-5 & novelty == "NEW"),
               cpgs_new_p5e8 = uniqueN(cpg_id[p_SMR < 5e-8 & novelty == "NEW"]))]
  cnt[, context := ct]
  # 2. sanity vs the instrument-based direction
  san <- S[in_instrument & p_SMR < 5e-8, .(n = .N, agree = fmt_ci(sum(direction == inst_dir), .N))]
  san_tier <- S[in_instrument & p_SMR < 5e-8, .(n = .N, agree = fmt_ci(sum(direction == inst_dir), .N)), by = inst_tier][order(inst_tier)]
  # 3. measured blood
  V <- merge(S, M, by = c("cpg_id", ".gkey"))
  v1 <- V[, .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N),
              baseline = sprintf("%.3f", max(mean(ref_dir == 1), mean(ref_dir == -1)))),
          by = .(novelty, p_bin)][order(novelty, p_bin)]
  v2 <- V[p_SMR < 5e-8, .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N)),
          by = .(novelty, heidi_status)][order(novelty, heidi_status)]
  v3 <- V[novelty == "NEW", .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N)),
          by = .(p_bin, heidi_status)][order(p_bin, heidi_status)]
  # nsnp / instrument strength inside the NEW p<5e-8 stratum
  v4 <- V[novelty == "NEW" & p_SMR < 5e-8, .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N)),
          by = .(expo_p = cut(p_Expo, c(0, 1e-30, 1e-15, 1e-10, 5e-8), labels = c("<1e-30", "1e-30..1e-15", "1e-15..1e-10", "1e-10..5e-8"), right = FALSE))][order(expo_p)]

  message(sprintf("\n== %s: %s tested pairs, %s at p<5e-8, %s NEW at p<5e-8 (%s HEIDI pass), %s NEW CpGs",
                  ct, big(cnt$tested), big(cnt$p5e8), big(cnt$new_p5e8), big(cnt$new_p5e8_heidi_pass), big(cnt$cpgs_new_p5e8)))
  message("   sign(b_SMR) vs instrument direction, p<5e-8: ", san$agree, " (n=", san$n, ")")
  print(v1); print(v2)

  report <- c(report, sprintf("## %s", ct), "",
              sprintf("Tested pairs %s; p_SMR < 5e-8: %s; NEW at p < 5e-8: %s (HEIDI pass %s), over %s CpGs.",
                      big(cnt$tested), big(cnt$p5e8), big(cnt$new_p5e8), big(cnt$new_p5e8_heidi_pass), big(cnt$cpgs_new_p5e8)),
              sprintf("Sanity: sign(b_SMR) agrees with the instrument-based direction on %s of %s shared p<5e-8 pairs.", san$agree, big(san$n)), "",
              "### sanity by instrument tier", "", md_table(san_tier), "",
              "### measured blood: by novelty x p_SMR", "", md_table(v1), "",
              "### measured blood: p<5e-8 by novelty x HEIDI", "", md_table(v2), "",
              "### measured blood: NEW pairs by p_SMR x HEIDI", "", md_table(v3), "",
              "### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p", "", md_table(v4), "")
  summary_rows[[ct]] <- cnt
  k <- V[novelty == "NEW" & p_SMR < 5e-8]
  for (hs in c("pass", "fail", "not_tested", "all")) {
    kk <- if (hs == "all") k else k[heidi_status == hs]
    if (nrow(kk)) strata_rows[[length(strata_rows) + 1L]] <-
      data.table(context = ct, stratum = paste0("NEW p<5e-8, HEIDI ", hs), n = nrow(kk),
                 agreement = sum(kk$direction == kk$ref_dir) / nrow(kk),
                 baseline = max(mean(kk$ref_dir == 1), mean(kk$ref_dir == -1)))
  }
  rm(S, V); invisible(gc())
}

if (length(summary_rows)) {
  Sm <- rbindlist(summary_rows); setcolorder(Sm, "context")
  report <- c(report, "## Coverage summary", "", md_table(Sm), "")
  print(Sm)
}
if (length(strata_rows)) {
  St <- rbindlist(strata_rows)
  St[, verdict := fcase(agreement >= 0.95 & n >= 200, "at smr_high level",
                        agreement >= 0.84 & n >= 200, "at smr_moderate level",
                        agreement >= 0.77 & n >= 200, "at catalogue_consensus level",
                        n < 200, "too few pairs",
                        default = "below ladder")]
  St[, `:=`(agreement = sprintf("%.3f", agreement), baseline = sprintf("%.3f", baseline))]
  report <- c(report, "## Verdict: NEW p<5e-8 pairs vs measured blood", "", md_table(St), "",
              "A NEW stratum earns the moderate rung only if its lower Wilson bound clears 0.84 on the measured blood reference;",
              "pairs also found by the instrument table or the eQTLGen layer are not new coverage and are scored only as sanity.", "")
  print(St)
}
dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
writeLines(report, out_md)
message("wrote ", out_md)
