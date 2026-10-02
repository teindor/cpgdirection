#!/usr/bin/env Rscript
# Build the `onco_eqtm_consensus` layer from Onco-eQTM per-cancer result tables.
#
# Source: Korra BT, Nishana M, Kumar R (2026) "Expression quantitative trait
# methylation across multiple cancer types with functional and therapeutic
# characterization using Onco-eQTM", NAR Genomics & Bioinformatics 8:lqag101,
# doi:10.1093/nargab/lqag101. Processed objects: doi:10.5281/zenodo.20813261;
# pipeline: https://github.com/CGnTLab/onco_eqtm; web DB:
# https://project.iith.ac.in/cgntlab/OncoeQTM/. Licence CC BY-NC 4.0.
#
# What the input looks like. Two layouts are accepted, detected per file:
#
#   (a) the authors' pipeline output `final_<CANCER>.csv` (scripts/analysis/
#       05_run_eqtm_pipeline.sh): consensus of Torch-eCpG and MatrixEQTL
#       (FDR < 0.05 in both, |r| > 0.3), with a `correlation` column;
#   (b) the web Data Download Center export
#       `OncoeQTM_<CANCER>_cis_results.csv` (POST cancerType=<CANCER>&module=cis
#       to download_module_data.php): the same inner join of both methods, but
#       WITHOUT the correlation column and without the |r| > 0.3 cut. For (b)
#       the correlation is recovered from the MatrixEQTL t statistic,
#       r = t / sqrt(df + t^2) with df = n_samples - 2 - n_covariates, using
#       the per-cancer sample sizes reported in the paper (Fig. 1A; 6,880
#       tumours) and the two covariates (age, sex) of the published model.
#       The |r| > --min-r cut (default 0.3, as in the paper) is then applied
#       here so both layouts yield the same consensus definition.
#
# Common columns: mt_id gt_id mt_chrom mt_chromStart mt_strand gt_chrom
# gt_chromStart gt_strand region mt_est mt_err mt_t mt_p beta t-stat p-value
# FDR [correlation]. mt_id is the 450K probe (cg...), gt_id the gene symbol
# (Xena HiSeqV2), region one of CIS / PROMOTER / DISTAL. Coordinates are
# TCGA/Xena, i.e. hg19 -- the same build as every other layer in cpgdirection,
# so no liftover is needed. 450K probes are a subset of EPIC v2 CpG identifiers
# (replicate suffixes are already collapsed to the cg ID in this package).
#
# The web site also lists composites (COADREAD, GBMLGG, LUNG) and six cancer
# types absent from the paper (CESC, LAML, OV, PRAD, UCEC, UCS). Composites are
# always skipped (they would double-count their members); the six extras are
# skipped unless --all-cancers is given, because no published sample size is
# available to recover r for them.
#
# What this script does NOT do by default: it does not append anything to
# `measured_eqtms`. In the pair ladder (R/direction_full.R) a measured record
# is the top rung with expected accuracy "~1.00 (measured, not predicted)" and
# is tissue-blind, so 5.25 million tumour pairs would silently outrank blood
# SMR S1 for a saliva study. Tumour evidence belongs in its own column set,
# reported alongside (`cpgd_onco_eqtm()`, `include_onco = TRUE` in
# cpg_gene_pairs), exactly as brain SMR is. `--append-solid` exists for the
# explicit decision to promote high-agreement pairs into the solid_tissue
# catalogue; it writes a SEPARATE file and labels every row it adds.
#
# Usage
#   Rscript tools/build_onco_eqtm.R --dir ~/onco_eqtm \
#       --out inst/extdata/onco_eqtm_consensus.csv.gz \
#       [--min-cancers 2] [--min-agreement 0.9] [--min-r 0.3] [--all-cancers] \
#       [--append-solid inst/extdata/measured_eqtms.csv.gz --append-out measured_eqtms_plus_onco.csv.gz]
#
# `--dir` is searched recursively for final_*.csv (the Zenodo archive keeps the
# per-cancer Results/ folders).

suppressPackageStartupMessages(library(data.table))

# write .rds when the target name says so (the Hub format), else gzipped csv
write_out <- function(d, path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) saveRDS(as.data.frame(d), path, compress = "xz")
  else fwrite(d, path, compress = "gzip")
}
read_any <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.table(readRDS(path)))
  fread(path, showProgress = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
in_dir        <- get_arg("--dir")
out_file      <- get_arg("--out", "onco_eqtm_consensus.csv.gz")
min_cancers   <- as.integer(get_arg("--min-cancers", "2"))
min_agreement <- as.numeric(get_arg("--min-agreement", "0.9"))
append_solid  <- get_arg("--append-solid")
append_out    <- get_arg("--append-out", "measured_eqtms_plus_onco.csv.gz")
min_r         <- as.numeric(get_arg("--min-r", "0.3"))
all_cancers   <- "--all-cancers" %in% args

# per-cancer primary-tumour sample sizes with complete multi-omics profiles
# (Korra et al. 2026, Fig. 1A; total 6,880) and the covariate count of the
# published model (age, sex). df = n - 2 - n_cov.
N_SAMPLES <- c(ACC = 79L, BLCA = 407L, BRCA = 783L, CHOL = 36L, COAD = 278L,
               DLBC = 48L, ESCA = 184L, GBM = 51L, HNSC = 520L, KICH = 66L,
               KIRC = 318L, KIRP = 274L, LGG = 516L, LIHC = 371L, LUAD = 453L,
               LUSC = 371L, MESO = 87L, PAAD = 178L, PCPG = 179L, READ = 92L,
               SARC = 259L, SKCM = 103L, STAD = 372L, TGCT = 150L, THCA = 505L,
               THYM = 120L, UVM = 80L)
N_COV <- 2L
COMPOSITES <- c("COADREAD", "GBMLGG", "LUNG")

if (is.null(in_dir) || !dir.exists(in_dir))
  stop("--dir must be a directory holding final_<CANCER>.csv tables.", call. = FALSE)
files <- list.files(in_dir,
                    pattern = "^(final_[A-Za-z0-9]+|OncoeQTM_[A-Za-z0-9]+_cis_results)\\.csv(\\.gz)?$",
                    recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("no final_<CANCER>.csv or OncoeQTM_<CANCER>_cis_results.csv files under ", in_dir)
cancer_of <- function(f) sub("^(final_|OncoeQTM_)([A-Za-z0-9]+?)(_cis_results)?\\.csv(\\.gz)?$", "\\2", basename(f))
cancers <- toupper(cancer_of(files))
drop <- cancers %in% COMPOSITES | (!all_cancers & !cancers %in% names(N_SAMPLES))
if (any(drop)) message("skipping: ", paste(unique(cancers[drop]), collapse = ", "),
                       " (composites, or no published sample size; use --all-cancers to keep the latter)")
files <- files[!drop]
if (!length(files)) stop("no usable cancer-type files left after skipping composites/extras")

read_one <- function(f) {
  d <- fread(f, showProgress = FALSE)
  cancer <- toupper(cancer_of(f))
  need <- c("mt_id", "gt_id")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(basename(f), " lacks columns: ", paste(miss, collapse = ", "))
  if ("correlation" %in% names(d)) {
    rr <- as.numeric(d$correlation)
    r_source <- "pipeline"
  } else {
    tcol <- if ("t-stat" %in% names(d)) "t-stat" else if ("mt_t" %in% names(d)) "mt_t" else NA
    if (is.na(tcol)) stop(basename(f), ": neither correlation nor a t statistic column")
    n <- N_SAMPLES[cancer]
    if (is.na(n)) {
      # extra cancer without a published n: keep the sign only, at a nominal
      # magnitude just above the cut so the row survives the |r| filter and
      # contributes to the sign vote but not to the r summaries (set to NA below)
      rr <- sign(as.numeric(d[[tcol]])) * (min_r + 1e-6)
      r_source <- "sign_only"
      message(sprintf("  %s: no sample size on record; direction kept, |r| unavailable", cancer))
    } else {
      tt <- as.numeric(d[[tcol]])
      df <- n - 2L - N_COV
      rr <- tt / sqrt(df + tt^2)
      r_source <- sprintf("t_to_r(n=%d)", n)
    }
  }
  m <- regexpr("cg[0-9]{6,}", d$mt_id, ignore.case = TRUE)
  d <- d[, list(
    cpg_id      = ifelse(m > 0, tolower(substr(mt_id, m, m + attr(m, "match.length") - 1L)), NA_character_),
    target_gene = toupper(trimws(as.character(gt_id))),
    cancer      = cancer,
    r           = rr,
    beta        = if ("beta" %in% names(d)) as.numeric(beta) else NA_real_,
    fdr         = if ("FDR" %in% names(d)) as.numeric(FDR) else NA_real_,
    region      = if ("region" %in% names(d)) toupper(as.character(region)) else NA_character_,
    cpg_pos     = if ("mt_chromStart" %in% names(d)) as.integer(mt_chromStart) else NA_integer_,
    gene_start  = if ("gt_chromStart" %in% names(d)) as.integer(gt_chromStart) else NA_integer_,
    gene_strand = if ("gt_strand" %in% names(d)) as.character(gt_strand) else NA_character_)]
  d <- d[!is.na(cpg_id) & nzchar(target_gene) & !grepl("^ENSG[0-9]+", target_gene) & !is.na(r)]
  n_before <- nrow(d)
  d <- d[abs(r) > min_r]
  message(sprintf("  %-6s %7d rows -> %7d with |r| > %.2f  (r from %s)",
                  cancer, n_before, nrow(d), min_r, r_source))
  unique(d, by = c("cpg_id", "target_gene"))
}
all <- rbindlist(lapply(files, read_one))
message(sprintf("read %d cancer types, %d CpG-gene-cancer rows, %d unique pairs",
                uniqueN(all$cancer), nrow(all), uniqueN(all[, paste(cpg_id, target_gene)])))

# signed TSS distance in the package's convention (negative = upstream of TSS
# on the gene's strand); gt_chromStart is the gene start, not strictly the TSS
# for minus-strand genes, so this is approximate and labelled as such.
all[, tss_dist_approx := ifelse(gene_strand %in% "-", gene_start - cpg_pos, cpg_pos - gene_start)]

cons <- all[, list(
  n_cancers      = .N,
  n_positive     = sum(r > 0),
  n_negative     = sum(r < 0),
  sign_agreement = max(mean(r > 0), mean(r < 0)),
  median_r       = median(r),
  max_abs_r      = max(abs(r)),
  region         = names(sort(table(region), decreasing = TRUE))[1],
  tss_dist_approx = as.integer(median(tss_dist_approx, na.rm = TRUE)),
  cancers        = paste(sort(unique(cancer)), collapse = ";")),
  by = c("cpg_id", "target_gene")]
cons[, direction := fifelse(sign_agreement >= min_agreement & n_cancers >= min_cancers,
                            fifelse(n_positive >= n_negative, 1, -1), NA_real_)]
cons[, onco_tier := fcase(
  n_cancers >= 5L & sign_agreement >= 0.9, "O1",   # many tumours, one sign
  n_cancers >= 2L & sign_agreement >= 0.9, "O2",   # few tumours, one sign
  n_cancers >= 2L,                         "O3",   # tumours disagree
  default = "O4")]                                  # single cancer type
cons[, source := "Onco-eQTM (Korra 2026); TCGA 450K hg19; CC BY-NC 4.0"]
setcolorder(cons, c("cpg_id", "target_gene", "direction", "onco_tier", "n_cancers",
                    "n_positive", "n_negative", "sign_agreement", "median_r", "max_abs_r",
                    "region", "tss_dist_approx", "cancers", "source"))
setorderv(cons, c("cpg_id", "target_gene"))
write_out(cons, out_file)

tt <- function(k) sum(cons$onco_tier == k)
message(sprintf("wrote %s: %d pairs, %d CpGs; tiers O1=%d O2=%d O3=%d O4=%d; %d pairs with a consensus sign",
                out_file, nrow(cons), uniqueN(cons$cpg_id),
                tt("O1"), tt("O2"), tt("O3"), tt("O4"), sum(!is.na(cons$direction))))
message(sprintf("sign split among consensus pairs: %.1f%% negative (higher methylation -> lower expression)",
                100 * mean(cons[!is.na(direction)]$direction < 0)))

# ---- optional, explicit: promote O1 pairs into the solid_tissue catalogue ----
if (!is.null(append_solid)) {
  if (!file.exists(append_solid)) stop("--append-solid file not found: ", append_solid)
  me <- read_any(append_solid)
  add <- cons[onco_tier == "O1" & !is.na(direction),
              list(cpg_id, target_gene, tissue = "solid_tissue", direction,
                   tss_dist = tss_dist_approx, source = "Onco-eQTM_O1")]
  if (!"source" %in% names(me)) me[, source := "packaged"]
  key_me <- me[, paste(cpg_id, toupper(target_gene), tissue)]
  add <- add[!paste(cpg_id, target_gene, tissue) %in% key_me]
  out <- rbindlist(list(me, add), use.names = TRUE, fill = TRUE)
  write_out(out, append_out)
  message(sprintf("appended %d O1 tumour pairs to solid_tissue -> %s (%d rows). ",
                  nrow(add), append_out, nrow(out)),
          "These rows will take the 'measured' rung with its ~1.00 accuracy label; ",
          "re-validate the solid-tissue accuracy strings in R/direction_full.R before shipping.")
}
