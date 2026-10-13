#!/usr/bin/env Rscript
# Build a cohort-resolved measured-eQTM table from the eQTM Atlas cohort
# downloads (Sriram ... Celedon 2026, https://shiny.crc.pitt.edu/eqtm_browser/,
# Downloads tab, one "<Cohort>_data.csv" per cohort; hg19 coordinates, EnsDb v75).
#
# What it does
#   1. reads every *_data.csv in --atlas-dir (columns: probeid, gene, pvalue, TSS,
#      dist, tissue, cohort, beta, pos, meth.chr, gene.chr, eQTM_label, num_cases,
#      num_controls, disease, statistic, FDR);
#   2. keeps cis rows with a usable sign (beta, else the test statistic) that pass
#      the cohort's own significance rule (FDR < 0.05 where an FDR column is
#      populated, else the p threshold the Atlas lists for that cohort);
#   3. maps the Atlas tissue to the package's tissue classes
#      (blood / nasal_epithelium / solid_tissue), keeping the raw tissue label;
#   4. reports, per cohort, how many of its signed cis pairs are ALREADY in the
#      packaged measured table and how often the signs agree -- this is how we
#      learn which cohorts the current build was made from;
#   5. writes a long table (one row per cohort x pair) and a collapsed table in the
#      packaged schema (cpg_id, target_gene, tissue, direction, tss_dist) plus
#      n_cohorts / cohort_agreement, with the sign taken by cohort majority and
#      ties dropped.
#
# Usage (from the cpgdirection package root)
#   Rscript inst/scripts/build_measured_eqtms.R --atlas-dir ~/Downloads/eqtm_atlas --out tools/out
#
# The merged table tools/out/measured_eqtms_v2.csv.gz is what ships as
# inst/extdata/measured_eqtms.csv.gz (copied by hand, see
# inst/scripts/make-extdata.md). Nothing in the package is changed by running
# this script; packaging is a separate decision.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
atlas_dir <- get_arg("--atlas-dir"); if (is.null(atlas_dir)) stop("--atlas-dir is required")
out_dir   <- get_arg("--out", "tools/out")
existing  <- get_arg("--existing", "inst/extdata/measured_eqtms.csv.gz")
# the packaged measured layer is a |CpG - TSS| <= 250 kb window (max 249,998 in
# inst/extdata/measured_eqtms.csv.gz); "cis" here means the same thing. The Atlas
# eQTM_label is not reliable for every cohort (FUSION rows at 41 kb are labelled
# "trans"), so cis is decided from the coordinates, not the label.
window    <- as.numeric(get_arg("--window", "250000"))
# What to do with packaged pairs the Atlas does not cover, per tissue class:
#   blood  -> kept as cohorts = "packaged_v1" (Framingham full set + HELIX, whose
#             Atlas exports are capped or unsigned)
#   solid  -> the 1,704 packaged solid pairs are TUMOUR-derived (the solid
#             catalogue is documented as tumour tissue, NEWS 1.8.x) and agree with
#             the normal-tissue cohorts on only 0.58 of shared pairs; tumour
#             evidence now lives in the Onco-eQTM layer, so by default they are
#             dropped from the measured rung. --packaged-solid keep retains them
#             under tissue "tumour"; --packaged-solid solid_tissue pools them.
packaged_solid <- get_arg("--packaged-solid", "drop")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# cohort p thresholds as listed on the Atlas Downloads page, used only where the
# file carries no FDR column worth using; FDR rules where the page lists a
# stricter FDR than 0.05
P_THRESH   <- c("Colonomics" = 1.7e-8, "EvoImmunoPop" = 4.99e-5, "HELIX" = 0.05,
                "Framingham" = 1e-7)
FDR_THRESH <- c("FUSION" = 0.01, "All Subjects" = 0.01)
rule_for <- function(tbl, cohort) {
  k <- names(tbl)[sapply(names(tbl), function(k) grepl(k, cohort, fixed = TRUE))]
  if (length(k)) tbl[[k[1]]] else NA_real_
}
tissue_class <- function(x) {
  x <- tolower(x)
  fifelse(grepl("blood|monocyte|wbc|leukocyte", x), "blood",
  fifelse(grepl("nasal", x), "nasal_epithelium", "solid_tissue"))
}

files <- list.files(atlas_dir, pattern = "_data( \\([0-9]+\\))?\\.csv$", full.names = TRUE)
if (!length(files)) stop("no *_data.csv in ", atlas_dir)
message(length(files), " cohort files")

read_cohort <- function(f) {
  d <- fread(f, select = c("probeid", "gene", "pvalue", "TSS", "dist", "tissue", "cohort",
                           "beta", "pos", "meth.chr", "gene.chr", "eQTM_label", "statistic",
                           "FDR", "num_cases", "num_controls"),
             colClasses = list(character = c("probeid", "gene", "tissue", "cohort",
                                             "meth.chr", "gene.chr", "eQTM_label")))
  n_all <- nrow(d); n_label_cis <- sum(d$eQTM_label == "cis", na.rm = TRUE)
  cohort <- d$cohort[1]
  d[, `:=`(mchr = sub("^chr", "", meth.chr), gchr = sub("^chr", "", gene.chr))]
  d[, dist_calc := suppressWarnings(as.numeric(TSS) - as.numeric(pos))]
  d[is.na(dist) | dist == 0, dist := dist_calc]
  d <- d[!is.na(dist) & abs(dist) <= window & mchr == gchr]
  n_cis <- nrow(d)
  d[, sgn := fifelse(!is.na(beta) & beta != 0, sign(beta),
             fifelse(!is.na(statistic) & statistic != 0, sign(statistic), NA_real_))]
  n_signed <- sum(!is.na(d$sgn))
  has_fdr <- any(!is.na(d$FDR))
  fth <- rule_for(FDR_THRESH, cohort); if (is.na(fth)) fth <- 0.05
  pth <- rule_for(P_THRESH, cohort);   if (is.na(pth)) pth <- 0.05
  d[, sig := if (has_fdr) !is.na(FDR) & FDR < fth else !is.na(pvalue) & pvalue < pth]
  d <- d[!is.na(sgn) & sig]
  d[, `:=`(cpg_id = probeid, target_gene = toupper(trimws(gene)),
           tissue_raw = tissue, tissue = tissue_class(tissue),
           chr = mchr, direction = sgn, tss_dist = dist,
           n = suppressWarnings(as.integer(num_cases)) + suppressWarnings(as.integer(num_controls)))]
  # one row per pair within a cohort: the most significant
  setorder(d, cpg_id, target_gene, pvalue)
  d <- unique(d, by = c("cpg_id", "target_gene"))
  list(meta = data.table(cohort = cohort, file = basename(f),
                         rows = n_all, label_cis = n_label_cis,
                         cis_window = n_cis, cis_signed = n_signed,
                         rule = if (has_fdr) sprintf("FDR<%g", fth) else sprintf("p<%g", pth),
                         kept = nrow(d), n_subjects = if (nrow(d)) d$n[1] else NA_integer_),
       data = d[, .(cohort, n_subjects = n, cpg_id, target_gene, tissue, tissue_raw, chr, pos,
                    tss = TSS, tss_dist, direction, beta, statistic, pvalue, FDR)])
}

res  <- lapply(files, read_cohort)
meta <- rbindlist(lapply(res, `[[`, "meta"))
L    <- rbindlist(lapply(res, `[[`, "data"))
message(sprintf("long table: %s cohort x pair rows", format(nrow(L), big.mark = ",")))

# ---- overlap with the packaged measured table -------------------------------------
E <- if (file.exists(existing)) fread(existing) else NULL
if (!is.null(E)) {
  E <- E[, .(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
         by = .(cpg_id, tissue, direction)]
  E <- unique(E[nzchar(target_gene)])
  J <- merge(L, E[, .(cpg_id, target_gene, tissue, e_dir = direction)],
             by = c("cpg_id", "target_gene", "tissue"), all.x = TRUE)
  ov <- J[, .(pairs = .N, in_package = sum(!is.na(e_dir)),
              share_in_package = round(mean(!is.na(e_dir)), 3),
              sign_agreement = round(mean(direction == e_dir, na.rm = TRUE), 3)),
          by = cohort][order(-share_in_package)]
  meta <- merge(meta, ov, by = "cohort", all.x = TRUE)
  # pairs in the package that no downloaded cohort explains
  Ex <- merge(E, unique(L[, .(cpg_id, target_gene, tissue, hit = TRUE)]),
              by = c("cpg_id", "target_gene", "tissue"), all.x = TRUE)
  unexplained <- Ex[is.na(hit), .N, by = tissue]
}

# ---- collapsed table in the packaged schema ----------------------------------------
C <- L[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) },
           n_cohorts = .N,
           cohort_agreement = round(max(mean(direction == 1), mean(direction == -1)), 3),
           tss_dist = tss_dist[which.min(pvalue)],
           cohorts = paste(sort(unique(cohort)), collapse = "|")),
       by = .(cpg_id, target_gene, tissue)]
ties <- C[is.na(direction)]
C <- C[!is.na(direction)]
if (!is.null(E)) {
  C <- merge(C, E[, .(cpg_id, target_gene, tissue, packaged_dir = direction)],
             by = c("cpg_id", "target_gene", "tissue"), all.x = TRUE)
  C[, status := fifelse(is.na(packaged_dir), "new",
                fifelse(packaged_dir == direction, "in_package_agree", "in_package_conflict"))]
  gain <- dcast(C[, .N, by = .(tissue, status)], tissue ~ status, value.var = "N", fill = 0L)
}

# ---- sign conventions: agreement between cohorts on the pairs they share -------------
# Two cohorts measuring the same CpG-gene pair should mostly agree; a cohort whose
# sign is flipped relative to every other shows up here as ~1 - the usual agreement.
P <- merge(L[, .(cohort, cpg_id, target_gene, direction)],
           L[, .(cohort2 = cohort, cpg_id, target_gene, direction2 = direction)],
           by = c("cpg_id", "target_gene"), allow.cartesian = TRUE)[cohort < cohort2]
xc <- P[, .(shared = .N, agreement = round(mean(direction == direction2), 3)),
        by = .(cohort, cohort2)][shared >= 20][order(-shared)]

# ---- merged table: Atlas pairs + packaged pairs the Atlas does not cover ------------
# Packaged pairs whose cohort is absent from the folder (or unsigned in the Atlas
# export) are kept with cohorts = "packaged_v1" so nothing measured is lost.
M <- copy(C)[, .(cpg_id, target_gene, tissue, direction, tss_dist, n_cohorts, cohort_agreement,
                 cohorts, source = "eqtm_atlas")]
if (!is.null(E)) {
  keep <- Ex[is.na(hit), .(cpg_id, target_gene, tissue, direction)]
  keep <- merge(keep, unique(fread(existing)[, .(cpg_id, tissue, tss_dist)]),
                by = c("cpg_id", "tissue"), all.x = TRUE)
  keep <- unique(keep, by = c("cpg_id", "target_gene", "tissue"))
  keep[, `:=`(n_cohorts = 1L, cohort_agreement = 1, cohorts = "packaged_v1", source = "packaged_v1")]
  if (packaged_solid == "drop") {
    keep <- keep[tissue != "solid_tissue"]
  } else if (packaged_solid == "keep") {
    keep[tissue == "solid_tissue", tissue := "tumour"]
  }
  M <- rbind(M, keep[, names(M), with = FALSE])
}
setorder(M, tissue, cpg_id, target_gene)

fwrite(L, file.path(out_dir, "measured_eqtms_atlas_long.csv.gz"))
fwrite(C, file.path(out_dir, "measured_eqtms_atlas.csv.gz"))
fwrite(M, file.path(out_dir, "measured_eqtms_v2.csv.gz"))

md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  h <- paste0("| ", paste(names(d), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  r <- apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(h, s, r), collapse = "\n")
}
rep <- c("# eQTM Atlas cohorts vs the packaged measured layer", "",
         sprintf("Atlas dir: `%s`; %d cohort files; cis window |CpG - TSS| <= %s bp, same chromosome; packaged table: `%s` (%s pairs).",
                 atlas_dir, length(files), format(window, big.mark = ","), existing,
                 if (!is.null(E)) format(nrow(E), big.mark = ",") else "absent"), "",
         "## Per cohort", "",
         "`label_cis` = rows the Atlas labels cis; `cis_window` = rows inside the window by coordinates; `cis_signed` = of those, rows with a sign (beta, else statistic); `kept` = signed, significant, one row per pair.",
         "`share_in_package` close to 1 means the current build already used this cohort.", "",
         md_table(meta), "")
if (!is.null(E)) rep <- c(rep, "## Packaged pairs not covered by any downloaded cohort", "",
                          md_table(unexplained), "",
                          "A large count here means a cohort that is not in this folder (e.g. Framingham, EVA-PR WBC) was used.", "")
rep <- c(rep, "## Collapsed table", "",
         sprintf("%s pairs across %s CpGs; %s dropped as sign ties between cohorts.",
                 format(nrow(C), big.mark = ","), format(uniqueN(C$cpg_id), big.mark = ","),
                 format(nrow(ties), big.mark = ",")), "",
         md_table(C[, .(pairs = .N, cpgs = uniqueN(cpg_id), neg_share = round(mean(direction == -1), 3),
                        multi_cohort = sum(n_cohorts > 1)), by = tissue]), "")
if (!is.null(E)) rep <- c(rep, "## What the Atlas adds to the packaged layer (collapsed pairs by status)", "",
                          md_table(gain), "")
rep <- c(rep, "## Sign agreement between cohorts on shared pairs (>= 20 shared)", "",
         "Agreement well below the others for one cohort against everyone = a flipped sign convention; near 0.5 against a different tissue = tissue-specific biology, not an error.", "",
         md_table(xc), "",
         "## Merged table (Atlas + packaged pairs the Atlas does not cover)", "",
         md_table(M[, .(pairs = .N, cpgs = uniqueN(cpg_id), from_atlas = sum(source == "eqtm_atlas"),
                        from_packaged_v1 = sum(source == "packaged_v1"),
                        neg_share = round(mean(direction == -1), 3)), by = tissue]), "")
rep <- c(rep,
         "## By cohort and distance to TSS", "",
         md_table(L[, .(n = .N, neg_share = round(mean(direction == -1), 3)),
                    by = .(cohort, bin = cut(abs(tss_dist), c(-Inf, 1500, 10000, 100000, Inf),
                                             labels = c("<=1.5kb", "1.5-10kb", "10-100kb", ">100kb")))][order(cohort, bin)]), "")
writeLines(rep, file.path(out_dir, "measured_eqtms_atlas.md"))
cat(paste(rep, collapse = "\n"), "\n")
message("wrote ", out_dir, "/measured_eqtms_atlas{.md,.csv.gz,_long.csv.gz} and measured_eqtms_v2.csv.gz")
