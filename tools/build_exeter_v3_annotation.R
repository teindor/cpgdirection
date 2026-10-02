#!/usr/bin/env Rscript
# Rebuild the packaged EPIC v2 probe-to-gene annotation union with the Exeter
# re-annotated manifest v3.0 (GENCODE release 49) in place of v2.0 (GENCODE 47).
#
# Why a patch and not a from-scratch rebuild: the shipped table is the UNION of
# five tracks and four of them (Illumina RefGene, Illumina GencodeV41, Zhou
# GENCODE v41, probe coordinates) are unchanged by this update. Replacing only
# the two Exeter-derived tracks keeps every other pair byte-identical, which is
# what makes the before/after diff below meaningful: anything that moved, moved
# because GENCODE 47 -> 49 moved it.
#
# Inputs
#   --exeter    Exeter re-annotated EPICv2 manifest v3.0 (CSV, may be .gz).
#               Zenodo record 20704849 (concept DOI 10.5281/zenodo.14933468),
#               "EPICv2 Re-annotated Manifest version 3.0 (Updated 15th June
#               2026)". Download by hand: the record is behind Zenodo's rate
#               limiter for automated fetches.
#   --current   The currently packaged epicv2_probe_gene_annotation table:
#               csv.gz (fat build) or the Hub .rds
#               (cpgdirectionData/hub_upload/epicv2_probe_gene_annotation.rds).
#               If omitted, it is pulled through the package's own
#               data backend (fat build extdata, or cpgdirectionData).
#   --out       Output path; a .rds name writes the Hub format, anything else
#               writes gzipped csv (default epicv2_probe_gene_annotation.csv.gz).
#   --report    Markdown diff report path (default exeter_v3_diff.md).
#
# Usage
#   Rscript tools/build_exeter_v3_annotation.R \
#     --exeter  ~/Downloads/EPICv2_reannotated_manifest_v3.0.csv \
#     --current inst/extdata/epicv2_probe_gene_annotation.csv.gz \
#     --out     inst/extdata/epicv2_probe_gene_annotation.csv.gz \
#     --report  tools/reports/exeter_v3_diff.md
#
# What the table looks like afterwards: identical column set, except that
# src_Exeter_GENCODEv47 becomes src_Exeter_GENCODEv49 and the annotation_source
# token "Exeter_GENCODEv47" becomes "Exeter_GENCODEv49". No package code keys
# on the old name (only roxygen docs mention it); update those with this file.

suppressPackageStartupMessages(library(data.table))

# ---- arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}
exeter_file  <- get_arg("--exeter")
current_file <- get_arg("--current")
out_file     <- get_arg("--out", "epicv2_probe_gene_annotation.csv.gz")
report_file  <- get_arg("--report", "exeter_v3_diff.md")

if (is.null(exeter_file) || !file.exists(exeter_file)) {
  stop("--exeter must point to the Exeter re-annotated manifest v3.0 CSV ",
       "(Zenodo 20704849).", call. = FALSE)
}

read_gz <- function(path, ...) {
  # .rds is what the cpgdirectionData Hub layers are stored as
  # (hub_upload/<name>.rds); csv/csv.gz is what the fat build ships.
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.table(readRDS(path)))
  # same base-R decompression the package uses, so no R.utils dependency
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".txt")
    inn <- gzfile(path, "rb"); outc <- file(tmp, "wb")
    while (length(chunk <- readBin(inn, "raw", 8L * 1024L * 1024L))) writeBin(chunk, outc)
    close(inn); close(outc)
    path <- tmp
  }
  fread(path, showProgress = FALSE, ...)
}

# write .rds when the target name says so (the Hub format), else gzipped csv
write_out <- function(d, path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) saveRDS(as.data.frame(d), path, compress = "xz")
  else fwrite(d, path, compress = "gzip")
}

# ---- current table ------------------------------------------------------------
if (!is.null(current_file)) {
  cur <- read_gz(current_file)
} else {
  if (!requireNamespace("cpgdirection", quietly = TRUE))
    stop("Either pass --current or install cpgdirection so the packaged table ",
         "can be loaded.", call. = FALSE)
  cur <- cpgdirection:::.cpgd_data("epicv2_probe_gene_annotation")
}
cur <- as.data.table(cur)
req <- c("cpg_id", "target_gene", "annotation_source", "n_annotation_sources",
         "src_Illumina_UCSC_RefGene", "src_Illumina_GencodeV41",
         "src_Zhou_GENCODEv41", "src_Exeter_regulatory")
miss <- setdiff(req, names(cur))
if (length(miss)) stop("current table lacks: ", paste(miss, collapse = ", "))
old_ex_col <- grep("^src_Exeter_GENCODEv\\d+$", names(cur), value = TRUE)
if (length(old_ex_col) != 1L) stop("expected exactly one src_Exeter_GENCODEv* column")
old_ex_tok <- sub("^src_", "", old_ex_col)
message(sprintf("current: %d pairs, %d CpGs, Exeter gene track = %s",
                nrow(cur), uniqueN(cur$cpg_id), old_ex_tok))

# ---- Exeter v3.0 ----------------------------------------------------------------
ex <- read_gz(exeter_file)
nm <- names(ex)
pick <- function(pattern, what, required = TRUE, prefer = NULL) {
  hits <- grep(pattern, nm, value = TRUE)
  if (!is.null(prefer)) hits <- c(intersect(prefer, hits), setdiff(hits, prefer))
  if (!length(hits)) {
    if (required) stop("Could not find the ", what, " column in the Exeter file. ",
                       "Columns present:\n  ", paste(nm, collapse = "\n  "))
    return(NA_character_)
  }
  hits[1]
}
id_col   <- pick("^(IlmnID|Name|Probe_ID|probeID)$", "probe identifier")
gene_col <- pick("^GENCODEv\\d+_Gene_Name$", "GENCODE gene name")
feat_col <- pick("^GENCODEv\\d+_Feature_Type$", "GENCODE feature type")
gh_col   <- pick("^In_GeneHancer$", "In_GeneHancer", required = FALSE)
# distance-based regulatory assignments carried over from v2.0; names vary
prom_col <- pick("Promoter.?2000", "Promoter_2000bp gene", required = FALSE)
enh_col  <- pick("Enhancer.?5000", "Enhancer_5000bp gene", required = FALSE)
gencode_release <- sub("^GENCODEv(\\d+)_.*$", "\\1", gene_col)
new_ex_tok <- paste0("Exeter_GENCODEv", gencode_release)
new_ex_col <- paste0("src_", new_ex_tok)
message(sprintf("Exeter file: %d rows; gene track = %s (GENCODE %s); regulatory cols: %s / %s",
                nrow(ex), gene_col, gencode_release,
                ifelse(is.na(prom_col), "none", prom_col),
                ifelse(is.na(enh_col), "none", enh_col)))

probe <- as.character(ex[[id_col]])
m <- regexpr("cg[0-9]{6,}", probe, ignore.case = TRUE)
cpg <- ifelse(m > 0, tolower(substr(probe, m, m + attr(m, "match.length") - 1L)), NA_character_)

# expand a semicolon-delimited gene field into (cpg, gene[, aux]) rows
expand <- function(cpg, genes, aux = NULL) {
  ok <- !is.na(cpg) & !is.na(genes) & nzchar(genes)
  cpg <- cpg[ok]; genes <- genes[ok]; if (!is.null(aux)) aux <- aux[ok]
  gl <- strsplit(genes, ";", fixed = TRUE)
  ln <- lengths(gl)
  o <- data.table(cpg_id = rep(cpg, ln), target_gene = toupper(trimws(unlist(gl))))
  if (!is.null(aux)) {
    al <- strsplit(aux, ";", fixed = TRUE)
    al <- mapply(function(a, n) if (length(a) == n) a else rep(NA_character_, n),
                 al, ln, SIMPLIFY = FALSE)
    o[, aux := unlist(al)]
  }
  o <- o[nzchar(target_gene) & !target_gene %in% c("NA", ".") &
           !grepl("^ENSG[0-9]+", target_gene)]
  unique(o)
}

ex_gene <- expand(cpg, as.character(ex[[gene_col]]), as.character(ex[[feat_col]]))
setnames(ex_gene, "aux", "exeter_feature")
ex_gene <- ex_gene[, list(exeter_feature = {
  f <- unique(exeter_feature[!is.na(exeter_feature) & nzchar(exeter_feature)])
  if (length(f)) paste(sort(f), collapse = ";") else NA_character_
}), by = c("cpg_id", "target_gene")]

ex_reg <- rbindlist(list(
  if (!is.na(prom_col)) expand(cpg, as.character(ex[[prom_col]]))[, reg := "Promoter_2000bp"],
  if (!is.na(enh_col))  expand(cpg, as.character(ex[[enh_col]]))[, reg := "Enhancer_5000bp"]
), use.names = TRUE)
if (nrow(ex_reg)) {
  ex_reg <- ex_reg[, list(exeter_regulatory = paste(sort(unique(reg)), collapse = ";")),
                   by = c("cpg_id", "target_gene")]
}
gh <- if (!is.na(gh_col)) {
  unique(data.table(cpg_id = cpg, in_genehancer = as.logical(ex[[gh_col]]))[!is.na(cpg_id)])[
    , list(in_genehancer = any(in_genehancer %in% TRUE)), by = "cpg_id"]
} else NULL

message(sprintf("Exeter v%s: %d gene-track pairs, %d regulatory pairs, %d CpGs with GeneHancer flag",
                gencode_release, nrow(ex_gene), nrow(ex_reg),
                if (is.null(gh)) 0L else sum(gh$in_genehancer)))

# ---- strip the OLD Exeter contribution -----------------------------------------
before <- copy(cur)
old_pairs <- cur[, paste(cpg_id, target_gene)]
cur[, (old_ex_col) := NULL]
cur[, src_Exeter_regulatory := FALSE]
cur[, exeter_feature := NA_character_]
if ("exeter_regulatory" %in% names(cur)) cur[, exeter_regulatory := NA_character_]
cur[, in_genehancer := NA]
non_ex <- c("src_Illumina_UCSC_RefGene", "src_Illumina_GencodeV41", "src_Zhou_GENCODEv41")
cur[, .keep := Reduce(`|`, lapply(.SD, function(x) x %in% TRUE)), .SDcols = non_ex]
dropped_only_exeter <- cur[.keep == FALSE]
cur <- cur[.keep == TRUE][, .keep := NULL]

# per-CpG coordinate/probe provenance to fill new rows
cpg_info <- unique(before[, c("cpg_id", "manifest_probe_id", "chr", "position", "strand"),
                          with = FALSE], by = "cpg_id")

# ---- add the NEW Exeter contribution ----------------------------------------------
cur[, (new_ex_col) := FALSE]
setkeyv(cur, c("cpg_id", "target_gene"))
setkeyv(ex_gene, c("cpg_id", "target_gene"))
cur[ex_gene, (new_ex_col) := TRUE]
cur[ex_gene, exeter_feature := i.exeter_feature]
new_gene_rows <- ex_gene[!cur, on = c("cpg_id", "target_gene")]
if (nrow(ex_reg)) {
  setkeyv(ex_reg, c("cpg_id", "target_gene"))
  cur[ex_reg, `:=`(src_Exeter_regulatory = TRUE, exeter_regulatory = i.exeter_regulatory)]
  new_reg_rows <- ex_reg[!cur, on = c("cpg_id", "target_gene")]
} else new_reg_rows <- ex_reg
new_rows <- merge(new_gene_rows, new_reg_rows, by = c("cpg_id", "target_gene"), all = TRUE)
if (nrow(new_rows)) {
  new_rows[, (new_ex_col) := !is.na(exeter_feature)]
  new_rows[, src_Exeter_regulatory := if ("exeter_regulatory" %in% names(new_rows))
    !is.na(exeter_regulatory) else FALSE]
  new_rows[, `:=`(array = "EPICv2", src_Illumina_UCSC_RefGene = FALSE,
                  src_Illumina_GencodeV41 = FALSE, src_Zhou_GENCODEv41 = FALSE,
                  refgene_group = NA_character_, zhou_dist_tss = NA_real_,
                  zhou_tx_types = NA_character_)]
  new_rows <- merge(new_rows, cpg_info, by = "cpg_id", all.x = TRUE)
  cur <- rbindlist(list(cur, new_rows), use.names = TRUE, fill = TRUE)
}
if (!is.null(gh)) {
  cur[, in_genehancer := NULL]
  cur <- merge(cur, gh, by = "cpg_id", all.x = TRUE)
}
for (cc in c(non_ex, new_ex_col, "src_Exeter_regulatory"))
  set(cur, which(is.na(cur[[cc]])), cc, FALSE)

# ---- recompute provenance ------------------------------------------------------
tracks <- c("Illumina_UCSC_RefGene", "Illumina_GencodeV41", "Zhou_GENCODEv41",
            new_ex_tok, "Exeter_regulatory")
flag_cols <- paste0("src_", tracks)
cur[, n_annotation_sources := Reduce(`+`, lapply(.SD, as.integer)), .SDcols = flag_cols]
cur[, annotation_source := apply(as.matrix(.SD), 1L, function(r) paste(tracks[r], collapse = ";")),
    .SDcols = flag_cols]
cur <- cur[n_annotation_sources > 0L]
setcolorder(cur, intersect(c("cpg_id", "target_gene", "array", "annotation_source",
                             "n_annotation_sources", flag_cols, "refgene_group",
                             "exeter_feature", "exeter_regulatory", "in_genehancer",
                             "zhou_dist_tss", "zhou_tx_types", "manifest_probe_id",
                             "chr", "position", "strand"), names(cur)))
setorderv(cur, c("cpg_id", "target_gene"))

# ---- diff report -----------------------------------------------------------------
new_pairs <- cur[, paste(cpg_id, target_gene)]
gained <- setdiff(new_pairs, old_pairs)
lost   <- setdiff(old_pairs, new_pairs)
old_feat <- before[, c("cpg_id", "target_gene", "exeter_feature"), with = FALSE]
cmp <- merge(old_feat, cur[, c("cpg_id", "target_gene", "exeter_feature"), with = FALSE],
             by = c("cpg_id", "target_gene"), suffixes = c("_v47", "_v49"))
feat_changed <- cmp[!is.na(exeter_feature_v47) & !is.na(exeter_feature_v49) &
                      exeter_feature_v47 != exeter_feature_v49]
tss_re <- "TSS(200|1500)"
tss_flip <- cmp[grepl(tss_re, exeter_feature_v47) != grepl(tss_re, exeter_feature_v49)]
agree_before <- round(100 * mean(before$n_annotation_sources >= 2), 1)
agree_after  <- round(100 * mean(cur$n_annotation_sources >= 2), 1)

rep <- c(
  sprintf("# Exeter re-annotation swap: %s -> %s", old_ex_tok, new_ex_tok),
  "", sprintf("Built %s from `%s`.", format(Sys.time(), "%Y-%m-%d %H:%M"), basename(exeter_file)), "",
  "| metric | before | after |", "|---|---:|---:|",
  sprintf("| CpG-gene pairs | %d | %d |", nrow(before), nrow(cur)),
  sprintf("| CpGs with >= 1 gene | %d | %d |", uniqueN(before$cpg_id), uniqueN(cur$cpg_id)),
  sprintf("| genes | %d | %d |", uniqueN(before$target_gene), uniqueN(cur$target_gene)),
  sprintf("| pairs supported by >= 2 tracks | %s%% | %s%% |", agree_before, agree_after),
  sprintf("| pairs on the Exeter gene track | %d | %d |",
          sum(before[[old_ex_col]] %in% TRUE), sum(cur[[new_ex_col]] %in% TRUE)),
  sprintf("| pairs on the Exeter regulatory track | %d | %d |",
          sum(before$src_Exeter_regulatory %in% TRUE), sum(cur$src_Exeter_regulatory %in% TRUE)),
  "",
  sprintf("* pairs gained: %d; pairs lost: %d (lost pairs were supported only by the old Exeter tracks: %d)",
          length(gained), length(lost), nrow(dropped_only_exeter)),
  sprintf("* pairs whose Exeter feature label changed: %d", nrow(feat_changed)),
  sprintf("* pairs that gained or lost a TSS200/TSS1500 label: %d", nrow(tss_flip)),
  "", "Examples of changed features (first 20):", "",
  "| cpg | gene | v47 feature | v49 feature |", "|---|---|---|---|",
  if (nrow(feat_changed)) feat_changed[1:min(20L, .N),
    sprintf("| %s | %s | %s | %s |", cpg_id, target_gene, exeter_feature_v47, exeter_feature_v49)] else "| (none) | | | |"
)
writeLines(rep, report_file)
write_out(cur, out_file)
message(sprintf("wrote %s (%d pairs, %d CpGs) and %s", out_file, nrow(cur),
                uniqueN(cur$cpg_id), report_file))
message("Remember: update the roxygen in R/manifest_annotation.R (", old_ex_tok,
        " -> ", new_ex_tok, "), NEWS.md, and the cpgdirectionData Hub resource.")
