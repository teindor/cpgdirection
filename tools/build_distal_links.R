#!/usr/bin/env Rscript
# Build the `distal_links` layer: EPIC v2 CpG -> gene links supported by
# regulatory-element evidence that is better than proximity, keyed by cg ID.
#
# Three sources, one table, one evidence class ("predicted or physical distal
# link"), kept apart from the five proximity annotation tracks and from every
# direction layer. A distal link proposes a TARGET; it carries no sign.
#
#   rE2G   ENCODE-rE2G enhancer-gene predictions (Gschwind et al. 2026, Nature,
#          doi:10.1038/s41586-026-10781-4; CC BY 4.0; GRCh38). One file per
#          biosample from the ENCODE portal, annotation type "element gene
#          regulatory interaction predictions" (e.g. ENCSR087OMI). Columns the
#          pipeline writes (EngreitzLab/ENCODE_rE2G): #chr start end name class
#          TargetGene TargetGeneTSS ... Score (or ENCODE-rE2G.Score / E2G.Score)
#          plus a thresholded indicator column in some releases.
#   scE2G  Single-cell enhancer-gene maps (Sheth et al. 2026, Nat Genet,
#          doi:10.1038/s41588-026-02695-8; CC BY 4.0; GRCh38). Same column
#          layout as rE2G, one file per cell cluster (PBMC subsets, BMMC, islets).
#   HCEA   Human Cell Epigenome Atlas loops (Zhou et al. 2026, Science,
#          doi:10.1126/science.adx0673; hg38). Cell-type loop calls as BEDPE:
#          chr1 start1 end1 chr2 start2 end2 [name score ...]. A CpG in one
#          anchor is linked to every gene whose TSS lies in the other anchor;
#          the gene TSS table is the package's own cpgd_gene_tss() lifted to
#          hg38, or any 4-column BED of TSSs you pass with --tss-hg38.
#
# Probe coordinates. The Illumina EPIC v2 manifest carries hg38 (CHR, MAPINFO).
# Intersect in hg38, then key on the cg ID; nothing here needs hg19 and nothing
# downstream needs hg38, because every consumer joins on cpg_id.
#
# Usage
#   Rscript tools/build_distal_links.R \
#     --manifest  ~/EPIC-8v2-0_A2.csv \
#     --re2g      ~/encode_re2g/              # dir of *.tsv(.gz), one per biosample
#     --sce2g     ~/sce2g/                     # dir of *.tsv(.gz), one per cell type
#     --loops     ~/hcea_loops/                # dir of *.bedpe(.gz), one per cell type
#     --tss-hg38  ~/gencode_v49_tss.bed        # chr start end gene (0-based) for --loops
#     --biosample-map tools/distal_biosamples.csv   # file,source,biosample,tissue_class
#     --min-score 0                            # keep links at or above this score
#     --keep-abc                               # also keep ABC-model files (see below)
#     --out inst/extdata/distal_links.csv.gz
#
# ENCODE rE2G annotation sets ship two prediction tables per biosample: the
# ENCODE-rE2G model (header carries ABC.Score.Feature or an rE2G score column;
# thresholded scores run ~0.2-1) and the older ABC model it was trained from
# (header carries ABC.Score and powerlaw.Score; thresholded scores start ~0.018).
# The two score scales must not share one column, so ABC files are skipped
# unless --keep-abc is given, in which case they enter as source = "ABC".
#
# --biosample-map is optional: a CSV with columns file (basename), source
# (rE2G|scE2G|HCEA), biosample (free text), tissue_class (one of blood, brain,
# other) so the package can filter by context. Without it the biosample is the
# file name and tissue_class is "other".
#
# Output columns
#   cpg_id target_gene source biosample tissue_class score element_chr
#   element_start element_end element_class

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
manifest  <- get_arg("--manifest")
re2g_dir  <- get_arg("--re2g")
sce2g_dir <- get_arg("--sce2g")
loop_dir  <- get_arg("--loops")
tss_bed   <- get_arg("--tss-hg38")
bs_map    <- get_arg("--biosample-map")
min_score <- as.numeric(get_arg("--min-score", "0"))
keep_abc  <- "--keep-abc" %in% args
out_file  <- get_arg("--out", "distal_links.csv.gz")

if (is.null(manifest) || !file.exists(manifest))
  stop("--manifest must point to the Illumina EPIC v2 manifest CSV (hg38 coordinates).")
if (is.null(re2g_dir) && is.null(sce2g_dir) && is.null(loop_dir))
  stop("give at least one of --re2g, --sce2g, --loops")

read_gz <- function(path, ...) {
  # .rds is what the cpgdirectionData Hub layers are stored as
  # (hub_upload/<name>.rds); csv/csv.gz is what the fat build ships.
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.table(readRDS(path)))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".txt")
    inn <- gzfile(path, "rb"); outc <- file(tmp, "wb")
    while (length(chunk <- readBin(inn, "raw", 8L * 1024L * 1024L))) writeBin(chunk, outc)
    close(inn); close(outc); path <- tmp
  }
  fread(path, showProgress = FALSE, ...)
}

# write .rds when the target name says so (the Hub format), else gzipped csv
write_out <- function(d, path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) saveRDS(as.data.frame(d), path, compress = "xz")
  else fwrite(d, path, compress = "gzip")
}

# ---- EPIC v2 probes in hg38 ------------------------------------------------------
first <- readLines(manifest, n = 200L, warn = FALSE)
assay_at <- grep("^\\[Assay\\]", first)[1]
m <- read_gz(manifest, skip = if (is.na(assay_at)) 0L else assay_at, fill = TRUE)
ctrl <- grep("^\\[Controls\\]", m[[1]])[1]
if (!is.na(ctrl)) m <- m[seq_len(ctrl - 1L)]
idc <- names(m)[names(m) %in% c("IlmnID", "Name", "Probe_ID")][1]
chrc <- names(m)[names(m) %in% c("CHR", "Chr", "chr")][1]
posc <- names(m)[names(m) %in% c("MAPINFO", "mapinfo", "Start_hg38", "pos")][1]
pid <- as.character(m[[idc]])
mm <- regexpr("cg[0-9]{6,}", pid, ignore.case = TRUE)
probes <- unique(data.table(
  cpg_id = ifelse(mm > 0, tolower(substr(pid, mm, mm + attr(mm, "match.length") - 1L)), NA_character_),
  chr = sub("^chr", "", as.character(m[[chrc]])),
  pos = suppressWarnings(as.integer(m[[posc]])))[!is.na(cpg_id) & !is.na(pos)], by = "cpg_id")
probes[, `:=`(start = pos, end = pos)]        # 1-bp interval, 1-based
setkeyv(probes, c("chr", "start", "end"))
message(sprintf("manifest: %d CpGs with hg38 coordinates", nrow(probes)))

# ---- biosample map ---------------------------------------------------------------
bsm <- if (!is.null(bs_map)) fread(bs_map) else NULL
lookup_bs <- function(f, src) {
  if (!is.null(bsm)) {
    r <- bsm[file == basename(f) & source == src]
    if (nrow(r)) return(list(biosample = r$biosample[1], tissue_class = r$tissue_class[1]))
  }
  list(biosample = sub("\\.(tsv|txt|bed|bedpe)(\\.gz)?$", "", basename(f)), tissue_class = "other")
}

# ---- E2G-style files (rE2G and scE2G share the layout) -------------------------
read_e2g <- function(f, src) {
  d <- read_gz(f)
  nm <- names(d)
  chr_c <- nm[nm %in% c("#chr", "chr", "chrom")][1]
  st_c  <- nm[nm %in% c("start", "chromStart")][1]
  en_c  <- nm[nm %in% c("end", "chromEnd")][1]
  g_c   <- nm[nm %in% c("TargetGene", "gene", "target_gene")][1]
  sc_c  <- nm[grepl("^(ENCODE-rE2G\\.Score|E2G\\.Score|Score|score|scE2G\\.Score)$", nm)][1]
  cl_c  <- nm[nm %in% c("class", "element_class")][1]
  if (anyNA(c(chr_c, st_c, en_c, g_c, sc_c)))
    stop(basename(f), ": cannot find chr/start/end/TargetGene/Score columns; have: ",
         paste(nm, collapse = ", "))
  # ENCODE sets carry an ABC-model table beside the rE2G one (see header note)
  is_abc <- src == "rE2G" && "ABC.Score" %in% nm &&
            !any(grepl("ABC\\.Score\\.Feature|rE2G|E2G\\.Score", nm))
  if (is_abc) {
    if (!keep_abc) { message("  skip ABC-model file: ", basename(f)); return(NULL) }
    src <- "ABC"
  }
  bs <- lookup_bs(f, if (src == "ABC") "rE2G" else src)
  e <- data.table(chr = sub("^chr", "", as.character(d[[chr_c]])),
                  start = as.integer(d[[st_c]]) + 1L,       # BED 0-based -> 1-based
                  end   = as.integer(d[[en_c]]),
                  target_gene = toupper(trimws(as.character(d[[g_c]]))),
                  score = as.numeric(d[[sc_c]]),
                  element_class = if (!is.na(cl_c)) as.character(d[[cl_c]]) else NA_character_)
  e <- e[!is.na(score) & score >= min_score & nzchar(target_gene)]
  if (!nrow(e)) return(NULL)
  setkeyv(e, c("chr", "start", "end"))
  ov <- foverlaps(probes, e, type = "within", nomatch = NULL)
  ov[, list(cpg_id, target_gene, source = src, biosample = bs$biosample,
            tissue_class = bs$tissue_class, score,
            element_chr = chr, element_start = start - 1L, element_end = end,
            element_class)]
}

# ---- loop files -------------------------------------------------------------------
read_loops <- function(f, tss) {
  d <- read_gz(f, header = FALSE)
  if (ncol(d) < 6L) stop(basename(f), ": BEDPE needs >= 6 columns")
  setnames(d, 1:6, c("c1", "s1", "e1", "c2", "s2", "e2"))
  d[, `:=`(c1 = sub("^chr", "", c1), c2 = sub("^chr", "", c2),
           s1 = as.integer(s1) + 1L, e1 = as.integer(e1),
           s2 = as.integer(s2) + 1L, e2 = as.integer(e2))]
  d[, loop_id := .I]
  sc <- if (ncol(d) >= 8L && is.numeric(d[[8L]])) as.numeric(d[[8L]]) else NA_real_
  d[, score := sc]
  bs <- lookup_bs(f, "HCEA")
  # anchor A holds the CpG, anchor B holds a TSS -- in both orientations
  link_side <- function(ac, as_, ae, bc, bs_, be) {
    A <- d[, list(loop_id, score, chr = get(ac), start = get(as_), end = get(ae))]
    B <- d[, list(loop_id, chr = get(bc), start = get(bs_), end = get(be))]
    setkeyv(A, c("chr", "start", "end")); setkeyv(B, c("chr", "start", "end"))
    pa <- foverlaps(probes, A, type = "within", nomatch = NULL)[, list(cpg_id, loop_id, score,
                                                                        element_chr = chr, element_start = start - 1L, element_end = end)]
    tb <- foverlaps(tss, B, type = "within", nomatch = NULL)[, list(loop_id, target_gene)]
    merge(pa, unique(tb), by = "loop_id", allow.cartesian = TRUE)
  }
  r <- rbindlist(list(link_side("c1", "s1", "e1", "c2", "s2", "e2"),
                      link_side("c2", "s2", "e2", "c1", "s1", "e1")))
  if (!nrow(r)) return(NULL)
  r[, list(cpg_id, target_gene, source = "HCEA", biosample = bs$biosample,
           tissue_class = bs$tissue_class, score, element_chr, element_start, element_end,
           element_class = "loop_anchor")]
}

parts <- list()
if (!is.null(re2g_dir)) {
  fs <- list.files(re2g_dir, pattern = "\\.(tsv|txt|bed)(\\.gz)?$", full.names = TRUE)
  message(sprintf("rE2G: %d files", length(fs)))
  parts$re2g <- rbindlist(lapply(fs, read_e2g, src = "rE2G"))
}
if (!is.null(sce2g_dir)) {
  fs <- list.files(sce2g_dir, pattern = "\\.(tsv|txt|bed)(\\.gz)?$", full.names = TRUE)
  message(sprintf("scE2G: %d files", length(fs)))
  parts$sce2g <- rbindlist(lapply(fs, read_e2g, src = "scE2G"))
}
if (!is.null(loop_dir)) {
  if (is.null(tss_bed) || !file.exists(tss_bed))
    stop("--loops needs --tss-hg38 (BED: chr start end gene, hg38)")
  tss <- read_gz(tss_bed, header = FALSE)[, 1:4]
  setnames(tss, c("chr", "start", "end", "target_gene"))
  tss[, `:=`(chr = sub("^chr", "", chr), start = as.integer(start) + 1L, end = as.integer(end),
             target_gene = toupper(trimws(target_gene)))]
  tss <- tss[!grepl("^ENSG[0-9]+", target_gene)]
  setkeyv(tss, c("chr", "start", "end"))
  fs <- list.files(loop_dir, pattern = "\\.(bedpe|bed|tsv)(\\.gz)?$", full.names = TRUE)
  message(sprintf("HCEA loops: %d files", length(fs)))
  parts$hcea <- rbindlist(lapply(fs, read_loops, tss = tss))
}
links <- rbindlist(parts, use.names = TRUE, fill = TRUE)
if (!nrow(links)) stop("no CpG-gene links produced")
links <- links[!grepl("^ENSG[0-9]+", target_gene)]
# one row per cpg x gene x source x biosample: keep the strongest element
setorderv(links, c("cpg_id", "target_gene", "source", "biosample", "score"),
          order = c(1L, 1L, 1L, 1L, -1L), na.last = TRUE)
links <- unique(links, by = c("cpg_id", "target_gene", "source", "biosample"))
setcolorder(links, c("cpg_id", "target_gene", "source", "biosample", "tissue_class", "score",
                     "element_chr", "element_start", "element_end", "element_class"))
write_out(links, out_file)
message(sprintf("wrote %s: %d link rows, %d CpG-gene pairs, %d CpGs, %d genes; by source: %s",
                out_file, nrow(links), uniqueN(links[, paste(cpg_id, target_gene)]),
                uniqueN(links$cpg_id), uniqueN(links$target_gene),
                paste(sprintf("%s=%d", names(table(links$source)), table(links$source)), collapse = ", ")))
