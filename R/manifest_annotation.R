# EPIC v2 probe -> gene annotation from the Illumina manifest.
# New in 2.4.0. This is the resource that makes a BARE CpG identifier
# sufficient for automatic target-gene discovery: cpg_gene_pairs("cg26261055")
# can recover CRHBP without the caller supplying a gene.
#
# The table is a MAPPING source, not proof of regulation. It answers "Illumina
# annotates this CpG to this gene/region", not "this CpG regulates this gene".
# Its provenance is therefore kept distinct from measured eQTMs and SMR
# throughout the pair workflow (mapping_source = "EPICv2_manifest",
# mapping_strength = "annotation").

#' EPIC v2 probe-to-gene annotation (long form, multi-source)
#'
#' The packaged long-form EPIC v2 probe-to-gene table: one row per CpG x gene,
#' with every annotated gene preserved and, from 2.5.0, every assignment
#' cross-checked across independent annotation tracks. A mis-assigned CpG sends
#' every downstream direction to the wrong gene, so no single source is
#' trusted alone; the shipped table is the UNION of five tracks with per-track
#' support flags and an agreement count:
#' \itemize{
#'   \item \code{Illumina_UCSC_RefGene} -- the official Illumina EPIC-8v2-0
#'     manifest RefGene annotation;
#'   \item \code{Illumina_GencodeV41} -- the manifest's GENCODE v41 track
#'     (symbol-named entries);
#'   \item \code{Zhou_GENCODEv41} -- the Zhou lab (SeSAMe) transcript-aware
#'     GENCODE v41 annotation, contributing \code{zhou_dist_tss} (signed bp to
#'     the nearest annotated TSS of that gene) and transcript types;
#'   \item \code{Exeter_GENCODEv47} -- the Exeter re-annotation (Mallabar-Rimmer
#'     et al. 2025, GENCODE release 47): gene-body and TSS features per
#'     transcript;
#'   \item \code{Exeter_regulatory} -- the Exeter distance-based regulatory
#'     assignments (\code{Promoter_2000bp} / \code{Enhancer_5000bp} element
#'     genes), plus the probe-level \code{in_genehancer} flag.
#' }
#' Semicolon-delimited multi-gene fields are always fully expanded -- never
#' reduced to the first gene -- and Ensembl-ID-only entries are excluded
#' (nothing else in the package keys on Ensembl IDs).
#'
#' \code{n_annotation_sources} is the number of independent tracks supporting
#' the pair (84.4\% of shipped pairs have 2 or more; 54.8\% have 3 or more).
#' Downstream consumers such as DMSA can weight or filter on it.
#'
#' @section What this table is and is not:
#' Annotation is positional: "this source annotates this CpG to this gene". It
#' does not establish that the CpG regulates the gene, and it is deliberately
#' kept distinct from measured eQTM and SMR evidence, which carry different
#' biological meaning. In \code{\link{cpg_gene_pairs}} it contributes candidate
#' targets; the direction evidence for each pair is resolved separately.
#' Probe RELIABILITY (does the probe read the right locus at all?) is likewise
#' a separate question, answered by \code{\link{cpgd_probe_qc}}.
#'
#' @section Coordinates:
#' \code{chr}, \code{position} and \code{strand} are hg38 manifest
#' coordinates, carried for provenance. The package's distance computations use
#' the hg19 positions from \code{\link{cpgd_cpg_positions}}, not these.
#'
#' @return A \code{data.table} with columns \code{cpg_id}, \code{target_gene},
#'   \code{array}, \code{annotation_source} (semicolon-joined supporting
#'   tracks), \code{n_annotation_sources}, the per-track logical flags
#'   \code{src_Illumina_UCSC_RefGene}, \code{src_Illumina_GencodeV41},
#'   \code{src_Zhou_GENCODEv41}, \code{src_Exeter_GENCODEv47},
#'   \code{src_Exeter_regulatory}, plus \code{refgene_group} (Illumina),
#'   \code{exeter_feature}, \code{in_genehancer}, \code{zhou_dist_tss},
#'   \code{zhou_tx_types}, \code{manifest_probe_id}, \code{chr},
#'   \code{position} and \code{strand}.
#' @examples
#' m <- cpgd_manifest_genes()
#' m[m$cpg_id == "cg26261055", ]   # CRHBP, 4 independent sources, TSS -899 bp
#' @export
cpgd_manifest_genes <- function() {
  if (!is.null(.cpgd_env$manifest_genes)) return(.cpgd_env$manifest_genes)
  M <- .cpgd_data("epicv2_probe_gene_annotation")
  req <- c("cpg_id", "target_gene", "array", "annotation_source", "refgene_group")
  miss <- setdiff(req, names(M))
  if (length(miss)) {
    stop("The EPIC v2 annotation table is missing expected columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  data.table::setkeyv(M, c("cpg_id", "target_gene"))
  .cpgd_env$manifest_genes <- M
  M
}


#' Rebuild the EPIC v2 probe-to-gene annotation from an Illumina manifest
#'
#' Regenerates the long-form table served by \code{\link{cpgd_manifest_genes}}
#' from the official Illumina EPIC v2 manifest CSV (e.g.
#' \code{EPIC-8v2-0_A2.csv}). Runtime users never need this; the packaged
#' resource ships with the package. Call it to refresh the table from a newer
#' manifest release or to customise the annotation.
#'
#' The builder:
#' \enumerate{
#'   \item reads the manifest, skipping the \code{[Assay]} header preamble and
#'     stopping before the \code{[Controls]} section;
#'   \item locates the probe identifier column (\code{IlmnID}/\code{Name}) and
#'     the \code{UCSC_RefGene_Name} / \code{UCSC_RefGene_Group} columns (plus
#'     \code{GencodeV41_Name} / \code{GencodeV41_Group} where present);
#'   \item splits every semicolon-delimited multi-gene field, pairing gene and
#'     group tokens per transcript -- a mismatched pair count yields
#'     \code{NA} groups rather than guessed ones, and the field is never
#'     reduced to its first gene;
#'   \item uppercases and trims gene symbols, dropping Ensembl-ID-only entries;
#'   \item collapses to one row per CpG x gene, pooling transcript groups and
#'     retaining the raw manifest fields for provenance;
#'   \item writes the compressed table.
#' }
#'
#' @param manifest_file Path to the Illumina EPIC v2 manifest CSV (may be
#'   gzipped).
#' @param out Path for the resulting \code{.csv.gz}. Defaults to a file named
#'   like the packaged resource, in the working directory. To make a rebuilt
#'   table the one \code{\link{cpgd_manifest_genes}} serves, overwrite the
#'   packaged file under \code{inst/extdata/} and reinstall.
#' @return Invisibly, the path written.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' writeLines(c("[Heading]", "Descriptor File Name,EPIC-8v2-0_TEST",
#'   "[Assay]",
#'   "IlmnID,Name,CHR,MAPINFO,Strand_FR,UCSC_RefGene_Name,UCSC_RefGene_Group",
#'   "cg10000001_TC21,cg10000001,chr1,1000,F,GENEA;GENEB,TSS200;TSS1500",
#'   "[Controls]", "ctl1,negative,0,0"), f)
#' out <- tempfile(fileext = ".csv.gz")
#' cpgd_build_manifest_genes(f, out = out)
#' read.csv(out)
#' @export
cpgd_build_manifest_genes <- function(manifest_file,
                                      out = "epicv2_probe_gene_annotation.csv.gz") {
  if (missing(manifest_file) || !file.exists(manifest_file)) {
    stop("Supply the path to an Illumina EPIC v2 manifest CSV ",
         "(e.g. EPIC-8v2-0_A2.csv).", call. = FALSE)
  }

  # The Illumina CSV has a [Heading] preamble before [Assay] and a [Controls]
  # section after the probes. Find the preamble without loading the whole file.
  con <- if (grepl("\\.gz$", manifest_file)) gzfile(manifest_file, "r")
         else file(manifest_file, "r")
  first <- readLines(con, n = 200L, warn = FALSE)
  close(con)
  assay_at <- grep("^\\[Assay\\]", first)[1]
  skip <- if (is.na(assay_at)) 0L else assay_at
  m <- .cpgd_fread(manifest_file, skip = skip, fill = TRUE)
  # rows past [Controls] have a bracketed first field; cut them
  ctrl <- grep("^\\[Controls\\]", m[[1]])[1]
  if (!is.na(ctrl)) m <- m[seq_len(ctrl - 1L)]

  nm <- names(m)
  idc <- nm[nm %in% c("IlmnID", "Name", "Probe_ID")][1]
  gnc <- nm[grepl("UCSC_RefGene_Name|^RefGene_Name$", nm)][1]
  ggc <- nm[grepl("UCSC_RefGene_Group|^RefGene_Group$", nm)][1]
  if (is.na(idc) || is.na(gnc)) {
    stop("Could not locate the probe identifier and UCSC_RefGene_Name columns ",
         "in the manifest.", call. = FALSE)
  }
  enc <- nm[grepl("GencodeV\\d+_Name", nm)][1]
  egc <- nm[grepl("GencodeV\\d+_Group", nm)][1]
  chrc <- nm[nm %in% c("CHR", "Chr", "chr")][1]
  posc <- nm[nm %in% c("MAPINFO", "mapinfo", "Start_hg38", "pos")][1]
  strc <- nm[nm %in% c("Strand_FR", "Strand", "strand")][1]

  probe_id <- as.character(m[[idc]])
  posn <- regexpr("cg[0-9]{6,}", probe_id, ignore.case = TRUE)
  lenn <- attr(posn, "match.length")
  cpg <- ifelse(posn > 0L, tolower(substr(probe_id, posn, posn + lenn - 1L)),
                NA_character_)

  base <- data.table::data.table(
    manifest_probe_id = probe_id,
    cpg_id   = cpg,
    chr      = if (!is.na(chrc)) as.character(m[[chrc]]) else NA_character_,
    position = if (!is.na(posc)) suppressWarnings(as.integer(m[[posc]])) else NA_integer_,
    strand   = if (!is.na(strc)) as.character(m[[strc]]) else NA_character_,
    ucsc_name  = as.character(m[[gnc]]),
    ucsc_group = if (!is.na(ggc)) as.character(m[[ggc]]) else "",
    gen_name   = if (!is.na(enc)) as.character(m[[enc]]) else "",
    gen_group  = if (!is.na(egc)) as.character(m[[egc]]) else "")
  base <- base[!is.na(base$cpg_id), ]

  ann <- .cpgd_manifest_long(base)
  data.table::fwrite(ann, out, compress = "gzip")
  message(sprintf("wrote %s: %d CpG-gene rows, %d CpGs, %d genes",
                  out, nrow(ann), length(unique(ann$cpg_id)),
                  length(unique(ann$target_gene))))
  invisible(out)
}


# Shared long-form expansion: every semicolon-delimited gene becomes its own
# row, gene and group tokens are paired per transcript, Gencode supplements
# UCSC, and the result is one row per CpG x gene with pooled groups.
.cpgd_manifest_long <- function(base) {
  expand_pairs <- function(d, name_col, group_col, source_label) {
    d <- d[nzchar(d[[name_col]]) & !is.na(d[[name_col]]), ]
    if (!nrow(d)) return(data.table::data.table())
    ns <- strsplit(d[[name_col]],  ";", fixed = TRUE)
    gs <- strsplit(d[[group_col]], ";", fixed = TRUE)
    ln <- lengths(ns)
    gs <- mapply(function(g, n) if (length(g) == n) g else rep(NA_character_, n),
                 gs, ln, SIMPLIFY = FALSE)
    o <- data.table::data.table(
      manifest_probe_id = rep(d$manifest_probe_id, ln),
      cpg_id   = rep(d$cpg_id, ln),
      chr      = rep(d$chr, ln),
      position = rep(d$position, ln),
      strand   = rep(d$strand, ln),
      target_gene    = toupper(trimws(unlist(ns))),
      refgene_group  = unlist(gs),
      refgene_source = source_label,
      raw_name  = rep(d[[name_col]], ln),
      raw_group = rep(d[[group_col]], ln))
    o <- o[nzchar(o$target_gene) & !o$target_gene %in% c("NA", "."), ]
    # Ensembl gene IDs are unusable as symbols anywhere else in the package
    # (lookup / measured / SMR tables key on symbols); where a symbol exists
    # the UCSC track already carries it.
    o[!grepl("^ENSG[0-9]+", o$target_gene), ]
  }

  long <- data.table::rbindlist(list(
    expand_pairs(base, "ucsc_name", "ucsc_group", "UCSC_RefGene"),
    expand_pairs(base, "gen_name",  "gen_group",  "GencodeV41")))
  if (!nrow(long)) {
    stop("No gene annotation rows could be extracted from the manifest.",
         call. = FALSE)
  }
  # UCSC provenance preferred where both tracks list the gene ("U" > "G")
  data.table::setorderv(long, c("cpg_id", "target_gene", "refgene_source"),
                        order = c(1L, 1L, -1L))
  ann <- long[, list(
    array             = "EPICv2",
    annotation_source = "Illumina_EPICv2",
    refgene_source    = get("refgene_source")[1L],
    refgene_group     = {
      g <- unique(get("refgene_group")[!is.na(get("refgene_group")) &
                                         nzchar(get("refgene_group"))])
      if (length(g)) paste(sort(g), collapse = ";") else NA_character_
    },
    manifest_probe_id = paste(sort(unique(get("manifest_probe_id"))), collapse = ";"),
    chr      = get("chr")[1L],
    position = get("position")[1L],
    strand   = get("strand")[1L],
    UCSC_RefGene_Name_raw  = get("raw_name")[1L],
    UCSC_RefGene_Group_raw = get("raw_group")[1L]),
    by = c("cpg_id", "target_gene")]
  data.table::setcolorder(ann, c(
    "cpg_id", "target_gene", "array", "annotation_source", "refgene_source",
    "refgene_group", "manifest_probe_id", "chr", "position", "strand",
    "UCSC_RefGene_Name_raw", "UCSC_RefGene_Group_raw"))
  data.table::setorderv(ann, c("cpg_id", "target_gene"))
  ann
}


#' EPIC v2 probe reliability (QC) table
#'
#' Per-CpG probe quality assembled from independent sources: the Zhou lab
#' (SeSAMe) EPIC v2 general mask -- degenerate or non-unique mapping, SNP
#' contamination at the extension base, and other design faults -- the Garvan
#' (Peters et al. 2024) experimental evidence -- WGBS-confirmed
#' cross-hybridization, BLAT-predicted off-targets, mapping mismatches -- and
#' a verification of the packaged hg19 position against the Zhou hg19
#' manifest. The two cross-hybridization angles corroborate each other (64\%
#' of WGBS-confirmed off-target probes are also sequence-masked by Zhou) and
#' each catches probes the other misses.
#'
#' A probe that reads the wrong locus assigns methylation to the wrong gene
#' entirely -- the one failure mode no downstream direction evidence can
#' repair. \code{\link{cpg_gene_pairs}} therefore consults this table and, by
#' default (\code{probe_qc = "exclude"}), drops CpGs whose every replicate
#' probe carries the general mask.
#'
#' Exclusion-grade flags follow the EVERY-replicate rule: \code{masked_general}
#' (Zhou) and \code{cross_hybridizing} (Garvan, WGBS-evidenced) are TRUE only
#' when every replicate probe of the CpG is flagged; a CpG with at least one
#' clean replicate is kept (\code{masked_partial}). \code{qc_exclude} is their
#' union and is the criterion \code{probe_qc = "exclude"} acts on.
#' Informational flags use the ANY-replicate rule: \code{ch_blat_predicted}
#' (computationally predicted off-targets) and \code{mapping_flagged}
#' (alignment mismatches) mark rows without excluding them.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{n_probes},
#'   \code{n_masked_general}, \code{mask_reasons}, \code{masked_general},
#'   \code{masked_partial}, \code{pos_hg19_verified} (\code{NA} where the CpG
#'   is absent from the verification source), \code{cross_hybridizing},
#'   \code{ch_blat_predicted}, \code{mapping_flagged},
#'   \code{n_offtargets_max}, \code{qc_exclude} and \code{qc_sources}.
#' @examples
#' qc <- cpgd_probe_qc()
#' head(qc[qc$masked_general == TRUE, ], 3)
#' @export
cpgd_probe_qc <- function() {
  if (!is.null(.cpgd_env$probe_qc)) return(.cpgd_env$probe_qc)
  Q <- .cpgd_data("epicv2_probe_qc")
  req <- c("cpg_id", "masked_general", "masked_partial", "mask_reasons",
           "pos_hg19_verified")
  miss <- setdiff(req, names(Q))
  if (length(miss)) {
    stop("The probe QC table is missing expected columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  data.table::setkeyv(Q, "cpg_id")
  .cpgd_env$probe_qc <- Q
  Q
}
