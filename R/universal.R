#' Build a gene TSS annotation for universal mode
#'
#' The catalogue-based prediction only covers genes some catalogue measured. To
#' predict for \emph{any} gene, the package needs transcription start sites,
#' which it does not ship — a gene annotation is large, versioned, and better
#' taken from Bioconductor than frozen into a package.
#'
#' Run this once. It writes a small table that
#' \code{\link{cpg_direction_universal}} then uses.
#'
#' @param out Path for the resulting CSV. Default \code{"gene_tss_hg19.csv"}.
#' @param source \code{"txdb"} uses \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   with \code{org.Hs.eg.db} for symbols. \code{"manifest"} reads an Illumina
#'   manifest instead, via \code{manifest_file}.
#' @param manifest_file Path to an Illumina EPIC manifest CSV, used when
#'   \code{source = "manifest"}.
#'
#' @return Invisibly, the path written. The table has \code{gene}, \code{chr},
#'   \code{tss} and \code{strand}.
#'
#' @details
#' TxDb is the more accurate source: real transcription start sites, strand
#' resolved, one canonical position per gene. The Illumina manifest's RefGene
#' fields are convenient but give a gene \emph{region} label per probe rather
#' than a TSS, so distances derived from it are cruder.
#'
#' @examples
#' \donttest{
#' # needs TxDb.Hsapiens.UCSC.hg19.knownGene and org.Hs.eg.db
#' cpgd_build_gene_tss(out = tempfile(fileext = ".csv"))
#' }
#' @export
cpgd_build_gene_tss <- function(out = "gene_tss_hg19.csv",
                                source = c("txdb", "manifest"),
                                manifest_file = NULL) {
  source <- match.arg(source)

  if (source == "txdb") {
    need <- c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db",
              "GenomicFeatures", "AnnotationDbi")
    miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) {
      stop("Install these first:\n  BiocManager::install(c(",
           paste(sprintf('"%s"', miss), collapse = ", "), "))", call. = FALSE)
    }
    txdb <- get("TxDb.Hsapiens.UCSC.hg19.knownGene",
                envir = asNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    g <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
    sym <- AnnotationDbi::mapIds(get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db")),
                                 keys = names(g), keytype = "ENTREZID",
                                 column = "SYMBOL", multiVals = "first")
    st <- as.character(GenomicRanges::strand(g))
    tss <- ifelse(st == "+", GenomicRanges::start(g), GenomicRanges::end(g))
    d <- data.table::data.table(
      gene   = toupper(as.character(sym)),
      chr    = as.character(GenomicRanges::seqnames(g)),
      tss    = as.integer(tss),
      strand = st)
    d <- d[!is.na(d$gene) & grepl("^chr[0-9XY]+$", d$chr), ]
    d <- unique(d, by = "gene")
  } else {
    if (is.null(manifest_file) || !file.exists(manifest_file))
      stop("Supply manifest_file when source = \"manifest\".", call. = FALSE)
    m <- .cpgd_fread(manifest_file)
    nm <- names(m)
    gcol <- nm[grepl("RefGene_Name|UCSC_RefGene_Name", nm)][1]
    ccol <- nm[grepl("^CHR$|^chr$|CHR_hg38|Chromosome", nm)][1]
    pcol <- nm[grepl("MAPINFO|^pos$|Start_hg38", nm)][1]
    if (any(is.na(c(gcol, ccol, pcol))))
      stop("Could not find gene, chromosome and position columns in the manifest.",
           call. = FALSE)
    g <- strsplit(as.character(m[[gcol]]), ";", fixed = TRUE)
    d <- data.table::data.table(
      gene = toupper(vapply(g, function(x) if (length(x)) x[1] else NA_character_, character(1))),
      chr  = paste0("chr", gsub("^chr", "", as.character(m[[ccol]]))),
      tss  = as.integer(m[[pcol]]), strand = NA_character_)
    d <- d[!is.na(d$gene) & d$gene != "", ]
    d <- d[, list(tss = as.integer(stats::median(get("tss"), na.rm = TRUE)),
                  strand = NA_character_), by = c("gene", "chr")]
    message("Manifest-derived positions approximate a gene region, not a TSS. ",
            "TxDb is more accurate.")
  }
  data.table::fwrite(d, out)
  message(sprintf("wrote %s (%d genes)", out, nrow(d)))
  invisible(out)
}


#' Direction for any CpG-gene pair, from distance alone
#'
#' Catalogue-based prediction only covers genes some catalogue measured, which
#' is why genes such as CRH, MC2R or HTR2C return
#' \code{NO_EXPRESSION_ANCHOR} — they are not transcribed in blood, so no blood
#' catalogue contains them. This function drops that restriction and predicts
#' from distance alone, for any pair whose coordinates you can supply.
#'
#' It is deliberately weaker and says so. The distance-only relationship reaches
#' AUC 0.616 within cohort against 0.793 for the full blood model, because it
#' uses none of the neighbourhood or gene-level evidence. Everything it returns
#' carries \code{evidence_tier = "U"}.
#'
#' @param cpgs Character vector of CpG identifiers.
#' @param genes Target genes, same length as \code{cpgs}. Optional: if omitted,
#'   the symbol is taken from the identifier where one is present, as in
#'   \code{"cg01961214_TC21_CRHR1"}.
#' @param gene_tss Optional. Defaults to \code{\link{cpgd_gene_tss}}, which
#'   builds the annotation from TxDb on first use and caches it. Pass a path or a
#'   data.frame with \code{gene}, \code{chr}, \code{tss} to override.
#' @param cpg_pos Optional data.frame with \code{cpg_id}, \code{chr}, \code{pos}
#'   in hg19. Defaults to the packaged EPIC v2 positions, 930,178 probes.
#' @param tissue Which tissue's curve to use, or \code{"all"} for one column per
#'   tissue plus a consensus.
#' @param max_dist Pairs further apart than this are not called. Default 1e6.
#' @param verbose Print a summary.
#'
#' @return A \code{data.table} with the CpG, gene, distance, and a direction and
#'   probability per tissue.
#'
#' @section What this is and is not:
#' The estimate is the empirical fraction of positive associations at that
#' distance in that tissue class, isotonically smoothed and reported with its
#' 95\% interval. It is a prior conditioned on distance, not a statement about
#' your specific CpG, and like everything else here it is conditional on the CpG
#' affecting expression at all. Use it to orient a hypothesis, not to settle one.
#'
#' @examples
#' \donttest{
#' # builds (and caches) the TxDb-derived TSS table on first use
#' cpg_direction_universal(c("cg03405789", "cg08215831"),
#'                         genes = c("CRH", "CRH"))
#' }
#' @export
cpg_direction_universal <- function(cpgs, genes, gene_tss = NULL,
                                    cpg_pos = NULL,
                                    tissue = "all",
                                    max_dist = 1e6, verbose = TRUE) {

  G <- if (is.null(gene_tss)) cpgd_gene_tss()
       else if (is.character(gene_tss) && length(gene_tss) == 1L)
         .cpgd_fread(gene_tss)
       else data.table::as.data.table(gene_tss)
  if (!all(c("gene", "chr", "tss") %in% names(G)))
    stop("gene_tss needs columns: gene, chr, tss.", call. = FALSE)
  G[, "gene" := toupper(get("gene"))]
  G <- unique(G, by = "gene")

  q <- .cpgd_parse_input(cpgs, genes)
  if (!nrow(q)) stop("No valid CpG identifiers found.", call. = FALSE)

  # Distance needs a target gene, but it does not have to be supplied
  # explicitly. Identifiers such as "Zcg01961214_TC21_LINC02210_CRHR1" already
  # carry it, and refusing to look would make this function the only one in the
  # package that ignores what the parser extracted. Prefer the full parsed
  # label, fall back to its last token, which is where the informative symbol
  # sits in Illumina-style compound annotations.
  known <- G$gene
  q[, "gene_use" := data.table::fifelse(
        toupper(get("given_gene")) %in% known, toupper(get("given_gene")),
        data.table::fifelse(toupper(get("given_gene2")) %in% known,
                            toupper(get("given_gene2")), NA_character_))]
  n_nogene <- sum(is.na(q$gene_use))
  if (n_nogene == nrow(q)) {
    stop("No target gene could be resolved for any input. Supply `genes`, or ",
         "use identifiers that carry the symbol (e.g. \"cg01961214_TC21_CRHR1\"). ",
         "cpgd_parse_check() shows what the parser extracted.", call. = FALSE)
  }
  if (n_nogene > 0L) {
    message(sprintf("%d of %d inputs had no gene symbol matching the annotation; dropped.",
                    n_nogene, nrow(q)))
  }
  q <- q[!is.na(get("gene_use")), ]

  if (is.null(cpg_pos)) cpg_pos <- cpgd_cpg_positions()
  P <- data.table::as.data.table(cpg_pos)
  names(P) <- tolower(names(P))
  if (!all(c("cpg_id", "chr", "pos") %in% names(P)))
    stop("cpg_pos needs columns: cpg_id, chr, pos.", call. = FALSE)
  P[, "cpg_id" := tolower(get("cpg_id"))]

  d <- merge(q[, c("cpg_id", "input", "gene_use"), with = FALSE], P, by = "cpg_id")
  d <- merge(d, G[, c("gene", "chr", "tss"), with = FALSE],
             by.x = "gene_use", by.y = "gene", suffixes = c("", "_gene"))
  if (!nrow(d)) {
    message("No CpG-gene pair could be resolved: check that gene symbols match ",
            "the annotation and that positions are hg19.")
    return(invisible(data.table::data.table()))
  }
  d <- d[get("chr") == get("chr_gene"), ]
  d[, "abs_dist" := abs(as.numeric(get("pos")) - as.numeric(get("tss")))]
  d <- d[get("abs_dist") <= max_dist, ]

  K <- data.table::fread(system.file("extdata", "distance_curves.csv",
                                     package = "cpgdirection"), showProgress = FALSE)
  tis <- if (identical(tissue, "all")) CPGD_TISSUES else match.arg(tissue, CPGD_TISSUES)
  short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")
  for (t in tis) {
    k <- K[get("tissue") == t, ]
    p <- stats::approx(log10(pmax(k$dist_mid, 1)), k$p_positive,
                       xout = log10(pmax(d$abs_dist, 1)), rule = 2)$y
    lo <- stats::approx(log10(pmax(k$dist_mid, 1)), k$ci_lo,
                        xout = log10(pmax(d$abs_dist, 1)), rule = 2)$y
    hi <- stats::approx(log10(pmax(k$dist_mid, 1)), k$ci_hi,
                        xout = log10(pmax(d$abs_dist, 1)), rule = 2)$y
    s <- short[[t]]
    # set(), not `[[<-`: the latter is base-R assignment, which copies the
    # data.table and marks it, so the `:=` calls further down then emit a
    # shallow-copy warning that has nothing to do with the caller's code.
    data.table::set(d, j = paste0("p_", s),   value = round(p, 4))
    data.table::set(d, j = paste0("dir_", s), value = ifelse(p >= .5, 1, -1))
    data.table::set(d, j = paste0("ci_", s),
                    value = paste0("[", round(lo, 2), ",", round(hi, 2), "]"))
  }
  if (length(tis) > 1) {
    D <- as.matrix(d[, paste0("dir_", short[tis]), with = FALSE])
    n_pos <- rowSums(D == 1); n_neg <- rowSums(D == -1)
    d[, "consensus_direction" := ifelse(n_pos == ncol(D), 1,
                                 ifelse(n_neg == ncol(D), -1, NA_real_))]
    d[, "tissue_agreement" := round(pmax(n_pos, n_neg) / ncol(D), 3)]
  }
  d[, "evidence_tier" := "U"]
  d[, "expected_accuracy" := "0.60-0.65 (distance only; weaker than tier A/B)"]
  data.table::setnames(d, "gene_use", "target_gene")
  keep <- c("input", "cpg_id", "target_gene", "abs_dist",
            unlist(lapply(short[tis], function(s) paste0(c("dir_", "p_", "ci_"), s))),
            "consensus_direction", "tissue_agreement", "evidence_tier",
            "expected_accuracy")
  d <- d[, intersect(keep, names(d)), with = FALSE]
  if (isTRUE(verbose)) {
    cat(sprintf("\nuniversal mode: %d CpG-gene pairs resolved from distance\n", nrow(d)))
    cat("  tier U - distance only, no neighbourhood or gene evidence.\n")
    cat("  AUC 0.616 within cohort, against 0.793 for the full blood model.\n")
    cat("  A prior conditioned on distance, not a statement about your CpG.\n\n")
    print(utils::head(as.data.frame(d), 10))
  }
  invisible(d)
}


#' Packaged CpG positions (EPIC v2, hg19)
#'
#' 930,178 probe positions, lifted to hg19 and verified at 99.8 percent exact
#' agreement against the HELIX catalogue. Used by
#' \code{\link{cpg_direction_universal}}; exported because it is generally
#' useful.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{chr}, \code{pos}.
#' @examples
#' p <- cpgd_cpg_positions()
#' p[p$cpg_id == "cg00000029", ]
#' @export
cpgd_cpg_positions <- function() {
  if (!is.null(.cpgd_env$positions)) return(.cpgd_env$positions)
  P <- .cpgd_data("cpg_positions_hg19")
  data.table::setkeyv(P, "cpg_id")
  .cpgd_env$positions <- P
  P
}


#' Gene TSS annotation, built once and cached
#'
#' Derives transcription start sites from
#' \code{TxDb.Hsapiens.UCSC.hg19.knownGene}, mapping Entrez identifiers to
#' symbols through \code{org.Hs.eg.db}. The result is cached under
#' \code{tools::R_user_dir("cpgdirection", "cache")}, so the cost is paid once
#' per machine rather than once per session.
#'
#' Called automatically by \code{\link{cpg_direction_universal}}. Call it
#' directly only to force a rebuild or to inspect the table.
#'
#' @param refresh Rebuild even if a cached copy exists.
#' @return A \code{data.table} with \code{gene}, \code{chr}, \code{tss},
#'   \code{strand}.
#' @examples \donttest{g <- cpgd_gene_tss(); g[gene == "CRH"]}
#' @export
cpgd_gene_tss <- function(refresh = FALSE) {
  if (!refresh && !is.null(.cpgd_env$gene_tss)) return(.cpgd_env$gene_tss)
  dir <- tools::R_user_dir("cpgdirection", "cache")
  f <- file.path(dir, "gene_tss_hg19.csv.gz")
  if (!refresh && file.exists(f)) {
    G <- .cpgd_fread(f)
    .cpgd_env$gene_tss <- G
    return(G)
  }
  message("Building the gene TSS annotation from TxDb (once per machine)...")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".csv")
  cpgd_build_gene_tss(out = tmp, source = "txdb")
  G <- data.table::fread(tmp, showProgress = FALSE)
  data.table::fwrite(G, f)
  message(sprintf("cached %d genes at %s", nrow(G), f))
  .cpgd_env$gene_tss <- G
  G
}
