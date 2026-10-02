# Distal CpG-to-gene links and tumour eQTM consensus. New in 2.99.6.
#
# Both tables are kept OUT of the default pair workflow. A distal link proposes
# a target gene on regulatory-element evidence (predicted enhancer-gene
# interaction, single-cell peak-gene link, chromatin loop) and carries no sign;
# the tumour consensus carries a sign measured in TCGA tumours, which is not
# the tissue any Project Alpha question is about. Each is reported alongside
# the peripheral record, as brain SMR is, and neither enters best_direction.

#' Distal CpG-to-gene links from enhancer-gene maps and chromatin loops
#'
#' The packaged \code{distal_links} layer: EPIC v2 CpGs linked to target genes
#' by evidence that is better than proximity, built by
#' \code{tools/build_distal_links.R} from three sources, each intersected with
#' the probe's hg38 coordinate and keyed on the cg identifier:
#' \itemize{
#'   \item \code{rE2G} -- ENCODE-rE2G enhancer-gene predictions (Gschwind et
#'     al. 2026, Nature; 1,458 biosamples; per-link logistic score, binarised
#'     by the authors at 70\% CRISPR recall);
#'   \item \code{scE2G} -- single-cell multiome enhancer-gene maps (Sheth et
#'     al. 2026, Nature Genetics; PBMC and bone-marrow cell clusters, islets);
#'   \item \code{HCEA} -- cell-type chromatin loops from the Human Cell
#'     Epigenome Atlas (Zhou et al. 2026, Science; snm3C-seq, 16 tissues
#'     including peripheral blood and cortex): the CpG sits in one anchor and
#'     the gene's TSS in the other.
#' }
#' Every row is one CpG x gene x source x biosample link with the source's own
#' score. The table is a TARGET-proposal resource. It says "an enhancer or loop
#' containing this CpG is predicted or observed to contact this gene in this
#' biosample"; it does not say in which direction methylation moves
#' expression, and it is deliberately not counted in
#' \code{n_annotation_sources}, which tallies positional annotation tracks.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{source}, \code{biosample}, \code{tissue_class} (\code{blood},
#'   \code{brain} or \code{other}), \code{score}, and the element coordinates
#'   \code{element_chr}, \code{element_start}, \code{element_end} (hg38,
#'   provenance only) and \code{element_class}. \code{NULL}, with a message,
#'   when the layer is not installed.
#' @examplesIf cpgd_has_data("distal_links")
#' d <- cpgd_distal_links()
#' d[d$cpg_id == "cg26261055", ]
#' @export
cpgd_distal_links <- function() {
  if (!is.null(.cpgd_env$distal_links)) return(.cpgd_env$distal_links)
  D <- .cpgd_data("distal_links", required = FALSE)
  if (is.null(D)) {
    message("distal_links layer not available in this installation.")
    return(NULL)
  }
  req <- c("cpg_id", "target_gene", "source", "biosample", "score")
  miss <- setdiff(req, names(D))
  if (length(miss)) {
    stop("The distal_links table is missing expected columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!"tissue_class" %in% names(D)) D[, "tissue_class" := "other"]
  data.table::setkeyv(D, c("cpg_id", "target_gene"))
  .cpgd_env$distal_links <- D
  D
}


#' Tumour eQTM consensus across 27 TCGA cancer types (Onco-eQTM)
#'
#' The packaged \code{onco_eqtm_consensus} layer, built by
#' \code{tools/build_onco_eqtm.R} from the Onco-eQTM per-cancer cis-eQTM tables
#' (Korra, Nishana and Kumar 2026, NAR Genomics and Bioinformatics; 6,880 TCGA
#' primary tumours, Illumina 450K, hg19; consensus of Torch-eCpG and
#' MatrixEQTL at FDR < 0.05 and |r| > 0.3). One row per CpG x gene, pooled
#' over the cancer types in which the pair was significant.
#'
#' \code{onco_tier} grades the cross-tumour agreement: \code{O1} = significant
#' in >= 5 cancer types with >= 90\% sign agreement; \code{O2} = 2-4 types,
#' >= 90\% agreement; \code{O3} = >= 2 types but the sign disagrees;
#' \code{O4} = a single cancer type. \code{direction} is the consensus sign
#' for O1/O2 and \code{NA} otherwise.
#'
#' This is tumour evidence. It never enters \code{best_direction}: the
#' peripheral ladder answers a blood, nasal or solid-tissue question, and a
#' sign measured across tumours is reported beside it (see
#' \code{include_onco} in \code{\link{cpg_gene_pairs}}) so that agreement or
#' disagreement with the peripheral call is visible. Cross-tumour agreement
#' (O1) is strong evidence that a CpG-gene coupling is constitutive rather
#' than context-specific, which is useful; it is not evidence about saliva.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{direction}, \code{onco_tier}, \code{n_cancers}, \code{n_positive},
#'   \code{n_negative}, \code{sign_agreement}, \code{median_r},
#'   \code{max_abs_r}, \code{region}, \code{tss_dist_approx}, \code{cancers}
#'   and \code{source}. \code{NULL}, with a message, when the layer is not
#'   installed.
#' @examplesIf cpgd_has_data("onco_eqtm_consensus")
#' o <- cpgd_onco_eqtm()
#' o[o$onco_tier == "O1", ][1:5, ]
#' @export
cpgd_onco_eqtm <- function() {
  if (!is.null(.cpgd_env$onco_eqtm)) return(.cpgd_env$onco_eqtm)
  O <- .cpgd_data("onco_eqtm_consensus", required = FALSE)
  if (is.null(O)) {
    message("onco_eqtm_consensus layer not available in this installation.")
    return(NULL)
  }
  req <- c("cpg_id", "target_gene", "direction", "onco_tier", "n_cancers",
           "sign_agreement")
  miss <- setdiff(req, names(O))
  if (length(miss)) {
    stop("The onco_eqtm_consensus table is missing expected columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  data.table::setkeyv(O, c("cpg_id", "target_gene"))
  .cpgd_env$onco_eqtm <- O
  O
}


# Collapse the long distal table to one row per CpG x gene with the summary
# columns the pair workflow reports. `tissue_class` restricts to biosamples of
# that class ("blood"/"brain"/"other"); NULL keeps all.
.cpgd_distal_summary <- function(D, cpgs, tissue_class = NULL) {
  # data.table resolves bare names as columns first: the argument is copied
  # to a name no column can have before it is used inside the frame.
  .want_cpgs <- cpgs
  .want_class <- tissue_class
  DD <- data.table::as.data.table(D)[get("cpg_id") %chin% .want_cpgs]
  if (!is.null(.want_class) && "tissue_class" %in% names(DD)) {
    DD <- DD[get("tissue_class") %chin% .want_class]
  }
  if (!nrow(DD)) {
    return(data.table::data.table(cpg_id = character(0), .gkey = character(0),
                                  target_gene = character(0),
                                  distal_sources = character(0),
                                  distal_n_biosamples = integer(0),
                                  distal_score_max = numeric(0),
                                  distal_biosamples = character(0)))
  }
  DD[, ".gkey" := gsub("_", "-", toupper(get("target_gene")))]
  DD[, list(
    target_gene         = get("target_gene")[1L],
    distal_sources      = paste(sort(unique(get("source"))), collapse = ";"),
    distal_n_biosamples = data.table::uniqueN(paste(get("source"), get("biosample"))),
    distal_score_max    = suppressWarnings(max(as.numeric(get("score")), na.rm = TRUE)),
    distal_biosamples   = paste(sort(unique(get("biosample"))), collapse = ";")),
    by = c("cpg_id", ".gkey")][
    , "distal_score_max" := data.table::fifelse(is.finite(get("distal_score_max")),
                                                 get("distal_score_max"), NA_real_)][]
}
