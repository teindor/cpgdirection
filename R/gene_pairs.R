# Pair-preserving API (2.4.0): automatic CpG -> gene target discovery and one
# direction/evidence record per CpG x gene pair.
#
# The package separates two questions that the one-row-per-CpG interfaces
# partly couple:
#   A. target discovery -- which gene or genes are plausible targets of this
#      CpG?
#   B. direction -- for a specific CpG x gene pair, what direction of
#      methylation-expression association is supported?
# cpg_gene_pairs() answers A by taking the UNION of every configured mapping
# source, then answers B separately for every resulting pair via
# .cpgd_resolve_pair_direction(). Nothing here collapses by cpg_id alone: the
# biological key is cpg_id + target_gene.

# priority used for choosing a pair's canonical spelling and mapping_primary;
# smaller is stronger biological meaning
.CPGD_MAPPING_PRIORITY <- c(
  measured_eQTM    = 1L,
  SMR              = 2L,
  tissue_lookup    = 3L,
  EPICv2_manifest  = 4L,
  distal_link      = 5L,
  brain_SMR        = 6L,
  input_annotation = 7L,
  requested_pair   = 8L)


#' All supported target genes and pair-specific directions for a set of CpGs
#'
#' The pair-preserving interface. Give it bare CpG identifiers and it first
#' discovers every plausible target gene automatically -- from the packaged
#' EPIC v2 manifest annotation, the tissue catalogue, directly measured eQTMs
#' and SMR (plus the brain layer on request, and any gene parsed from the
#' input) -- and then resolves direction and evidence separately for every
#' resulting CpG x gene pair. One row per supported pair.
#'
#' \strong{You do not have to supply genes.} \code{genes = NULL} means
#' auto-discover all supported targets; it never means "cannot calculate
#' direction" and it never forces a nearest-gene fallback. Supplying
#' \code{genes} is an optional restriction (\code{gene_mode = "filter"}) or an
#' explicit pair request (\code{gene_mode = "pairwise"}), never a requirement.
#'
#' @param cpgs Character vector of CpG identifiers, with or without suffixes
#'   (\code{"cg26261055"}, \code{"cg00645497_BC21_T1_CD38"}), a data.frame, or
#'   a file path -- the same input forms as
#'   \code{\link{cpg_expression_direction}}. Every input row is preserved:
#'   nothing is deduplicated by \code{cpg_id} alone.
#' @param genes Optional character vector of gene symbols. With the default
#'   \code{gene_mode = "filter"} it is a set: targets are auto-discovered
#'   first and then restricted to genes in this set (any length). With
#'   \code{gene_mode = "pairwise"} it must be the same length as \code{cpgs}
#'   and \code{cpgs[i]} x \code{genes[i]} are evaluated as explicit pairs.
#' @param tissue Tissue context for the catalogue mapping layer and the
#'   \code{lookup_*} audit columns. Direction resolution consults all three
#'   packaged catalogues for consensus, exactly as
#'   \code{cpg_expression_direction()} does.
#' @param gene_mode \code{"filter"} (default) or \code{"pairwise"}; see
#'   \code{genes}. The mode is never inferred from vector lengths -- state
#'   \code{"pairwise"} explicitly when you mean it.
#' @param annotation_mode How a gene parsed from an input suffix is used.
#'   \code{"augment"} (default): it is one more candidate target alongside the
#'   auto-discovered ones. \code{"strict"}: inputs that carry a suffix gene
#'   are restricted to their annotated gene(s); bare inputs still
#'   auto-discover.
#' @param include Which mapping sources contribute candidate targets. Any
#'   subset of \code{c("manifest", "lookup", "measured", "smr", "distal")};
#'   default the first four. \code{"distal"} is opt-in: it adds targets from
#'   the \code{\link{cpgd_distal_links}} layer (ENCODE-rE2G / scE2G
#'   enhancer-gene predictions and Human Cell Epigenome Atlas loops), which
#'   propose a gene on regulatory-element evidence rather than proximity. A
#'   distal link is a TARGET proposal with no sign; the pair's direction is
#'   still resolved by the evidence ladder, and the link's provenance is
#'   reported in \code{distal_sources}, \code{distal_n_biosamples},
#'   \code{distal_score_max} and \code{distal_biosamples}. The biosamples
#'   consulted follow \code{tissue}: blood-class biosamples for
#'   \code{tissue = "blood"}, all biosamples otherwise.
#' @param probe_qc How to treat CpGs whose probes are flagged as unreliable by
#'   the packaged QC table (\code{\link{cpgd_probe_qc}}: cross-hybridizing or
#'   degenerate mapping, SNP-contaminated extension base, unverified hg19
#'   position). \code{"exclude"} (default) drops CpGs whose \emph{every}
#'   replicate probe carries the general mask -- a probe that reads the wrong
#'   locus assigns methylation to the wrong gene, and no direction evidence
#'   can repair that. \code{"flag"} keeps them, marked in
#'   \code{probe_masked}/\code{probe_mask_reasons}. \code{"ignore"} skips QC
#'   entirely. Exclusions are reported and listed in
#'   \code{attr(result, "qc_excluded_cpgs")}.
#' @param include_brain Add brain SMR targets (\code{\link{cpgd_brain_directions}})
#'   to candidate discovery and report \code{brain_direction}/\code{brain_tier}
#'   per pair. Default \code{FALSE}: brain evidence is context-specific and is
#'   never mixed into a default peripheral run. Even when included, brain
#'   directions are reported alongside and never enter \code{best_direction},
#'   which is peripheral; use \code{\link{cpg_brain_bridge}} for the brain
#'   question.
#' @param include_onco Add the tumour eQTM consensus
#'   (\code{\link{cpgd_onco_eqtm}}; Onco-eQTM, 27 TCGA cancer types) as a
#'   candidate-discovery source and report \code{onco_direction},
#'   \code{onco_tier}, \code{onco_n_cancers} and \code{onco_agreement} per
#'   pair. Default \code{FALSE}. As with brain, tumour directions are reported
#'   alongside and never enter \code{best_direction}; \code{onco_agreement}
#'   says whether the tumour sign matches the peripheral call.
#' @param direction_policy \code{"best"} (default) returns the pair-level
#'   best-evidence record with the standard audit columns.
#'   \code{"all_evidence"} additionally keeps every per-tissue catalogue and
#'   distance column.
#' @param universal Include the distance layer in pair-level resolution.
#'   Default \code{TRUE}; requires the Bioconductor annotation packages on
#'   first use, and degrades gracefully when unavailable.
#' @param sources Advanced/testing: a named list overriding packaged
#'   resources -- any of \code{manifest}, \code{lookup} (named list by tissue,
#'   or one table for \code{tissue}), \code{measured}, \code{smr},
#'   \code{brain}, \code{cpg_pos}, \code{gene_tss}, \code{distance_curves}.
#' @param verbose Print the discovery summary. Default \code{TRUE}.
#'
#' @return A \code{data.table} of class \code{cpgd_pairs}, one row per
#'   CpG x gene pair, with (among others):
#'   \describe{
#'     \item{input_id, input}{which submitted row(s) the CpG came from
#'       (\code{input_ids}/\code{inputs} collapse several)}
#'     \item{cpg_id, target_gene, pair_id}{the biological key}
#'     \item{mapping_sources, mapping_primary, mapping_strength}{how the gene
#'       entered the candidate set: \code{EPICv2_manifest},
#'       \code{tissue_lookup}, \code{measured_eQTM}, \code{SMR},
#'       \code{distal_link}, \code{brain_SMR}, \code{tumour_eQTM},
#'       \code{input_annotation}, \code{requested_pair} --
#'       plus \code{has_*} flags. Mapping provenance is deliberately distinct
#'       from direction evidence: manifest annotation says "Illumina annotates
#'       this CpG to this gene", not "this CpG regulates this gene"}
#'     \item{best_direction, best_evidence, best_confidence, direction_tier,
#'       probability_plus1}{the pair's direction record, resolved by the same
#'       evidence ladder as \code{cpg_expression_direction()}, strictly within
#'       \code{cpg_id + target_gene}}
#'     \item{usable, abstain_reason}{pair-specific: one gene of a CpG can be
#'       usable while another abstains}
#'     \item{n_targets_for_cpg, is_coeffect}{co-effect accounting}
#'   }
#'
#' @section Pair invariance:
#' For a fixed pair the direction record does not depend on how the gene
#' entered the candidate set. \code{cpg_gene_pairs("cg26261055")} (manifest
#' discovers CRHBP), \code{cpg_gene_pairs("cg26261055_CRHBP")} (suffix
#' supplies it) and \code{cpg_gene_pairs("cg26261055", genes = "CRHBP",
#' gene_mode = "pairwise")} return the same \code{best_direction},
#' \code{best_evidence} and confidence for cg26261055 x CRHBP; only the
#' mapping provenance differs.
#'
#' @section Co-effects:
#' A CpG may legitimately have opposite directions for different target genes;
#' both rows are retained and there is no single CpG-level direction.
#' Abstention is likewise pair-specific.
#'
#' @examples
#' # hermetic demonstration with injected fixture sources (fast, no data):
#' fx <- data.frame(cpg_id = "cg10000001", target_gene = "GENEA",
#'                  array = "EPICv2", annotation_source = "Illumina_EPICv2",
#'                  refgene_group = "TSS200")
#' p <- cpg_gene_pairs("cg10000001", universal = FALSE, verbose = FALSE,
#'                     sources = list(manifest = fx, lookup = NULL,
#'                                    measured = NULL, smr = NULL,
#'                                    probe_qc = NULL))
#' p[, c("cpg_id", "target_gene", "mapping_primary", "best_evidence"),
#'   with = FALSE]
#'
#' if (cpgd_has_data("lookup_blood_hg19")) {
#' # real resources: bare CpGs, targets discovered automatically
#' p <- cpg_gene_pairs(c("cg26261055", "cg06495038", "cg26196496"))
#'
#' # DMSA-style use: discovered targets restricted to a gene set
#' p <- cpg_gene_pairs(c("cg26261055", "cg06495038"),
#'                     genes = c("CRHBP", "F2RL1"),
#'                     gene_mode = "filter")
#' }
#' @export
cpg_gene_pairs <- function(cpgs,
                           genes = NULL,
                           tissue = c("blood", "nasal_epithelium", "solid_tissue"),
                           gene_mode = c("filter", "pairwise"),
                           annotation_mode = c("augment", "strict"),
                           include = c("manifest", "lookup", "measured", "smr"),
                           include_brain = FALSE,
                           include_onco = FALSE,
                           probe_qc = c("exclude", "flag", "ignore"),
                           direction_policy = c("best", "all_evidence"),
                           universal = TRUE,
                           sources = NULL,
                           verbose = TRUE) {

  tissue <- match.arg(tissue, CPGD_TISSUES)
  gene_mode <- match.arg(gene_mode)
  annotation_mode <- match.arg(annotation_mode)
  probe_qc <- match.arg(probe_qc)
  direction_policy <- match.arg(direction_policy)
  include <- match.arg(include,
                       c("manifest", "lookup", "measured", "smr", "distal"),
                       several.ok = TRUE)

  if (gene_mode == "pairwise" && is.null(genes)) {
    stop("gene_mode = \"pairwise\" needs `genes`, one per element of `cpgs`.",
         call. = FALSE)
  }

  # Every parsed row survives: the pair workflow never deduplicates by cpg_id
  # alone (that would fold cg123_GENE_A and cg123_GENE_B into one row).
  q <- .cpgd_parse_input(cpgs,
                         genes = if (gene_mode == "pairwise") genes else NULL,
                         dedupe = FALSE)
  if (!nrow(q)) {
    stop("No valid CpG identifiers found in the input. Identifiers must ",
         "contain 'cg' followed by at least six digits. ",
         "Use cpgd_parse_check() to see what the parser extracted.",
         call. = FALSE)
  }
  n_submitted <- attr(q, "n_raw")
  u <- unique(q$cpg_id)

  # ---- probe QC: a probe that reads the wrong locus assigns methylation to
  # the wrong gene entirely, a failure no downstream direction evidence can
  # repair. By default CpGs whose EVERY replicate probe carries the general
  # mask are excluded; probe_qc = "flag" keeps them marked instead. ---------
  QC <- NULL
  qc_excluded <- character(0)
  if (probe_qc != "ignore") {
    QC <- .cpgd_source_table(sources, "probe_qc", cpgd_probe_qc)
    if (!is.null(QC) && nrow(QC)) {
      # qc_exclude (present from 2.5.0 tables that integrate the Garvan
      # cross-hybridization evidence) is the combined exclusion criterion;
      # older tables fall back to the Zhou mask alone
      bad <- if ("qc_exclude" %in% names(QC)) {
        unique(QC$cpg_id[QC$qc_exclude %in% TRUE])
      } else {
        unique(QC$cpg_id[QC$masked_general %in% TRUE])
      }
      if (probe_qc == "exclude") {
        qc_excluded <- intersect(u, bad)
        if (length(qc_excluded)) {
          u <- setdiff(u, qc_excluded)
          q <- q[!q$cpg_id %chin% qc_excluded, ]
          message(length(qc_excluded), " CpG(s) excluded as unreliable ",
                  "probes (all replicates flagged by the probe QC mask: ",
                  "cross-hybridizing/degenerate mapping/SNP-contaminated). ",
                  "Use probe_qc = \"flag\" to keep them, marked.")
        }
        if (!nrow(q)) {
          warning("Every submitted CpG was excluded by probe QC. ",
                  "Re-run with probe_qc = \"flag\" to inspect them.",
                  call. = FALSE)
          empty <- data.table::data.table(cpg_id = character(0),
                                          target_gene = character(0))
          data.table::setattr(empty, "qc_excluded_cpgs", qc_excluded)
          data.table::setattr(empty, "class", c("cpgd_pairs", class(empty)))
          return(invisible(empty))
        }
      }
    } else if (verbose) {
      message("probe QC table unavailable; probe reliability not checked.")
    }
  }

  # ---- candidate discovery: the union of every configured mapping source --
  cand <- list()
  src_counts <- c(manifest = NA_integer_, lookup = NA_integer_,
                  measured = NA_integer_, smr = NA_integer_,
                  distal = NA_integer_, brain = NA_integer_,
                  onco = NA_integer_, input = NA_integer_)
  # gene symbols known to any loaded source, dash-normalised; used to decide
  # whether a compound parsed tail ("LINC02210_CRHR1") is itself a gene or a
  # parsing artefact
  known_gkeys <- character(0)

  T_man <- if ("manifest" %in% include)
    .cpgd_source_table(sources, "manifest", cpgd_manifest_genes) else NULL
  if ("manifest" %in% include) {
    if (!is.null(T_man) && nrow(T_man)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(T_man$target_gene))))
      mm <- T_man[get("cpg_id") %chin% u,
                  c("cpg_id", "target_gene"), with = FALSE]
      src_counts["manifest"] <- nrow(mm)
      if (nrow(mm)) cand$manifest <- mm[, "mapping_source" := "EPICv2_manifest"]
    } else if (verbose) {
      message("EPIC v2 manifest annotation unavailable; ",
              "manifest targets not discovered.")
    }
  }
  if ("lookup" %in% include) {
    L <- .cpgd_source_lookup(sources, tissue, selected = tissue)
    if (!is.null(L) && nrow(L)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(L$target_gene))))
      ll <- L[get("cpg_id") %chin% u & !is.na(get("target_gene")),
              c("cpg_id", "target_gene"), with = FALSE]
      src_counts["lookup"] <- nrow(unique(ll))
      if (nrow(ll)) cand$lookup <- ll[, "mapping_source" := "tissue_lookup"]
    }
  }
  if ("measured" %in% include) {
    Me <- .cpgd_source_table(sources, "measured", cpgd_measured_eqtms)
    if (!is.null(Me) && nrow(Me)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(Me$target_gene))))
      mm <- Me[get("cpg_id") %chin% u,
               c("cpg_id", "target_gene"), with = FALSE]
      src_counts["measured"] <- nrow(unique(mm))
      if (nrow(mm)) cand$measured <- mm[, "mapping_source" := "measured_eQTM"]
    }
  }
  if ("smr" %in% include) {
    S <- .cpgd_source_table(sources, "smr", cpgd_smr_directions)
    if (!is.null(S) && nrow(S)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(S$target_gene))))
      ss <- S[get("cpg_id") %chin% u,
              c("cpg_id", "target_gene"), with = FALSE]
      src_counts["smr"] <- nrow(unique(ss))
      if (nrow(ss)) cand$smr <- ss[, "mapping_source" := "SMR"]
    }
  }
  if (isTRUE(include_brain)) {
    B <- .cpgd_source_table(sources, "brain", cpgd_brain_directions)
    if (!is.null(B) && nrow(B)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(B$target_gene))))
      bb <- B[get("cpg_id") %chin% u,
              c("cpg_id", "target_gene"), with = FALSE]
      src_counts["brain"] <- nrow(unique(bb))
      if (nrow(bb)) cand$brain <- bb[, "mapping_source" := "brain_SMR"]
    }
  }
  # Distal links: a regulatory element containing the CpG is predicted or
  # observed to contact the gene. Opt-in, and consulted per tissue class so a
  # blood question is not answered with an islet enhancer map. These are
  # target proposals, not annotation: they never raise n_annotation_sources.
  dist_summary <- NULL
  if ("distal" %in% include) {
    D <- .cpgd_source_table(sources, "distal", cpgd_distal_links)
    if (!is.null(D) && nrow(D)) {
      tc <- if (identical(tissue, "blood")) "blood" else NULL
      dist_summary <- .cpgd_distal_summary(D, u, tissue_class = tc)
      if (!nrow(dist_summary) && !is.null(tc)) {
        # no blood-class biosample in the layer: fall back to everything,
        # and say so in the biosample column rather than silently
        dist_summary <- .cpgd_distal_summary(D, u, tissue_class = NULL)
      }
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(D$target_gene))))
      dd <- dist_summary[, c("cpg_id", "target_gene"), with = FALSE]
      src_counts["distal"] <- nrow(unique(dd))
      if (nrow(dd)) cand$distal <- dd[, "mapping_source" := "distal_link"]
    } else if (verbose) {
      message("distal_links layer unavailable; distal targets not discovered.")
    }
  }
  if (isTRUE(include_onco)) {
    O <- .cpgd_source_table(sources, "onco", cpgd_onco_eqtm)
    if (!is.null(O) && nrow(O)) {
      known_gkeys <- c(known_gkeys,
                       gsub("_", "-", toupper(unique(O$target_gene))))
      oo <- O[get("cpg_id") %chin% u,
              c("cpg_id", "target_gene"), with = FALSE]
      src_counts["onco"] <- nrow(unique(oo))
      if (nrow(oo)) cand$onco <- oo[, "mapping_source" := "tumour_eQTM"]
    }
  }
  known_gkeys <- unique(known_gkeys)

  # input/parsed annotation: provenance, not an exclusive target
  ig <- .cpgd_input_gene_pairs(q, known_gkeys = known_gkeys,
                               pairwise = gene_mode == "pairwise")
  if (!is.null(ig) && nrow(ig)) {
    src_counts["input"] <- nrow(unique(ig[, c("cpg_id", "target_gene"), with = FALSE]))
    cand$input <- ig[, "mapping_source" :=
                       if (gene_mode == "pairwise") "requested_pair"
                       else "input_annotation"]
  }

  cand <- data.table::rbindlist(cand, use.names = TRUE, fill = TRUE)
  if (!nrow(cand)) {
    warning("No target gene was found from the EPIC-v2 manifest, packaged ",
            "catalogue, measured eQTM, SMR, or other enabled mapping ",
            "resources for any submitted CpG. You may supply an explicit ",
            "gene hypothesis with gene_mode = \"pairwise\".", call. = FALSE)
    empty <- data.table::data.table(cpg_id = character(0),
                                    target_gene = character(0))
    data.table::setattr(empty, "class", c("cpgd_pairs", class(empty)))
    return(invisible(empty))
  }
  cand[, "target_gene" := toupper(trimws(as.character(get("target_gene"))))]
  cand <- cand[!is.na(get("target_gene")) & nzchar(get("target_gene"))]
  cand[, ".gkey" := gsub("_", "-", get("target_gene"))]

  # ---- pair table: canonical spelling + full mapping provenance -----------
  cand[, ".pri" := .CPGD_MAPPING_PRIORITY[get("mapping_source")]]
  data.table::setorderv(cand, c("cpg_id", ".gkey", ".pri"))
  pairs <- cand[, list(
    target_gene = get("target_gene")[1L],
    mapping_sources = paste(unique(get("mapping_source")), collapse = ";"),
    mapping_primary = get("mapping_source")[1L],
    has_manifest         = any(get("mapping_source") == "EPICv2_manifest"),
    has_lookup           = any(get("mapping_source") == "tissue_lookup"),
    has_measured_eqtm    = any(get("mapping_source") == "measured_eQTM"),
    has_smr              = any(get("mapping_source") == "SMR"),
    has_brain            = any(get("mapping_source") == "brain_SMR"),
    has_distal           = any(get("mapping_source") == "distal_link"),
    has_onco             = any(get("mapping_source") == "tumour_eQTM"),
    has_input_annotation = any(get("mapping_source") %in%
                                 c("input_annotation", "requested_pair"))),
    by = c("cpg_id", ".gkey")]
  pairs[, "mapping_strength" := data.table::fcase(
    get("has_measured_eqtm"), "measured",
    get("has_smr"), "genetic_instrument",
    get("has_lookup"), "catalogue_model",
    default = "annotation")]

  # Discovery is judged BEFORE any gene restriction: a CpG whose targets are
  # later removed by `genes=` was discovered, not unmapped, and must not be
  # reported as having no target.
  discovered_cpgs <- unique(pairs$cpg_id)

  # ---- annotation_mode = "strict": an input that names a gene is restricted
  # to its annotated gene(s); bare inputs still auto-discover ---------------
  if (annotation_mode == "strict" && !is.null(ig) && nrow(ig)) {
    allowed <- unique(ig[, c("cpg_id", ".gkey"), with = FALSE])
    annotated_cpgs <- unique(allowed$cpg_id)
    in_allowed <- paste0(pairs$cpg_id, "\r", pairs$.gkey) %chin%
      paste0(allowed$cpg_id, "\r", allowed$.gkey)
    pairs <- pairs[!(get("cpg_id") %chin% annotated_cpgs) | in_allowed]
  }

  # ---- gene_mode ----------------------------------------------------------
  if (gene_mode == "pairwise") {
    # exactly the requested pairs, in request order; provenance from the union
    reqp <- unique(ig[, c("cpg_id", ".gkey", "target_gene"), with = FALSE])
    pairs <- pairs[reqp[, c("cpg_id", ".gkey"), with = FALSE],
                   on = c("cpg_id", ".gkey"), nomatch = NULL]
  } else if (!is.null(genes)) {
    fset <- unique(gsub("_", "-", toupper(trimws(as.character(genes)))))
    pairs <- pairs[get(".gkey") %chin% fset]
  }

  if (!nrow(pairs)) {
    warning("Target discovery succeeded but no pair survived the gene ",
            "restriction. Check the supplied `genes` against the discovered ",
            "targets (run again with genes = NULL to see them).",
            call. = FALSE)
    empty <- data.table::data.table(cpg_id = character(0),
                                    target_gene = character(0))
    data.table::setattr(empty, "class", c("cpgd_pairs", class(empty)))
    return(invisible(empty))
  }

  # ---- direction, resolved WITHIN each pair -------------------------------
  res <- .cpgd_resolve_pair_direction(
    pairs[, c("cpg_id", "target_gene"), with = FALSE],
    tissue = tissue, sources = sources, universal = universal)
  res[, ".gkey" := gsub("_", "-", toupper(get("target_gene")))]
  out <- merge(pairs, res[, setdiff(names(res), "target_gene"), with = FALSE],
               by = c("cpg_id", ".gkey"), all.x = TRUE, sort = FALSE)

  # brain audit columns, only when explicitly requested
  if (isTRUE(include_brain)) {
    B <- .cpgd_source_table(sources, "brain", cpgd_brain_directions)
    if (!is.null(B) && nrow(B)) {
      BB <- data.table::as.data.table(B)
      BB <- BB[, list(cpg_id = get("cpg_id"),
                      .gkey = gsub("_", "-", toupper(get("target_gene"))),
                      brain_direction = as.numeric(get("direction")),
                      brain_tier = as.character(get("tier")))]
      BB <- unique(BB, by = c("cpg_id", ".gkey"))
      out <- merge(out, BB, by = c("cpg_id", ".gkey"), all.x = TRUE, sort = FALSE)
    } else {
      out[, c("brain_direction", "brain_tier") := list(NA_real_, NA_character_)]
    }
  }

  # distal provenance columns, only when distal discovery was on
  if ("distal" %in% include) {
    if (!is.null(dist_summary) && nrow(dist_summary)) {
      DS <- dist_summary[, setdiff(names(dist_summary), "target_gene"), with = FALSE]
      out <- merge(out, DS, by = c("cpg_id", ".gkey"), all.x = TRUE, sort = FALSE)
    } else {
      out[, c("distal_sources", "distal_n_biosamples", "distal_score_max",
              "distal_biosamples") :=
            list(NA_character_, NA_integer_, NA_real_, NA_character_)]
    }
  }

  # tumour audit columns, only when explicitly requested. onco_agreement
  # compares the tumour sign with the peripheral best_direction of the SAME
  # pair; NA where either is missing.
  if (isTRUE(include_onco)) {
    O <- .cpgd_source_table(sources, "onco", cpgd_onco_eqtm)
    if (!is.null(O) && nrow(O)) {
      OO <- data.table::as.data.table(O)
      OO <- OO[, list(cpg_id = get("cpg_id"),
                      .gkey = gsub("_", "-", toupper(get("target_gene"))),
                      onco_direction = as.numeric(get("direction")),
                      onco_tier = as.character(get("onco_tier")),
                      onco_n_cancers = as.integer(get("n_cancers")),
                      onco_sign_agreement = as.numeric(get("sign_agreement")))]
      OO <- unique(OO, by = c("cpg_id", ".gkey"))
      out <- merge(out, OO, by = c("cpg_id", ".gkey"), all.x = TRUE, sort = FALSE)
    } else {
      out[, c("onco_direction", "onco_tier", "onco_n_cancers",
              "onco_sign_agreement") :=
            list(NA_real_, NA_character_, NA_integer_, NA_real_)]
    }
    out[, "onco_agreement" := data.table::fifelse(
          is.na(get("onco_direction")) | is.na(get("best_direction")), NA,
          get("onco_direction") == get("best_direction"))]
  }

  # ---- input provenance: trace every pair back to submitted rows ----------
  inp <- q[, list(
    input_id  = get("input_id")[1L],
    input     = get("input")[1L],
    input_ids = paste(get("input_id"), collapse = ";"),
    inputs    = paste(get("input"), collapse = ";")),
    by = "cpg_id"]
  out <- merge(out, inp, by = "cpg_id", all.x = TRUE, sort = FALSE)

  out[, "pair_id" := paste0(get("cpg_id"), "|", get("target_gene"))]
  out[, "n_targets_for_cpg" := .N, by = "cpg_id"]
  out[, "is_coeffect" := get("n_targets_for_cpg") > 1L]
  out[, ".gkey" := NULL]

  # ---- probe QC columns on the surviving rows -----------------------------
  if (probe_qc != "ignore" && !is.null(QC) && nrow(QC)) {
    qcj <- data.table::as.data.table(QC)[, intersect(
      c("cpg_id", "masked_general", "masked_partial", "mask_reasons",
        "cross_hybridizing", "mapping_flagged", "pos_hg19_verified"),
      names(QC)), with = FALSE]
    qcj <- unique(qcj, by = "cpg_id")
    data.table::setnames(
      qcj,
      old = c("masked_general", "masked_partial", "mask_reasons"),
      new = c("probe_masked", "probe_masked_partial", "probe_mask_reasons"),
      skip_absent = TRUE)
    out <- merge(out, qcj, by = "cpg_id", all.x = TRUE, sort = FALSE)
    if ("probe_masked" %in% names(out)) {
      data.table::set(out, which(is.na(out$probe_masked)), "probe_masked", FALSE)
    }
    if ("probe_masked_partial" %in% names(out)) {
      data.table::set(out, which(is.na(out$probe_masked_partial)),
                      "probe_masked_partial", FALSE)
    }
  }

  # ---- CpGs with no discovered target: visible, never silent --------------
  unmapped <- setdiff(u, discovered_cpgs)
  if (length(unmapped)) {
    warning(length(unmapped), " CpG(s) had no target gene from the EPIC-v2 ",
            "manifest, packaged catalogue, measured eQTM, SMR, or other ",
            "enabled mapping resources: ",
            paste(utils::head(unmapped, 5), collapse = ", "),
            if (length(unmapped) > 5) ", ..." else "",
            ". You may supply an explicit gene hypothesis with ",
            "gene_mode = \"pairwise\".", call. = FALSE)
  }

  # ---- column order and policy -------------------------------------------
  core <- c("input_id", "input", "cpg_id", "target_gene", "pair_id",
            "mapping_sources", "mapping_primary", "mapping_strength",
            "best_direction", "best_evidence", "best_confidence",
            "direction_tier", "probability_plus1", "best_expected_accuracy",
            "usable", "abstain_reason",
            "n_targets_for_cpg", "is_coeffect",
            "has_manifest", "has_lookup", "has_measured_eqtm", "has_smr",
            "has_brain", "has_distal", "has_onco", "has_input_annotation",
            "measured_direction", "measured_tissues",
            "lookup_direction", "lookup_probability_plus1",
            "lookup_confidence", "lookup_evidence_tier",
            "n_tissues_calling", "tissue_agreement",
            "smr_direction", "smr_tier", "p_SMR", "n_instruments",
            "instrument_agreement", "p_HEIDI", "heidi_status",
            "abs_dist", "p_universal", "dir_universal", "dist_unanimous",
            "brain_direction", "brain_tier",
            "distal_sources", "distal_n_biosamples", "distal_score_max",
            "distal_biosamples",
            "onco_direction", "onco_tier", "onco_n_cancers",
            "onco_sign_agreement", "onco_agreement",
            "probe_masked", "probe_masked_partial", "probe_mask_reasons",
            "cross_hybridizing", "mapping_flagged", "pos_hg19_verified",
            "input_ids", "inputs")
  extra <- c(paste0("cat_dir_",  c("blood", "nasal", "solid")),
             paste0("cat_conf_", c("blood", "nasal", "solid")),
             paste0("cat_tier_", c("blood", "nasal", "solid")),
             paste0("cat_prob_", c("blood", "nasal", "solid")),
             paste0("cat_status_", c("blood", "nasal", "solid")),
             "dist_dir_blood", "dist_dir_nasal", "dist_dir_solid")
  keep <- if (direction_policy == "all_evidence") c(core, extra) else core
  drop <- setdiff(names(out), keep)
  if (direction_policy != "all_evidence" && length(drop)) {
    out[, (drop) := NULL]
  }
  data.table::setcolorder(out, intersect(keep, names(out)))
  data.table::setorderv(out, c("input_id", "cpg_id", "target_gene"),
                        na.last = TRUE)

  data.table::setattr(out, "n_submitted", n_submitted)
  data.table::setattr(out, "n_unique_cpgs", length(u))
  data.table::setattr(out, "n_unmapped", length(unmapped))
  data.table::setattr(out, "unmapped_cpgs", unmapped)
  data.table::setattr(out, "source_counts", src_counts)
  data.table::setattr(out, "tissue", tissue)
  data.table::setattr(out, "gene_mode", gene_mode)
  data.table::setattr(out, "probe_qc", probe_qc)
  data.table::setattr(out, "qc_excluded_cpgs", qc_excluded)
  data.table::setattr(out, "class", c("cpgd_pairs", class(out)))
  if (isTRUE(verbose)) print(out)
  invisible(out)
}


# The genes the input itself carries (suffix-parsed or supplied pairwise), as
# candidate pairs.
#
# The parser extracts two candidates per input: the full tail after the cg
# token ("LINC02210_CRHR1", but also artefacts like "T1_CD38" from panel
# column names) and its last token ("CRHR1", "CD38"). The last token is always
# kept -- an unknown symbol there is a user hypothesis and costs at most a
# no-evidence row, while dropping it would cost a co-effect. The compound full
# form is kept only when it differs from the last token AND matches a gene
# some loaded source actually knows (dash/underscore-insensitively), because
# an unknown compound is overwhelmingly a parsing artefact of the panel
# naming, not a gene. Pairwise-supplied genes are the caller's explicit
# request and are always kept verbatim.
.cpgd_input_gene_pairs <- function(q, known_gkeys = character(0),
                                   pairwise = FALSE) {
  d <- q[!is.na(q$given_gene) | !is.na(q$given_gene2), ]
  if (!nrow(d)) return(NULL)

  g1 <- toupper(trimws(as.character(d$given_gene)))
  g2 <- toupper(trimws(as.character(d$given_gene2)))
  keep1 <- !is.na(g1) & nzchar(g1) &
    (isTRUE(pairwise) |
       is.na(g2) | g1 == g2 |
       gsub("_", "-", g1) %chin% known_gkeys)
  keep2 <- !is.na(g2) & nzchar(g2) & (is.na(g1) | g2 != g1 | !keep1)

  long <- data.table::rbindlist(list(
    data.table::data.table(cpg_id = d$cpg_id[keep1], target_gene = g1[keep1]),
    data.table::data.table(cpg_id = d$cpg_id[keep2], target_gene = g2[keep2])))
  long <- unique(long[!is.na(get("target_gene")) & nzchar(get("target_gene"))])
  if (!nrow(long)) return(NULL)
  long[, ".gkey" := gsub("_", "-", get("target_gene"))]
  long
}


#' @export
print.cpgd_pairs <- function(x, ...) {
  # A column subset or copy loses the discovery attributes; the summary would
  # then print misleading zeros. Show those as plain tables instead.
  if (is.null(attr(x, "n_submitted"))) {
    return(print(data.table::as.data.table(unclass(x)), ...))
  }
  n_sub <- attr(x, "n_submitted")
  n_cpg <- attr(x, "n_unique_cpgs")
  sc <- attr(x, "source_counts")
  cat("\ncpgdirection automatic CpG-gene discovery\n\n")
  cat(sprintf("  submitted inputs:                %6d\n",
              if (is.null(n_sub)) NA_integer_ else n_sub))
  cat(sprintf("  unique canonical CpGs:           %6d\n",
              if (is.null(n_cpg)) length(unique(x$cpg_id)) else n_cpg))
  if (!is.null(sc)) {
    cat("\n  target discovery (candidate pairs per source):\n")
    lab <- c(manifest = "EPIC-v2 manifest targets:",
             lookup   = "tissue lookup targets:",
             measured = "measured eQTM targets:",
             smr      = "SMR targets:",
             distal   = "distal-link targets:",
             brain    = "brain targets:",
             onco     = "tumour eQTM targets:",
             input    = "input-annotation targets:")
    for (k in names(lab)) {
      if (!is.na(sc[[k]])) cat(sprintf("    %-30s %6d\n", lab[[k]], sc[[k]]))
    }
  }
  if (nrow(x)) {
    ntar <- x[, list(n = data.table::uniqueN(get("target_gene"))), by = "cpg_id"]
    cat(sprintf("\n  CpGs with >=1 target:            %6d\n", nrow(ntar)))
    cat(sprintf("  CpGs with >1 target/co-effect:   %6d\n", sum(ntar$n > 1L)))
    cat(sprintf("  total unique CpG-gene pairs:     %6d\n",
                data.table::uniqueN(x$pair_id)))
    if ("usable" %in% names(x)) {
      cat("\n  pair-level direction:\n")
      cat(sprintf("    usable:                        %6d\n",
                  sum(x$usable %in% TRUE)))
      cat(sprintf("    abstain/no direction:          %6d\n",
                  sum(!(x$usable %in% TRUE))))
    }
    if ("best_evidence" %in% names(x)) {
      ev <- sort(table(x$best_evidence), decreasing = TRUE)
      cat("\n  best_evidence:\n")
      for (nm in names(ev)) cat(sprintf("    %-26s %6d\n", nm, ev[[nm]]))
    }
  }
  nunm <- attr(x, "n_unmapped")
  if (!is.null(nunm) && nunm > 0) {
    cat(sprintf("\n  CpGs with no discovered target:  %6d (see warning)\n", nunm))
  }
  qcm <- attr(x, "probe_qc")
  qce <- attr(x, "qc_excluded_cpgs")
  if (!is.null(qcm) && qcm != "ignore") {
    cat("\n  probe QC (mask / position verification):\n")
    if (qcm == "exclude") {
      cat(sprintf("    excluded as unreliable:        %6d CpGs\n",
                  if (is.null(qce)) 0L else length(qce)))
    } else if ("probe_masked" %in% names(x)) {
      cat(sprintf("    masked (kept, flagged):        %6d CpGs\n",
                  data.table::uniqueN(x$cpg_id[x$probe_masked %in% TRUE])))
    }
    if ("probe_masked_partial" %in% names(x)) {
      cat(sprintf("    partially masked replicates:   %6d CpGs\n",
                  data.table::uniqueN(x$cpg_id[x$probe_masked_partial %in% TRUE])))
    }
    if ("pos_hg19_verified" %in% names(x)) {
      cat(sprintf("    hg19 position not verified:    %6d CpGs\n",
                  data.table::uniqueN(x$cpg_id[!(x$pos_hg19_verified %in% TRUE)])))
    }
  }
  cat("\n  Mapping provenance (mapping_sources) says how a gene entered the\n",
      "  candidate set; best_evidence says what the direction rests on. They\n",
      "  are different facts and are kept apart. A CpG can carry opposite\n",
      "  directions for different genes; read pairs, not CpGs.\n", sep = "")
  cat("\n")
  cols <- intersect(c("cpg_id", "target_gene", "mapping_primary",
                      "best_direction", "best_evidence", "usable"), names(x))
  if (nrow(x)) {
    print(utils::head(as.data.frame(x[, cols, with = FALSE]), 10))
    if (nrow(x) > 10) cat(sprintf("... %d more rows\n", nrow(x) - 10))
  }
  invisible(x)
}
