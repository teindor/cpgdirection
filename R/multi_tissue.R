#' Direction across all tissues, with a cross-tissue consensus
#'
#' Runs the prediction in every packaged tissue and returns one row per CpG with
#' a column per tissue plus an aggregate. This is usually what you want when the
#' sampled tissue is a proxy for something else: saliva, blood and buccal
#' epithelium are rarely the tissue of interest, and a direction that holds
#' across several tissues is more useful than one tied to whichever proxy
#' happened to be collected.
#'
#' The aggregate is justified empirically rather than assumed. Among CpG-gene
#' pairs measured in at least two of the packaged catalogues, direction agrees
#' in 90.5\% of cases (1,267 pairs), and excess concordance on shared pairs runs
#' +0.41 to +0.46 within and between blood and epithelium. Where tissues
#' disagree, \code{consensus_direction} is \code{NA} and \code{tissue_agreement}
#' records how badly.
#'
#' @param cpgs Character vector of CpG identifiers, as for
#'   \code{\link{cpg_expression_direction}}.
#' @param genes Optional target genes, same length as \code{cpgs}.
#' @param min_tissues Minimum number of tissues that must return a call before a
#'   consensus is reported. Default 2.
#' @param allow_nearest_gene Passed through; see
#'   \code{\link{cpg_expression_direction}}.
#' @param verbose Print a summary. Default \code{TRUE}.
#'
#' @return A \code{data.table} with, per CpG: a direction, confidence and target
#'   gene for each tissue; \code{consensus_direction},
#'   \code{consensus_confidence}, \code{n_tissues_calling} and
#'   \code{tissue_agreement}; and the measured-evidence columns described in
#'   \code{\link{cpgd_measured_eqtms}}.
#'
#' @section Interpretation:
#' \code{consensus_direction} is the direction agreed by the calling tissues. It
#' is \strong{not} a claim about any one tissue, and it is not a substitute for
#' measured RNA. It is the best available answer to "if this CpG affects
#' expression of this gene anywhere we have looked, which way does it go".
#'
#' A CpG can legitimately have opposite directions for different genes, and
#' 20.3\% of CpGs with measured evidence are associated with more than one gene.
#' Read the gene columns, not just the direction.
#'
#' @examplesIf cpgd_has_data("lookup_blood_hg19")
#' res <- cpg_direction_all_tissues(c("cg00050692", "cg00000029"),
#'                                  verbose = FALSE)
#' res[, c("cpg_id", "best_direction", "best_evidence"), with = FALSE]
#' @export
cpg_direction_all_tissues <- function(cpgs, genes = NULL, min_tissues = 2L,
                                      allow_nearest_gene = FALSE, verbose = TRUE) {
  cpg_expression_direction(cpgs, genes = genes, tissue = "all",
                           min_tissues = min_tissues,
                           allow_nearest_gene = allow_nearest_gene,
                           universal = FALSE, verbose = verbose)
}


#' Measured CpG-gene associations across the packaged catalogues
#'
#' The CpG-gene pairs that were directly measured in at least one eQTM
#' cohort, rather than predicted. From 2.99.8 the table is built by
#' \code{inst/scripts/build_measured_eqtms.R} from the eQTM Atlas cohort exports
#' (Sriram et al. 2026): EVA-PR nasal epithelium and white blood cells
#' (Kim et al. 2020, 2023), GTEx whole blood and seven organs (Oliva et al.
#' 2023), FUSION skeletal muscle (Taylor et al. 2019), Colonomics colon
#' (Diez-Villanueva et al. 2021) and Framingham whole blood (Keshawarz et
#' al. 2023), plus the Framingham and HELIX (Ruiz-Arenas et al. 2022) blood
#' pairs of the earlier build that the Atlas exports do not cover. Pairs are
#' cis (|CpG - TSS| <= 250 kb), significant by the cohort's own rule, with
#' the sign taken by cohort majority; sign ties are dropped.
#' \code{tissue} is \code{blood}, \code{nasal_epithelium} or
#' \code{solid_tissue} (normal muscle, lung, ovary, colon, breast, prostate,
#' testis pooled). Tumour-derived pairs are no longer in this table; tumour
#' evidence is the Onco-eQTM layer.
#'
#' These matter for a specific reason. Nearest-gene annotation is wrong often
#' enough to be a problem: a substantial share of measured pairs involve a
#' gene that is \emph{not} the nearest gene to the CpG, and many CpGs have
#' measured associations with more than one gene. A CpG annotated to gene X
#' on the array may have its measured relationship with gene Y.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{tissue}, \code{direction}, \code{tss_dist}, and (2.99.8)
#'   \code{n_cohorts}, \code{cohort_agreement}, \code{cohorts} (pipe-separated
#'   cohort names, or \code{packaged_v1}) and \code{source}.
#' @examples
#' m <- cpgd_measured_eqtms()
#' m[m$cpg_id == "cg00000165", ]
#' @export
cpgd_measured_eqtms <- function() {
  if (!is.null(.cpgd_env$measured)) return(.cpgd_env$measured)
  M <- .cpgd_data("measured_eqtms")
  data.table::setkeyv(M, "cpg_id")
  .cpgd_env$measured <- M
  M
}


# Attach measured-evidence columns and, crucially, flag where a measured gene
# differs from the gene the caller asked about. Ruiz-Arenas et al. found
# nearest-gene annotation recovers only about half of true target genes, so a
# measured association with a DIFFERENT gene is a finding, not a footnote.
.cpgd_add_measured <- function(out) {
  M <- tryCatch(cpgd_measured_eqtms(), error = function(e) NULL)
  if (is.null(M)) {
    out[, c("measured_genes", "measured_tissues", "measured_direction",
            "annotation_mismatch", "mismatch_note") :=
          list(NA_character_, NA_character_, NA_character_, NA, NA_character_)]
    return(out)
  }
  agg <- M[, list(
      measured_genes    = paste(sort(unique(get("target_gene"))), collapse = ";"),
      measured_tissues  = paste(sort(unique(get("tissue"))), collapse = ";"),
      measured_direction = paste(sort(unique(paste0(get("target_gene"), ":",
                                 ifelse(get("direction") > 0, "+1", "-1")))), collapse = ";")),
      by = "cpg_id"]
  out <- merge(out, agg, by = "cpg_id", all.x = TRUE)
  req <- toupper(out$requested_gene)
  mg  <- strsplit(out$measured_genes, ";", fixed = TRUE)
  mm <- mapply(function(r, g) {
      # A CpG with NO measured record is not a mismatch, it is an absence of
      # evidence. strsplit(NA) returns a length-one element holding NA rather
      # than an empty vector, so testing length() alone lets every unmeasured
      # CpG through and flags it. Test the value.
      if (is.na(r) || is.null(g) || !length(g) || all(is.na(g))) return(NA)
      # Panel labels are often compound ("LINC02210_CRHR1",
      # "ENSG00000253356_STAR"); the informative symbol is the last token.
      # Accept the whole label, both separator forms, and that last token.
      cand <- unique(c(r, gsub("_", "-", r), gsub("-", "_", r),
                       utils::tail(strsplit(r, "[_-]")[[1]], 1L)))
      !any(cand %in% g)
    }, req, mg)
  out[, "annotation_mismatch" := as.logical(mm)]
  out[, "mismatch_note" := ""]
  out[which(get("annotation_mismatch")), "mismatch_note" :=
        paste0("measured association is with ", get("measured_genes"),
               ", not the requested gene - check which gene this CpG actually tracks")]
  out
}


#' Causal CpG-gene directions from two-sample SMR
#'
#' 413,982 CpG-gene pairs over 103,285 CpGs and 14,056 genes, derived by
#' summary-data Mendelian randomisation with DNA methylation as the exposure
#' (GoDMC blood mQTL meta-analysis, ~32,000 individuals) and gene expression as
#' the outcome (eQTLGen blood cis-eQTL, 31,684 individuals).
#'
#' This layer answers a question the catalogues cannot: it does not ask whether
#' anyone measured this CpG against this gene, but whether a genetic instrument
#' for methylation at the CpG moves expression of the gene. Coverage is
#' therefore not limited to pairs some study happened to test, and the median
#' CpG-gene distance is 135 kb, well beyond the range where the catalogue models
#' have support.
#'
#' @section What the tiers are worth:
#' Validated against the directly measured eQTMs shipped with this package,
#' on the pairs where both exist:
#' \tabular{lrrl}{
#'   tier \tab pairs \tab validated on \tab concordance \cr
#'   S1   \tab 29,456  \tab 2,141 \tab 95.9\% (+/- 0.8) \cr
#'   S2   \tab 366,838 \tab 6,008 \tab 84.9\% (+/- 0.9) \cr
#'   S3   \tab 17,688  \tab   456 \tab 70.4\% (+/- 4.2) \cr
#' }
#' against a majority-class baseline of 59.8\% on the same overlap. S1 requires
#' more than one instrument, all agreeing in sign, at p_SMR < 5e-8. S2 is a
#' single instrument at the same threshold. S3 has multiple instruments that
#' \emph{disagree}, and scores below S2 accordingly: instrument disagreement is
#' a warning, not extra evidence.
#'
#' @section What SMR does and does not establish:
#' A significant SMR result is consistent with methylation affecting expression,
#' and equally with \strong{linkage} -- two distinct causal variants in LD, one
#' acting on each trait. Agreement among multiple independent instruments
#' (tier S1) argues against linkage but does not exclude it.
#'
#' Both datasets are whole blood, so the causal claim is about blood whatever
#' tissue the gene matters in. Palindromic SNPs and unresolvable allele
#' mismatches were dropped rather than assumed; 21,698 instrument-gene
#' comparisons were refused on those grounds.
#'
#' @section HEIDI, and why it is reported but not applied:
#' Distinguishing mediation from linkage is what HEIDI is for. It was run on all
#' 413,982 pairs with SMR v1.03 against the 1000 Genomes European reference:
#' 22,421 pairs passed (5.4\%), 323,464 failed (78.1\%), and 68,097 (16.4\%) had
#' fewer than three instruments after LD pruning and could not be tested.
#' \strong{Not tested is not passed} and the two are never merged.
#'
#' The tiering does \strong{not} condition on HEIDI, because on these data it
#' does not discriminate in the direction that would justify a filter. Against
#' the measured eQTMs, passing pairs were concordant 80.8\% of the time and
#' failing pairs 87.3\% (p = 0.002) -- the wrong way round. The cause is power:
#' HEIDI's ability to reject rises with instrument count, and so does the
#' precision of the Wald ratio. Concordance by instruments surviving pruning was
#' 69.9\% (<3), 74.8\% (3-5), 74.0\% (6-10) and 88.0\% (11-20), while mean
#' instrument count was 17.3 among failing pairs against 12.0 among passing
#' ones. Filtering on HEIDI would discard the best-instrumented and most
#' accurate pairs first.
#'
#' The columns are shipped so you can condition on them if you wish. Do so with
#' the above in view.
#'
#' @section Replication across outcome studies (1.1.0):
#' The same GoDMC instruments were run against a second, independent outcome
#' study: the FinnGen immune-multiome single-nucleus cis-eQTLs (PBMC and eight
#' L1 cell types), SMR 1.03 per chromosome with the dense GoDMC cis mQTL set
#' as exposure (\code{tools/build_smr_replication.R}; validation in
#' \code{tools/reports/validate_celltype_smr.md} and
#' \code{tools/reports/build_smr_replication.md}). Each row now carries
#' \code{outcome}, \code{replication}, the FinnGen side
#' (\code{finngen_direction}, \code{finngen_tier}, \code{finngen_p_SMR},
#' \code{finngen_heidi_status}, \code{finngen_contexts},
#' \code{finngen_evidence}: \code{"smr_top_snp"} when the sign is the SMR
#' estimate at the top cis-mQTL, \code{"instrument_vote"} when it is the vote
#' of the pair's genome-wide-significant mQTL instruments, \code{"both"} when
#' the two exist and agree in source) and \code{smr_rung}, which is what the
#' ladder acts on:
#' \tabular{lll}{
#'   rung \tab rule \tab validated against measured blood \cr
#'   replicated \tab eQTLGen S1/S2, FinnGen PBMC same sign \tab S1 0.986 (n=1,795), S2 0.971 (n=3,387) \cr
#'   high / moderate \tab eQTLGen S1 / S2, no FinnGen PBMC pair \tab 0.954 (n=827) / 0.889 (n=2,191) \cr
#'   moderate \tab FinnGen PBMC only, SMR p < 5e-8 at the top cis-mQTL \tab 0.904 (n=1,798; Wilson 0.89-0.92) \cr
#'   weak \tab eQTLGen S3; or S1/S2 with FinnGen PBMC OPPOSITE sign (eQTLGen sign kept) \tab 0.70 (n=227); discordant S1 0.92 (n=111), S2 0.59 (n=302) \cr
#'   reported \tab FinnGen PBMC pairs below that threshold, and every cell-type-only pair \tab 0.74 (n=519) -- below the ladder \cr
#' }
#' Two facts drove the design. Concordant pairs are more accurate than either
#' study alone, because replication on an independent outcome is the one check
#' that catches a wrong-signed Wald ratio; and where the studies disagree the
#' larger study is right 77-92\% of the time, so its sign is kept but
#' demoted. The instruments are shared, so the two outcomes have correlated
#' errors and the agreement is less independent than it looks; the rung is
#' licensed by the validation against measured eQTMs, not by the agreement
#' itself. HEIDI conditions nothing, here as elsewhere in the package: against
#' eQTLGen it rejects the best-instrumented pairs (pass 0.808 vs fail 0.873
#' against measured), and against FinnGen its apparent separation rested on
#' 110 pairs; the top-cis-mQTL SMR threshold reaches 0.90 on 1,798 and needs no
#' heterogeneity filter (HEIDI-passing admitted pairs run 0.95-0.97 on n=150,
#' which is noted, and below the 200-pair bar for a rung of its own). Cell-type contexts agree with PBMC on 99\% of shared pairs and PBMC is right
#' 81-94\% of the time when they differ, so a cell type never overrides PBMC;
#' they add coverage (\code{finngen_contexts}), not sign.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{direction}, \code{smr_tier}, \code{p_SMR}, \code{n_instruments},
#'   \code{instrument_agreement}, \code{cpg_gene_dist},
#'   \code{top_instrument}, from 2.2.5 \code{p_HEIDI}, \code{nsnp_HEIDI}
#'   and \code{heidi_status} (\code{"pass"}, \code{"fail"}, \code{"not_tested"}),
#'   and from 1.1.0 \code{outcome}, \code{replication}, \code{finngen_direction},
#'   \code{finngen_tier}, \code{finngen_p_SMR}, \code{finngen_heidi_status},
#'   \code{finngen_contexts}, \code{finngen_evidence} and \code{smr_rung}.
#' @examplesIf cpgd_has_data("smr_directions")
#' s <- cpgd_smr_directions()
#' s[s$cpg_id == "cg18998365", ]
#' @export
cpgd_smr_directions <- function() {
  if (!is.null(.cpgd_env$smr)) return(.cpgd_env$smr)
  S <- .cpgd_data("smr_directions")
  data.table::setkeyv(S, c("cpg_id", "target_gene"))
  .cpgd_env$smr <- S
  S
}
