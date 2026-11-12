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
#' @param age_group,require_detected,accuracy Passed through; see
#'   \code{\link{cpg_expression_direction}}, \code{\link{cpgd_smr_directions}}
#'   and \code{\link{cpgd_accuracy_values}}.
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
                                      allow_nearest_gene = FALSE, verbose = TRUE,
                                      age_group = c("adult", "child", "any"),
                                      require_detected = FALSE,
                                      accuracy = getOption("cpgdirection.accuracy", "lower")) {
  age_group <- match.arg(age_group)
  cpg_expression_direction(cpgs, genes = genes, tissue = "all",
                           min_tissues = min_tissues,
                           allow_nearest_gene = allow_nearest_gene,
                           universal = FALSE, verbose = verbose,
                           age_group = age_group,
                           require_detected = require_detected,
                           accuracy = accuracy)
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
#' al. 2023), plus the remaining blood pairs of the earlier build
#' (\code{packaged_v1}: Framingham and HELIX Atlas exports) that no newer
#' source covers. Pairs are cis (|CpG - TSS| <= 250 kb),
#' significant by the cohort's own rule, with the sign taken by cohort
#' majority (one vote per cohort); sign ties are dropped.
#' \code{tissue} is \code{blood}, \code{nasal_epithelium} or
#' \code{solid_tissue} (normal muscle, lung, ovary, colon, breast, prostate,
#' testis pooled; from 1.1.0 also osteoarthritic cartilage and synovium). Tumour-derived pairs are no longer in this table; tumour
#' evidence is the Onco-eQTM layer.
#'
#' @section Blood cohorts added in 1.1.0 (November 2026 build):
#' The blood class grows from 26,555 to 45,311 pairs (6,997 measured in more
#' than one cohort; \code{tools/build_nov2026_extdata.R}, report in
#' \code{tools/reports/build_nov2026_extdata.md}):
#' \itemize{
#'   \item \strong{HELIX child blood} (Ruiz-Arenas et al. 2022; Dryad
#'     doi:10.5061/dryad.fxpnvx0t0; ages 6-11, n = 823): the signed,
#'     cell-type-adjusted significant pairs, cohort
#'     \code{HELIX_2022_Dryad_adj}, source \code{helix_dryad_2022}. They agree
#'     with the blood pairs of the other cohorts at 0.988 (n = 4,091, Wilson
#'     0.984-0.991) and replace the 8,432 unsigned \code{packaged_v1} HELIX
#'     rows of the earlier build.
#'   \item \strong{MESA sorted cells} (Wang et al. 2021, Table S2): CD14+
#'     monocytes (n = 1,201) and CD4+ T cells (n = 214), cohorts
#'     \code{MESA_monocytes} and \code{MESA_Tcells}, source
#'     \code{wang2021_S2}; 0.93 against the shipped blood cohorts.
#' }
#' The 250 kb window is enforced on every row; 10 duplicated rows were removed,
#' 118 pairs dropped as cohort ties and 13 existing pairs changed sign. Liver
#' (Wang et al. 2021, 0.841 against shipped solid tissue, interval straddling
#' the 0.85 bar) is held out. Nasal epithelium is unchanged.
#'
#' @section Solid-tissue source added in 1.1.0 (osteoarthritis joint tissues):
#' The solid class grows from 37,462 to 49,832 pairs
#' (\code{tools/build_kreitmaier2026_extdata.R}, report in
#' \code{tools/reports/build_kreitmaier2026_extdata.md}) with the cartilage
#' and synovium eQTMs of Kreitmaier et al. 2026 (Ann Rheum Dis,
#' doi:10.1016/j.ard.2026.01.020; low-grade cartilage n = 227, high-grade
#' n = 127, synovium n = 226; full statistics from the Musculoskeletal
#' Knowledge Portal, read by \code{tools/ingest_kreitmaier2026_oa.R}): the
#' pairs significant by the study's own rule (gene at 5\% FDR, pair at the
#' gene's Bonferroni threshold), within 250 kb, 15,596 pairs over 13,101
#' CpGs, cohorts \code{Kreitmaier2026_cartilage} and
#' \code{Kreitmaier2026_synovium}, source \code{kreitmaier2026_msk}. All
#' tissues come from the same patients, so the two cohorts are one donor
#' group and never count as a replication of each other. Against the solid
#' pairs already shipped the signs agree at 0.926 (n = 3,042); against GTEx
#' solid tissues 0.969 (n = 1,500) and against FUSION muscle 0.909
#' (n = 2,405). 12,554 pairs are new, 184 are dropped as cohort ties and 30
#' existing pairs change sign; pairs seen in two or more donor groups with
#' every cohort agreeing rise from 1,530 to 3,730. Sub-threshold pairs and
#' the fat-pad set (n = 45, no significant pair) are not included. The
#' tissue is osteoarthritic; the class label stays \code{solid_tissue}, and
#' \code{cohorts} says where a sign was measured.
#'
#' @section Sources added in 1.1.0 (kidney, placenta, adipose, GTP whole blood):
#' Four published eQTM tables, each one cohort and one donor group, each
#' admitted because its signs agree with the pairs already measured at 0.85
#' or better (\code{tools/ingest_nov2026_measured_sources.R},
#' \code{tools/build_nov2026_sources_extdata.R}, report in
#' \code{tools/reports/build_nov2026_sources_extdata.md}):
#' \itemize{
#'   \item \strong{Kidney tubule} (Liu et al. 2022, Nat Genet,
#'     doi:10.1038/s41588-022-01097-w, Supplementary Table 12; n = 414, EPIC,
#'     RNA-seq, global FDR below 0.05): 29,988 pairs within 250 kb, cohort
#'     \code{Liu2022_kidney}, source \code{liu2022_kidney}, solid tissue;
#'     0.897 against all shared measured pairs (n = 6,156), 0.908 against
#'     solid tissue (n = 4,685). The published sheet lists its CpG and
#'     Ensembl identifiers in a different row order from its statistics, so
#'     each CpG is recovered from the row's chromosome and position; the
#'     recovery reproduces the 16,136 pairs of the study's Supplementary
#'     Table 11 it can be checked on, with identical effect sizes.
#'   \item \strong{Placenta} (Delahaye et al. 2018, PLoS Genet,
#'     doi:10.1371/journal.pgen.1007785, S20 Table; n = 74, 450K, RNA-seq):
#'     2,314 pairs, cohort \code{Delahaye2018_placenta}, source
#'     \code{delahaye2018_placenta}, solid tissue; 0.901 (n = 434).
#'   \item \strong{Adipose tissue} (Voisin et al. 2016, PLoS One,
#'     doi:10.1371/journal.pone.0157776, S20 Table; n = 118, 450K): 427
#'     pairs, cohort \code{Voisin2016_adipose}, source
#'     \code{voisin2016_adipose}, solid tissue; 0.913 (n = 277).
#'   \item \strong{Whole blood, Grady Trauma Project} (Kennedy et al. 2018,
#'     BMC Genomics, doi:10.1186/s12864-018-4842-3, Additional file 1;
#'     n = 333, all African American, 450K, expression array): the
#'     same-chromosome pairs of the file (p below 1e-5) within 250 kb, 3,735
#'     pairs, cohort \code{Kennedy2018_GTP_blood}, source
#'     \code{kennedy2018_gtp}, blood; 0.932 (n = 2,259), and 0.939 against
#'     the blood pairs (n = 2,156).
#' }
#' Solid tissue grows from 49,832 to 76,958 pairs and blood from 45,311 to
#' 46,790; 412 solid and 100 blood pairs are dropped as cohort ties and 2
#' existing solid pairs change sign. Pairs seen in two or more donor groups
#' with every cohort agreeing rise from 3,730 to 6,891 in solid tissue and
#' from 6,333 to 7,159 in blood. CpGs of the 450K array that the package's
#' position table does not hold are placed by the Illumina 450K manifest
#' (build 37). The class labels stay \code{solid_tissue} and \code{blood};
#' \code{cohorts} says where a sign was measured. Retina, liver, pancreatic
#' islets and fibroblasts were appraised with the same test, fell below the
#' bar (0.73 to 0.82), and are not part of this table.
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
#'   cohort names, or \code{packaged_v1}) and \code{source}
#'   (\code{eqtm_atlas}, \code{helix_dryad_2022}, \code{wang2021_S2},
#'   \code{kreitmaier2026_msk}, \code{liu2022_kidney},
#'   \code{delahaye2018_placenta}, \code{voisin2016_adipose},
#'   \code{kennedy2018_gtp} or \code{packaged_v1}).
#' @seealso \code{\link{cpgd_brain_measured}} for measured pairs in brain
#'   cortex, which are kept out of this table and used only by
#'   \code{\link{cpg_brain_bridge}}.
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
#'   moderate \tab eQTLGen S3 whose majority sign FinnGen PBMC reproduces \tab 0.855 (n=214; Wilson 0.80-0.90) \cr
#'   weak \tab eQTLGen S3 without a PBMC pair; or S1 with FinnGen PBMC OPPOSITE sign (eQTLGen sign kept) \tab 0.70 (n=227); discordant S1 0.92 (n=111) \cr
#'   reported \tab FinnGen PBMC pairs below that threshold; every cell-type-only pair; eQTLGen S2/S3 pairs FinnGen PBMC contradicts (evidence conflicts, no sign called) \tab 0.74 (n=519); discordant S2 0.59 (n=302), S3 0.38 (n=68) \cr
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
#' @section A third outcome study and a per-pair accuracy (1.1.0, INTERVAL rebuild):
#' The GoDMC instruments were also run against the INTERVAL whole-blood
#' RNA-seq cis-eQTLs (SMR and HEIDI on the full summary statistics;
#' \code{tools/build_smr_layer_interval.py}, report in
#' \code{tools/reports/build_smr_layer_interval.md}). The table grows from
#' 549,316 to 768,862 pairs (219,546 seen in INTERVAL alone) and each row
#' carries the INTERVAL side: \code{interval_direction},
#' \code{interval_tier}, \code{interval_p_SMR}, \code{interval_p_HEIDI},
#' \code{interval_evidence} (\code{"both"}: the top-SNP SMR test passes
#' 5e-8 and the clumped instruments vote the same way; \code{"smr_top_snp"};
#' \code{"instrument_vote"}: only the vote exists; \code{"conflict"}: the
#' two disagree, the top-SNP sign is kept) and \code{interval_status}
#' (\code{agrees}, \code{disagrees}, \code{absent}, \code{interval_only}).
#'
#' \strong{Rung rule.} With three bulk outcomes (eQTLGen, FinnGen PBMC,
#' INTERVAL) a pair is \code{replicated} when at least two give the same
#' sign, none dissents, and either the instruments are concordant (tier S1)
#' or one outcome reaches p < 1e-40. An absent outcome is neutral. A dissent
#' demotes: two against one is \code{weak}, one against one is
#' \code{reported}. Only an INTERVAL record whose top-SNP test passes and
#' does not conflict with its instrument vote can promote a pair; any
#' INTERVAL dissent demotes; a pair seen in INTERVAL alone enters by its own
#' tier when its record is of that kind and as \code{weak} otherwise.
#' Against the sign measured in blood (pairs with a measured blood sign):
#' \tabular{lll}{
#'   rung \tab pairs \tab SMR sign = measured blood sign \cr
#'   replicated \tab 108,135 \tab 0.963 (n = 13,283) \cr
#'   high \tab 14,284 \tab 0.949 (n = 942) \cr
#'   moderate \tab 528,204 \tab 0.859 (n = 8,608) \cr
#'   weak \tab 52,878 \tab 0.755 (n = 1,377) \cr
#'   reported \tab 65,361 \tab 0.412 (n = 1,667) \cr
#' }
#'
#' \strong{Per-pair accuracy.} \code{smr_accuracy} is the probability that
#' the pair's SMR sign equals the sign measured in blood, from a
#' gradient-boosted model of the pair's SMR evidence (strength and agreement
#' of each outcome, tiers, instruments, HEIDI, distance, the called sign),
#' fitted leaving the pair's chromosome out, so no figure comes from a model
#' that saw the pair or its neighbors (25,875 validated pairs; AUC 0.85;
#' mean calibration gap 0.009 over deciles). It is a figure for blood, and
#' like every SMR figure it is conditional on the pair being coupled.
#'
#' \strong{What the ladder uses.} The direction ladder no longer ranks SMR
#' pairs by rung. It groups them into accuracy bands (0.95 or higher, 0.90
#' to 0.95, 0.80 to 0.90, 0.65 to 0.80, 0.55 to 0.65), validates each band
#' in every tissue against that tissue's independent yardstick like any
#' other source, and ranks it on the lower bound of that figure
#' (\code{\link{cpgd_accuracy_values}}). A pair below 0.55, and any pair on
#' the \code{reported} rung, supplies no direction; a pair below 0.50 is
#' not flipped. In blood (pairs that can supply a direction; sign against
#' measured blood): 0.95 or higher 0.975 (n = 13,587), 0.90 to 0.95 0.917
#' (n = 4,269), 0.80 to 0.90 0.866 (n = 3,352), 0.65 to 0.80 0.743
#' (n = 1,970), 0.55 to 0.65 0.643 (n = 502). The rung stays in
#' \code{smr_rung} and in \code{best_evidence} as the row's provenance.
#' The sections below describe the eQTLGen and FinnGen analyses that the
#' rebuild extends; their figures by rung predate INTERVAL.
#'
#' @section replicated and high are close, and ranked by the data:
#' Scored against five measured references -- adult blood, child blood, nasal
#' epithelium, solid tissue and adult cortex
#' (\code{tools/validate_ladder_transport.R}, 4 Nov 2026) -- the two rungs
#' differ by 0.00-0.03 with overlapping intervals in every one. From 1.1.0 the
#' direction ladder ranks every SMR stratum (replicated S1 and S2, high, the
#' moderate and weak strata) by its agreement with each tissue's independent
#' yardstick (\code{tools/validate_rung_accuracy.R};
#' \code{\link{cpgd_accuracy_values}}): in adult blood replicated S1 (0.93)
#' and high (0.91) sit above replicated S2 (0.85).
#'
#' @section Every figure is conditional on coupling:
#' An SMR accuracy is the probability that the sign is right \emph{given that
#' the pair is coupled in the tissue}. An SMR direction says which way a CpG
#' would act on the gene if it acts on it, not that it does. In child blood
#' (HELIX, n = 823, cell-adjusted, within 250 kb) most SMR pairs show no
#' detectable coupling: of the 241,044 pairs HELIX tested, the share with
#' nominal p < 0.05 is 0.40 (replicated), 0.40 (high), 0.20 (moderate), 0.28
#' (weak) and 0.21 (reported), against 0.05 under the null; 60\% of
#' replicated and 80\% of moderate pairs are not detected. Sign agreement
#' among the nominal hits, corrected for the null hits they contain,
#' reproduces the significant-pair figures below, so the two analyses agree
#' (helix_full_calibration_ruling.md, 4 Nov 2026). \code{smr_moderate} is
#' large and diluted -- 160k tested pairs, 4.1\% HELIX-significant, close to
#' the reported rung's 4.6\%.
#'
#' @section Child blood and age_group:
#' The ladder's figures for the SMR accuracy bands
#' (\code{\link{cpgd_accuracy_values}}; adult blood against MESA and
#' BLUEPRINT, child blood against HELIX pairs outside the catalogue's
#' training set; agreement, n pairs):
#' \tabular{lll}{
#'   SMR band (per-pair accuracy) \tab adult blood \tab child blood (HELIX 6-11 y) \cr
#'   0.95 or higher \tab 0.961 (11,549) \tab 0.926 (917) \cr
#'   0.90 to 0.95 \tab 0.904 (4,424) \tab 0.928 (567) \cr
#'   0.80 to 0.90 \tab 0.846 (3,719) \tab 0.893 (438) \cr
#'   0.65 to 0.80 \tab 0.730 (2,368) \tab 0.764 (288) \cr
#'   0.55 to 0.65 \tab 0.618 (607) \tab 26 pairs: the adult figure is quoted \cr
#'   below 0.55, or the reported rung \tab never answers \tab never answers \cr
#' }
#' Age moves no direction and no rung: where HELIX scores a band on 100 or
#' more pairs its figure is within 0.05 of the adult one (higher in three
#' bands, lower in the top one). The share of pairs HELIX detects at all, by
#' rung, is given in the section above.
#' \code{age_group} therefore changes only (1) which validation figures the
#' ladders quote (every blood source HELIX scores on 100+ pairs, not only
#' SMR), (2) which age-specific evidence columns are returned and (3) an
#' optional filter. \code{"adult"} is the default and gives the same output as
#' leaving the argument out. \code{"child"} adds \code{child_blood_helix} --
#' \code{"sig"} (HELIX significant pair), \code{"nominal"} (p < 0.05, not
#' significant), \code{"null"} (tested, p >= 0.05) or \code{"untested"} --
#' and quotes the child-blood figures; \code{"any"} quotes both. A
#' \code{"null"} means "not detected at n = 823", not "absent", and never
#' changes a rung. \code{require_detected = TRUE} (child only) keeps the pairs
#' whose \code{child_blood_helix} is \code{"sig"} or \code{"nominal"}.
#' \strong{Project Alpha users (children and adolescents) should set
#' \code{age_group = "child"}.}
#'
#' @section What the age references cover:
#' Blood: children aged 6-11 (HELIX) and adults. Cortex: 1-25 years (Jaffe
#' DLPFC, age-adjusted) and adults (ROSMAP). There is no adolescent blood
#' reference, no fetal or infant cortex, and no age information for nasal
#' epithelium or solid tissue.
#'
#' @section Cross-tissue transport:
#' Sign agreement against each independent measured reference (Wilson
#' intervals in \code{tools/reports/validate_ladder_transport.md}; measured
#' blood excludes the reference's own rows):
#' \tabular{llllll}{
#'   source \tab adult blood \tab child blood \tab nasal \tab solid \tab adult cortex \cr
#'   SMR replicated \tab 0.914 (1,712) \tab 0.970 (6,444) \tab 0.777 (728) \tab 0.803 (2,016) \tab 0.772 (1,055) \cr
#'   SMR high \tab 0.927 (589) \tab 0.970 (1,233) \tab 0.780 (236) \tab 0.815 (406) \tab 0.739 (218) \cr
#'   SMR moderate \tab 0.803 (2,878) \tab 0.895 (6,576) \tab 0.623 (1,692) \tab 0.698 (3,112) \tab 0.666 (1,693) \cr
#'   SMR weak \tab 0.710 (238) \tab 0.711 (454) \tab 0.635 (104) \tab 0.663 (169) \tab 0.725 (91) \cr
#'   SMR reported \tab 0.616 (521) \tab 0.668 (1,309) \tab 0.541 (362) \tab 0.561 (651) \tab 0.533 (396) \cr
#'   measured blood \tab -- \tab 0.988 \tab 0.911 (1,813) \tab 0.899 (4,149) \tab 0.840 (1,654) \cr
#' }
#' This is the 3 Nov 2026 transport analysis. The rung order is monotone in
#' all five references. Age does not move the ladder; tissue does. Measured
#' signs transport across tissue far better than blood SMR signs: measured
#' blood beats SMR moderate by 0.17 in cortex, 0.20 in solid tissue and 0.29
#' in nasal epithelium. The ladder's own order and figures
#' (\code{\link{cpgd_accuracy_values}}) come from a later validation that
#' scores every source, SMR and catalogue alike, against one yardstick per
#' tissue that is independent of it; for adult blood that is MESA monocytes
#' and T cells, which no catalogue was trained on. There, in adult blood, SMR
#' moderate (0.80 [0.77, 0.82]) ranks below every measured blood sign and
#' just above the catalogue consensus (0.80 [0.76, 0.83], formed from
#' predicted calls); the ladder ranks by the lower bound.
#'
#' @param age_group \code{"adult"} (default), \code{"child"} or \code{"any"}.
#'   Changes which validation figures are quoted and which age-specific
#'   columns are returned, never a direction or a rung. \code{"adult"}
#'   returns the table without \code{child_blood_helix}, as before 1.1.0.
#' @param require_detected With \code{age_group = "child"} only: keep the
#'   pairs whose \code{child_blood_helix} is \code{"sig"} or
#'   \code{"nominal"}; a message reports how many were dropped. With any
#'   other \code{age_group} it is ignored with a warning.
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{direction}, \code{smr_tier}, \code{p_SMR}, \code{n_instruments},
#'   \code{instrument_agreement}, \code{cpg_gene_dist},
#'   \code{top_instrument}, from 2.2.5 \code{p_HEIDI}, \code{nsnp_HEIDI}
#'   and \code{heidi_status} (\code{"pass"}, \code{"fail"}, \code{"not_tested"}),
#'   and from 1.1.0 \code{outcome}, \code{replication}, \code{finngen_direction},
#'   \code{finngen_tier}, \code{finngen_p_SMR}, \code{finngen_heidi_status},
#'   \code{finngen_contexts}, \code{finngen_evidence} and \code{smr_rung};
#'   from the INTERVAL rebuild \code{interval_direction},
#'   \code{interval_tier}, \code{interval_p_SMR}, \code{interval_evidence},
#'   \code{interval_p_HEIDI}, \code{interval_status} and
#'   \code{smr_accuracy};
#'   with \code{age_group = "child"} or \code{"any"} also
#'   \code{child_blood_helix} (\code{NA} where the resource predates it).
#' @examplesIf cpgd_has_data("smr_directions")
#' s <- cpgd_smr_directions()
#' s[s$cpg_id == "cg18998365", ]
#'
#' # children: the HELIX child-blood detection column
#' k <- cpgd_smr_directions(age_group = "child")
#' table(k$child_blood_helix, useNA = "ifany")
#' @export
cpgd_smr_directions <- function(age_group = c("adult", "child", "any"),
                                require_detected = FALSE) {
  age_group <- match.arg(age_group)
  require_detected <- .cpgd_require_detected(age_group, require_detected)
  S <- .cpgd_smr_table()
  has_child <- "child_blood_helix" %in% names(S)
  if (identical(age_group, "adult")) {
    if (!has_child) return(S)
    return(S[, setdiff(names(S), "child_blood_helix"), with = FALSE])
  }
  require_detected <- .cpgd_check_child_column(S, age_group, require_detected)
  # a copy, never the cache: a user's := on the result must not reach the
  # table every ladder reads
  S <- data.table::copy(S)
  if (!has_child) S[, "child_blood_helix" := NA_character_]
  if (isTRUE(require_detected)) {
    keep <- S$child_blood_helix %in% CPGD_CHILD_DETECTED
    message(sprintf(paste0("require_detected: %s of %s SMR pairs dropped ",
                           "(child_blood_helix not sig/nominal)."),
                    format(sum(!keep), big.mark = ","),
                    format(length(keep), big.mark = ",")))
    S <- S[keep]
  }
  S
}

# The full cached SMR table, every column, by reference. Internal callers use
# this rather than the exported accessor, which drops child_blood_helix for
# the adult default and would copy a 549k-row table on every call.
.cpgd_smr_table <- function() {
  if (!is.null(.cpgd_env$smr)) return(.cpgd_env$smr)
  S <- .cpgd_data("smr_directions")
  data.table::setkeyv(S, c("cpg_id", "target_gene"))
  .cpgd_env$smr <- S
  S
}
