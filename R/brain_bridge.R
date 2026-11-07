# The brain bridge: what a peripheral measurement licenses about the brain.
# New in 2.3.0. Everything here is built from open resources (IMAGE-CpG, AMAZE,
# GSE59685, Sommerer et al., Gibbs GSE15745, Brain-mMeta x BrainMeta v2);
# components derived from CC-BY-NC sources are deliberately NOT shipped.

# Brain-layer accuracy, with the same-person overlap removed (4 Nov 2026,
# tools/validate_brain_overlap.R). Brain-mMeta, the exposure arm,
# meta-analyses the ROSMAP, Hannon et al. and Jaffe et al. (GSE74193) cortex
# mQTLs (Qi et al. 2018, Nat Commun 9:2282); BrainMeta v2, the outcome arm,
# includes ROSMAP (Qi et al. 2022, Nat Genet 54:1355). Agreement with ROSMAP is
# therefore in-sample (T1 0.995, T2 0.939). The quoted figure is the agreement
# expected against an adult cortex cohort that shares no donor with either
# arm, estimated without comparing any two donor-sharing datasets: each tier's
# agreement with donor-independent comparators (blood without EVA-PR/GTEx,
# EVA-PR nasal, solid tissue without GTEx), as a ratio to ROSMAP's own on the
# same pairs, times ROSMAP's sign accuracy (anchored on Jaffe at p < 0.01);
# cluster bootstrap by CpG, 1,000 resamples. The in-sample ROSMAP figures
# showed no dose-response by SMR strength, and T2 shows no inflation at all.
CPGD_BRAIN_ACCURACY <- c(
  T1 = paste0("0.93-1.00 against an independent adult cortex cohort (0.964, ",
              "same-person overlap removed; in-sample ROSMAP DLPFC 0.995, n=1,521; ",
              "GSE15745 27K 0.92-0.96)"),
  T2 = paste0("0.90-0.98 against an independent adult cortex cohort (0.943, ",
              "same-person overlap removed; in-sample ROSMAP DLPFC 0.939, n=2,904; ",
              "GSE15745 27K 0.67-0.77)"))

# What each rung of the brain-target ladder is worth in brain (ruling 3,
# tools/validate_ladder_transport.R, 4 Nov 2026). The blood figures are
# transport to adult cortex (ROSMAP), which shares no sample with any blood
# source. Order: measured brain > brain T1 > brain T2 > measured blood > blood
# SMR; blood evidence is a labelled fallback for a brain question, never above
# brain SMR.
CPGD_BRAIN_LADDER <- c("measured_brain", "brain_T1", "brain_T2",
                       "measured_blood", "blood_SMR")
# The blood SMR fallback, ranked by its accuracy in adult cortex.
CPGD_BRAIN_SMR_ORDER <- c("replicated", "high", "weak", "moderate")
CPGD_BRAIN_TARGET_ACCURACY <- c(
  measured_brain_replicated = "~1.00 (measured in cortex, same sign in ROSMAP and Jaffe)",
  measured_brain_single = paste0("0.93-0.98 (measured in one cortex cohort; ROSMAP and ",
                                 "Jaffe agree on 174 of 180 shared pairs)"),
  brain_T1 = CPGD_BRAIN_ACCURACY[["T1"]],
  brain_T2 = CPGD_BRAIN_ACCURACY[["T2"]],
  measured_blood = paste0("0.84 against adult cortex (measured blood eQTMs vs ROSMAP, ",
                          "n=1,654; below the 0.85 bar: fallback only)"),
  blood_SMR_replicated = "0.77 against adult cortex (blood SMR replicated, n=1,055; fallback only)",
  blood_SMR_high       = "0.74 against adult cortex (blood SMR high, n=218; fallback only)",
  blood_SMR_moderate   = "0.67 against adult cortex (blood SMR moderate, n=1,693; fallback only)",
  blood_SMR_weak       = "0.73 against adult cortex (blood SMR weak, n=91; fallback only)")

# The age caveat the bridge repeats, by age_group (lead18b_result_ruling.md,
# 4 Nov 2026: with age removed, 1-25 y cortex eQTM signs match adult ROSMAP
# cortex at 0.975, the same as Jaffe adults > 25 y at 0.980; the earlier 0.73
# developmental figure was age confounding).
CPGD_BRAIN_CAVEAT <- c(
  adult = paste0("Brain directions are validated in adult cortex (ROSMAP aged DLPFC; ",
                 "Jaffe DLPFC > 25 y, age-adjusted, 0.980 vs ROSMAP) and hold in ",
                 "postnatal cortex from 1 year (Jaffe 1-25 y, age-adjusted, 0.975); ",
                 "fetal and infant (< 1 y) cortex is untested."),
  child = paste0("Brain directions are validated in cortex aged 1-25 y (Jaffe, ",
                 "age-adjusted, 0.975 vs adult cortex); fetal/< 1 y untested."),
  any   = paste0("Brain directions are validated in postnatal cortex from 1 year to ",
                 "old age (ROSMAP adult; Jaffe 1-25 y, age-adjusted, 0.975 vs ROSMAP); ",
                 "fetal and infant (< 1 y) cortex is untested."))

#' Brain causal directions (tier T1/T2)
#'
#' The genome-wide brain layer: 220,383 CpG-gene pairs with an SMR causal
#' direction (Brain-mMeta mQTLs as exposure, BrainMeta v2 cis-eQTLs as outcome),
#' graded by a score validated out-of-sample in blood and against the
#' independent Gibbs brain catalogue.
#'
#' @section What the tiers are worth in brain:
#' Against measured brain eQTMs:
#' \tabular{llll}{
#'   reference \tab overlap with the SMR arms \tab T1 \tab T2 \cr
#'   GSE15745 (Gibbs 2010, 27K) \tab none (verified absent from both) \tab 0.92-0.96 \tab 0.67-0.77 \cr
#'   ROSMAP aged DLPFC (450K) \tab in both arms: in-sample \tab 0.995 (n=1,521) \tab 0.939 (n=2,904) \cr
#'   Jaffe DLPFC 25-78 y (450K) \tab methylation donors in the exposure arm \tab 0.969 (n=65) \tab 0.918 (n=61) \cr
#'   independent adult cortex, overlap removed \tab none (estimated) \tab 0.964 [0.926, 0.998] \tab 0.943 [0.899, 0.984] \cr
#' }
#' Brain-mMeta meta-analyses the ROSMAP, Hannon et al. and Jaffe et al.
#' (GSE74193) cortex mQTLs and BrainMeta v2 includes ROSMAP, so the 450K
#' figures are not independent validations: the ROSMAP figures are in-sample
#' on both arms and the Jaffe figures share donors with the exposure arm. The
#' last row removes the same-person overlap without comparing any two
#' donor-sharing datasets (\code{tools/validate_brain_overlap.R}, 4 Nov 2026):
#' each tier's agreement with donor-independent comparators, as a ratio to
#' ROSMAP's own agreement with them on the same pairs, scaled by ROSMAP's sign
#' accuracy; cluster bootstrap by CpG. It is the agreement expected against
#' an adult cortex cohort that shares no donor with either arm, and the
#' figure \code{\link{cpg_brain_bridge}} quotes. The in-sample inflation is
#' T1 0.995 to 0.964; T2 shows none (0.939 to 0.943). Restricting to pairs within 100 kb changes nothing (T1
#' 0.974, T2 0.794 against the Wang et al. 2021 brain pairs, the same as at
#' any distance), so the earlier distance caveat is dropped. Every figure is
#' conditional on an association existing.
#'
#' Blood evidence transports to adult cortex less well than this layer:
#' measured blood eQTMs agree with measured cortex 0.84 (n = 1,654) and the
#' blood SMR rungs 0.67-0.77 (\code{tools/validate_ladder_transport.R}), so a
#' peripheral direction is never substituted for a brain one;
#' \code{\link{cpg_brain_bridge}} uses blood only as a labelled fallback.
#'
#' @section Age:
#' Validated in postnatal cortex from 1 year to old age: with age regressed
#' out, cis eQTM signs in 1-25-year cortex (Jaffe, n = 94) match adult ROSMAP
#' cortex at 0.975 and Jaffe adults over 25 at 0.980. Fetal and infant
#' (< 1 year) cortex is untested. The earlier 0.73 "developmental" figure was
#' age confounding, not biology.
#'
#' @return A \code{data.table}: \code{cpg_id}, \code{target_gene},
#'   \code{direction}, \code{tier}, \code{score}, \code{n_instruments},
#'   \code{instrument_agreement}, \code{cpg_gene_dist}, \code{hla},
#'   \code{inv17q21}, \code{p_SMR}, \code{p_HEIDI}, \code{heidi_status}.
#' @examplesIf cpgd_has_data("brain_directions")
#' b <- cpgd_brain_directions()
#' head(b)
#' @export
cpgd_brain_directions <- function() {
  if (!is.null(.cpgd_env$brain)) return(.cpgd_env$brain)
  B <- .cpgd_data("brain_directions")
  data.table::setkeyv(B, "cpg_id")
  .cpgd_env$brain <- B
  B
}

#' Measured CpG-gene associations in brain cortex
#'
#' Directly measured cis-eQTM pairs in dorsolateral prefrontal cortex, the
#' measured rung of \code{\link{cpg_brain_bridge}}. Kept apart from
#' \code{\link{cpgd_measured_eqtms}}, whose classes are peripheral, so that
#' brain evidence never enters a peripheral \code{best_direction}.
#'
#' Two cohorts, shipped as \code{inst/extdata/brain_measured_eqtms.csv.gz}
#' with one row per pair per cohort:
#' \itemize{
#'   \item \code{ROSMAP_DLPFC_aged} (Wang et al. 2021, Table S2; n = 502,
#'     covariates including age regressed; source \code{wang2021_S2}): 8,815
#'     pairs over 8,150 CpGs.
#'   \item \code{Jaffe_DLPFC_ageadj} (lead 18b; GSE74193 450K x GSE30272
#'     expression, donors paired by age-sex-race fingerprint, model
#'     \code{expr ~ meth + ns(age, 4) + sex + race}, postnatal 1-25 y and
#'     adult > 25 y sets at FDR < 0.05, one cohort vote, majority sign;
#'     source \code{lead18b_GSE74193_GSE30272}): 383 pairs.
#' }
#' The two agree on 174 of 180 shared pairs (0.967; 0.85 admission bar).
#' The Wang developmental (fetal to 25 y, age not regressed) pairs are not
#' used anywhere: most were age-trend artefacts.
#'
#' @param collapse \code{TRUE} (default) returns one row per pair: one vote
#'   per cohort (a cohort whose own rows disagree abstains), the majority
#'   sign, ties dropped -- the rule \code{\link{cpgd_measured_eqtms}} uses.
#'   \code{FALSE} returns the per-cohort rows as shipped.
#' @return With \code{collapse = TRUE}, a \code{data.table} keyed on
#'   \code{cpg_id}, \code{target_gene}, with \code{direction} (+1/-1),
#'   \code{n_cohorts}, \code{cohort_agreement}, \code{cohorts}
#'   (pipe-separated), \code{source} (pipe-separated) and \code{tss_dist}.
#'   With \code{collapse = FALSE}, the shipped columns \code{cpg_id},
#'   \code{target_gene}, \code{direction}, \code{corr}, \code{FDR},
#'   \code{tss_dist}, \code{cohort}, \code{source}.
#' @examples
#' bm <- cpgd_brain_measured()
#' table(bm$n_cohorts)
#' bm[bm$cpg_id == "cg06846259", ]
#' @export
cpgd_brain_measured <- function(collapse = TRUE) {
  if (!isTRUE(collapse) && !identical(collapse, FALSE)) {
    stop("`collapse` must be TRUE or FALSE.", call. = FALSE)
  }
  # copies, never the cache the bridge reads
  data.table::copy(.cpgd_brain_measured_table(collapse))
}

# The cached measured-brain tables, by reference, for internal use.
.cpgd_brain_measured_table <- function(collapse = TRUE) {
  if (is.null(.cpgd_env$brain_measured_raw)) {
    R <- .cpgd_data("brain_measured_eqtms")
    if (!"cohort" %in% names(R)) {
      R[, "cohort" := if ("source" %in% names(R)) as.character(get("source"))
                      else "brain_measured"]
    }
    R[, "target_gene" := toupper(trimws(as.character(get("target_gene"))))]
    .cpgd_env$brain_measured_raw <- R
  }
  if (!isTRUE(collapse)) return(.cpgd_env$brain_measured_raw)
  if (!is.null(.cpgd_env$brain_measured)) return(.cpgd_env$brain_measured)
  M <- .cpgd_cohort_vote(.cpgd_env$brain_measured_raw)
  .cpgd_env$brain_measured <- M
  M
}

# One vote per cohort per pair, majority sign, ties dropped. A cohort whose
# own rows for a pair disagree casts no vote. `D` needs cpg_id, target_gene,
# direction, cohort; source and tss_dist are carried when present.
.cpgd_cohort_vote <- function(D) {
  D <- data.table::as.data.table(D)
  V <- data.table::data.table(
    cpg_id = as.character(D$cpg_id),
    target_gene = toupper(trimws(as.character(D$target_gene))),
    d = sign(suppressWarnings(as.numeric(D$direction))),
    cohort = as.character(D$cohort),
    source = if ("source" %in% names(D)) as.character(D$source) else NA_character_,
    tss = if ("tss_dist" %in% names(D)) suppressWarnings(as.numeric(D$tss_dist))
          else NA_real_)
  V <- V[!is.na(get("d")) & get("d") != 0]
  V <- V[, list(v = sign(sum(get("d"))),
                source = get("source")[1L],
                tss = get("tss")[!is.na(get("tss"))][1L]),
         by = c("cpg_id", "target_gene", "cohort")]
  V <- V[get("v") != 0]
  M <- V[, list(plus  = sum(get("v") > 0),
                minus = sum(get("v") < 0),
                cohorts = paste(sort(unique(get("cohort"))), collapse = "|"),
                source  = paste(sort(unique(stats::na.omit(get("source")))), collapse = "|"),
                tss_dist = get("tss")[!is.na(get("tss"))][1L]),
         by = c("cpg_id", "target_gene")]
  M <- M[get("plus") != get("minus")]
  M[, "direction" := as.integer(sign(get("plus") - get("minus")))]
  M[, "n_cohorts" := as.integer(get("plus") + get("minus"))]
  M[, "cohort_agreement" := round(pmax(get("plus"), get("minus")) / get("n_cohorts"), 3)]
  M <- M[, c("cpg_id", "target_gene", "direction", "n_cohorts",
             "cohort_agreement", "cohorts", "source", "tss_dist"), with = FALSE]
  data.table::setkeyv(M, c("cpg_id", "target_gene"))
  M
}

#' The graded periphery-to-brain deliverable
#'
#' 728 CpGs where a tier-T1 brain causal direction coexists with robust
#' peripheral-brain methylation concordance (rho >= 0.5 on BOTH IMAGE-CpG
#' arrays, in at least two peripheral tissues, outside the extended MHC), with
#' per-cohort support flags and grades. Grade A++ (n = 205; 184 saliva-usable)
#' is supported by four independent cohorts.
#' @return A \code{data.table}, one row per CpG; see the paper for column
#'   definitions. \code{strict_qc == FALSE} rows sit on probes failing the
#'   stringent probe filter and should be treated as provisional.
#' @examples
#' d <- cpgd_bridge_deliverable()
#' head(d)
#' @export
cpgd_bridge_deliverable <- function() {
  if (!is.null(.cpgd_env$bridge)) return(.cpgd_env$bridge)
  D <- .cpgd_data("bridge_deliverable")
  data.table::setkeyv(D, "cpg_id")
  .cpgd_env$bridge <- D
  D
}

#' What does a peripheral CpG measurement license about the brain?
#'
#' For each input CpG: the measured cortex associations and the brain causal
#' direction and target gene(s) where they exist, the peripheral-brain
#' concordance evidence for the tissue you measured, the deliverable grade
#' where the CpG reaches the graded set, and otherwise a genome-wide
#' synthesised saliva-bridge score.
#'
#' @section The brain-target ladder:
#' Where the input names a gene, \code{requested_gene_best_direction} answers
#' for that gene from the strongest evidence held, in this order
#' (\code{tools/validate_ladder_transport.R}, 4 Nov 2026):
#' \tabular{ll}{
#'   \code{requested_gene_best_evidence} \tab worth in brain \cr
#'   \code{measured_brain} \tab measured in cortex (\code{\link{cpgd_brain_measured}}); ~1.00 in both cohorts, 0.93-0.98 in one \cr
#'   \code{brain_T1} \tab 0.964 [0.926, 0.998] against independent adult cortex (0.995 in-sample on ROSMAP) \cr
#'   \code{brain_T2} \tab 0.943 [0.899, 0.984] against independent adult cortex (0.939 in-sample on ROSMAP) \cr
#'   \code{measured_blood} \tab 0.84 against adult cortex (below the 0.85 bar) \cr
#'   \code{blood_SMR} \tab 0.67-0.77 against adult cortex, by rung \cr
#' }
#' Blood evidence is a labelled fallback for a brain question and never ranks
#' above brain SMR. \code{requested_gene_direction} and
#' \code{requested_gene_tier} keep their earlier meaning: the brain SMR layer
#' alone, never substituted. See \code{\link{cpgd_brain_directions}} for why
#' the 450K figures are in-sample.
#'
#' @section The caveats this function will keep repeating:
#' \enumerate{
#'   \item Directions are BRAIN-derived first. Blood evidence transports to
#'     adult cortex at 0.84 (measured) and 0.67-0.77 (SMR) and enters only as
#'     the labelled fallback of the ladder above.
#'   \item Age: validated in postnatal cortex from 1 year to old age; fetal
#'     and infant (< 1 year) cortex is untested. The wording follows
#'     \code{age_group} and is stored in \code{attr(result, "brain_caveat")}.
#'   \item Peripheral-brain concordance couplings rest on n = 27 living-brain
#'     donors; per-CpG estimates have a cross-resource reliability of r ~ 0.2-0.4.
#'     Grades exist precisely because single estimates cannot be trusted alone.
#'   \item Sex-chromosome and discontinued probes produce replicable artefactual
#'     concordance; \code{strict_qc} is the guard. This function never returns
#'     an X/Y probe as bridge-usable.
#' }
#'
#' @param cpgs CpG identifiers in any of the formats
#'   \code{cpg_expression_direction()} accepts.
#' @param tissue Which peripheral tissue you measured: "saliva" (default),
#'   "blood" or "buccal".
#' @param age_group \code{"adult"} (default), \code{"child"} or \code{"any"}:
#'   which age the caveat and the printed population refer to. It changes no
#'   direction and no rung of the ladder. \code{"child"} (for Project Alpha
#'   children and adolescents) and \code{"any"} also return
#'   \code{requested_gene_child_blood_helix} for blood-SMR fallback answers.
#' @param require_detected With \code{age_group = "child"} only: the blood-SMR
#'   fallback uses only pairs HELIX child blood detected (\code{"sig"} or
#'   \code{"nominal"}); measured evidence is never filtered. Ignored, with a
#'   warning, for any other \code{age_group}.
#' @param accuracy \code{"lower"} (default) or \code{"median"}: which single
#'   number \code{requested_gene_accuracy} carries for
#'   \code{requested_gene_best_accuracy}, the lower 95\% bound of the rung's
#'   validation in cortex or its point estimate (for a pair measured in both
#'   cortex cohorts: the floor set by their agreement, and the single-cohort
#'   accuracy it implies). Defaults to
#'   \code{getOption("cpgdirection.accuracy", "lower")}. The brain figures are
#'   the same-person-overlap-adjusted ones. See
#'   \code{\link{cpgd_accuracy_values}}.
#' @return A \code{data.table} of class \code{cpgd_bridge}: per CpG,
#'   \code{bridge_grade}, \code{bridge_usable} (concordance evidence in YOUR
#'   tissue), \code{rho_epic}/\code{rho_450k} (IMAGE-CpG, n = 27),
#'   \code{amaze_replicated}, \code{hannon_cortical_rho} (blood, n = 75),
#'   \code{sommerer_buccal_R} (n = 120), \code{ensemble_pctile} (genome-wide
#'   synthesised saliva score), \code{strict_qc}, the brain layer's
#'   \code{t1_genes}/\code{t1_directions}, and the measured cortex pairs:
#'   \code{measured_brain_genes}, \code{measured_brain_directions} and
#'   \code{measured_brain_n_cohorts} (\code{;}-separated, aligned) and
#'   \code{measured_brain_cohorts}. Where a gene was named in the input,
#'   \code{requested_gene}, \code{requested_gene_direction} and
#'   \code{requested_gene_tier} give the brain SMR answer for that gene, and
#'   \code{requested_gene_best_direction}, \code{requested_gene_best_evidence}
#'   and \code{requested_gene_best_accuracy} the answer of the brain-target
#'   ladder, with \code{requested_gene_accuracy} the same figure as one
#'   number (for a weighting engine:
#'   \code{requested_gene_best_direction * (2 * requested_gene_accuracy - 1)}).
#'   \code{attr(result, "brain_caveat")} holds the age caveat and
#'   \code{attr(result, "accuracy")} the value chosen.
#'
#' @section Co-regulation columns:
#'   \code{co_up} and \code{co_down} list T1 partner genes this CpG's
#'   methylation moves up and down respectively - genes within one column rise
#'   and fall together (causal co-direction, from the brain layer itself); a
#'   CpG with both columns filled is bivalent and its partners move in opposite
#'   directions. \code{hic_genes} lists genes contacted by the same
#'   adult-cortex enhancer the CpG sits in (PsychENCODE): candidate
#'   co-regulation by physical contact, with no direction attached - the
#'   evidence class that corroborated causal targets 2.5x over nearest genes.
#' @examplesIf cpgd_has_data("saliva_bridge_scores")
#' cpg_brain_bridge("cg06846259", tissue = "saliva")   # POMC, grade A++
#'
#' # a named gene, answered by the brain-target ladder, for a child cohort
#' r <- cpg_brain_bridge("cg06846259_POMC", age_group = "child")
#' as.data.frame(r)[, c("requested_gene", "requested_gene_best_direction",
#'                      "requested_gene_best_evidence")]
#' @export
cpg_brain_bridge <- function(cpgs, tissue = c("saliva", "blood", "buccal"),
                             age_group = c("adult", "child", "any"),
                             require_detected = FALSE,
                             accuracy = getOption("cpgdirection.accuracy", "lower")) {
  tissue <- match.arg(tissue)
  age_group <- match.arg(age_group)
  require_detected <- .cpgd_require_detected(age_group, require_detected)
  accuracy <- .cpgd_accuracy_arg(accuracy)
  q <- .cpgd_parse_input(cpgs)
  if (!nrow(q)) stop("No CpG identifiers recognised.", call. = FALSE)

  out <- data.table::data.table(cpg_id = unique(q$cpg_id))
  # Genes carried in the input (panel names like cg..._TC21_MC2R, or a gene
  # column in a data.frame) are honoured, not silently discarded: the
  # requested gene is echoed back with its own T1 answer, while t1_genes
  # still lists every partner so nothing is hidden.
  gq <- q[!duplicated(q$cpg_id), c("cpg_id", "given_gene", "given_gene2")]
  out <- merge(out, gq, by = "cpg_id", all.x = TRUE)

  ## brain layer: T1/T2 genes and directions per CpG
  B <- cpgd_brain_directions()
  bt <- B[B$tier == "T1"]
  agg <- bt[, list(t1_genes = paste(get("target_gene"), collapse = ";"),
                   t1_directions = paste(get("direction"), collapse = ";")),
            by = "cpg_id"]
  out <- merge(out, agg, by = "cpg_id", all.x = TRUE)
  t2 <- unique(B[B$tier == "T2"]$cpg_id)
  out[, "has_t2_direction" := out$cpg_id %in% t2]

  ## measured cortex pairs (the top rung for a brain question): every gene the
  ## CpG was measured against, signs and cohort counts aligned with the genes
  MB <- tryCatch(.cpgd_brain_measured_table(), error = function(e) NULL)
  mcols <- c("measured_brain_genes", "measured_brain_directions",
             "measured_brain_n_cohorts", "measured_brain_cohorts")
  if (!is.null(MB) && nrow(MB)) {
    mb <- MB[MB$cpg_id %chin% out$cpg_id]
    magg <- mb[, list(
      measured_brain_genes      = paste(get("target_gene"), collapse = ";"),
      measured_brain_directions = paste(get("direction"), collapse = ";"),
      measured_brain_n_cohorts  = paste(get("n_cohorts"), collapse = ";"),
      measured_brain_cohorts    = paste(sort(unique(unlist(strsplit(
        get("cohorts"), "|", fixed = TRUE)))), collapse = "|")),
      by = "cpg_id"]
    out <- merge(out, magg, by = "cpg_id", all.x = TRUE)
  } else {
    out[, (mcols) := NA_character_]
  }

  ## deliverable
  D <- cpgd_bridge_deliverable()
  ok_col   <- paste0(tissue, "_ok")
  am_raw   <- paste0(tissue, "_amaze_raw")
  am_adj   <- paste0(tissue, "_amaze_adj")
  rho_e    <- paste0("rho_", tissue, "_epic")
  rho_4    <- paste0("rho_", tissue, "_450k")   # buccal has no 450K column
  # NOT the deliverable's `genes`/`directions`: those duplicate t1_genes /
  # t1_directions and duplicated columns invite people to trust whichever they
  # noticed first.
  dcols <- intersect(c("cpg_id", "grade_v4", "strict_qc", ok_col, am_raw, am_adj,
                       rho_e, rho_4, "hannon_cortical_raw", "sommerer_buccal_R"),
                     names(D))
  out <- merge(out, D[, dcols, with = FALSE], by = "cpg_id", all.x = TRUE)
  data.table::setnames(out,
    old = intersect(c("grade_v4", ok_col, am_raw, am_adj, rho_e, rho_4,
                      "hannon_cortical_raw"), names(out)),
    new = intersect(c("bridge_grade", "bridge_usable", "amaze_replicated_raw",
                      "amaze_replicated_adj", "rho_epic", "rho_450k",
                      "hannon_cortical_rho"),
                    c("bridge_grade", "bridge_usable", "amaze_replicated_raw",
                      "amaze_replicated_adj", "rho_epic", "rho_450k",
                      "hannon_cortical_rho")[
                        c("grade_v4", ok_col, am_raw, am_adj, rho_e, rho_4,
                          "hannon_cortical_raw") %in% names(out)]),
    skip_absent = TRUE)

  ## genome-wide synthesised score (saliva-bridge; reported for every tissue,
  ## labelled for what it is)
  S <- .cpgd_data("saliva_bridge_scores", required = FALSE)
  if (!is.null(S)) {
    S <- S[, c("cpg_id", "ensemble_pctile"), with = FALSE]
    out <- merge(out, S, by = "cpg_id", all.x = TRUE)
  }

  ## co-regulation, where knowledge exists - two evidence classes, kept apart:
  ##
  ## (1) causal co-direction, from our own brain layer: genes moved by THIS
  ##     CpG's methylation in the same direction rise and fall together.
  ##     co_up = T1 partners at +1, co_down = T1 partners at -1. A CpG with
  ##     entries in both is bivalent - its partners move in OPPOSITE directions,
  ##     which is co-regulation but not co-activation, and listing them in one
  ##     column would assert the wrong thing.
  cu <- bt[bt$direction ==  1L,
           list(co_up   = paste(get("target_gene"), collapse = ";")), by = "cpg_id"]
  cd <- bt[bt$direction == -1L,
           list(co_down = paste(get("target_gene"), collapse = ";")), by = "cpg_id"]
  out <- merge(out, cu, by = "cpg_id", all.x = TRUE)
  out <- merge(out, cd, by = "cpg_id", all.x = TRUE)

  ## (2) physical co-targets, from adult-cortex Hi-C: genes contacted by the
  ##     same enhancer the CpG sits in. Candidate co-regulation only - contact
  ##     is not direction - but this evidence class corroborated the causal
  ##     target 2.5x over the nearest gene, so it earns a labelled column.
  HC <- .cpgd_data("hic_cotargets", required = FALSE)
  if (!is.null(HC)) {
    out <- merge(out, HC, by = "cpg_id", all.x = TRUE)
  }

  ## if the input named a gene, answer for that gene specifically
  has_req <- !all(is.na(out$given_gene)) || !all(is.na(out$given_gene2))
  if (has_req) {
    B1 <- cpgd_brain_directions()
    B1 <- B1[B1$tier %in% c("T1", "T2")]
    key <- data.table::data.table(
      cpg_id = out$cpg_id,
      g1 = toupper(out$given_gene), g2 = toupper(out$given_gene2))
    hit <- function(gcol) {
      m <- merge(key[, c("cpg_id", gcol), with = FALSE],
                 B1[, list(cpg_id = get("cpg_id"),
                           g = toupper(get("target_gene")),
                           d = get("direction"), tr = get("tier"))],
                 by.x = c("cpg_id", gcol), by.y = c("cpg_id", "g"))
      m[!duplicated(m$cpg_id)]
    }
    h <- hit("g1"); h2 <- hit("g2")
    out[, "requested_gene" := data.table::fifelse(!is.na(out$given_gene),
                                                  out$given_gene, out$given_gene2)]
    idx1 <- match(out$cpg_id, h$cpg_id);  idx2 <- match(out$cpg_id, h2$cpg_id)
    out[, "requested_gene_direction" := data.table::fifelse(
          !is.na(idx1), h$d[idx1],
          data.table::fifelse(!is.na(idx2), h2$d[idx2], NA_integer_))]
    out[, "requested_gene_tier" := data.table::fifelse(
          !is.na(idx1), h$tr[idx1],
          data.table::fifelse(!is.na(idx2), h2$tr[idx2], NA_character_))]

    ## the brain-target ladder for the named gene: measured brain > brain T1
    ## > brain T2 > measured blood > blood SMR. The brain SMR rungs are taken
    ## from requested_gene_tier above, so the two answers can never disagree
    ## about what the brain layer holds.
    lad <- .cpgd_brain_target_ladder(
      out$cpg_id, toupper(out$given_gene), toupper(out$given_gene2),
      brain_tier = out$requested_gene_tier,
      brain_dir = out$requested_gene_direction,
      MB = MB, age_group = age_group, require_detected = require_detected)
    out[, "requested_gene_best_direction" := lad$d]
    out[, "requested_gene_best_evidence"  := lad$ev]
    out[, "requested_gene_best_accuracy"  := lad$acc]
    data.table::set(out, j = "requested_gene_accuracy",
                    value = .cpgd_accuracy_number(lad$acc, accuracy))
    if (!identical(age_group, "adult"))
      out[, "requested_gene_child_blood_helix" := lad$cb]
  }
  out[, c("given_gene", "given_gene2") := NULL]

  data.table::setattr(out, "tissue", tissue)
  data.table::setattr(out, "age_group", age_group)
  data.table::setattr(out, "accuracy", accuracy)
  data.table::setattr(out, "brain_caveat", CPGD_BRAIN_CAVEAT[[age_group]])
  data.table::setattr(out, "class", c("cpgd_bridge", class(out)))
  out[]
}

# The brain-target ladder for one named gene per row. `g1`/`g2` are the two
# spellings the parser offers for the requested gene (the whole tail and its
# last token); each source is searched for g1, then g2, before the next,
# weaker source is consulted. Genes match case-insensitively and with
# dash/underscore equivalence, as in the pair API.
.cpgd_brain_target_ladder <- function(cpg_id, g1, g2, brain_tier, brain_dir,
                                      MB = NULL, age_group = "adult",
                                      require_detected = FALSE) {
  n <- length(cpg_id)
  gk <- function(x) gsub("_", "-", toupper(x))
  k1 <- gk(g1); k2 <- gk(g2)
  d  <- rep(NA_integer_, n); ev <- rep(NA_character_, n)
  acc <- rep(NA_character_, n); cb <- rep(NA_character_, n)

  # look up a keyed source table for whichever spelling hits first
  look <- function(Tb, need) {
    res <- list(d = rep(NA_real_, n), a = rep(NA_character_, n),
                cb = rep(NA_character_, n))
    if (is.null(Tb) || !nrow(Tb)) return(res)
    for (k in list(k1, k2)) {
      want <- need & is.na(res$d) & !is.na(k)
      if (!any(want)) next
      # built outside `[`, where `cpg_id` would resolve to Tb's own column
      qk <- data.table::data.table(cpg_id = cpg_id, g = k)
      h <- Tb[qk, on = c("cpg_id", "g")]
      take <- want & !is.na(h$d)
      res$d[take] <- h$d[take]; res$a[take] <- h$a[take]; res$cb[take] <- h$cb[take]
    }
    res
  }
  keyed <- function(cp, g, dd, a, cbv = NA_character_) {
    Tb <- data.table::data.table(cpg_id = as.character(cp), g = gk(g),
                                 d = as.numeric(dd), a = a, cb = cbv)
    Tb <- Tb[!is.na(get("d"))]
    Tb <- unique(Tb, by = c("cpg_id", "g"))
    data.table::setkeyv(Tb, c("cpg_id", "g"))
    Tb
  }

  # 1. measured in cortex
  if (!is.null(MB) && nrow(MB)) {
    rep_ <- MB$n_cohorts >= 2L & MB$cohort_agreement >= 0.99
    a <- ifelse(rep_, CPGD_BRAIN_TARGET_ACCURACY[["measured_brain_replicated"]],
                CPGD_BRAIN_TARGET_ACCURACY[["measured_brain_single"]])
    r <- look(keyed(MB$cpg_id, MB$target_gene, MB$direction, a), rep(TRUE, n))
    i <- which(!is.na(r$d))
    d[i] <- as.integer(r$d[i]); ev[i] <- "measured_brain"; acc[i] <- r$a[i]
  }
  # 2-3. brain SMR, exactly as requested_gene_tier reports it
  for (tr in c("T1", "T2")) {
    i <- which(is.na(ev) & brain_tier %in% tr & !is.na(brain_dir))
    d[i] <- as.integer(brain_dir[i]); ev[i] <- paste0("brain_", tr)
    acc[i] <- CPGD_BRAIN_TARGET_ACCURACY[[paste0("brain_", tr)]]
  }
  # 4. measured blood (fallback)
  if (any(is.na(ev))) {
    M <- tryCatch(cpgd_measured_eqtms(), error = function(e) NULL)
    if (!is.null(M) && nrow(M)) {
      # computed outside `[`, where `cpg_id` would resolve to M's own column
      keep <- M$tissue == "blood" & M$cpg_id %chin% cpg_id
      M <- M[keep]
      r <- look(keyed(M$cpg_id, M$target_gene, sign(as.numeric(M$direction)),
                      CPGD_BRAIN_TARGET_ACCURACY[["measured_blood"]]), is.na(ev))
      i <- which(is.na(ev) & !is.na(r$d))
      d[i] <- as.integer(r$d[i]); ev[i] <- "measured_blood"; acc[i] <- r$a[i]
    }
  }
  # 5. blood SMR (fallback; the reported rung never supplies a sign)
  if (any(is.na(ev))) {
    S <- tryCatch(.cpgd_smr_table(), error = function(e) NULL)
    if (!is.null(S) && nrow(S)) {
      require_detected <- .cpgd_check_child_column(S, age_group, require_detected)
      rung <- .cpgd_smr_rung(S)
      cbs <- .cpgd_smr_child(S)
      ok <- rung %in% c("replicated", "high", "moderate", "weak") &
            S$cpg_id %chin% cpg_id
      if (isTRUE(require_detected)) {
        det <- cbs %in% CPGD_CHILD_DETECTED
        message(sprintf(paste0("require_detected: %s blood-SMR pair(s) for these CpGs ",
                               "dropped from the fallback (child_blood_helix not ",
                               "sig/nominal); measured evidence unaffected."),
                        format(sum(ok & !det), big.mark = ",")))
        ok <- ok & det
      }
      key <- paste0("blood_SMR_", rung[ok])
      # most accurate in cortex first, so unique() keeps it when spellings
      # collide: replicated 0.772 > high 0.739 > weak 0.725 > moderate 0.666
      # against ROSMAP (tools/reports/validate_ladder_transport.md)
      o <- order(match(rung[ok], CPGD_BRAIN_SMR_ORDER))
      a <- unname(CPGD_BRAIN_TARGET_ACCURACY[key])
      Tb <- keyed(S$cpg_id[ok][o], S$target_gene[ok][o],
                  as.numeric(S$direction[ok][o]), a[o], cbs[ok][o])
      r <- look(Tb, is.na(ev))
      i <- which(is.na(ev) & !is.na(r$d))
      d[i] <- as.integer(r$d[i]); ev[i] <- "blood_SMR"; acc[i] <- r$a[i]
      cb[i] <- r$cb[i]
    }
  }
  data.table::data.table(d = d, ev = ev, acc = acc, cb = cb)
}

#' @export
print.cpgd_bridge <- function(x, ...) {
  # A column subset keeps the class but not the columns the card reads; show
  # it as the plain table it now is.
  core <- c("cpg_id", "bridge_grade", "bridge_usable", "strict_qc", "t1_genes",
            "has_t2_direction", "rho_epic", "hannon_cortical_rho", "co_up", "co_down")
  if (!all(core %in% names(x))) {
    return(print(data.table::as.data.table(unclass(x)), ...))
  }
  tis <- attr(x, "tissue")
  n <- nrow(x)
  arrow <- function(d) if (is.na(d)) "" else if (as.integer(d) > 0)
    "raises" else "lowers"

  if (n <= 10) {
    ## card layout: the answer first, evidence beneath it, one CpG at a time
    for (i in seq_len(n)) {
      r <- x[i]
      head_bits <- c(r$cpg_id,
        if (!is.na(r$bridge_grade)) sprintf("grade %s", r$bridge_grade) else "not in the graded set",
        if (isTRUE(r$strict_qc %in% c(TRUE, "TRUE"))) "strict QC"
        else if (!is.na(r$bridge_grade)) "PROVISIONAL: discontinued probe" else NULL,
        if (isTRUE(r$bridge_usable %in% c(TRUE, "TRUE"))) sprintf("usable from %s", tis)
        else if (!is.na(r$bridge_grade)) sprintf("NOT robust from %s", tis) else NULL)
      cat(paste(head_bits, collapse = "  |  "), "\n")

      if ("requested_gene" %in% names(x) && !is.na(r$requested_gene)) {
        bev <- if ("requested_gene_best_evidence" %in% names(x))
          r$requested_gene_best_evidence else NA_character_
        if (!is.na(bev) && bev == "measured_brain")
          cat(sprintf("  %s: measured in cortex (best brain evidence): higher methylation here %s %s expression\n",
                      r$requested_gene, arrow(r$requested_gene_best_direction),
                      r$requested_gene))
        if (!is.na(r$requested_gene_direction))
          cat(sprintf("  %s: in brain, higher methylation here %s %s expression (%s)\n",
                      r$requested_gene, arrow(r$requested_gene_direction),
                      r$requested_gene, r$requested_gene_tier))
        else if (!is.na(bev) && bev %in% c("measured_blood", "blood_SMR"))
          cat(sprintf(paste0("  %s: no brain evidence; BLOOD fallback (%s): higher ",
                             "methylation here %s %s expression - %s\n"),
                      r$requested_gene, bev, arrow(r$requested_gene_best_direction),
                      r$requested_gene, r$requested_gene_best_accuracy))
        else if (is.na(bev) || bev != "measured_brain")
          cat(sprintf("  %s: not a brain partner of this CpG (no substitution made)\n",
                      r$requested_gene))
      }
      if (!is.na(r$t1_genes)) {
        # One plain statement per direction, never a sign convention the reader
        # has to decode: "-POMC" cost a question; a sentence would not have.
        g <- strsplit(r$t1_genes, ";", fixed = TRUE)[[1]]
        d <- strsplit(r$t1_directions, ";", fixed = TRUE)[[1]]
        lo <- g[d != "1"]; hi <- g[d == "1"]
        parts <- c(
          if (length(lo)) sprintf("lowers %s expression", paste(lo, collapse = ", ")),
          if (length(hi)) sprintf("raises %s expression", paste(hi, collapse = ", ")))
        cat(sprintf("  brain (T1): higher methylation here %s\n",
                    paste(parts, collapse = " and ")))
      } else if (isTRUE(r$has_t2_direction)) {
        cat("  brain: T2 evidence only (see cpgd_brain_directions())\n")
      } else cat("  brain: no tiered causal direction\n")
      if ("measured_brain_genes" %in% names(x) && !is.na(r$measured_brain_genes)) {
        g <- strsplit(r$measured_brain_genes, ";", fixed = TRUE)[[1]]
        d <- strsplit(r$measured_brain_directions, ";", fixed = TRUE)[[1]]
        lo <- g[d != "1"]; hi <- g[d == "1"]
        parts <- c(
          if (length(lo)) sprintf("lowers %s expression", paste(lo, collapse = ", ")),
          if (length(hi)) sprintf("raises %s expression", paste(hi, collapse = ", ")))
        cat(sprintf("  brain (measured, %s): higher methylation here %s\n",
                    gsub("|", ", ", r$measured_brain_cohorts, fixed = TRUE),
                    paste(parts, collapse = " and ")))
      }

      cw <- character(0)
      if (!is.na(r$rho_epic))
        cw <- c(cw, sprintf("%s rho %.2f/%.2f (two arrays, n=27)", tis,
                            r$rho_epic,
                            if ("rho_450k" %in% names(x) && !is.na(r$rho_450k)) r$rho_450k else NA))
      if ("amaze_replicated_raw" %in% names(x) && !is.na(r$bridge_grade))
        cw <- c(cw, sprintf("AMAZE %s",
                 if (isTRUE(r$amaze_replicated_raw %in% c(TRUE,"TRUE")) ||
                     isTRUE(r$amaze_replicated_adj %in% c(TRUE,"TRUE"))) "replicated" else "not replicated"))
      if (!is.na(r$hannon_cortical_rho))
        cw <- c(cw, sprintf("blood-brain %.2f (n=75)", r$hannon_cortical_rho))
      if ("sommerer_buccal_R" %in% names(x) && !is.na(r$sommerer_buccal_R))
        cw <- c(cw, sprintf("buccal-brain %.2f (n=120)", r$sommerer_buccal_R))
      if (length(cw)) cat("  concordance: ", paste(cw, collapse = " | "), "\n", sep = "")
      if ("ensemble_pctile" %in% names(x) && !is.na(r$ensemble_pctile))
        cat(sprintf("  synthesised bridge score: %.1f percentile\n", r$ensemble_pctile))
      co <- character(0)
      if (!is.na(r$co_up))   co <- c(co, paste0("up-together: ", r$co_up))
      if (!is.na(r$co_down)) co <- c(co, paste0("down-together: ", r$co_down))
      if (length(co) == 2) co <- c(co, "(bivalent CpG: the two sets move oppositely)")
      if (length(co)) cat("  co-regulated: ", paste(co, collapse = "  "), "\n", sep = "")
      if ("hic_genes" %in% names(x) && !is.na(r$hic_genes))
        cat("  same-enhancer contacts (Hi-C, no direction): ", r$hic_genes, "\n", sep = "")
      cat("\n")
    }
    cat("Caveats: brain evidence first (measured cortex > T1 > T2); blood is a\n")
    cat("labelled fallback (0.84 measured, 0.67-0.77 SMR against adult cortex);\n")
    cat("couplings rest on n = 27 donors - trust grades over any single rho.\n")
    .cpgd_bridge_age_lines(x)
    cat("Full columns: as.data.frame(x).\n")
    return(invisible(x))
  }

  ## many CpGs: summary + compact table, no 20-column dump
  cat(sprintf("cpg_brain_bridge: %d CpGs, measured tissue = %s\n", n, tis))
  ng <- sum(!is.na(x$bridge_grade))
  if (ng) {
    tb <- table(x$bridge_grade)
    cat(sprintf("  graded deliverable: %d  (%s)\n", ng,
                paste(sprintf("%s %d", names(tb), tb), collapse = ", ")))
  } else cat("  graded deliverable: 0\n")
  cat(sprintf("  usable from %s: %d   T1 brain direction: %d   T2 only: %d\n",
              tis, sum(x$bridge_usable %in% c(TRUE, "TRUE")),
              sum(!is.na(x$t1_genes)),
              sum(x$has_t2_direction & is.na(x$t1_genes))))
  if ("measured_brain_genes" %in% names(x))
    cat(sprintf("  measured in cortex: %d\n", sum(!is.na(x$measured_brain_genes))))
  if ("requested_gene_best_evidence" %in% names(x)) {
    tb <- table(factor(x$requested_gene_best_evidence, levels = CPGD_BRAIN_LADDER))
    cat(sprintf("  named genes answered by: %s\n",
                paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")))
  }
  if ("ensemble_pctile" %in% names(x)) {
    k <- sum(is.na(x$bridge_grade) & x$ensemble_pctile >= 95, na.rm = TRUE)
    cat(sprintf("  ungraded but top-5%% synthesised bridge score: %d\n", k))
  }
  show <- intersect(c("cpg_id", "bridge_grade", "bridge_usable", "t1_genes",
                      "t1_directions", "ensemble_pctile", "strict_qc"), names(x))
  cat("\n")
  print(data.table::as.data.table(x)[, show, with = FALSE], nrows = 20)
  cat("\n")
  .cpgd_bridge_age_lines(x)
  cat("Full columns: as.data.frame(x). Caveats: see ?cpg_brain_bridge.\n")
  invisible(x)
}

# The population the bridge's figures describe, and its age caveat.
.cpgd_bridge_age_lines <- function(x) {
  ag <- attr(x, "age_group")
  cav <- attr(x, "brain_caveat")
  if (!is.null(cav)) {
    lab <- if (is.null(ag)) "" else sprintf("Age group %s: ", ag)
    cat(paste(strwrap(paste0(lab, cav), width = 72), collapse = "\n"), "\n", sep = "")
  }
  invisible(NULL)
}
