#' Predict the direction of a CpG's effect on gene expression
#'
#' The main entry point. By default it runs every method the package has and
#' returns one row per CpG carrying all of them, plus a single
#' \code{best_direction} column holding the answer from the strongest evidence
#' available for that CpG.
#'
#' Every source that holds a sign for the CpG and the requested gene is a
#' candidate: a measured eQTM (by the tissue it was measured in, one study or
#' replicated), the SMR causal layer (by rung and stratum), the catalogue
#' models (their consensus, or one tissue's model alone, by tier), the tumour
#' consensus, and the distance curves. The most accurate candidate in the
#' call's tissue answers.
#'
#' @section The ladder is ranked by accuracy, per tissue:
#' The order is not fixed. \code{tools/validate_rung_accuracy.R} scores every
#' source in each tissue against a measured yardstick independent of it (adult
#' blood: MESA monocytes and T cells; nasal epithelium: EVA-PR; solid tissue:
#' FUSION and GTEx, each held out in turn), and the sources are ranked there by
#' the lower bound of the 95\% interval of that agreement
#' (\code{\link{cpgd_accuracy_values}}) -- the number a call returns by
#' default, so the answering source always carries the highest weight
#' available and a small sample cannot outrank a large one on a noisy point.
#' Nothing else enters the order: no tier, stratum or parent placement. The
#' tissue is \code{target_tissue}, else blood. A source whose point estimate
#' there is below 0.5, or that no yardstick scores at all, never answers in
#' that tissue. A catalogue
#' consensus is between predicted calls on the same gene (a DIRECT_eQTM row is
#' a measured sign and answers through the measured rung of its tissue);
#' models naming different genes are single calls about their own genes.
#' Where predicted calls on the same gene disagree, no source ranked at or
#' below the distance curves settles it (\code{tissue_conflict}). Age never
#' changes the order. In adult blood replicated SMR S1 ranks above a single
#' measured blood study, SMR moderate just above the catalogue consensus, and
#' single calls of the nasal and solid-tissue models fall below the distance
#' curves.
#' \code{best_rung} names the source that answered.
#'
#' Passing a specific \code{tissue} instead returns the single-tissue result,
#' with the column set documented under \emph{Value} in earlier versions.
#'
#' \strong{One row per CpG, by design.} \code{cpg_expression_direction()} is
#' the convenience one-result-per-CpG interface: it returns one requested/best
#' target per CpG. For automatic discovery of every supported target/co-effect
#' and pair-specific direction evidence, use \code{\link{cpg_gene_pairs}} --
#' a CpG can have several target genes, with legitimately opposite directions,
#' and this interface by construction shows only one of them. DMSA and any
#' other many-to-many consumer must use the pair API.
#'
#' @param cpgs Character vector of CpG identifiers, with or without suffixes.
#'   A single file path to a text or CSV file is also accepted.
#' @param genes Optional character vector of target genes, the same length as
#'   \code{cpgs}. Where identifiers already carry the symbol, as in
#'   \code{"cg02079741_TC21_POMC"}, it is parsed automatically.
#' @param tissue \code{"all"} (default) runs everything and adds the
#'   \code{best_*} columns. \code{"blood"}, \code{"nasal_epithelium"} or
#'   \code{"solid_tissue"} returns that tissue alone.
#' @param min_tissues Tissues that must agree before
#'   \code{consensus_direction} is reported. Default 2. In the ladder a
#'   consensus always needs two models at least: with \code{min_tissues = 1} a
#'   lone model's call is still ranked and quoted as a single call.
#' @param min_confidence Optional numeric in \code{[0, 1]} overriding the
#'   tier-specific abstention thresholds.
#' @param allow_nearest_gene If a gene is requested but no pairing for it is
#'   held, the nearest catalogued gene is a \emph{different} gene. By default
#'   those rows abstain. See \code{\link{cpgd_accuracy_table}}.
#' @param min_distance_info Floor on the pooled distance probability, as
#'   \code{|p - 0.5| * 2}, below which a unanimous distance call is withheld.
#'   Default 0.02. It is a floor, not the main gate: whether a distance call is
#'   made at all is decided by whether the three tissue curves agree in sign.
#' @param universal Include the distance-only fallback. Default \code{TRUE}.
#'   Requires the Bioconductor annotation packages, and builds a TSS table on
#'   first use. Set \code{FALSE} to skip it.
#' @param target_tissue Optional, one of \code{"blood"},
#'   \code{"nasal_epithelium"} or \code{"solid_tissue"}: the tissue the
#'   question is about. It sets the order of the ladder and the figures each
#'   call carries (blood when not given). Where the tissues
#'   disagree, it also resolves using this tissue instead of abstaining. Declare it on
#'   biological grounds \strong{before} looking at the votes; the resulting rows
#'   are worth that tissue's accuracy, not the consensus figure, and are still
#'   flagged in \code{direction_uncertain}. Setting it also fills any row that
#'   would otherwise abstain with that tissue's distance-curve sign, labelled
#'   \code{targeted_last_resort} in \code{best_evidence} and carrying that
#'   curve's validated figure (\code{best_rung} = \code{distance_curve_<tissue>}),
#'   so \code{best_direction} comes back complete.
#' @param html_out Optional. \code{TRUE} writes
#'   \code{cpgdirection_results.html} into \code{tempdir()} and reports the
#'   path; a string is used as the file name (a relative name also lands in
#'   \code{tempdir()}; give an absolute path to keep the report). See
#'   \code{\link{cpgd_report}}.
#' @param verbose Emit a summary with \code{message()}. Default \code{TRUE}.
#' @param age_group \code{"adult"} (default), \code{"child"} or \code{"any"}:
#'   which population the blood accuracy figures in
#'   \code{best_expected_accuracy} refer to. \code{"child"} quotes HELIX child
#'   blood (6-11 y) for every source HELIX scores on 100+ pairs (measured,
#'   SMR, tumour, the nasal model's calls; the adult figure otherwise) and adds
#'   \code{smr_child_blood_helix}; \code{"any"} quotes both and adds the
#'   column. Directions and rungs never change: age did not
#'   move the ladder in validation (see \code{\link{cpgd_smr_directions}}).
#'   \code{"adult"} equals the default. Project Alpha users
#'   (children and adolescents) should set \code{"child"}. Only
#'   \code{tissue = "all"} carries the SMR layer; with a single tissue
#'   \code{age_group} and \code{require_detected} have no effect.
#' @param require_detected With \code{age_group = "child"} only: drop SMR
#'   pairs whose \code{child_blood_helix} is not \code{"sig"} or
#'   \code{"nominal"} before the ladder runs; measured, catalogue and distance
#'   evidence are untouched and a message reports how many pairs were dropped.
#'   Ignored, with a warning, for any other \code{age_group}.
#' @param accuracy Which single number \code{best_accuracy} (and, with one
#'   tissue, \code{accuracy}) carries for each accuracy string:
#'   \code{"lower"} (default), the lower bound of the 95\% interval of the
#'   source's validation in the tissue, never below 0.5, or
#'   \code{"median"}, its point estimate. The default
#'   comes from \code{getOption("cpgdirection.accuracy", "lower")}, so a
#'   session can switch it for every call, including those an engine such as
#'   DMSA makes. See \code{\link{cpgd_accuracy_values}}.
#'
#' @section One number for weighting engines:
#'   \code{best_expected_accuracy} quotes a validation in words.
#'   \code{best_accuracy} carries the same validation as one number, the
#'   probability that \code{best_direction} is right given that the pair is
#'   coupled: by default the lower bound of its 95\% interval, never below
#'   0.5, with \code{accuracy = "median"} the point estimate. An engine that
#'   weights by
#'   sign certainty can use \eqn{P(d = +1)} = \code{best_accuracy} where
#'   \code{best_direction} is \code{+1} and \code{1 - best_accuracy} where
#'   it is \code{-1}, i.e. the multiplier
#'   \code{best_direction * (2 * best_accuracy - 1)}. Every row with a
#'   direction has one; it is \code{NA} only where the call abstains. A
#'   source no yardstick scores carries the neutral 0.5. With
#'   \code{age_group = "any"} it
#'   is the smaller of the adult and the child figure. Every value and its
#'   source: \code{\link{cpgd_accuracy_values}}.
#'
#' @return With \code{tissue = "all"}, a \code{data.table} of class
#'   \code{cpgd_full}, one row per unique CpG:
#'   \describe{
#'     \item{best_direction}{+1, -1 or \code{NA}: the answer from the strongest
#'       evidence available}
#'     \item{best_evidence}{which kind of source it came from — \code{measured},
#'       \code{smr_replicated}, \code{smr_high}, \code{catalogue_consensus},
#'       \code{smr_moderate},
#'       \code{catalogue_single}, \code{onco_consensus}, \code{smr_weak},
#'       \code{distance_only}, \code{distance_tissue_conflict},
#'       \code{distance_uninformative}, \code{tissue_conflict} or
#'       \code{no_evidence}}
#'     \item{best_rung}{the exact source that answered, as ranked (e.g.
#'       \code{measured_blood_single}, \code{smr_replicated_S2},
#'       \code{catalogue_single_nasal_A}); the state for an abstention;
#'       \code{attr(result, "rank_tissue")} names the tissue whose order was
#'       used}
#'     \item{best_confidence, best_expected_accuracy}{what that answer is worth}
#'     \item{best_accuracy}{the same figure as one number (see the section
#'       One number for weighting engines); \code{attr(result, "accuracy")}
#'       records which value was chosen}
#'     \item{consensus_direction, n_tissues_calling, tissue_agreement}{the
#'       catalogue layer}
#'     \item{gene_*, dir_*, conf_*, tier_*, status_*}{per tissue}
#'     \item{gene_universal, dist_universal, dir_universal, p_universal}{the
#'       distance layer}
#'     \item{smr_direction, smr_tier, smr_p, smr_n_instruments, smr_gene_dist,
#'       smr_agreement}{the Mendelian randomisation layer}
#'     \item{smr_gene}{the gene the SMR direction refers to, never left to
#'       inference}
#'     \item{smr_gene_match}{whether \code{smr_gene} is the same gene the rest of
#'       the row concerns. Only matching rows may supply \code{best_direction};
#'       non-matching rows are reported in full but describe a different pair.
#'       Only 58.6\% of CpGs present in both the SMR table and a catalogue share
#'       even one gene between them, so this is the common case rather than an
#'       edge case}
#'     \item{smr_heidi_status, smr_p_heidi}{the HEIDI heterogeneity test:
#'       \code{"pass"}, \code{"fail"} or \code{"not_tested"} (fewer than three
#'       instruments survived LD pruning). Reported, never acted on, for the
#'       eQTLGen and the FinnGen rows alike: against eQTLGen HEIDI did not
#'       separate concordant from discordant directions, and FinnGen-only pairs
#'       are admitted on the SMR test at the top cis-mQTL instead. See
#'       \code{?cpgd_smr_directions}}
#'     \item{smr_rung, smr_replication, smr_outcome}{how the SMR row enters the
#'       ladder: \code{smr_rung} is \code{replicated} (same sign in eQTLGen and
#'       FinnGen PBMC), \code{high}, \code{moderate}, \code{weak} or
#'       \code{reported} (never supplies \code{best_direction});
#'       \code{smr_replication} says why; \code{smr_outcome} names the
#'       outcome study behind the sign}
#'     \item{smr_child_blood_helix}{with \code{age_group = "child"} or
#'       \code{"any"} only: whether HELIX child blood detected the SMR pair
#'       (\code{sig}, \code{nominal}, \code{null}, \code{untested}).
#'       Informational; \code{null} means not detected at n = 823, not absent}
#'     \item{smr_in_table}{whether the CpG has any SMR evidence at all, for any
#'       gene. Distinguishes "no instrument exists" from "the instrument points
#'       at a different gene", which are different facts with different remedies}
#'     \item{onco_direction, onco_tier, onco_n_cancers, onco_median_r,
#'       onco_gene, onco_gene_match, onco_agreement}{the tumour consensus layer
#'       (Onco-eQTM, 27 TCGA cancer types; \code{\link{cpgd_onco_eqtm}}),
#'       keyed and gated exactly like SMR: only a gene match may reach
#'       \code{best_direction}, and only tier O1 does so (as
#'       \code{onco_consensus}, ranked by its accuracy like every source).
#'       \code{onco_agreement} compares the tumour sign with a peripheral
#'       \code{best_direction} for the same pair; \code{NA} when the tumour
#'       sign is itself the best direction}
#'     \item{measured_genes, annotation_mismatch, mismatch_note}{measured
#'       evidence, and whether it concerns a gene other than the one requested}
#'     \item{measured_n_cohorts, measured_cohort_agreement}{for a measured
#'       call, how many cohorts measured the pair (strata and cell types of one
#'       study counted separately) and the share agreeing on its sign (2.99.8,
#'       eQTM Atlas build). The ladder ranks a measured sign by the tissue it
#'       was measured in and whether 2+ independent studies (donor groups)
#'       there agree}
#'   }
#'   With a single \code{tissue}, the per-tissue result instead: \code{status},
#'   \code{direction}, \code{probability_plus1}, \code{confidence},
#'   \code{evidence_tier}, \code{expected_accuracy} and \code{accuracy} (the
#'   same figure as one number), \code{call} and \code{note}.
#'
#' @section The two alternative codings:
#' Where the evidence leaves a sign nearly free to flip — a distance prior at
#' \code{p = 0.5009}, say — \code{best_direction} abstains. That is honest but
#' leaves nothing to analyse, so two further columns bracket the ambiguity:
#' \code{best_direction_filled} takes the point estimate for those rows and
#' \code{best_direction_flipped} takes its opposite. \code{direction_uncertain}
#' flags which rows differ.
#'
#' \strong{These are a sensitivity analysis, not a menu.} The question they
#' answer is "does my downstream result survive both codings", not "which coding
#' gives me a result". Picking whichever column fits the data better is
#' selecting on the outcome: it doubles the effective number of tests and biases
#' towards significance. If a finding holds under both, the ambiguity does not
#' matter. If it does not, the finding rests on coin flips and should be
#' reported as unresolved.
#'
#' @section Why the distance layer reports conflicts:
#' The three tissue curves are not parallel. Blood asymptotes at 0.449 and never
#' crosses 0.5 at any distance; nasal epithelium reaches 0.594 and solid tissue
#' 0.861. Beyond roughly 20 kb they therefore disagree, often sharply.
#'
#' Averaging them produces a pooled probability near 0.50 that looks like
#' ignorance but is nothing of the kind: blood is saying \code{-1} and solid
#' tissue \code{+1}, both with confidence. Those rows are reported as
#' \code{distance_tissue_conflict}, with \code{dist_dir_blood},
#' \code{dist_dir_nasal} and \code{dist_dir_solid} exposed so the tissue
#' relevant to your question can be read directly. A distance call enters
#' \code{best_direction} only when all three curves agree.
#'
#' \strong{A structural consequence worth knowing.} The blood curve peaks at
#' 0.449. Because a distance call requires unanimity, and blood never reaches
#' 0.5, \code{distance_only} can never return \code{+1} for any CpG at any
#' distance. Every direction that layer contributes is \code{-1} by
#' construction. This is a property of the fitted curves, not of your data, and
#' it is the reason a panel dominated by distance calls will look uniformly
#' negative. Treat \code{distance_only} as a one-way instrument and read
#' \code{dist_dir_solid} where the positive signal would be.
#'
#' @section Declaring a target tissue:
#' Blood holds an effective veto over distance calls, because unanimity is
#' required and blood never votes \code{+1}. That is correct when the question
#' is about blood and wrong when it is not. A saliva-derived panel of
#' brain-expressed genes is asking about brain, and blood's opinion should not
#' silence it.
#'
#' \code{target_tissue} resolves conflicts using a named tissue. It is tissue
#' \strong{selection}, not voting: name the tissue first, on the biology, then
#' read what it says. Promoting whichever tissue happens to break the tie the way
#' you hoped is selecting on the outcome, and a majority rule tuned after seeing
#' the votes is the same thing wearing a better hat.
#'
#' The price is explicit. Rows resolved this way carry that tissue's own
#' figure (\code{catalogue_model_<tissue>_<tier>} or
#' \code{distance_curve_<tissue>} in \code{best_rung};
#' \code{\link{cpgd_accuracy_values}}). Tier-A calls of the tissue's model:
#' blood 0.75 against MESA, solid tissue 0.71 against FUSION and GTEx; the
#' nasal model has no independent nasal yardstick, so it quotes its
#' cross-tissue figure (0.69). The tissue's distance curve: blood 0.64, solid
#' tissue 0.54 (a tumour curve), nasal 0.60 (cross-tissue). The solid-tissue
#' model is trained on tumour, provisional, and its brain and kidney arms did
#' not exceed their own majority baselines under leave-one-cohort-out. If solid
#' tissue is your brain proxy, that last clause is the one to quote in your
#' limitations.
#'
#' @section How to read best_direction:
#' It is the direction supported by the best evidence this package holds for
#' that CpG, not a claim about any one tissue and not a substitute for measured
#' RNA in your own samples. Always read \code{best_rung} and
#' \code{best_accuracy} alongside it: in adult blood a replicated measured
#' sign agrees with an independent cohort 0.97 of the time and a
#' \code{distance_only} call 0.67.
#'
#' \code{tissue_conflict} means the tissues returned opposite directions. That
#' is reported as \code{NA} rather than resolved, because a weaker method should
#' not break a tie between two stronger ones.
#'
#' \strong{Everything here is conditional on the CpG being an eQTM at all.} The
#' catalogues behind the models contain only significant associations, so no
#' output can tell you whether an effect exists. \code{no_evidence} is a correct
#' answer.
#'
#' @examplesIf cpgd_has_data("lookup_blood_hg19")
#' res <- cpg_expression_direction(c("cg02079741_TC21_POMC", "cg00000029"),
#'                                 universal = FALSE, verbose = FALSE)
#' res[, c("cpg_id", "best_direction", "best_evidence",
#'         "best_expected_accuracy"), with = FALSE]
#'
#' # a single tissue returns the per-tissue column set instead
#' one <- cpg_expression_direction("cg00000029", tissue = "blood",
#'                                 verbose = FALSE)
#' one[, c("cpg_id", "target_gene", "direction", "evidence_tier", "call"),
#'     with = FALSE]
#' @export
cpg_expression_direction <- function(cpgs,
                                     genes = NULL,
                                     tissue = "all",
                                     min_tissues = 2L,
                                     min_confidence = NULL,
                                     allow_nearest_gene = FALSE,
                                     universal = TRUE,
                                     min_distance_info = 0.02,
                                     target_tissue = NULL,
                                     html_out = NULL,
                                     verbose = TRUE,
                                     age_group = c("adult", "child", "any"),
                                     require_detected = FALSE,
                                     accuracy = getOption("cpgdirection.accuracy", "lower")) {

  age_group <- match.arg(age_group)
  require_detected <- .cpgd_require_detected(age_group, require_detected)
  accuracy <- .cpgd_accuracy_arg(accuracy)

  if (!identical(tissue, "all")) {
    r1 <- .cpgd_one_tissue(cpgs, genes = genes, tissue = tissue,
                           min_confidence = min_confidence,
                           allow_nearest_gene = allow_nearest_gene,
                           verbose = verbose, accuracy = accuracy)
    if (!is.null(html_out)) .cpgd_maybe_report(r1, html_out)
    return(invisible(r1))
  }

  short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")

  # ---- layer 2: the catalogue models, one pass per tissue -----------------
  per <- list()
  for (t in CPGD_TISSUES) {
    r <- .cpgd_one_tissue(cpgs, genes = genes, tissue = t,
                          min_confidence = min_confidence,
                          allow_nearest_gene = allow_nearest_gene, verbose = FALSE)
    s <- short[[t]]
    keep <- data.table::data.table(
      cpg_id = r$cpg_id,
      input  = r$input,
      requested_gene = if ("requested_gene" %in% names(r)) r$requested_gene else NA_character_)
    # data.table::set(), not `[[<-`. The latter is base-R assignment, which
    # copies the table and marks it, and every later `:=` then emits a
    # shallow-copy warning that has nothing to do with the user's code.
    data.table::set(keep, j = paste0("gene_", s), value = r$target_gene)
    data.table::set(keep, j = paste0("dir_", s),
                    value = suppressWarnings(as.numeric(
                      ifelse(r$call %in% c("1", "-1"), r$call, NA))))
    data.table::set(keep, j = paste0("conf_", s), value = r$confidence)
    data.table::set(keep, j = paste0("tier_", s),
                    value = if ("evidence_tier" %in% names(r)) r$evidence_tier
                            else NA_character_)
    data.table::set(keep, j = paste0("status_", s), value = r$status)
    per[[t]] <- keep
  }
  out <- Reduce(function(a, b)
    merge(a, b[, setdiff(names(b), c("input", "requested_gene")), with = FALSE],
          by = "cpg_id", all = TRUE), per)

  n <- nrow(out)
  if (!n) {
    message("No CpG identifiers resolved.")
    return(invisible(out))
  }

  dcols <- paste0("dir_",  short)
  ccols <- paste0("conf_", short)
  tcols <- paste0("tier_", short)
  D  <- as.matrix(out[, dcols, with = FALSE])
  Cf <- as.matrix(out[, ccols, with = FALSE])
  TT <- as.matrix(out[, tcols, with = FALSE])

  n_call <- rowSums(!is.na(D))
  n_pos  <- rowSums(D ==  1, na.rm = TRUE)
  n_neg  <- rowSums(D == -1, na.rm = TRUE)
  agree  <- ifelse(n_call > 0, pmax(n_pos, n_neg) / n_call, NA_real_)
  cons   <- ifelse(n_call >= min_tissues & agree == 1,
                   ifelse(n_pos > n_neg, 1, -1), NA_real_)
  cons   <- ifelse(is.na(cons) & n_call == 1L & min_tissues <= 1L,
                   ifelse(n_pos > 0, 1, -1), cons)

  out[, "n_tissues_calling"  := n_call]
  out[, "tissue_agreement"   := round(agree, 3)]
  out[, "consensus_direction" := cons]

  # the direction when exactly one tissue spoke, kept separate from consensus,
  # and which tissue's model it was (its accuracy depends on that)
  single <- vapply(seq_len(n), function(i) {
    v <- D[i, ][!is.na(D[i, ])]
    if (length(v) == 1L) v[1] else NA_real_
  }, numeric(1))
  single_tis <- vapply(seq_len(n), function(i) {
    w <- which(!is.na(D[i, ]))
    if (length(w) == 1L) unname(short[w]) else NA_character_
  }, character(1))

  # worst tier among the tissues that actually called: a consensus is only as
  # good as its weakest contributor
  worst_tier <- vapply(seq_len(n), function(i) {
    v <- TT[i, ][!is.na(D[i, ])]
    v <- v[!is.na(v)]
    if (!length(v)) NA_character_ else max(v)
  }, character(1))

  # Average confidence over the tissues that actually CALLED, not over every
  # tissue holding a number. Including abstainers dragged a single confident
  # call down towards the abstention threshold and made catalogue_single rows
  # look far weaker than they are.
  conf_cat <- vapply(seq_len(n), function(i) {
    v <- Cf[i, ][!is.na(D[i, ])]
    v <- v[!is.na(v)]
    if (!length(v)) NA_real_ else mean(v)
  }, numeric(1))

  # Park these on the table rather than keeping them as loose vectors. The
  # merges below re-sort by cpg_id, and a parallel vector would silently end up
  # attached to the wrong rows.
  out[, ".single"     := single]
  out[, ".single_tis" := single_tis]
  out[, ".worst_tier" := worst_tier]
  out[, ".conf_cat"   := conf_cat]

  # ---- layer 1: directly measured associations ----------------------------
  out <- .cpgd_add_measured(out)

  req  <- toupper(out$requested_gene)
  last <- vapply(strsplit(ifelse(is.na(req), "", req), "[_-]"),
                 function(p) if (length(p)) p[length(p)] else NA_character_,
                 character(1))

  md <- rep(NA_real_, nrow(out)); mn <- rep(NA_integer_, nrow(out)); ma <- rep(NA_real_, nrow(out))
  M <- tryCatch(cpgd_measured_eqtms(), error = function(e) NULL)
  if (!is.null(M)) {
    # Match the measured record to the gene actually requested, not to whatever
    # the CpG happens to have been measured against: a measured association with
    # a DIFFERENT gene is reported separately, in annotation_mismatch.
    nco <- if ("n_cohorts" %in% names(M)) suppressWarnings(as.integer(M$n_cohorts)) else rep(1L, nrow(M))
    agr <- if ("cohort_agreement" %in% names(M)) suppressWarnings(as.numeric(M$cohort_agreement)) else rep(1, nrow(M))
    MM <- data.table::data.table(cpg_id = M$cpg_id,
                                 g = toupper(M$target_gene),
                                 d = sign(as.numeric(M$direction)),
                                 nco = ifelse(is.na(nco), 1L, nco),
                                 agr = ifelse(is.na(agr), 1, agr))
    MM <- MM[, list(d = if (length(unique(get("d"))) == 1L) get("d")[1] else NA_real_,
                    nco = sum(get("nco")), agr = min(get("agr"))),
             by = c("cpg_id", "g")]
    data.table::setkeyv(MM, c("cpg_id", "g"))
    h1 <- MM[data.table::data.table(cpg_id = out$cpg_id, g = req), on = c("cpg_id", "g")]
    h2 <- MM[data.table::data.table(cpg_id = out$cpg_id, g = last), on = c("cpg_id", "g")]
    d1 <- h1$d; d2 <- h2$d
    md <- ifelse(!is.na(d1), d1, d2)
    mn <- ifelse(!is.na(d1), h1$nco, ifelse(!is.na(d2), h2$nco, NA_integer_))
    ma <- ifelse(!is.na(d1), h1$agr, ifelse(!is.na(d2), h2$agr, NA_real_))
  }
  out[, "measured_direction_requested" := md]
  out[, "measured_n_cohorts" := as.integer(mn)]
  out[, "measured_cohort_agreement" := as.numeric(ma)]
  # the same, per tissue the pair was measured in: the ladder ranks a measured
  # sign by where it was measured
  MT <- .cpgd_measured_by_tissue(M, out$cpg_id, req, last)
  for (s in names(MT)) {
    data.table::set(out, j = paste0(".meas_d_", s), value = MT[[s]]$d)
    data.table::set(out, j = paste0(".meas_rep_", s), value = MT[[s]]$rep)
  }

  # ---- layer 1b: causal directions from SMR --------------------------------
  # GoDMC mQTL as exposure, eQTLGen cis-eQTL as outcome. This asks whether a
  # genetic instrument for methylation moves expression, rather than whether
  # some study happened to measure the pair, so its coverage is not bounded by
  # what the catalogues tested. Median CpG-gene distance 135 kb: it reaches
  # exactly the range where the distance curves go flat.
  scols <- c("smr_direction", "smr_tier", "smr_p", "smr_n_instruments",
             "smr_gene_dist", "smr_gene", "smr_gene_match", "smr_in_table",
             "smr_heidi_status", "smr_p_heidi",
             "smr_rung", "smr_replication", "smr_outcome")
  with_child <- !identical(age_group, "adult")
  S <- tryCatch(.cpgd_smr_table(), error = function(e) NULL)
  if (!is.null(S) && nrow(S)) {
    # heidi_status is absent from tables built before 2.2.5, so read it
    # defensively rather than assuming the shipped file is current: a user with
    # a stale cached copy should get NA, not an error.
    hs <- if ("heidi_status" %in% names(S)) as.character(S$heidi_status)
          else rep(NA_character_, nrow(S))
    ph <- if ("p_HEIDI" %in% names(S)) suppressWarnings(as.numeric(S$p_HEIDI))
          else rep(NA_real_, nrow(S))
    SS <- data.table::data.table(
      cpg_id = S$cpg_id, g = toupper(S$target_gene),
      d = as.numeric(S$direction), tier = as.character(S$smr_tier),
      p = as.numeric(S$p_SMR), ni = as.integer(S$n_instruments),
      gd = suppressWarnings(as.numeric(S$cpg_gene_dist)),
      hstat = hs, pheidi = ph,
      rung = .cpgd_smr_rung(S),
      rep  = .cpgd_smr_col(S, "replication", "eqtlgen_only"),
      outc = .cpgd_smr_col(S, "outcome", "eQTLGen"),
      cb   = .cpgd_smr_child(S))
    require_detected <- .cpgd_check_child_column(S, age_group, require_detected)
    # Whether the CpG has SMR evidence at all, judged BEFORE any filter, so a
    # CpG whose pairs require_detected removed is not reported as having no
    # instrument.
    in_table <- unique(SS$cpg_id)
    # age_group = "child", require_detected = TRUE: only SMR pairs HELIX child
    # blood detected (sig or nominal) stay in the layer. Measured, catalogue
    # and distance evidence are not filtered.
    if (isTRUE(require_detected)) {
      drop <- !(SS$cb %in% CPGD_CHILD_DETECTED)
      message(sprintf(paste0("require_detected: %s SMR pair(s) for these CpGs ",
                             "dropped (child_blood_helix not sig/nominal); ",
                             "measured and catalogue evidence unaffected."),
                      format(sum(drop & SS$cpg_id %chin% out$cpg_id), big.mark = ",")))
      SS <- SS[!drop]
    }
    # one row per pair: strongest tier, then smallest p
    data.table::setorderv(SS, c("cpg_id", "g", "tier", "p"))
    SS <- SS[, .SD[1L], by = c("cpg_id", "g")]
    data.table::setkeyv(SS, c("cpg_id", "g"))

    # Whether the CpG appears in the SMR table AT ALL, independent of gene. The
    # print method needs this to avoid asserting a coverage explanation for a
    # CpG that is present and simply had no row for the gene we asked about.
    out[, "smr_in_table" := out$cpg_id %chin% in_table]

    # Candidate gene keys, in priority order: what the caller asked for, then
    # the target genes the catalogue layers resolved on their own. A hit on any
    # of these is a GENE MATCH -- the causal estimate and the catalogue
    # prediction are talking about the same pair.
    #
    # The catalogue-derived candidates apply only where the caller named nothing.
    # If you asked for CREBBP, a direction for ADCY9 is not an answer to your
    # question, and substituting one would be worse than returning nothing.
    gcols <- intersect(paste0("gene_", short), names(out))
    cand  <- c(list(req, last),
               if (length(gcols)) lapply(gcols, function(cc) toupper(out[[cc]])))
    asked <- !is.na(req)

    d <- rep(NA_real_, n); ti <- rep(NA_character_, n); pp <- rep(NA_real_, n)
    ni <- rep(NA_integer_, n); gd <- rep(NA_real_, n); gg <- rep(NA_character_, n)
    hh <- rep(NA_character_, n); hp <- rep(NA_real_, n)
    rg <- rep(NA_character_, n); rp <- rep(NA_character_, n); oc <- rep(NA_character_, n)
    cbv <- rep(NA_character_, n)
    for (k in seq_along(cand)) {
      key <- cand[[k]]
      if (k > 2L) key[asked] <- NA_character_
      need <- is.na(d) & !is.na(key)
      if (!any(need)) next
      hit <- SS[data.table::data.table(cpg_id = out$cpg_id, g = key),
                on = c("cpg_id", "g")]
      take <- need & !is.na(hit$d)
      if (!any(take)) next
      d[take]  <- hit$d[take];    ti[take] <- hit$tier[take]
      pp[take] <- hit$p[take];    ni[take] <- hit$ni[take]
      gd[take] <- hit$gd[take];   gg[take] <- key[take]
      hh[take] <- hit$hstat[take]; hp[take] <- hit$pheidi[take]
      rg[take] <- hit$rung[take];  rp[take] <- hit$rep[take]; oc[take] <- hit$outc[take]
      cbv[take] <- hit$cb[take]
    }
    matched <- !is.na(d)

    # Where no gene matched, still REPORT the CpG's strongest SMR pair rather
    # than nothing, naming the gene it concerns.
    #
    # The two layers largely disagree about which gene a CpG regulates: only
    # 58.6% of CpGs present in both share even one gene between them, and since
    # the catalogue commits to a single target the realised match rate is lower
    # still. Gating the layer on gene identity therefore hid causal evidence for
    # the majority of the 103,285 CpGs it covers. Withholding an S1 estimate at
    # p = 1e-30 because it names RBL2 while the catalogue named SCP2 tells the
    # user nothing they can act on; showing it, labelled, tells them a great deal.
    #
    # What such a row must NOT do is supply best_direction, because that column
    # answers a question about a specific pair. `smr_gene_match` carries the
    # distinction and every downstream use of the layer is gated on it.
    SB <- SS[, .SD[1L], by = "cpg_id"]        # SS is already ordered tier, then p
    data.table::setkeyv(SB, "cpg_id")
    hb <- SB[data.table::data.table(cpg_id = out$cpg_id), on = "cpg_id"]
    fill <- !matched & !is.na(hb$d)
    d[fill]  <- hb$d[fill];    ti[fill] <- hb$tier[fill]
    pp[fill] <- hb$p[fill];    ni[fill] <- hb$ni[fill]
    gd[fill] <- hb$gd[fill];   gg[fill] <- hb$g[fill]
    hh[fill] <- hb$hstat[fill]; hp[fill] <- hb$pheidi[fill]
    rg[fill] <- hb$rung[fill];  rp[fill] <- hb$rep[fill]; oc[fill] <- hb$outc[fill]
    cbv[fill] <- hb$cb[fill]

    out[, "smr_direction"     := d]
    out[, "smr_tier"          := ti]
    out[, "smr_p"             := pp]
    out[, "smr_n_instruments" := ni]
    out[, "smr_gene_dist"     := gd]
    # Which gene the SMR direction refers to. Never left to inference: without
    # this column a reported direction is ambiguous as to its subject.
    out[, "smr_gene"          := gg]
    out[, "smr_gene_match"    := matched]
    # Reported, never acted on. HEIDI did not discriminate on these data -- see
    # ?cpgd_smr_directions -- so nothing in this package conditions on it, and a
    # user who wants to must do so deliberately.
    out[, "smr_heidi_status"  := hh]
    out[, "smr_p_heidi"       := hp]
    # How the row enters the ladder (replicated / high / moderate / weak /
    # reported), why, and which outcome study the sign comes from.
    out[, "smr_rung"          := rg]
    out[, "smr_replication"   := rp]
    out[, "smr_outcome"       := oc]
    # Child-blood detection of the reported SMR pair, only when asked for: the
    # adult default carries no age-specific column.
    if (with_child) out[, "smr_child_blood_helix" := cbv]
  } else {
    out[, (scols) := list(NA_real_, NA_character_, NA_real_, NA_integer_,
                          NA_real_, NA_character_, FALSE, FALSE,
                          NA_character_, NA_real_,
                          NA_character_, NA_character_, NA_character_)]
    if (with_child) out[, "smr_child_blood_helix" := NA_character_]
  }

  # ---- layer 1c: tumour consensus (Onco-eQTM) -----------------------------
  # Keyed on the same gene candidates as SMR (requested gene, its last token,
  # then the catalogue-assigned genes where nothing was requested), and gated
  # the same way: only a gene MATCH may reach best_direction, and only tier O1
  # (see ?cpgd_onco_eqtm). Reported for every CpG the layer covers.
  ocols <- c("onco_direction", "onco_tier", "onco_n_cancers", "onco_median_r",
             "onco_gene", "onco_gene_match")
  OO <- .cpgd_onco_keyed()
  if (!is.null(OO)) {
    gcols <- intersect(paste0("gene_", short), names(out))
    cand  <- c(list(req, last),
               if (length(gcols)) lapply(gcols, function(cc) toupper(out[[cc]])))
    asked <- !is.na(req)
    od <- rep(NA_real_, n); ot <- rep(NA_character_, n); on_ <- rep(NA_integer_, n)
    om <- rep(NA_real_, n); og <- rep(NA_character_, n)
    for (k in seq_along(cand)) {
      key <- cand[[k]]
      if (k > 2L) key[asked] <- NA_character_
      need <- is.na(ot) & !is.na(key)
      if (!any(need)) next
      hit <- OO[data.table::data.table(cpg_id = out$cpg_id, .gkey = gsub("_", "-", key)),
                on = c("cpg_id", ".gkey")]
      take <- need & !is.na(hit$tier)
      if (!any(take)) next
      od[take] <- hit$d[take];  ot[take] <- hit$tier[take]
      on_[take] <- hit$nc[take]; om[take] <- hit$mr[take]; og[take] <- key[take]
    }
    matched <- !is.na(ot)
    # no gene matched: still report the CpG's strongest tumour pair, named,
    # exactly as the SMR layer does -- visible, but barred from best_direction
    OB <- OO[, .SD[1L], by = "cpg_id"]
    hb <- OB[data.table::data.table(cpg_id = out$cpg_id), on = "cpg_id"]
    fill <- !matched & !is.na(hb$tier)
    od[fill] <- hb$d[fill];  ot[fill] <- hb$tier[fill]
    on_[fill] <- hb$nc[fill]; om[fill] <- hb$mr[fill]; og[fill] <- hb$g[fill]
    out[, "onco_direction"  := od]
    out[, "onco_tier"       := ot]
    out[, "onco_n_cancers"  := on_]
    out[, "onco_median_r"   := om]
    out[, "onco_gene"       := og]
    out[, "onco_gene_match" := matched]
  } else {
    out[, (ocols) := list(NA_real_, NA_character_, NA_integer_, NA_real_,
                          NA_character_, FALSE)]
  }

  # ---- layer 3: distance only --------------------------------------------
  ucols <- c("gene_universal", "dist_universal", "dir_universal", "p_universal",
             "agree_universal", "dist_unanimous",
             "dist_dir_blood", "dist_dir_nasal", "dist_dir_solid")
  uni <- NULL
  if (isTRUE(universal)) {
    uni <- tryCatch(
      cpg_direction_universal(cpgs, genes = genes, tissue = "all", verbose = FALSE),
      error = function(e) {
        if (isTRUE(verbose))
          message("distance layer unavailable (", conditionMessage(e),
                  "); catalogue layers still reported.")
        NULL
      })
  }
  if (!is.null(uni) && nrow(uni) && length(intersect(paste0("p_", short), names(uni)))) {
    pc <- intersect(paste0("p_", short), names(uni))
    P3 <- as.matrix(uni[, pc, with = FALSE])
    # Take the pooled probability rather than demanding the three curves agree
    # on a side of 0.5. The blood curve asymptotes at 0.449 and never crosses,
    # while nasal and solid rise past it, so beyond ~30 kb they disagree by
    # construction. Requiring consensus there threw the estimate away and
    # reported "no evidence" for pairs whose distance was perfectly well known.
    pm   <- rowMeans(P3, na.rm = TRUE)
    nval <- rowSums(!is.na(P3))
    npos <- rowSums(P3 >= 0.5, na.rm = TRUE)
    # Unanimity across the three tissue curves, NOT the pooled mean, decides
    # whether a distance call is made. Averaging curves that point opposite ways
    # produces p ~ 0.50 and looks like ignorance, when in fact blood says -1 and
    # solid tissue says +1 and both say it clearly. That is a conflict to report,
    # not a coin flip to abstain on, and it is the same rule the catalogue layer
    # already uses.
    unanimous <- nval > 0L & (npos == nval | npos == 0L)
    u <- data.table::data.table(
      cpg_id          = uni$cpg_id,
      gene_universal  = uni$target_gene,
      dist_universal  = as.numeric(uni$abs_dist),
      dir_universal   = ifelse(!unanimous | is.nan(pm), NA_real_,
                               ifelse(npos == nval, 1, -1)),
      p_universal     = round(pm, 4),
      agree_universal = round(pmax(npos, nval - npos) / pmax(nval, 1L), 3),
      dist_unanimous  = unanimous)
    for (k in seq_along(pc)) {
      s <- sub("^p_", "", pc[k])
      data.table::set(u, j = paste0("dist_dir_", s),
                      value = ifelse(is.na(P3[, k]), NA_real_,
                                     ifelse(P3[, k] >= 0.5, 1, -1)))
    }
    u <- unique(u, by = "cpg_id")
    out <- merge(out, u, by = "cpg_id", all.x = TRUE)
  } else {
    out[, (ucols) := list(NA_character_, NA_real_, NA_real_, NA_real_, NA_real_,
                          NA, NA_real_, NA_real_, NA_real_)]
  }

  # ---- the ladder ---------------------------------------------------------
  # Every source that has a sign for the row is a candidate; the most accurate
  # one in the call's tissue answers (target_tissue, else blood). The order is
  # not written here: .cpgd_rung_order() reads it from the per-tissue
  # validation (R/rung_accuracy.R, tools/validate_rung_accuracy.R). which()
  # throughout: these conditions carry NA, and a logical NA in a subscripted
  # assignment is an error rather than a skipped row.
  rank_tissue <- if (is.null(target_tissue)) "blood" else match.arg(target_tissue, CPGD_TISSUES)
  # the catalogue calls, rebuilt from the table (merges re-sorted it): a
  # DIRECT_eQTM call is a measured sign; a consensus or a conflict is between
  # predicted calls on the same gene (.cpgd_catalogue_cand())
  mcol <- function(pre) as.matrix(out[, paste0(pre, short), with = FALSE])
  CC <- .cpgd_catalogue_cand(mcol("dir_"), mcol("tier_"), mcol("status_"), mcol("gene_"),
                             min_tissues)
  colnames(CC$single) <- colnames(CC$direct) <- unname(short)
  TTm <- mcol("tier_"); colnames(TTm) <- unname(short)
  cand <- list()
  # measured, by the tissue it was measured in, for the requested gene; where
  # the measured layer has nothing in a tissue, that tissue's DIRECT_eQTM call
  for (s in c("blood", "nasal", "solid")) {
    dv <- out[[paste0(".meas_d_", s)]]; rp <- out[[paste0(".meas_rep_", s)]]
    if (is.null(dv)) { dv <- rep(NA_real_, nrow(out)); rp <- rep(FALSE, nrow(out)) }
    cand[[paste0("measured_", s, "_replicated")]] <- ifelse(rp %in% TRUE, dv, NA_real_)
    cand[[paste0("measured_", s, "_single")]] <-
      ifelse(rp %in% TRUE, NA_real_, ifelse(is.na(dv), CC$direct[, s], dv))
  }
  # Only SMR rows that concern the SAME gene as the rest of the row may supply
  # best_direction. Where smr_gene_match is FALSE the estimate is still reported
  # in the smr_* columns, but it answers a question about a different pair. The
  # reported rung never supplies a sign.
  smr_ok <- out$smr_gene_match %in% TRUE & !is.na(out$smr_direction)
  rung <- out$smr_rung
  smr_src <- data.table::fcase(rung %in% "replicated", "smr_replicated",
                               rung %in% "high", "smr_high",
                               rung %in% "moderate", "smr_moderate",
                               rung %in% "weak", "smr_weak", default = NA_character_)
  skey <- .cpgd_smr_key(smr_src, out$smr_tier, out$smr_outcome, out$smr_replication)
  for (k in unique(stats::na.omit(skey)))
    cand[[paste0("smr_", k)]] <- ifelse(smr_ok & skey %in% k, out$smr_direction, NA_real_)
  # consensus needs two models at least, whatever min_tissues says: a lone
  # model's call is a single call, ranked and quoted as that model's
  for (tr in c("A", "B"))
    cand[[paste0("catalogue_consensus_", tr)]] <-
      ifelse(!is.na(CC$cons) & CC$worst %in% tr, CC$cons, NA_real_)
  for (s in c("blood", "nasal", "solid")) for (tr in c("A", "B"))
    cand[[paste0("catalogue_single_", s, "_", tr)]] <-
      ifelse(!is.na(CC$single[, s]) & TTm[, s] %in% tr, CC$single[, s], NA_real_)
  onco_ok <- out$onco_gene_match %in% TRUE & !is.na(out$onco_direction) &
             out$onco_tier %in% "O1"
  cand[["onco_consensus"]] <- ifelse(onco_ok, out$onco_direction, NA_real_)
  # distance: all three tissue curves pointing the same way is the evidence
  dist_ok <- out$dist_unanimous %in% TRUE & !is.na(out$dir_universal) &
             abs(out$p_universal - 0.5) * 2 >= min_distance_info
  cand[["distance_only"]] <- ifelse(dist_ok %in% TRUE, out$dir_universal, NA_real_)
  # Tissues pointing opposite ways is a finding, not a gap: from the distance
  # rung down, no weaker source settles a disagreement between catalogue models.
  conflict <- CC$conflict
  lad  <- .cpgd_apply_ladder(cand, .cpgd_rung_order(rank_tissue), conflict, nrow(out))
  best <- lad$best; rkey <- lad$key
  # Distance states where nothing answered: unanimous but sitting on 0.50 is
  # genuinely nothing; curves pointing opposite ways is a disagreement to
  # report, not average away (per-tissue signs in dist_dir_*).
  i <- which(is.na(rkey) & out$dist_unanimous %in% TRUE)
  rkey[i] <- "distance_uninformative"
  i <- which(is.na(rkey) & !is.na(out$dist_universal))
  rkey[i] <- "distance_tissue_conflict"
  rkey[is.na(rkey)] <- "no_evidence"
  src <- .cpgd_rung_evidence(rkey)

  # ---- declared target tissue ---------------------------------------------
  # Where the tissues disagree, a caller who has declared in advance which
  # tissue their question is about can have that tissue's answer instead of an
  # abstention. This is tissue SELECTION, not voting: the tissue is named before
  # the votes are seen, and its own accuracy is what the row is then worth.
  if (!is.null(target_tissue)) {
    tt <- match.arg(target_tissue, CPGD_TISSUES)
    ts <- short[[tt]]
    dcol <- paste0("dist_dir_", ts)
    ccol <- paste0("dir_", ts)
    tcol <- paste0("tier_", ts)
    if (dcol %in% names(out)) {
      i <- which(src == "distance_tissue_conflict" & !is.na(out[[dcol]]))
      best[i] <- out[[dcol]][i];  src[i] <- "distance_targeted"
      rkey[i] <- paste0("distance_curve_", ts)
    }
    if (ccol %in% names(out)) {
      i <- which(src == "tissue_conflict" & !is.na(out[[ccol]]))
      best[i] <- out[[ccol]][i];  src[i] <- "catalogue_targeted"
      # the declared tissue's call: its model's figure, or, for a DIRECT_eQTM
      # row, the figure of a measured sign
      rkey[i] <- ifelse(out[[paste0("status_", ts)]][i] %in% "DIRECT_eQTM",
                        paste0("measured_", ts, "_single"),
                        paste0("catalogue_model_", ts, "_", out[[tcol]][i]))
    }
    # Last resort: anything that would still abstain takes the target tissue's
    # distance-curve sign, if that tissue has one. The rows are labelled, and
    # carry that curve's validated figure: the validation scores the same
    # sign (the curve's p >= 0.5) on every yardstick pair.
    if (dcol %in% names(out)) {
      i <- which(is.na(best) & !is.na(out[[dcol]]))
      best[i] <- out[[dcol]][i];  src[i] <- "targeted_last_resort"
      rkey[i] <- paste0("distance_curve_", ts)
    }
  }

  bconf <- rep(NA_real_, nrow(out))
  bconf[which(src == "measured")] <- 1
  i <- which(src %in% c("catalogue_consensus", "catalogue_single"))
  bconf[i] <- round(out$.conf_cat[i], 4)
  i <- which(src == "distance_only")
  bconf[i] <- round(abs(out$p_universal[i] - 0.5) * 2, 4)

  # what the answer is worth: the validation of its source in the ranked tissue
  bacc <- .cpgd_rung_label(ifelse(is.na(best), NA_character_, rkey), rank_tissue, age_group)

  # the tier of the catalogue source that answered (its own model's tier for a
  # single call)
  wt <- out$.worst_tier
  i <- which(grepl("^catalogue_(consensus|single|model)_.*_?[AB]$", rkey))
  wt[i] <- substring(rkey[i], nchar(rkey[i]))
  wt[which(src %in% c("distance_only", "distance_uninformative",
                      "distance_tissue_conflict", "distance_targeted",
                      "targeted_last_resort"))] <- "U"
  wt[which(src == "measured")] <- "M"
  i <- which(src %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak"))
  wt[i] <- out$smr_tier[i]

  # ---- the two alternative codings ---------------------------------------
  # A row is "uncertain" when its sign would flip under a perturbation smaller
  # than the evidence: the distance prior sitting at p = 0.5009 is the clear
  # case. For those rows `best_direction` abstains, which is honest but leaves
  # nothing to analyse. These two columns bracket the ambiguity instead:
  # _filled takes the point estimate, _flipped takes its opposite.
  #
  # They are a SENSITIVITY ANALYSIS, not a menu. The question they answer is
  # "does my downstream result survive both codings", not "which coding gives
  # me a result". Choosing whichever column fits the data better is selecting
  # on the outcome, and doubles the effective number of tests.
  pu  <- out$p_universal
  pt  <- ifelse(!is.na(pu), ifelse(pu >= 0.5, 1, -1), NA_real_)
  unc <- (src %in% c("distance_uninformative", "distance_tissue_conflict",
                     "distance_targeted", "catalogue_targeted",
                     "targeted_last_resort")) |
         (!is.na(bconf) & bconf < min_distance_info)
  unc[is.na(unc)] <- FALSE

  filled <- best
  i <- which(unc & is.na(filled) & !is.na(pt))
  filled[i] <- pt[i]
  flipped <- filled
  j <- which(unc & !is.na(flipped))
  flipped[j] <- -flipped[j]

  out[, "direction_uncertain"     := unc]
  out[, "best_direction_filled"   := filled]
  out[, "best_direction_flipped"  := flipped]

  # SMR frequently loses the ladder to a catalogue consensus and then vanishes
  # from the output entirely, even where it disagrees. A correlational
  # prediction and a causal estimate pointing opposite ways is a finding, not
  # noise to be resolved by precedence, so record it.
  # Agreement is only meaningful between two statements about the same pair.
  # Where the SMR row concerns a different gene there is nothing to agree or
  # disagree with, and reporting FALSE there would manufacture a contradiction
  # out of two compatible facts.
  out[, "smr_agreement" := data.table::fifelse(
        is.na(get("smr_direction")) | is.na(best) |
          !(get("smr_gene_match") %in% TRUE), NA,
        get("smr_direction") == best)]

  out[, "onco_agreement" := data.table::fifelse(
        is.na(get("onco_direction")) | is.na(best) |
          !(get("onco_gene_match") %in% TRUE) | src == "onco_consensus", NA,
        get("onco_direction") == best)]

  out[, "best_direction"  := best]
  out[, "best_evidence"   := src]
  data.table::set(out, j = "best_rung", value = rkey)
  out[, "best_confidence" := bconf]
  out[, "best_expected_accuracy" := bacc]
  # the same figure as one number, for engines that weight by it
  data.table::set(out, j = "best_accuracy",
                  value = .cpgd_accuracy_number(bacc, accuracy))
  out[, "best_tier" := worst_tier]
  out[, "best_label" := data.table::fcase(
        get("best_direction") ==  1, "higher methylation -> HIGHER expression",
        get("best_direction") == -1, "higher methylation -> LOWER expression",
        default = NA_character_)]

  out[, "note" := ""]
  out[get("best_evidence") == "tissue_conflict", "note" :=
        "tissues returned opposite directions; not resolved"]
  out[get("best_evidence") == "no_evidence", "note" :=
        "no catalogue holds this CpG-gene pair and no distance could be computed"]
  out[get("best_evidence") == "distance_uninformative", "note" :=
        "all three distance curves agree in sign but sit on 0.5; no usable signal"]
  out[get("best_evidence") == "targeted_last_resort", "note" :=
        "LAST RESORT: no measured, SMR, catalogue or unanimous distance evidence; this is the target tissue's distance-curve sign, worth that curve's figure"]
  out[get("best_evidence") %in% c("distance_targeted", "catalogue_targeted"),
      "note" := "tissues disagreed; resolved by the target tissue you declared - worth that tissue's accuracy, not the consensus figure"]
  out[get("best_evidence") == "distance_tissue_conflict", "note" :=
        "distance curves disagree by tissue - see dist_dir_blood / _nasal / _solid; NOT a coin flip"]
  out[get("best_evidence") == "distance_only", "note" :=
        "distance-only prior; weakest tier, treat as orienting a hypothesis"]
  out[which(get("smr_agreement") == FALSE), "note" := paste0(
        get("note"), " | SMR (causal, blood) says ",
        ifelse(get("smr_direction") > 0, "+1", "-1"),
        " at tier ", get("smr_tier"),
        " but a stronger layer supplied the opposite sign - inspect before using")]
  out[get("best_evidence") %in% c("smr_high", "smr_moderate", "smr_weak"), "note" :=
        "causal direction from SMR (blood); consistent with methylation affecting expression AND with linkage - see smr_heidi_status and smr_replication"]
  out[get("best_evidence") == "smr_replicated", "note" :=
        "causal direction from SMR (blood), same sign from two independent outcome studies (eQTLGen and FinnGen PBMC) on the GoDMC instruments"]
  out[get("best_evidence") == "smr_weak" & get("smr_replication") %in% "discordant", "note" :=
        "SMR (blood): eQTLGen and FinnGen PBMC give OPPOSITE signs for this pair; the eQTLGen sign is kept at the weak rung - inspect before using"]
  out[get("best_evidence") == "smr_moderate" & !(get("smr_outcome") %in% c(NA, "eQTLGen")), "note" :=
        "causal direction from SMR (blood), FinnGen PBMC outcome only (pair absent from eQTLGen; SMR p<5e-8 at the top cis-mQTL), see ?cpgd_smr_directions"]
  out[get("best_evidence") == "catalogue_single", "note" :=
        "only one tissue called; no cross-tissue corroboration"]
  out[get("best_evidence") == "onco_consensus", "note" :=
        "direction from the tumour consensus (Onco-eQTM O1, >=5 cancer types); no measured, SMR S1/S2 or catalogue call for this pair"]
  out[which(get("onco_agreement") == FALSE), "note" := paste0(
        get("note"), " | tumour consensus (Onco-eQTM O1/O2) says ",
        ifelse(get("onco_direction") > 0, "+1", "-1"),
        " but a peripheral layer supplied the opposite sign")]
  out[which(get("annotation_mismatch")), "note" := trimws(paste0(
        get("note"), " | measured association is with ", get("measured_genes"),
        ", not the requested gene"), which = "left")]

  ord <- c("input", "cpg_id", "requested_gene",
           "best_direction", "best_label", "best_confidence", "best_evidence",
           "best_rung", "best_tier", "best_expected_accuracy", "best_accuracy",
           "direction_uncertain", "best_direction_filled", "best_direction_flipped",
           "smr_direction", "smr_gene", "smr_gene_match", "smr_tier", "smr_p",
           "smr_n_instruments", "smr_gene_dist", "smr_agreement", "smr_in_table",
           "smr_heidi_status", "smr_p_heidi",
           "smr_rung", "smr_replication", "smr_outcome", "smr_child_blood_helix",
           "onco_direction", "onco_tier", "onco_n_cancers", "onco_median_r",
           "onco_gene", "onco_gene_match", "onco_agreement",
           "consensus_direction", "n_tissues_calling", "tissue_agreement",
           as.vector(rbind(paste0("gene_", short), paste0("dir_", short),
                           paste0("conf_", short), paste0("tier_", short),
                           paste0("status_", short))),
           ucols,
           "measured_genes", "measured_tissues", "measured_direction",
           "measured_direction_requested", "measured_n_cohorts",
           "measured_cohort_agreement",
           "annotation_mismatch", "mismatch_note", "note")
  drop <- intersect(c(".single", ".single_tis", ".worst_tier", ".conf_cat",
                     paste0(".meas_", rep(c("d_", "rep_"), each = 3),
                            c("blood", "nasal", "solid"))), names(out))
  if (length(drop)) out[, (drop) := NULL]
  out <- out[, intersect(ord, names(out)), with = FALSE]
  data.table::setattr(out, "target_tissue",
                      if (is.null(target_tissue)) NULL else
                        match.arg(target_tissue, CPGD_TISSUES))
  # The population the SMR figures describe; set only off the adult default,
  # so the default object is unchanged.
  if (with_child) data.table::setattr(out, "age_group", age_group)
  if (isTRUE(require_detected)) data.table::setattr(out, "require_detected", TRUE)
  data.table::setattr(out, "accuracy", accuracy)
  data.table::setattr(out, "rank_tissue", rank_tissue)
  data.table::setattr(out, "class", c("cpgd_full", class(out)))
  if (isTRUE(verbose)) .cpgd_message_print(out)
  if (!is.null(html_out)) .cpgd_maybe_report(out, html_out)
  invisible(out)
}


# html_out = TRUE writes a default name into tempdir() (cpgd_report's default
# directory); a string is used as the file name. Reporting must never take down
# a completed analysis, so failures here warn rather than stop.
.cpgd_maybe_report <- function(res, html_out) {
  f <- if (isTRUE(html_out)) "cpgdirection_results.html" else as.character(html_out)[1]
  tryCatch(cpgd_report(res, file = f),
           error = function(e)
             warning("results computed but the HTML report failed: ",
                     conditionMessage(e), call. = FALSE))
}


#' @export
print.cpgd_full <- function(x, ...) {
  cat(sprintf("\ncpgdirection: %d CpGs, all methods\n", nrow(x)))
  # Any level the ladder can emit but the constant does not know about would be
  # silently dropped from the tally, so surface it rather than lose it.
  lev <- union(CPGD_EVIDENCE, unique(stats::na.omit(x$best_evidence)))
  ev <- table(factor(x$best_evidence, levels = lev))
  cat("\n  best_direction came from:\n")
  for (nm in names(ev)) if (ev[[nm]] > 0) cat(sprintf("    %-20s %4d\n", nm, ev[[nm]]))
  called <- sum(!is.na(x$best_direction))
  cat(sprintf("\n  direction assigned : %d of %d (%.0f%%)\n",
              called, nrow(x), 100 * called / max(nrow(x), 1)))
  if (called) {
    tb <- table(x$best_direction)
    cat(sprintf("    +1 %d    -1 %d\n",
                if ("1" %in% names(tb)) tb[["1"]] else 0,
                if ("-1" %in% names(tb)) tb[["-1"]] else 0))
  }
  if ("annotation_mismatch" %in% names(x)) {
    k <- sum(x$annotation_mismatch %in% TRUE)
    if (k > 0) cat(sprintf(
      "\n  *** %d CpGs have a MEASURED association with a gene other than the\n      one requested. See measured_genes.\n", k))
  }
  nd <- sum(x$best_evidence == "distance_only")
  if (nd > 0 && nd / max(sum(!is.na(x$best_direction)), 1) > 0.5) {
    cat(sprintf(
      "\n  NOTE: %d of %d directions rest on distance alone. Inverse associations\n",
      nd, called),
      "  are the majority class in every catalogue, so these sit close to the\n",
      "  base rate. Expect them to be mostly -1, and read them as a prior.\n", sep = "")
  }
  if ("smr_direction" %in% names(x)) {
    ns  <- sum(!is.na(x$smr_direction))
    nsu <- sum(x$best_evidence %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak"))
    nd  <- sum(x$smr_agreement %in% FALSE)
    cat("\n  Mendelian randomisation layer:\n")
    if (ns == 0) {
      # Two different reasons produce no SMR direction, and they send you to
      # different places. Absent from the table means no instrument exists.
      # Present but unmatched means the instrument exists and points at some
      # other gene. Blaming array coverage in the second case is simply wrong,
      # and it cost an afternoon of looking in the wrong file.
      nin <- if ("smr_in_table" %in% names(x)) sum(x$smr_in_table %in% TRUE) else 0L
      if (nin > 0L && isTRUE(attr(x, "require_detected"))) {
        cat(sprintf("    no SMR direction reported: %d of these CpGs have SMR pairs, but\n", nin),
            "    require_detected removed them (not detected in child blood).\n", sep = "")
      } else if (nin == 0L) {
        cat("    no SMR evidence for any CpG here: none of them appears in the\n",
            "    table at all. GoDMC is a 450K meta-analysis, so EPIC-only probes\n",
            "    were never assayed and cannot have an mQTL.\n", sep = "")
      } else {
        cat(sprintf("    no SMR direction reported, but %d of these CpGs DO have SMR\n", nin),
            "    evidence - for genes other than the one matched here. The\n",
            "    instrument exists; it points somewhere else. See\n",
            "    cpgd_smr_directions() for the genes it does reach.\n", sep = "")
      }
    } else {
      nm_ <- if ("smr_gene_match" %in% names(x)) sum(x$smr_gene_match %in% TRUE) else ns
      cat(sprintf("    evidence for %d CpGs; %d concern the same gene as the rest of the row\n",
                  ns, nm_))
      cat(sprintf("    of those %d: supplied best_direction for %d, outranked for %d\n",
                  nm_, nsu, max(nm_ - nsu, 0)))
      if (ns - nm_ > 0)
        cat(sprintf("    %d report a causal direction for a DIFFERENT gene - see smr_gene.\n",
                    ns - nm_),
            "        Not a contradiction and not promoted to best_direction:\n",
            "        the catalogue and the instrument disagree about which gene\n",
            "        this CpG regulates, which is worth knowing on its own.\n", sep = "")
      if (nd > 0)
        cat(sprintf("    *** %d CpGs where SMR DISAGREES with the direction reported.\n", nd),
            "        A causal estimate and a correlational prediction pointing\n",
            "        opposite ways. See smr_agreement and smr_direction.\n", sep = "")
    }
  }
  if ("onco_direction" %in% names(x) && any(!is.na(x$onco_direction))) {
    no  <- sum(!is.na(x$onco_direction))
    nm_ <- if ("onco_gene_match" %in% names(x)) sum(x$onco_gene_match %in% TRUE) else no
    nou <- sum(x$best_evidence %in% "onco_consensus")
    nod <- sum(x$onco_agreement %in% FALSE)
    cat("\n  tumour consensus layer (Onco-eQTM, 27 TCGA cancer types):\n")
    cat(sprintf("    evidence for %d CpGs; %d concern the same gene as the rest of the row;\n", no, nm_))
    cat(sprintf("    tier O1 supplied best_direction for %d (onco_consensus; 0.70-0.93 validated)\n", nou))
    if (nod > 0)
      cat(sprintf("    %d CpGs where the tumour sign disagrees with the peripheral direction - see onco_agreement\n", nod))
  }
  nna <- sum(is.na(x$best_direction))
  if (nna > 0 && is.null(attr(x, "target_tissue"))) {
    cat(sprintf("\n  %d CpGs have no direction. To fill them from a tissue you name\n", nna),
        "  in advance (on biology, before seeing the answers):\n",
        "    cpg_expression_direction(my_cpgs, target_tissue = \"solid_tissue\")\n", sep = "")
  }
  nt <- sum(x$best_evidence %in% c("distance_targeted", "catalogue_targeted"))
  tt <- attr(x, "target_tissue")
  if (nt > 0 && !is.null(tt)) {
    cat(sprintf("\n  %d conflicts resolved using the target tissue you declared: %s\n", nt, tt),
        "    These rows are worth THAT TISSUE's accuracy, not the consensus figure,\n",
        "    and remain flagged in direction_uncertain. Declaring the tissue after\n",
        "    seeing which way it votes is selecting on the outcome.\n", sep = "")
  }
  nc <- sum(x$best_evidence == "distance_tissue_conflict")
  if (nc > 0 && all(c("dist_dir_blood", "dist_dir_solid") %in% names(x))) {
    k <- x[x$best_evidence == "distance_tissue_conflict", ]
    cat(sprintf("\n  %d CpGs where the distance curves DISAGREE BY TISSUE:\n", nc),
        sprintf("    blood +1 %-3d   nasal +1 %-3d   solid +1 %-3d  (of %d)\n",
                sum(k$dist_dir_blood %in% 1), sum(k$dist_dir_nasal %in% 1),
                sum(k$dist_dir_solid %in% 1), nc),
        "    These are not coin flips. Read the per-tissue columns and pick the\n",
        "    tissue your question is about.\n", sep = "")
  }
  nu <- sum(x$direction_uncertain %in% TRUE)
  if (nu > 0) {
    cat(sprintf("\n  %d CpGs have a sign that is nearly free to flip. For these,\n", nu),
        "  best_direction_filled takes the point estimate and\n",
        "  best_direction_flipped takes its opposite.\n",
        "  Run your analysis with BOTH. If the conclusion survives both, the\n",
        "  ambiguity does not matter. If it does not, the conclusion rests on\n",
        "  coin flips. Do NOT pick the column that fits better - that is\n",
        "  selecting on the outcome.\n", sep = "")
  }
  cat("\n  Read best_rung and best_accuracy alongside best_direction: in adult\n",
      "  blood a replicated measured sign is worth 0.97, distance_only 0.67\n",
      "  (cpgd_accuracy_values()). All of it is conditional on the CpG being an eQTM.\n",
      sep = "")
  ag <- attr(x, "age_group")
  if (!is.null(ag) && ag %in% names(CPGD_AGE_POPULATION)) {
    cat(sprintf("\n  Blood accuracy figures refer to: %s,\n", CPGD_AGE_POPULATION[[ag]]),
        "  conditional on the pair being coupled. Directions and rungs do not\n",
        "  change with age; see ?cpgd_smr_directions.\n", sep = "")
    if ("smr_child_blood_helix" %in% names(x)) {
      tb <- table(factor(x$smr_child_blood_helix, levels = CPGD_CHILD_BLOOD_LEVELS))
      cat(sprintf("  child_blood_helix of the SMR pairs shown: %s\n",
                  paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")))
    }
    if (isTRUE(attr(x, "require_detected")))
      cat("  require_detected: SMR pairs HELIX did not detect were removed before the ladder ran.\n")
  }
  cols <- intersect(c("cpg_id", "requested_gene", "best_direction", "best_evidence",
                      "best_confidence", "n_tissues_calling"), names(x))
  cat("\n")
  print(utils::head(as.data.frame(x[, cols, with = FALSE]), 10))
  if (nrow(x) > 10) cat(sprintf("... %d more rows\n", nrow(x) - 10))
  invisible(x)
}


# ---------------------------------------------------------------------------
# Pair-level direction resolution (2.4.0).
#
# The evidence ladder above decides ONE direction per CpG, choosing a target
# gene along the way. The pair workflow asks a different question: for a fixed
# cpg_id + target_gene pair, what direction does the evidence support? This
# resolver applies the SAME precedence as cpg_expression_direction() --
# measured > smr_high (S1) > catalogue_consensus > smr_moderate (S2) >
# catalogue_single > smr_weak (S3) > tissue_conflict > distance -- but keyed
# strictly on cpg_id + target_gene: no source may promote evidence that
# concerns a different target gene into a pair's best_direction. The result
# for a fixed pair is therefore invariant to HOW the gene entered the
# candidate set (manifest discovery, an input suffix, or an explicit genes=).
# ---------------------------------------------------------------------------

# `pairs`: data.table with cpg_id + target_gene (one row per pair).
# `sources`: optional named list overriding the packaged resources, for tests
#   and advanced use: measured, smr, lookup (a named list by tissue, or a
#   single table applied to `tissue`), cpg_pos, gene_tss, distance_curves.
# Returns `pairs` with the direction/evidence columns appended.
.cpgd_resolve_pair_direction <- function(pairs,
                                         tissue = "blood",
                                         sources = NULL,
                                         universal = TRUE,
                                         min_tissues = 2L,
                                         min_confidence = NULL,
                                         min_distance_info = 0.02,
                                         age_group = "adult",
                                         accuracy = "lower") {
  P <- data.table::copy(data.table::as.data.table(pairs))
  stopifnot(all(c("cpg_id", "target_gene") %in% names(P)))
  P <- unique(P, by = c("cpg_id", "target_gene"))
  # Genes are matched case-insensitively and with dash/underscore equivalence,
  # because catalogue spellings use dashes where panel names use underscores.
  P[, ".gkey" := gsub("_", "-", toupper(get("target_gene")))]
  n <- nrow(P)

  thrA <- if (is.null(min_confidence)) 0.30 else min_confidence
  thrB <- if (is.null(min_confidence)) 0.50 else min_confidence

  # ---- measured eQTMs, this pair only ------------------------------------
  M <- .cpgd_source_table(sources, "measured", cpgd_measured_eqtms)
  if (!is.null(M) && nrow(M)) {
    # 2.99.8: the table carries n_cohorts / cohort_agreement per tissue row
    # (eQTM Atlas build); older tables do not, and read as one cohort.
    nco <- if ("n_cohorts" %in% names(M)) suppressWarnings(as.integer(M$n_cohorts)) else rep(1L, nrow(M))
    agr <- if ("cohort_agreement" %in% names(M)) suppressWarnings(as.numeric(M$cohort_agreement)) else rep(1, nrow(M))
    MM <- data.table::data.table(
      cpg_id = M$cpg_id,
      .gkey  = gsub("_", "-", toupper(M$target_gene)),
      d      = sign(as.numeric(M$direction)),
      tis    = as.character(M$tissue),
      nco    = ifelse(is.na(nco), 1L, nco),
      agr    = ifelse(is.na(agr), 1, agr))
    MM <- MM[, list(
      d   = if (length(unique(get("d"))) == 1L) get("d")[1] else NA_real_,
      tis = paste(sort(unique(get("tis"))), collapse = ";"),
      nco = sum(get("nco")),
      agr = min(get("agr"))),
      by = c("cpg_id", ".gkey")]
    data.table::setkeyv(MM, c("cpg_id", ".gkey"))
    hit <- MM[P[, c("cpg_id", ".gkey"), with = FALSE], on = c("cpg_id", ".gkey")]
    P[, "measured_direction" := hit$d]
    P[, "measured_tissues"   := hit$tis]
    P[, "measured_n_cohorts" := ifelse(is.na(hit$d), NA_integer_, as.integer(hit$nco))]
    P[, "measured_cohort_agreement" := ifelse(is.na(hit$d), NA_real_, hit$agr)]
  } else {
    P[, c("measured_direction", "measured_tissues",
          "measured_n_cohorts", "measured_cohort_agreement") :=
        list(NA_real_, NA_character_, NA_integer_, NA_real_)]
  }

  # ---- SMR, this pair only ------------------------------------------------
  # require_detected is applied by the caller (cpg_gene_pairs), which hands the
  # filtered table in through `sources$smr`, so discovery and resolution see
  # the same SMR evidence.
  with_child <- !identical(age_group, "adult")
  S <- .cpgd_source_table(sources, "smr", .cpgd_smr_table)
  if (!is.null(S) && nrow(S)) {
    hs <- if ("heidi_status" %in% names(S)) as.character(S$heidi_status)
          else rep(NA_character_, nrow(S))
    ph <- if ("p_HEIDI" %in% names(S)) suppressWarnings(as.numeric(S$p_HEIDI))
          else rep(NA_real_, nrow(S))
    ia <- if ("instrument_agreement" %in% names(S))
            suppressWarnings(as.numeric(S$instrument_agreement))
          else rep(NA_real_, nrow(S))
    SS <- data.table::data.table(
      cpg_id = S$cpg_id,
      .gkey  = gsub("_", "-", toupper(S$target_gene)),
      d = as.numeric(S$direction), tier = as.character(S$smr_tier),
      p = as.numeric(S$p_SMR), ni = as.integer(S$n_instruments),
      ia = ia, hstat = hs, pheidi = ph,
      rung = .cpgd_smr_rung(S),
      rep  = .cpgd_smr_col(S, "replication", "eqtlgen_only"),
      outc = .cpgd_smr_col(S, "outcome", "eQTLGen"),
      cb   = .cpgd_smr_child(S))
    # one row per pair: strongest tier, then smallest p -- same rule as the
    # CpG-level ladder uses
    data.table::setorderv(SS, c("cpg_id", ".gkey", "tier", "p"))
    SS <- SS[, .SD[1L], by = c("cpg_id", ".gkey")]
    data.table::setkeyv(SS, c("cpg_id", ".gkey"))
    hit <- SS[P[, c("cpg_id", ".gkey"), with = FALSE], on = c("cpg_id", ".gkey")]
    P[, "smr_direction"        := hit$d]
    P[, "smr_tier"             := hit$tier]
    P[, "p_SMR"                := hit$p]
    P[, "n_instruments"        := hit$ni]
    P[, "instrument_agreement" := hit$ia]
    P[, "p_HEIDI"              := hit$pheidi]
    P[, "heidi_status"         := hit$hstat]
    P[, "smr_rung"             := hit$rung]
    P[, "smr_replication"      := hit$rep]
    P[, "smr_outcome"          := hit$outc]
    if (with_child) P[, "smr_child_blood_helix" := hit$cb]
  } else {
    P[, c("smr_direction", "smr_tier", "p_SMR", "n_instruments",
          "instrument_agreement", "p_HEIDI", "heidi_status",
          "smr_rung", "smr_replication", "smr_outcome") :=
        list(NA_real_, NA_character_, NA_real_, NA_integer_, NA_real_,
             NA_real_, NA_character_,
             NA_character_, NA_character_, NA_character_)]
    if (with_child) P[, "smr_child_blood_helix" := NA_character_]
  }

  # ---- tumour consensus (Onco-eQTM), this pair only ----------------------
  # Always reported; tier O1 alone may supply best_direction (see the
  # onco_consensus rung below and ?cpgd_onco_eqtm for its validation).
  OO <- .cpgd_onco_keyed(sources)
  if (!is.null(OO)) {
    hit <- OO[P[, c("cpg_id", ".gkey"), with = FALSE], on = c("cpg_id", ".gkey")]
    P[, "onco_direction"      := hit$d]
    P[, "onco_tier"           := hit$tier]
    P[, "onco_n_cancers"      := hit$nc]
    P[, "onco_sign_agreement" := hit$agr]
    P[, "onco_median_r"       := hit$mr]
  } else {
    P[, c("onco_direction", "onco_tier", "onco_n_cancers",
          "onco_sign_agreement", "onco_median_r") :=
        list(NA_real_, NA_character_, NA_integer_, NA_real_, NA_real_)]
  }

  # ---- catalogue models, per tissue, this pair only ----------------------
  short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")
  for (t in CPGD_TISSUES) {
    s <- short[[t]]
    L <- .cpgd_source_lookup(sources, t, selected = tissue)
    if (is.null(L) || !nrow(L)) {
      P[, (paste0("cat_dir_",  s)) := NA_real_]
      P[, (paste0("cat_conf_", s)) := NA_real_]
      P[, (paste0("cat_tier_", s)) := NA_character_]
      P[, (paste0("cat_prob_", s)) := NA_real_]
      P[, (paste0("cat_status_", s)) := NA_character_]
      next
    }
    LL <- data.table::data.table(
      cpg_id = L$cpg_id,
      .gkey  = gsub("_", "-", toupper(L$target_gene)),
      d = suppressWarnings(as.numeric(L$direction)),
      pr = suppressWarnings(as.numeric(L$probability_plus1)),
      cf = suppressWarnings(as.numeric(L$confidence)),
      ti = as.character(L$evidence_tier),
      st = as.character(L$status))
    # a catalogue can hold several rows for one pair; keep the most confident
    data.table::setorderv(LL, c("cpg_id", ".gkey", "cf"), order = c(1L, 1L, -1L))
    LL <- LL[, .SD[1L], by = c("cpg_id", ".gkey")]
    data.table::setkeyv(LL, c("cpg_id", ".gkey"))
    hit <- LL[P[, c("cpg_id", ".gkey"), with = FALSE], on = c("cpg_id", ".gkey")]
    # the same call rule as .cpgd_one_tissue(): tier A calls at conf >= thrA,
    # tier B at conf >= thrB, tier C abstains, DIRECT_eQTM rows alike (their
    # measured sign reaches the ladder through the measured rungs)
    call <- data.table::fcase(
      is.na(hit$d), NA_real_,
      hit$ti == "A" & !is.na(hit$cf) & hit$cf >= thrA, hit$d,
      hit$ti == "B" & !is.na(hit$cf) & hit$cf >= thrB, hit$d,
      default = NA_real_)
    P[, (paste0("cat_dir_",  s)) := call]
    P[, (paste0("cat_conf_", s)) := hit$cf]
    P[, (paste0("cat_tier_", s)) := hit$ti]
    P[, (paste0("cat_prob_", s)) := hit$pr]
    P[, (paste0("cat_status_", s)) := hit$st]
  }

  D  <- as.matrix(P[, paste0("cat_dir_",  short), with = FALSE])
  Cf <- as.matrix(P[, paste0("cat_conf_", short), with = FALSE])
  Pr <- as.matrix(P[, paste0("cat_prob_", short), with = FALSE])
  TT <- as.matrix(P[, paste0("cat_tier_", short), with = FALSE])
  n_call <- rowSums(!is.na(D))
  n_pos  <- rowSums(D ==  1, na.rm = TRUE)
  n_neg  <- rowSums(D == -1, na.rm = TRUE)
  agree  <- ifelse(n_call > 0, pmax(n_pos, n_neg) / n_call, NA_real_)
  cons   <- ifelse(n_call >= min_tissues & agree == 1,
                   ifelse(n_pos > n_neg, 1, -1), NA_real_)
  single <- vapply(seq_len(n), function(i) {
    v <- D[i, ][!is.na(D[i, ])]
    if (length(v) == 1L) v[1] else NA_real_
  }, numeric(1))
  conf_cat <- vapply(seq_len(n), function(i) {
    v <- Cf[i, ][!is.na(D[i, ])]
    v <- v[!is.na(v)]
    if (!length(v)) NA_real_ else mean(v)
  }, numeric(1))
  prob_cat <- vapply(seq_len(n), function(i) {
    v <- Pr[i, ][!is.na(D[i, ])]
    v <- v[!is.na(v)]
    if (!length(v)) NA_real_ else mean(v)
  }, numeric(1))
  worst_tier <- vapply(seq_len(n), function(i) {
    v <- TT[i, ][!is.na(D[i, ])]
    v <- v[!is.na(v)]
    if (!length(v)) NA_character_ else max(v)
  }, character(1))
  P[, "n_tissues_calling" := n_call]
  P[, "tissue_agreement"  := round(agree, 3)]

  # selected-tissue audit columns, named as the spec's pair schema expects
  ssel <- short[[tissue]]
  P[, "lookup_direction"         := P[[paste0("cat_dir_",  ssel)]]]
  P[, "lookup_probability_plus1" := P[[paste0("cat_prob_", ssel)]]]
  P[, "lookup_confidence"        := P[[paste0("cat_conf_", ssel)]]]
  P[, "lookup_evidence_tier"     := P[[paste0("cat_tier_", ssel)]]]

  # ---- distance, this pair only ------------------------------------------
  dist_p <- .cpgd_pair_distance(P, sources = sources, universal = universal)
  P[, "abs_dist"       := dist_p$abs_dist]
  P[, "p_universal"    := dist_p$p_universal]
  P[, "dir_universal"  := dist_p$dir_universal]
  P[, "dist_unanimous" := dist_p$dist_unanimous]
  P[, "dist_dir_blood" := dist_p$dist_dir_blood]
  P[, "dist_dir_nasal" := dist_p$dist_dir_nasal]
  P[, "dist_dir_solid" := dist_p$dist_dir_solid]

  # ---- the ladder, ranked by accuracy in `tissue`, within the pair --------
  # The same candidates and the same order as the CpG-level ladder
  # (.cpgd_rung_order(); R/rung_accuracy.R), resolved strictly within the pair.
  cand <- list()
  # the catalogue calls: a DIRECT_eQTM call is a measured sign; a consensus or
  # a conflict is between predicted calls (.cpgd_catalogue_cand())
  pcol <- function(pre) as.matrix(P[, paste0(pre, short), with = FALSE])
  CC <- .cpgd_catalogue_cand(pcol("cat_dir_"), pcol("cat_tier_"), pcol("cat_status_"),
                             NULL, min_tissues)
  colnames(CC$single) <- colnames(CC$direct) <- unname(short)
  TTm <- pcol("cat_tier_"); colnames(TTm) <- unname(short)
  MT <- .cpgd_measured_by_tissue(M, P$cpg_id, P$.gkey,
                                 norm = function(x) gsub("_", "-", toupper(x)))
  for (s in names(MT)) {
    dv <- MT[[s]]$d
    cand[[paste0("measured_", s, "_replicated")]] <- ifelse(MT[[s]]$rep %in% TRUE, dv, NA_real_)
    cand[[paste0("measured_", s, "_single")]]     <-
      ifelse(MT[[s]]$rep %in% TRUE, NA_real_, ifelse(is.na(dv), CC$direct[, s], dv))
  }
  smr_ok <- !is.na(P$smr_direction)
  rung <- P$smr_rung
  smr_src <- data.table::fcase(rung %in% "replicated", "smr_replicated",
                               rung %in% "high", "smr_high",
                               rung %in% "moderate", "smr_moderate",
                               rung %in% "weak", "smr_weak", default = NA_character_)
  skey <- .cpgd_smr_key(smr_src, P$smr_tier, P$smr_outcome, P$smr_replication)
  for (k in unique(stats::na.omit(skey)))
    cand[[paste0("smr_", k)]] <- ifelse(smr_ok & skey %in% k, P$smr_direction, NA_real_)
  # consensus needs two models at least, whatever min_tissues says
  for (tr in c("A", "B"))
    cand[[paste0("catalogue_consensus_", tr)]] <-
      ifelse(!is.na(CC$cons) & CC$worst %in% tr, CC$cons, NA_real_)
  for (s in c("blood", "nasal", "solid")) for (tr in c("A", "B"))
    cand[[paste0("catalogue_single_", s, "_", tr)]] <-
      ifelse(!is.na(CC$single[, s]) & TTm[, s] %in% tr, CC$single[, s], NA_real_)
  onco_ok <- !is.na(P$onco_direction) & P$onco_tier %in% "O1"
  cand[["onco_consensus"]] <- ifelse(onco_ok, P$onco_direction, NA_real_)
  dist_ok <- P$dist_unanimous %in% TRUE & !is.na(P$dir_universal) &
             abs(P$p_universal - 0.5) * 2 >= min_distance_info
  cand[["distance_only"]] <- ifelse(dist_ok %in% TRUE, P$dir_universal, NA_real_)
  conflict <- CC$conflict
  lad  <- .cpgd_apply_ladder(cand, .cpgd_rung_order(tissue), conflict, n)
  best <- lad$best; rkey <- lad$key
  i <- which(is.na(rkey) & P$dist_unanimous %in% TRUE)
  rkey[i] <- "distance_uninformative"
  i <- which(is.na(rkey) & !is.na(P$abs_dist))
  rkey[i] <- "distance_tissue_conflict"
  rkey[is.na(rkey)] <- "no_evidence"
  src <- .cpgd_rung_evidence(rkey)

  bconf <- rep(NA_real_, n)
  bconf[which(src == "measured")] <- 1
  i <- which(src %in% c("catalogue_consensus", "catalogue_single"))
  bconf[i] <- round(conf_cat[i], 4)
  i <- which(src == "distance_only")
  bconf[i] <- round(abs(P$p_universal[i] - 0.5) * 2, 4)

  prob <- rep(NA_real_, n)
  i <- which(src %in% c("catalogue_consensus", "catalogue_single"))
  prob[i] <- round(prob_cat[i], 4)
  i <- which(src %in% c("distance_only", "distance_uninformative"))
  prob[i] <- P$p_universal[i]
  i <- which(src == "measured")
  prob[i] <- ifelse(best[i] > 0, 1, 0)

  tier <- worst_tier
  i <- which(grepl("^catalogue_(consensus|single)_.*[AB]$", rkey))
  tier[i] <- substring(rkey[i], nchar(rkey[i]))
  tier[which(src == "measured")] <- "M"
  i <- which(src %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak"))
  tier[i] <- P$smr_tier[i]
  tier[which(src == "onco_consensus")] <- "O1"
  tier[which(src %in% c("distance_only", "distance_uninformative",
                        "distance_tissue_conflict"))] <- "U"

  bacc <- .cpgd_rung_label(ifelse(is.na(best), NA_character_, rkey), tissue, age_group)

  # agreement is a statement about two INDEPENDENT signs for the same pair;
  # where the tumour sign IS best_direction there is nothing to agree with
  P[, "onco_agreement" := data.table::fifelse(
        is.na(get("onco_direction")) | is.na(best) | src == "onco_consensus", NA,
        get("onco_direction") == best)]

  P[, "best_direction"  := best]
  P[, "best_evidence"   := src]
  data.table::set(P, j = "best_rung", value = rkey)
  P[, "best_confidence" := bconf]
  P[, "direction_tier"  := tier]
  P[, "probability_plus1" := prob]
  P[, "best_expected_accuracy" := bacc]
  data.table::set(P, j = "best_accuracy",
                  value = .cpgd_accuracy_number(bacc, accuracy))
  P[, "usable" := !is.na(best)]
  P[, "abstain_reason" := data.table::fcase(
      !is.na(best), "",
      src == "tissue_conflict", "catalogue tissues returned opposite directions; not resolved",
      src == "distance_uninformative", "distance curves agree in sign but sit on 0.5; no usable signal",
      src == "distance_tissue_conflict", "distance curves disagree by tissue - see dist_dir_blood/_nasal/_solid",
      src == "no_evidence", "no direction evidence held for this CpG-gene pair",
      default = "no direction evidence held for this CpG-gene pair")]
  P[, ".gkey" := NULL]
  P
}


# Resolve a source table override, falling back to the packaged resource.
# Failure to load a packaged resource degrades to NULL: a missing optional
# layer must weaken the answer, not kill the call.
.cpgd_source_table <- function(sources, name, default_fn) {
  if (!is.null(sources) && name %in% names(sources)) {
    v <- sources[[name]]
    if (is.null(v)) return(NULL)
    return(data.table::as.data.table(v))
  }
  tryCatch(default_fn(), error = function(e) NULL)
}

# The lookup override may be a named list by tissue, or a single table taken
# to be the SELECTED tissue's catalogue (the other tissues then contribute
# nothing, which is what a single-tissue fixture intends).
.cpgd_source_lookup <- function(sources, t, selected = "blood") {
  if (!is.null(sources) && "lookup" %in% names(sources)) {
    lk <- sources$lookup
    if (is.null(lk)) return(NULL)
    if (is.data.frame(lk)) {
      if (identical(t, selected)) return(data.table::as.data.table(lk))
      return(NULL)
    }
    if (!is.null(lk[[t]])) return(data.table::as.data.table(lk[[t]]))
    return(NULL)
  }
  tryCatch(.cpgd_lookup(t), error = function(e) NULL)
}


# Distance layer for explicit pairs: hg19 CpG position against the gene TSS,
# scored on the three packaged distance curves with the same unanimity rule as
# the CpG-level ladder. Anything unavailable (annotation packages, positions,
# the gene itself) degrades to NA columns rather than an error.
.cpgd_pair_distance <- function(P, sources = NULL, universal = TRUE,
                                max_dist = 1e6) {
  blank <- data.table::data.table(
    abs_dist = rep(NA_real_, nrow(P)), p_universal = NA_real_,
    dir_universal = NA_real_, dist_unanimous = NA,
    dist_dir_blood = NA_real_, dist_dir_nasal = NA_real_,
    dist_dir_solid = NA_real_)
  if (!isTRUE(universal)) return(blank)

  G <- .cpgd_source_table(sources, "gene_tss", cpgd_gene_tss)
  Cp <- .cpgd_source_table(sources, "cpg_pos", cpgd_cpg_positions)
  K <- .cpgd_source_table(sources, "distance_curves", function() {
    data.table::fread(system.file("extdata", "distance_curves.csv",
                                  package = "cpgdirection"), showProgress = FALSE)
  })
  if (is.null(G) || is.null(Cp) || is.null(K)) return(blank)
  if (!all(c("gene", "chr", "tss") %in% names(G))) return(blank)

  G <- data.table::as.data.table(G)
  G[, ".gkey" := gsub("_", "-", toupper(get("gene")))]
  G <- unique(G, by = ".gkey")
  Cp <- data.table::as.data.table(Cp)
  names(Cp) <- tolower(names(Cp))
  if (!all(c("cpg_id", "chr", "pos") %in% names(Cp))) return(blank)

  q <- data.table::data.table(cpg_id = P$cpg_id,
                              .gkey = gsub("_", "-", toupper(P$target_gene)))
  q[, ".row" := .I]
  q <- merge(q, Cp[, c("cpg_id", "chr", "pos"), with = FALSE],
             by = "cpg_id", all.x = TRUE, sort = FALSE)
  q <- merge(q, G[, c(".gkey", "chr", "tss"), with = FALSE],
             by = ".gkey", all.x = TRUE, sort = FALSE,
             suffixes = c("", "_gene"))
  q[, "abs_dist" := ifelse(!is.na(get("chr")) & !is.na(get("chr_gene")) &
                             get("chr") == get("chr_gene"),
                           abs(as.numeric(get("pos")) - as.numeric(get("tss"))),
                           NA_real_)]
  q[which(get("abs_dist") > max_dist), "abs_dist" := NA_real_]
  data.table::setorderv(q, ".row")

  out <- blank
  out[, "abs_dist" := q$abs_dist]
  short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")
  pm <- matrix(NA_real_, nrow = nrow(q), ncol = 3,
               dimnames = list(NULL, unname(short)))
  ok <- !is.na(q$abs_dist)
  if (any(ok)) {
    for (t in names(short)) {
      k <- K[K$tissue == t, ]
      if (!nrow(k)) next
      pm[ok, short[[t]]] <- stats::approx(
        log10(pmax(k$dist_mid, 1)), k$p_positive,
        xout = log10(pmax(q$abs_dist[ok], 1)), rule = 2)$y
    }
  }
  nval <- rowSums(!is.na(pm))
  npos <- rowSums(pm >= 0.5, na.rm = TRUE)
  pmean <- ifelse(nval > 0, rowMeans(pm, na.rm = TRUE), NA_real_)
  unan <- nval > 0L & (npos == nval | npos == 0L)
  out[, "p_universal"    := round(pmean, 4)]
  out[, "dist_unanimous" := ifelse(nval > 0L, unan, NA)]
  out[, "dir_universal"  := ifelse(unan & !is.na(pmean),
                                   ifelse(npos == nval & nval > 0L, 1, -1),
                                   NA_real_)]
  out[, "dist_dir_blood" := ifelse(is.na(pm[, "blood"]), NA_real_,
                                   ifelse(pm[, "blood"] >= 0.5, 1, -1))]
  out[, "dist_dir_nasal" := ifelse(is.na(pm[, "nasal"]), NA_real_,
                                   ifelse(pm[, "nasal"] >= 0.5, 1, -1))]
  out[, "dist_dir_solid" := ifelse(is.na(pm[, "solid"]), NA_real_,
                                   ifelse(pm[, "solid"] >= 0.5, 1, -1))]
  out
}
