#' Predict the direction of a CpG's effect on gene expression
#'
#' The main entry point. By default it runs every method the package has and
#' returns one row per CpG carrying all of them, plus a single
#' \code{best_direction} column holding the answer from the strongest evidence
#' available for that CpG.
#'
#' Three sources are consulted, in this order:
#' \enumerate{
#'   \item a directly measured eQTM for the requested gene, where one exists;
#'   \item the catalogue models, run in all three packaged tissues, with a
#'     consensus where they agree;
#'   \item distance to the transcription start site, which applies to any
#'     CpG-gene pair and is much weaker.
#' }
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
#'   \code{consensus_direction} is reported. Default 2.
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
#'   \code{"nasal_epithelium"} or \code{"solid_tissue"}. Where the tissues
#'   disagree, resolve using this tissue instead of abstaining. Declare it on
#'   biological grounds \strong{before} looking at the votes; the resulting rows
#'   are worth that tissue's accuracy, not the consensus figure, and are still
#'   flagged in \code{direction_uncertain}. Setting it also fills any row that
#'   would otherwise abstain with that tissue's sign, labelled
#'   \code{targeted_last_resort}, so \code{best_direction} comes back complete.
#' @param html_out Optional. \code{TRUE} writes
#'   \code{cpgdirection_results.html} into the working directory; a string is
#'   used as the file name. See \code{\link{cpgd_report}}.
#' @param verbose Print a summary. Default \code{TRUE}.
#'
#' @return With \code{tissue = "all"}, a \code{data.table} of class
#'   \code{cpgd_full}, one row per unique CpG:
#'   \describe{
#'     \item{best_direction}{+1, -1 or \code{NA}: the answer from the strongest
#'       evidence available}
#'     \item{best_evidence}{which source it came from — \code{measured},
#'       \code{smr_high}, \code{catalogue_consensus}, \code{smr_moderate},
#'       \code{catalogue_single}, \code{onco_consensus}, \code{smr_weak},
#'       \code{distance_only}, \code{distance_tissue_conflict},
#'       \code{distance_uninformative}, \code{tissue_conflict} or
#'       \code{no_evidence}, in ladder order}
#'     \item{best_confidence, best_expected_accuracy}{what that answer is worth}
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
#'       instruments survived LD pruning). Reported, never acted on: on these
#'       data HEIDI did not separate concordant from discordant directions, and
#'       filtering on it would discard the best-instrumented pairs first. See
#'       \code{?cpgd_smr_directions}}
#'     \item{smr_in_table}{whether the CpG has any SMR evidence at all, for any
#'       gene. Distinguishes "no instrument exists" from "the instrument points
#'       at a different gene", which are different facts with different remedies}
#'     \item{onco_direction, onco_tier, onco_n_cancers, onco_median_r,
#'       onco_gene, onco_gene_match, onco_agreement}{the tumour consensus layer
#'       (Onco-eQTM, 27 TCGA cancer types; \code{\link{cpgd_onco_eqtm}}),
#'       keyed and gated exactly like SMR: only a gene match may reach
#'       \code{best_direction}, and only tier O1 does so (as
#'       \code{onco_consensus}, below every single-tissue peripheral call).
#'       \code{onco_agreement} compares the tumour sign with a peripheral
#'       \code{best_direction} for the same pair; \code{NA} when the tumour
#'       sign is itself the best direction}
#'     \item{measured_genes, annotation_mismatch, mismatch_note}{measured
#'       evidence, and whether it concerns a gene other than the one requested}
#'   }
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
#' The price is explicit. Rows resolved this way carry that tissue's accuracy:
#' 0.77-0.85 for blood, 0.70-0.87 for nasal epithelium (not externally
#' validated), and 0.55-0.70 for solid tissue, which is tumour, provisional, and
#' whose brain and kidney arms did not exceed their own majority baselines under
#' leave-one-cohort-out. If solid tissue is your brain proxy, that last clause is
#' the one to quote in your limitations.
#'
#' @section How to read best_direction:
#' It is the direction supported by the best evidence this package holds for
#' that CpG, not a claim about any one tissue and not a substitute for measured
#' RNA in your own samples. Always read \code{best_evidence} alongside it: a
#' \code{measured} row and a \code{distance_only} row are not comparable, and
#' the difference in what they are worth is roughly 1.00 against 0.62.
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
                                     verbose = TRUE) {

  if (!identical(tissue, "all")) {
    r1 <- .cpgd_one_tissue(cpgs, genes = genes, tissue = tissue,
                           min_confidence = min_confidence,
                           allow_nearest_gene = allow_nearest_gene,
                           verbose = verbose)
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

  # the direction when exactly one tissue spoke, kept separate from consensus
  single <- vapply(seq_len(n), function(i) {
    v <- D[i, ][!is.na(D[i, ])]
    if (length(v) == 1L) v[1] else NA_real_
  }, numeric(1))

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
  out[, ".worst_tier" := worst_tier]
  out[, ".conf_cat"   := conf_cat]

  # ---- layer 1: directly measured associations ----------------------------
  out <- .cpgd_add_measured(out)

  req  <- toupper(out$requested_gene)
  last <- vapply(strsplit(ifelse(is.na(req), "", req), "[_-]"),
                 function(p) if (length(p)) p[length(p)] else NA_character_,
                 character(1))

  md <- rep(NA_real_, nrow(out))
  M <- tryCatch(cpgd_measured_eqtms(), error = function(e) NULL)
  if (!is.null(M)) {
    # Match the measured record to the gene actually requested, not to whatever
    # the CpG happens to have been measured against: a measured association with
    # a DIFFERENT gene is reported separately, in annotation_mismatch.
    MM <- data.table::data.table(cpg_id = M$cpg_id,
                                 g = toupper(M$target_gene),
                                 d = sign(as.numeric(M$direction)))
    MM <- MM[, list(d = if (length(unique(get("d"))) == 1L) get("d")[1] else NA_real_),
             by = c("cpg_id", "g")]
    data.table::setkeyv(MM, c("cpg_id", "g"))
    d1 <- MM[data.table::data.table(cpg_id = out$cpg_id, g = req),
             on = c("cpg_id", "g")]$d
    d2 <- MM[data.table::data.table(cpg_id = out$cpg_id, g = last),
             on = c("cpg_id", "g")]$d
    md <- ifelse(!is.na(d1), d1, d2)
  }
  out[, "measured_direction_requested" := md]

  # ---- layer 1b: causal directions from SMR --------------------------------
  # GoDMC mQTL as exposure, eQTLGen cis-eQTL as outcome. This asks whether a
  # genetic instrument for methylation moves expression, rather than whether
  # some study happened to measure the pair, so its coverage is not bounded by
  # what the catalogues tested. Median CpG-gene distance 135 kb: it reaches
  # exactly the range where the distance curves go flat.
  scols <- c("smr_direction", "smr_tier", "smr_p", "smr_n_instruments",
             "smr_gene_dist", "smr_gene", "smr_gene_match", "smr_in_table",
             "smr_heidi_status", "smr_p_heidi")
  S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
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
      hstat = hs, pheidi = ph)
    # one row per pair: strongest tier, then smallest p
    data.table::setorderv(SS, c("cpg_id", "g", "tier", "p"))
    SS <- SS[, .SD[1L], by = c("cpg_id", "g")]
    data.table::setkeyv(SS, c("cpg_id", "g"))

    # Whether the CpG appears in the SMR table AT ALL, independent of gene. The
    # print method needs this to avoid asserting a coverage explanation for a
    # CpG that is present and simply had no row for the gene we asked about.
    out[, "smr_in_table" := out$cpg_id %chin% unique(SS$cpg_id)]

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
  } else {
    out[, (scols) := list(NA_real_, NA_character_, NA_real_, NA_integer_,
                          NA_real_, NA_character_, FALSE, FALSE,
                          NA_character_, NA_real_)]
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
  # Strongest evidence wins. Order matters and is the whole point of the column.
  # which() throughout: these conditions carry NA, and a logical NA in a
  # subscripted assignment is an error rather than a skipped row.
  best <- rep(NA_real_, nrow(out))
  src  <- rep(NA_character_, nrow(out))
  wt   <- out$.worst_tier

  i <- which(!is.na(out$measured_direction_requested))
  best[i] <- out$measured_direction_requested[i];  src[i] <- "measured"

  # Only SMR rows that concern the SAME gene as the rest of the row may supply
  # best_direction. Where smr_gene_match is FALSE the estimate is still reported
  # in the smr_* columns, but it answers a question about a different pair and
  # must not be promoted into a column that answers this one.
  smr_ok <- out$smr_gene_match %in% TRUE & !is.na(out$smr_direction)

  # S1 (95.9%) outranks every predicted layer; it is beaten only by a measured
  # eQTM. Slotted by validated accuracy, not by which layer is newest.
  i <- which(is.na(src) & smr_ok & out$smr_tier == "S1")
  best[i] <- out$smr_direction[i];                 src[i] <- "smr_high"

  i <- which(is.na(src) & !is.na(out$consensus_direction) &
             out$n_tissues_calling >= min_tissues)
  best[i] <- out$consensus_direction[i];           src[i] <- "catalogue_consensus"

  # S2 (84.9%) sits alongside the catalogue models, below their consensus.
  i <- which(is.na(src) & smr_ok & out$smr_tier == "S2")
  best[i] <- out$smr_direction[i];                 src[i] <- "smr_moderate"

  i <- which(is.na(src) & out$n_tissues_calling == 1L & !is.na(out$.single))
  best[i] <- out$.single[i];                       src[i] <- "catalogue_single"

  # Onco-eQTM O1 (0.82 blood / 0.93 nasal / 0.70 SMR-S1 against the peripheral
  # references; tools/validate_onco.R) sits below every single-tissue
  # peripheral call and above S3 and the distance curves.
  onco_ok <- out$onco_gene_match %in% TRUE & !is.na(out$onco_direction) &
             out$onco_tier %in% "O1"
  i <- which(is.na(src) & onco_ok)
  best[i] <- out$onco_direction[i];                src[i] <- "onco_consensus"

  # S3 (70.4%) is weaker than a single-tissue catalogue call: its instruments
  # disagree, which is a warning rather than extra evidence.
  i <- which(is.na(src) & smr_ok & out$smr_tier == "S3")
  best[i] <- out$smr_direction[i];                 src[i] <- "smr_weak"

  # Tissues pointing opposite ways is a finding, not a gap. Falling through to
  # the distance curve here would let the weakest method settle a disagreement
  # between two stronger ones.
  i <- which(is.na(src) & out$n_tissues_calling >= 2L &
             !is.na(out$tissue_agreement) & out$tissue_agreement < 1)
  src[i] <- "tissue_conflict"

  # Distance layer, three outcomes rather than two.
  #
  # Where all three tissue curves point the same way, that agreement is the
  # evidence, and it is reported even when the pooled probability is modest.
  i <- which(is.na(src) & out$dist_unanimous %in% TRUE &
             !is.na(out$dir_universal) &
             abs(out$p_universal - 0.5) * 2 >= min_distance_info)
  best[i] <- out$dir_universal[i];                 src[i] <- "distance_only"

  # Unanimous but sitting on 0.50 in every tissue: genuinely nothing.
  i <- which(is.na(src) & out$dist_unanimous %in% TRUE)
  src[i] <- "distance_uninformative"

  # Tissues pointing opposite ways. This is the case the old pooled mean hid:
  # blood at 0.44 and solid tissue at 0.60 average to 0.52 and look like
  # ignorance, when both are speaking clearly and disagreeing. Report the
  # disagreement and expose the per-tissue signs; do not average it away.
  i <- which(is.na(src) & !is.na(out$dist_universal))
  src[i] <- "distance_tissue_conflict"

  src[is.na(src)] <- "no_evidence"

  # ---- declared target tissue ---------------------------------------------
  # Where the tissues disagree, a caller who has declared in advance which
  # tissue their question is about can have that tissue's answer instead of an
  # abstention. This is tissue SELECTION, not voting: the tissue is named before
  # the votes are seen, and its own accuracy is what the row is then worth.
  #
  # It is deliberately not a majority rule. Promoting whichever tissue happens
  # to break a tie the way you hoped is selecting on the outcome; naming the
  # tissue on biological grounds first is not.
  if (!is.null(target_tissue)) {
    tt <- match.arg(target_tissue, CPGD_TISSUES)
    ts <- short[[tt]]
    dcol <- paste0("dist_dir_", ts)
    ccol <- paste0("dir_", ts)
    if (dcol %in% names(out)) {
      i <- which(src == "distance_tissue_conflict" & !is.na(out[[dcol]]))
      best[i] <- out[[dcol]][i];  src[i] <- "distance_targeted"
    }
    if (ccol %in% names(out)) {
      i <- which(src == "tissue_conflict" & !is.na(out[[ccol]]))
      best[i] <- out[[ccol]][i];  src[i] <- "catalogue_targeted"
    }
    # Last resort: anything that would still abstain takes the target tissue's
    # sign, if that tissue has one. Naming a tissue is a statement that its
    # answer is the one you want wherever a better one is unavailable, so
    # stopping at conflicts alone would have been half a decision. Rows filled
    # here are labelled and carry no accuracy figure, because they have none:
    # they are a stated prior, not a measurement.
    if (dcol %in% names(out)) {
      i <- which(is.na(best) & !is.na(out[[dcol]]))
      best[i] <- out[[dcol]][i];  src[i] <- "targeted_last_resort"
    }
  }

  bconf <- rep(NA_real_, nrow(out))
  bconf[which(src == "measured")] <- 1
  i <- which(src %in% c("catalogue_consensus", "catalogue_single"))
  bconf[i] <- round(out$.conf_cat[i], 4)
  i <- which(src == "distance_only")
  bconf[i] <- round(abs(out$p_universal[i] - 0.5) * 2, 4)

  bacc <- rep(NA_character_, nrow(out))
  bacc[which(src == "measured")] <- "~1.00 (measured, not predicted)"
  bacc[which(src == "smr_high")]     <- "0.95-0.97 (SMR tier S1, concordant instruments; validated n=2,141)"
  bacc[which(src == "smr_moderate")] <- "0.84-0.86 (SMR tier S2, single instrument; validated n=6,008)"
  bacc[which(src == "smr_weak")]     <- "0.66-0.75 (SMR tier S3, instruments disagree; validated n=456)"
  bacc[which(src == "onco_consensus")] <- CPGD_ONCO_ACCURACY
  bacc[which(src == "catalogue_consensus" & wt == "A")] <- "0.77-0.87 (tier A, tissues agree)"
  bacc[which(src == "catalogue_consensus" & wt == "B")] <- "0.64-0.84 (tier B, tissues agree)"
  bacc[which(src == "catalogue_single"    & wt == "A")] <- "0.62-0.87 (tier A, one tissue only)"
  bacc[which(src == "catalogue_single"    & wt == "B")] <- "0.55-0.84 (tier B, one tissue only)"
  bacc[which(src == "distance_only")] <- "0.60-0.65 (distance only, tier U)"
  if (!is.null(target_tissue)) {
    tt <- match.arg(target_tissue, CPGD_TISSUES)
    tacc <- switch(tt,
      blood            = "0.77-0.85 (blood declared as target tissue)",
      nasal_epithelium = "0.70-0.87 (nasal declared as target; NOT externally validated)",
      solid_tissue     = "0.55-0.70 (solid tissue declared as target; tumour, provisional, brain and kidney below baseline under LOCO)")
    bacc[which(src %in% c("distance_targeted", "catalogue_targeted"))] <- tacc
    bacc[which(src == "targeted_last_resort")] <-
      "unquantified (last resort: the target tissue's sign where nothing better exists)"
  }

  wt[which(src %in% c("distance_only", "distance_uninformative",
                      "distance_tissue_conflict", "distance_targeted",
                      "targeted_last_resort"))] <- "U"
  wt[which(src == "measured")] <- "M"
  i <- which(src %in% c("smr_high", "smr_moderate", "smr_weak"))
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
  out[, "best_confidence" := bconf]
  out[, "best_expected_accuracy" := bacc]
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
        "LAST RESORT: no measured, SMR, catalogue or unanimous distance evidence; this is the target tissue's distance sign and nothing more"]
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
        "causal direction from SMR (blood); consistent with methylation affecting expression AND with linkage - HEIDI not run"]
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
           "best_tier", "best_expected_accuracy",
           "direction_uncertain", "best_direction_filled", "best_direction_flipped",
           "smr_direction", "smr_gene", "smr_gene_match", "smr_tier", "smr_p",
           "smr_n_instruments", "smr_gene_dist", "smr_agreement", "smr_in_table",
           "smr_heidi_status", "smr_p_heidi",
           "onco_direction", "onco_tier", "onco_n_cancers", "onco_median_r",
           "onco_gene", "onco_gene_match", "onco_agreement",
           "consensus_direction", "n_tissues_calling", "tissue_agreement",
           as.vector(rbind(paste0("gene_", short), paste0("dir_", short),
                           paste0("conf_", short), paste0("tier_", short),
                           paste0("status_", short))),
           ucols,
           "measured_genes", "measured_tissues", "measured_direction",
           "measured_direction_requested",
           "annotation_mismatch", "mismatch_note", "note")
  drop <- intersect(c(".single", ".worst_tier", ".conf_cat"), names(out))
  if (length(drop)) out[, (drop) := NULL]
  out <- out[, intersect(ord, names(out)), with = FALSE]
  data.table::setattr(out, "target_tissue",
                      if (is.null(target_tissue)) NULL else
                        match.arg(target_tissue, CPGD_TISSUES))
  data.table::setattr(out, "class", c("cpgd_full", class(out)))
  if (isTRUE(verbose)) print(out)
  if (!is.null(html_out)) .cpgd_maybe_report(out, html_out)
  invisible(out)
}


# html_out = TRUE writes a default name into getwd(); a string is used as the
# file name. Reporting must never take down a completed analysis, so failures
# here warn rather than stop.
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
    nsu <- sum(x$best_evidence %in% c("smr_high", "smr_moderate", "smr_weak"))
    nd  <- sum(x$smr_agreement %in% FALSE)
    cat("\n  Mendelian randomisation layer:\n")
    if (ns == 0) {
      # Two different reasons produce no SMR direction, and they send you to
      # different places. Absent from the table means no instrument exists.
      # Present but unmatched means the instrument exists and points at some
      # other gene. Blaming array coverage in the second case is simply wrong,
      # and it cost an afternoon of looking in the wrong file.
      nin <- if ("smr_in_table" %in% names(x)) sum(x$smr_in_table %in% TRUE) else 0L
      if (nin == 0L) {
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
  cat("\n  Read best_evidence alongside best_direction: measured is worth ~1.00,\n",
      "  distance_only ~0.62. All of it is conditional on the CpG being an eQTM.\n", sep = "")
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
                                         min_distance_info = 0.02) {
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
    MM <- data.table::data.table(
      cpg_id = M$cpg_id,
      .gkey  = gsub("_", "-", toupper(M$target_gene)),
      d      = sign(as.numeric(M$direction)),
      tis    = as.character(M$tissue))
    MM <- MM[, list(
      d   = if (length(unique(get("d"))) == 1L) get("d")[1] else NA_real_,
      tis = paste(sort(unique(get("tis"))), collapse = ";")),
      by = c("cpg_id", ".gkey")]
    data.table::setkeyv(MM, c("cpg_id", ".gkey"))
    hit <- MM[P[, c("cpg_id", ".gkey"), with = FALSE], on = c("cpg_id", ".gkey")]
    P[, "measured_direction" := hit$d]
    P[, "measured_tissues"   := hit$tis]
  } else {
    P[, c("measured_direction", "measured_tissues") :=
        list(NA_real_, NA_character_)]
  }

  # ---- SMR, this pair only ------------------------------------------------
  S <- .cpgd_source_table(sources, "smr", cpgd_smr_directions)
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
      ia = ia, hstat = hs, pheidi = ph)
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
  } else {
    P[, c("smr_direction", "smr_tier", "p_SMR", "n_instruments",
          "instrument_agreement", "p_HEIDI", "heidi_status") :=
        list(NA_real_, NA_character_, NA_real_, NA_integer_, NA_real_,
             NA_real_, NA_character_)]
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
    # tier B at conf >= thrB, tier C abstains; DIRECT_eQTM rows always call
    call <- data.table::fcase(
      is.na(hit$d), NA_real_,
      hit$st == "DIRECT_eQTM", hit$d,
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

  # ---- the ladder, unchanged in order, applied within the pair -----------
  best <- rep(NA_real_, n)
  src  <- rep(NA_character_, n)

  i <- which(!is.na(P$measured_direction))
  best[i] <- P$measured_direction[i];  src[i] <- "measured"

  smr_ok <- !is.na(P$smr_direction)
  i <- which(is.na(src) & smr_ok & P$smr_tier == "S1")
  best[i] <- P$smr_direction[i];       src[i] <- "smr_high"

  i <- which(is.na(src) & !is.na(cons))
  best[i] <- cons[i];                  src[i] <- "catalogue_consensus"

  i <- which(is.na(src) & smr_ok & P$smr_tier == "S2")
  best[i] <- P$smr_direction[i];       src[i] <- "smr_moderate"

  i <- which(is.na(src) & n_call == 1L & !is.na(single))
  best[i] <- single[i];                src[i] <- "catalogue_single"

  # Onco-eQTM O1 (>= 5 cancer types, one sign): 0.82 blood / 0.93 nasal / 0.70
  # SMR-S1 against the peripheral references (tools/validate_onco.R). Below
  # every single-tissue peripheral call, above S3 and the distance curves.
  onco_ok <- !is.na(P$onco_direction) & P$onco_tier %in% "O1"
  i <- which(is.na(src) & onco_ok)
  best[i] <- P$onco_direction[i];      src[i] <- "onco_consensus"

  i <- which(is.na(src) & smr_ok & P$smr_tier == "S3")
  best[i] <- P$smr_direction[i];       src[i] <- "smr_weak"

  i <- which(is.na(src) & n_call >= 2L & !is.na(agree) & agree < 1)
  src[i] <- "tissue_conflict"

  i <- which(is.na(src) & P$dist_unanimous %in% TRUE &
               !is.na(P$dir_universal) &
               abs(P$p_universal - 0.5) * 2 >= min_distance_info)
  best[i] <- P$dir_universal[i];       src[i] <- "distance_only"

  i <- which(is.na(src) & P$dist_unanimous %in% TRUE)
  src[i] <- "distance_uninformative"

  i <- which(is.na(src) & !is.na(P$abs_dist))
  src[i] <- "distance_tissue_conflict"

  src[is.na(src)] <- "no_evidence"

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
  tier[which(src == "measured")] <- "M"
  i <- which(src %in% c("smr_high", "smr_moderate", "smr_weak"))
  tier[i] <- P$smr_tier[i]
  tier[which(src == "onco_consensus")] <- "O1"
  tier[which(src %in% c("distance_only", "distance_uninformative",
                        "distance_tissue_conflict"))] <- "U"

  bacc <- rep(NA_character_, n)
  bacc[which(src == "measured")]     <- "~1.00 (measured, not predicted)"
  bacc[which(src == "smr_high")]     <- "0.95-0.97 (SMR tier S1, concordant instruments; validated n=2,141)"
  bacc[which(src == "smr_moderate")] <- "0.84-0.86 (SMR tier S2, single instrument; validated n=6,008)"
  bacc[which(src == "smr_weak")]     <- "0.66-0.75 (SMR tier S3, instruments disagree; validated n=456)"
  bacc[which(src == "catalogue_consensus" & worst_tier == "A")] <- "0.77-0.87 (tier A, tissues agree)"
  bacc[which(src == "catalogue_consensus" & worst_tier == "B")] <- "0.64-0.84 (tier B, tissues agree)"
  bacc[which(src == "catalogue_single"    & worst_tier == "A")] <- "0.62-0.87 (tier A, one tissue only)"
  bacc[which(src == "catalogue_single"    & worst_tier == "B")] <- "0.55-0.84 (tier B, one tissue only)"
  bacc[which(src == "onco_consensus")] <- CPGD_ONCO_ACCURACY
  bacc[which(src == "distance_only")] <- "0.60-0.65 (distance only, tier U)"

  # agreement is a statement about two INDEPENDENT signs for the same pair;
  # where the tumour sign IS best_direction there is nothing to agree with
  P[, "onco_agreement" := data.table::fifelse(
        is.na(get("onco_direction")) | is.na(best) | src == "onco_consensus", NA,
        get("onco_direction") == best)]

  P[, "best_direction"  := best]
  P[, "best_evidence"   := src]
  P[, "best_confidence" := bconf]
  P[, "direction_tier"  := tier]
  P[, "probability_plus1" := prob]
  P[, "best_expected_accuracy" := bacc]
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
