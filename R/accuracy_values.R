# One number per accuracy string, for engines that weight by it.
#
# Every best_expected_accuracy / expected_accuracy / requested_gene_best_accuracy
# string the package emits has one row here. The ladder strings are built from
# the per-tissue validation (CPGD_RUNG_ACCURACY, R/rung_accuracy.R) and the
# brain strings from the cortex validation (CPGD_BRAIN_TARGET_ACCURACY); the
# numeric columns (best_accuracy, accuracy, requested_gene_accuracy) are looked
# up here by that exact string, so number and text cannot describe different
# validations.
#
#   lower  - the lower bound of the Wilson 95% interval of the source's
#            agreement with the tissue's independent yardstick (for the brain
#            T1/T2 rows, the cluster bootstrap interval). The default.
#   median - the point estimate of the same agreement.
# Built at run time, so it does not depend on the order R collates the files.

.cpgd_acc_row <- function(label, rung, population, median, lower, upper, n,
                          interval, basis) {
  data.frame(rung = rung, population = population,
             median = median, lower = lower, upper = upper,
             n = n, interval = interval, basis = basis, label = label,
             stringsAsFactors = FALSE)
}

.cpgd_accuracy_values_table <- function() {
  r <- .cpgd_acc_row
  W <- "Wilson 95%"
  rows <- list()
  # ---- the ladder sources, every tissue and age they are quoted for -------
  combos <- rbind(data.frame(target = "blood", age = c("adult", "child")),
                  data.frame(target = c("nasal_epithelium", "solid_tissue"), age = "adult"))
  for (s in names(CPGD_RUNG_DESC)) for (j in seq_len(nrow(combos))) {
    f <- .cpgd_rung_figure1(s, combos$target[j], combos$age[j])
    if (is.null(f)) {
      rows[[length(rows) + 1L]] <- r(.cpgd_rung_text(s, NULL), s, "any", 0.5, 0.5,
        NA_real_, "0", "none (no validation pair)",
        "No yardstick pair: the neutral weight 0.5, a sign with no evidence it beats chance")
      next
    }
    rows[[length(rows) + 1L]] <- r(
      .cpgd_rung_text(s, f), s, f$population, f$agree, f$lo, f$hi,
      format(f$n, big.mark = ","),
      if (identical(f$basis, "own")) W else "Wilson 95%, all yardsticks pooled",
      if (identical(f$basis, "own")) sprintf("%s vs %s", CPGD_RUNG_DESC[[s]], f$yardstick)
      else sprintf("%s; only %s pair(s) in %s, so the figure pools every yardstick",
                   CPGD_RUNG_DESC[[s]], f$n_here, f$population))
  }
  # ---- the brain-target ladder (cpg_brain_bridge), adult cortex ----------
  rows <- c(rows, list(
    r(CPGD_BRAIN_TARGET_ACCURACY[["measured_brain_replicated"]], "measured_brain_replicated",
      "adult cortex", 0.983, 0.929, NA_real_, "180", "derived from a Wilson 95% count",
      paste("Bounded below by one cortex cohort's agreement with the other (ROSMAP vs",
            "Jaffe, 174 of 180, 0.967 [0.929, 0.985]); median: the single-cohort",
            "accuracy that agreement implies, (1 + sqrt(2p - 1)) / 2")),
    r(CPGD_BRAIN_TARGET_ACCURACY[["measured_brain_single"]], "measured_brain_single",
      "adult cortex", 0.967, 0.929, 0.985, "180", W,
      "ROSMAP and Jaffe agree on 174 of 180 shared pairs"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["brain_T1"]], "brain_T1", "adult cortex",
      0.964, 0.926, 0.998, NA_character_, "cluster bootstrap 95% (by CpG)",
      "Agreement expected against a donor-independent adult cortex cohort (tools/validate_brain_overlap.R)"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["brain_T2"]], "brain_T2", "adult cortex",
      0.943, 0.899, 0.984, NA_character_, "cluster bootstrap 95% (by CpG)",
      "Agreement expected against a donor-independent adult cortex cohort (tools/validate_brain_overlap.R)"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["measured_blood"]], "brain_fallback_measured_blood",
      "adult cortex", 0.840, 0.822, 0.857, "1,654", W,
      "measured blood eQTMs vs ROSMAP cortex (tools/reports/validate_ladder_transport.md)"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["blood_SMR_replicated"]], "brain_fallback_blood_SMR_replicated",
      "adult cortex", 0.772, 0.745, 0.796, "1,055", W, "blood SMR replicated vs ROSMAP cortex"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["blood_SMR_high"]], "brain_fallback_blood_SMR_high",
      "adult cortex", 0.739, 0.676, 0.792, "218", W, "blood SMR high vs ROSMAP cortex"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["blood_SMR_weak"]], "brain_fallback_blood_SMR_weak",
      "adult cortex", 0.725, 0.626, 0.806, "91", W, "blood SMR weak vs ROSMAP cortex"),
    r(CPGD_BRAIN_TARGET_ACCURACY[["blood_SMR_moderate"]], "brain_fallback_blood_SMR_moderate",
      "adult cortex", 0.666, 0.643, 0.688, "1,693", W, "blood SMR moderate vs ROSMAP cortex")))
  V <- do.call(rbind, rows)
  V <- V[!duplicated(V$label), , drop = FALSE]
  rownames(V) <- NULL
  V
}

#' One accuracy number per evidence rung, for weighting engines
#'
#' The table behind the numeric accuracy columns: \code{best_accuracy} in
#' \code{\link{cpg_expression_direction}} and \code{\link{cpg_gene_pairs}},
#' \code{accuracy} in the single-tissue result and in
#' \code{\link{cpg_direction_universal}}, and
#' \code{requested_gene_accuracy} in \code{\link{cpg_brain_bridge}}. Each row
#' is one accuracy string the package quotes (\code{label}, exactly as it
#' appears in \code{best_expected_accuracy}, \code{expected_accuracy} or
#' \code{requested_gene_best_accuracy}), with the validation it comes from.
#'
#' @section Which value is returned:
#' \describe{
#'   \item{\code{"lower"} (default)}{the lower bound of the Wilson 95\%
#'     interval of the source's agreement with the tissue's independent
#'     yardstick (the cluster-bootstrap interval for brain T1/T2). The number
#'     a call returns is never below 0.5: where a small sample's bound dips
#'     under it, the conservative number is "uninformative", not a weight that
#'     would reverse the sign. The \code{lower} column of this table shows the
#'     bound as computed.}
#'   \item{\code{"median"}}{the point estimate of the same agreement. For a
#'     cortex pair measured in both cohorts, which no third cohort scores, it
#'     is the single-cohort accuracy implied by their agreement \eqn{p},
#'     \eqn{(1 + \sqrt{2p - 1}) / 2}, not 1, so no weight becomes certainty.}
#' }
#' A source with fewer than 30 yardstick pairs in a tissue carries its
#' cross-tissue figure (every yardstick pooled, whatever its n; the interval
#' shows it) instead of its own. A source no yardstick scores at all is quoted
#' as "unquantified" and carries the neutral weight 0.5 (lower and median
#' alike): a sign, but no evidence it beats chance. Child
#' blood carries the HELIX figure where HELIX scores the source on 100 or more
#' pairs, the adult (MESA) figure otherwise. Choose the value with the
#' \code{accuracy} argument of the direction functions, or for a whole session
#' (including calls an engine such as DMSA makes on your behalf) with
#' \code{options(cpgdirection.accuracy = "median")}. With
#' \code{age_group = "any"} the adult and the child blood figure are both
#' quoted and the number is the smaller of the two.
#'
#' @section Where the figures come from:
#' \code{tools/validate_rung_accuracy.R} scores every evidence source in each
#' tissue against a measured yardstick independent of it: MESA monocytes and T
#' cells for adult blood (not a catalogue training cohort, not an SMR arm);
#' HELIX child blood pairs outside the blood catalogue's training set; EVA-PR
#' nasal epithelium pairs outside the nasal catalogue's training set (the nasal
#' catalogue is trained on EVA-PR, so its nasal figure is same-study); FUSION
#' muscle and GTEx solid tissues, each held out in turn. The same figures order
#' the ladder in each tissue (see \code{\link{cpg_expression_direction}}).
#' Agreement with a yardstick that itself carries error sits a little below
#' the probability of being right; every source in a tissue faces the same
#' yardstick, so their order is unaffected.
#'
#' @section What the number means:
#' The probability that the called sign is right, given that the CpG-gene pair
#' is coupled in the tissue (an SMR direction says which way a CpG would act if
#' it acts on the gene, not that it does). For an engine that weights by sign
#' certainty: \eqn{P(d = +1)} is \code{best_accuracy} where
#' \code{best_direction} is \code{+1} and \code{1 - best_accuracy} where it is
#' \code{-1}, and the expected-sign multiplier \eqn{2P(d = +1) - 1} is
#' \code{best_direction * (2 * best_accuracy - 1)}. Every row with a
#' direction has a number; only an abstention is \code{NA}.
#'
#' @return A data.frame: \code{rung} (the source, as in \code{best_rung}),
#'   \code{population}, \code{median}, \code{lower}, \code{upper},
#'   \code{n} (yardstick pairs), \code{interval}, \code{basis} and
#'   \code{label}.
#' @examples
#' v <- cpgd_accuracy_values()
#' v[v$population == "adult blood", c("rung", "median", "lower", "n")]
#' @export
cpgd_accuracy_values <- function() .cpgd_accuracy_values_table()

# The accuracy argument: "lower" (default) or "median"; the session option
# cpgdirection.accuracy supplies the default.
.cpgd_accuracy_arg <- function(accuracy) {
  if (is.null(accuracy)) accuracy <- "lower"
  if (!is.character(accuracy) || length(accuracy) != 1L || is.na(accuracy) ||
      !accuracy %in% c("lower", "median")) {
    stop('`accuracy` must be "lower" or "median" (or set ',
         'options(cpgdirection.accuracy = ...) to one of them).', call. = FALSE)
  }
  accuracy
}

# Accuracy strings -> one number each. With age_group = "any" a string holds
# the adult and the child text joined by " | "; the number is the smaller.
# A string with no row in the table gives NA (the tests forbid that for
# every string the package itself emits).
.cpgd_accuracy_number <- function(labels, accuracy = "lower") {
  out <- rep(NA_real_, length(labels))
  ok <- which(!is.na(labels))
  if (!length(ok)) return(out)
  V <- .cpgd_accuracy_values_table()
  val <- V[[accuracy]]
  u <- unique(labels[ok])
  num <- vapply(u, function(s) {
    v <- val[match(strsplit(s, " | ", fixed = TRUE)[[1]], V$label)]
    if (!length(v) || anyNA(v)) NA_real_ else min(v)
  }, numeric(1), USE.NAMES = FALSE)
  out[ok] <- num[match(labels[ok], u)]
  # A sign the ladder supplies has a point estimate of at least 0.5 in its
  # tissue (sources below that never answer there); where the lower bound of a
  # small sample dips under 0.5, the conservative number is "uninformative",
  # 0.5, never a weight that would reverse the sign.
  pmax(out, 0.5)
}
