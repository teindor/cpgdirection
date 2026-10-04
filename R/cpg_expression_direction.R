#' Predict the direction of a CpG's effect on gene expression
#'
#' Takes CpG identifiers and returns, for each, the predicted direction of the
#' association between methylation at that site and expression of a nearby gene:
#' \code{+1} means higher methylation goes with higher expression, \code{-1}
#' means higher methylation goes with lower expression.
#'
#' Identifiers are stripped to their canonical form before lookup, so EPIC v2
#' suffixes are handled automatically: \code{cg00000029_TC21} and
#' \code{cg00050692_TC21_DNMT3A} both resolve to their \code{cg} identifier. In
#' the second case the gene is parsed from the name and used as the target,
#' which is usually what you want.
#'
#' No model is fitted at run time. The function is a lookup against a
#' precomputed table, so results are identical on every call.
#'
#' @param cpgs Character vector of CpG identifiers, with or without suffixes.
#'   A single file path to a text or CSV file is also accepted.
#' @param genes Optional character vector of target genes, the same length as
#'   \code{cpgs}. Supply this when you know which gene you care about. Without
#'   it the function falls back to the nearest transcription start site, which
#'   is frequently not the intended gene; the \code{gene_source} column records
#'   which path was taken.
#' @param tissue Tissue context. Only \code{"blood"} is packaged in this
#'   version; the table is derived from blood catalogues and validated in blood
#'   and epithelium. Passing another value raises an error rather than silently
#'   returning blood predictions.
#' @param min_confidence Optional numeric in \code{[0, 1]}. Rows below this are
#'   returned with \code{call = "ABSTAIN"}. Defaults to the tier-specific
#'   thresholds used in the paper (0.30 for tier A, 0.50 for tier B; tier C
#'   always abstains).
#' @param verbose Print a short summary. Default \code{TRUE}.
#'
#' @return A \code{data.table} of class \code{cpgd_result}, one row per input
#'   CpG, with columns:
#'   \describe{
#'     \item{input}{the identifier as supplied}
#'     \item{cpg_id}{canonical CpG identifier used for lookup}
#'     \item{target_gene}{gene the prediction refers to}
#'     \item{gene_source}{\code{user_supplied}, \code{parsed_from_name},
#'       \code{catalogue} or \code{nearest_TSS}}
#'     \item{tss_dist}{distance from CpG to the gene's transcription start site}
#'     \item{status}{\code{DIRECT_eQTM}, \code{PREDICTED}, \code{ABSTAIN},
#'       \code{NO_EXPRESSION_ANCHOR} or \code{UNMAPPED}}
#'     \item{direction}{\code{+1}, \code{-1} or \code{NA}}
#'     \item{direction_label}{human-readable form of \code{direction}}
#'     \item{probability_plus1}{calibrated probability of a positive association}
#'     \item{confidence}{0-1; distance of the probability from 0.5, doubled}
#'     \item{evidence_tier}{\code{A}, \code{B} or \code{C}; see Details}
#'     \item{expected_accuracy}{what a direction at this tier is worth}
#'     \item{call}{\code{+1}, \code{-1} or \code{ABSTAIN} after thresholds}
#'     \item{note}{caveat in plain language, where one applies}
#'   }
#'
#' @details
#' \strong{Evidence tiers.} Tier A: catalogued CpGs lie within 50 kb and the
#' target gene is itself catalogued. Tier B: one of those holds. Tier C:
#' neither, and the function always abstains.
#'
#' \strong{What a direction is worth.} Tier A predictions were 76.6\% accurate
#' on a pre-registered locked test set (baseline 59.2\%) and 85.4\% in an
#' independent white-blood-cell cohort. Tier B is weaker. See
#' \code{\link{cpgd_accuracy_table}}.
#'
#' \strong{The limitation that governs everything.} Every figure above is
#' conditional on the CpG being an eQTM. The model gives the direction of an
#' effect; it cannot tell you whether an effect exists, because the catalogues
#' it learned from contain only significant associations. \code{ABSTAIN} and
#' \code{NO_EXPRESSION_ANCHOR} are correct answers, not failures.
#'
#' \strong{Genome build.} The lookup table is hg19. EPIC v2 manifests ship
#' GRCh38, so if you are matching by position rather than by probe identifier,
#' lift over first.
#'
#' @noRd
.cpgd_one_tissue <- function(cpgs,
                             genes = NULL,
                             tissue = c("blood", "nasal_epithelium", "solid_tissue"),
                             min_confidence = NULL,
                             allow_nearest_gene = FALSE,
                             verbose = TRUE) {

  tissue <- match.arg(tissue, CPGD_TISSUES)

  q <- .cpgd_parse_input(cpgs, genes)
  if (!nrow(q)) {
    ex <- utils::head(as.character(if (is.data.frame(cpgs)) cpgs[[1]] else cpgs), 3)
    stop("No valid CpG identifiers found in the input. Identifiers must contain ",
         "a 'cg' followed by at least six digits, anywhere in the string. ",
         "First values seen: ", paste(sQuote(ex), collapse = ", "),
         ". Use cpgd_parse_check() to see what the parser extracted.", call. = FALSE)
  }
  n_raw <- attr(q, "n_raw")
  n_in <- if (is.null(n_raw)) nrow(q) else n_raw
  n_dup <- n_in - nrow(q)

  L <- .cpgd_lookup(tissue)
  hit <- L[q, on = "cpg_id", allow.cartesian = TRUE, nomatch = NA]

  # A right join fills non-matching rows with NA, so an identifier absent from
  # the table looks exactly like one present but with no gene nearby. They are
  # different facts and a user needs to tell them apart: the first usually means
  # a typo or a probe from another array, the second is a real biological
  # answer. Capture it now, before `status` is overwritten below.
  hit[, "cpgd_absent" := is.na(get("status"))]

  # A supplied gene may appear with underscores where the catalogue uses dashes,
  # and panel names often carry a compound tail whose informative symbol is the
  # last token. Accept any of these forms.
  hit[, "g_dash"  := gsub("_", "-", get("given_gene"))]
  hit[, "g2_dash" := gsub("_", "-", get("given_gene2"))]
  hit[, "match_user" := !is.na(get("target_gene")) & (
        (!is.na(get("given_gene"))  & (get("target_gene") == get("given_gene")  |
                                       get("target_gene") == get("g_dash"))) |
        (!is.na(get("given_gene2")) & (get("target_gene") == get("given_gene2") |
                                       get("target_gene") == get("g2_dash"))))]
  hit[is.na(get("match_user")), "match_user" := FALSE]

  # one row per CpG: the requested gene if we hold it, else the nearest TSS
  hit[, "gene_rank" := fifelse(get("match_user"), 0L, 1L)]
  hit[, "adist" := abs(as.numeric(get("tss_dist")))]
  hit[is.na(get("adist")), "adist" := Inf]
  setorderv(hit, c("cpg_id", "gene_rank", "adist"))
  out <- hit[, .SD[1L], by = "cpg_id"]

  out[, "gene_source" := fifelse(
        is.na(get("target_gene")), NA_character_,
        fifelse(get("match_user"),
                fifelse(get("gene_origin") == "parsed", "parsed_from_name", "user_supplied"),
                fifelse(get("status") == "DIRECT_eQTM", "catalogue", "nearest_TSS")))]

  out[is.na(get("target_gene")), "status" := "NO_EXPRESSION_ANCHOR"]
  out[get("cpgd_absent") %in% TRUE, "status" := "UNMAPPED"]

  # If a gene was requested but the table holds no pairing for it, the nearest
  # catalogued gene is a DIFFERENT gene. Returning a direction for it would
  # answer a question nobody asked: a CpG near POMC pairs with DNMT3A because
  # they are 50 kb apart, and a reader skimming the output would take the
  # direction as being about POMC. Refuse by default.
  out[, "requested_gene" := get("given_gene")]
  out[!is.na(get("requested_gene")) & !get("match_user") &
        !is.na(get("target_gene")), "status" := "REQUESTED_GENE_UNAVAILABLE"]

  thrA <- if (is.null(min_confidence)) 0.30 else min_confidence
  thrB <- if (is.null(min_confidence)) 0.50 else min_confidence
  out[, "call" := fcase(
        get("status") == "REQUESTED_GENE_UNAVAILABLE" & !isTRUE(allow_nearest_gene), "ABSTAIN",
        is.na(get("direction")), "ABSTAIN",
        get("evidence_tier") == "A" & get("confidence") >= thrA, as.character(get("direction")),
        get("evidence_tier") == "B" & get("confidence") >= thrB, as.character(get("direction")),
        default = "ABSTAIN")]
  out[get("call") == "ABSTAIN" & get("status") == "PREDICTED", "status" := "ABSTAIN"]

  out[, "direction_label" := fcase(
        get("direction") == 1,  "higher methylation -> HIGHER expression",
        get("direction") == -1, "higher methylation -> LOWER expression",
        default = NA_character_)]

  acc <- switch(tissue,
    blood            = c(A = "0.77-0.85", B = "0.64-0.85"),
    nasal_epithelium = c(A = "0.84-0.87", B = "0.70-0.84"),
    solid_tissue     = c(A = "0.62-0.70", B = "0.55-0.66"))
  out[, "expected_accuracy" := fcase(
        get("status") == "DIRECT_eQTM", "~1.00 (measured, not predicted)",
        get("status") == "PREDICTED" & get("evidence_tier") == "A", acc[["A"]],
        get("status") == "PREDICTED" & get("evidence_tier") == "B", acc[["B"]],
        default = NA_character_)]

  out[, "note" := ""]
  out[get("status") == "NO_EXPRESSION_ANCHOR", "note" :=
        "no gene with an expression anchor within 250 kb"]
  out[get("status") == "UNMAPPED", "note" :=
        "identifier not present in the lookup table - check the probe ID and the array version"]
  out[get("status") == "ABSTAIN", "note" := "below the confidence threshold; no call made"]
  out[get("gene_source") == "nearest_TSS" & get("status") != "NO_EXPRESSION_ANCHOR",
      "note" := "gene assigned by nearest TSS, not a measured association - verify"]
  out[get("evidence_tier") == "C", "note" := "uncharted locus; tier C always abstains"]
  out[get("status") == "REQUESTED_GENE_UNAVAILABLE", "note" := paste0(
        "no prediction available for the requested gene; nearest catalogued gene is ",
        get("target_gene"), ". Set allow_nearest_gene = TRUE to use it, but it is a ",
        "DIFFERENT gene")]

  keep <- c("input", "cpg_id", "tissue", "requested_gene", "target_gene", "gene_source", "tss_dist", "status",
            "direction", "direction_label", "probability_plus1", "confidence",
            "evidence_tier", "expected_accuracy", "call", "note")
  unmapped <- setdiff(q$cpg_id, out$cpg_id)
  out <- out[, intersect(keep, names(out)), with = FALSE]
  if (length(unmapped)) {
    u <- q[q$cpg_id %in% unmapped, c("input", "cpg_id"), with = FALSE]
    u[, "status" := "UNMAPPED"]
    u[, "note" := "identifier not present in the lookup table"]
    out <- data.table::rbindlist(list(out, u), fill = TRUE)
  }

  data.table::setattr(out, "class", c("cpgd_result", class(out)))
  data.table::setattr(out, "n_submitted", n_in)
  data.table::setattr(out, "n_duplicate", n_dup)
  data.table::setattr(out, "tissue", tissue)
  out[, "tissue" := tissue]
  if (isTRUE(verbose)) print(out)
  invisible(out)
}


#' Expected accuracy by evidence tier
#'
#' The validation figures behind the \code{expected_accuracy} column, so they can
#' be inspected rather than taken on trust.
#'
#' @return A \code{data.frame}.
#' @examples
#' cpgd_accuracy_table()
#' @export
cpgd_accuracy_table <- function() {
  data.frame(
    tissue = c(rep("blood", 3), rep("nasal_epithelium", 2), rep("solid_tissue", 2),
               "any (tumour consensus)", "any"),
    status = c("DIRECT_eQTM", "PREDICTED (tier A)", "PREDICTED (tier B)",
               "PREDICTED (tier A)", "PREDICTED (tier B)",
               "PREDICTED (tier A)", "PREDICTED (tier B)",
               "onco_consensus (Onco-eQTM O1)", "majority baseline"),
    expected_accuracy = c("~1.00", "0.77-0.85", "0.64-0.85",
                          "0.84-0.87", "0.70-0.84",
                          "0.62-0.70", "0.55-0.66", "0.70-0.93", "0.55-0.63"),
    basis = c("read from a catalogue, not predicted",
              "HELIX locked test 0.766; EVA-PR white blood cell 0.854",
              "HELIX locked test 0.642; EVA-PR white blood cell 0.854",
              "grouped CV AUC 0.865; NOT independently validated (both cohorts, one study)",
              "grouped CV AUC 0.865; NOT independently validated",
              "grouped CV AUC 0.664; tumour tissue; brain and kidney below baseline",
              "grouped CV AUC 0.664; provisional",
              "agreement with measured blood 0.82 (n=1,455), nasal 0.93 (n=1,221), SMR S1 0.70 (n=1,033); tools/validate_onco.R",
              "always predicting an inverse association"),
    stringsAsFactors = FALSE)
}


#' Information about the bundled lookup table
#'
#' @param tissue Which tissue's table to describe: \code{"blood"},
#'   \code{"nasal_epithelium"} or \code{"solid_tissue"}.
#' @return A named list: tissue, number of pairs, CpGs and genes, genome build,
#'   file path and package version.
#' @examplesIf cpgd_has_data("lookup_blood_hg19")
#' cpgd_lookup_info("blood")
#' @export
cpgd_lookup_info <- function(tissue = c("blood", "nasal_epithelium", "solid_tissue")) {
  tissue <- match.arg(tissue)
  L <- .cpgd_lookup(tissue)
  list(tissue = tissue,
       pairs = nrow(L),
       cpgs = length(unique(L$cpg_id)),
       genes = length(unique(L$target_gene)),
       genome_build = "hg19",
       file = {
         p <- .cpgd_lookup_path(tissue)
         if (is.na(p)) "cpgdirectionData (ExperimentHub)" else p
       },
       package_version = as.character(utils::packageVersion("cpgdirection")))
}


#' @export
print.cpgd_result <- function(x, ...) {
  n_in <- attr(x, "n_submitted")
  cat(sprintf("\ncpgdirection: %d CpGs submitted, %d returned\n",
              if (is.null(n_in)) nrow(x) else n_in, nrow(x)))
  n_dup <- attr(x, "n_duplicate")
  if (!is.null(n_dup) && n_dup > 0)
    cat(sprintf("  (%d duplicate probe identifiers collapsed)\n", n_dup))
  st <- table(x$status)
  for (nm in names(st)) cat(sprintf("  %-22s %d\n", nm, st[[nm]]))
  cat("\n+1 = higher methylation -> HIGHER expression",
      "\n-1 = higher methylation -> LOWER expression\n")
  cat("\nAll accuracies are CONDITIONAL on the CpG being an eQTM. This gives the\n",
      "direction of an effect, not evidence that one exists. ABSTAIN and\n",
      "NO_EXPRESSION_ANCHOR are correct answers. See cpgd_accuracy_table().\n", sep = "")
  cat("\n")
  print(utils::head(as.data.frame(x[, intersect(c("cpg_id", "target_gene",
        "direction", "confidence", "evidence_tier", "call"), names(x)), with = FALSE]), 10))
  if (nrow(x) > 10) cat(sprintf("... %d more rows\n", nrow(x) - 10))
  invisible(x)
}


#' @export
summary.cpgd_result <- function(object, ...) {
  called <- object[object$call %in% c("1", "-1"), , drop = FALSE]
  list(n = nrow(object),
       by_status = table(object$status),
       by_tier = table(object$evidence_tier, useNA = "ifany"),
       n_called = nrow(called),
       coverage = if (nrow(object)) nrow(called) / nrow(object) else NA_real_,
       direction_of_calls = table(called$call))
}
