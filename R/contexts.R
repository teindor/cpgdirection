# Reported tissue contexts and per-CpG region flags (1.1.0; PI rulings of 11-12
# Nov 2026). Both are information: neither votes in the measured table and
# neither can supply or change best_direction.

# What each context is, and its sign agreement with the measured table (the
# admission test of a voting source, bar 0.85; tools/reports/build_reported_contexts.md).
CPGD_REPORTED_CONTEXTS <- data.frame(
  context = c("fibroblasts", "liver", "pancreatic_islets", "retina", "saliva"),
  evidence = c("measured_eqtm", "measured_eqtm", "measured_eqtm", "measured_eqtm", "smr"),
  study = c("Wagner et al. 2014 (Genome Biol)",
            "Bonder et al. 2014 (BMC Genomics); Horvath liver of Wang et al. 2021",
            "Olsson et al. 2014 (PLoS Genet)",
            "Advani et al. 2024",
            "GoDMC cis-mQTL instruments against saliva cis-eQTLs (SMR)"),
  agreement_with_measured = c(0.721, 0.847, 0.805, 0.842, 0.762),
  agreement_lo = c(0.674, 0.826, 0.788, 0.792, 0.705),
  agreement_hi = c(0.764, 0.865, 0.822, 0.882, 0.811),
  n_shared = c(380L, 1306L, 2149L, 253L, 248L),
  stringsAsFactors = FALSE)

#' Reported tissue contexts: signs measured or estimated outside the voting tissues
#'
#' CpG-gene signs from tissues the package has no voting class for: retina,
#' liver, pancreatic islets, fibroblasts (published eQTM tables) and saliva
#' (SMR of the GoDMC cis-mQTL instruments against saliva cis-eQTLs).
#' They are \strong{reported only}: they are not part of
#' \code{\link{cpgd_measured_eqtms}}, cast no vote, and never supply or change
#' \code{best_direction} in \code{\link{cpg_gene_pairs}} or
#' \code{\link{cpg_expression_direction}}.
#'
#' @section Why these are reported and not voting:
#' A measured source votes when its signs agree with the pairs the package
#' already measures at 0.85 or better. Within 250 kb these sources agree at
#' 0.72 (fibroblasts, n = 380), 0.81 (pancreatic islets, n = 2,149), 0.84
#' (retina, n = 253) and 0.85 with an interval that straddles the bar (liver,
#' 0.847, 95\% interval 0.826-0.865, n = 1,306); the saliva SMR context agrees
#' at 0.76 (n = 248). Lower agreement can reflect tissue-specific regulation
#' as much as error, so the signs are kept, labelled with their tissue, and
#' left out of the ladder until more datasets per tissue allow a validated
#' figure (\code{attr(result, "contexts")} holds this table). Liver comes from
#' two independent cohorts; the 337 pairs both hold have the same sign in
#' every case (\code{n_cohorts} = 2).
#'
#' @param context Optional: one or more of \code{"fibroblasts"},
#'   \code{"liver"}, \code{"pancreatic_islets"}, \code{"retina"},
#'   \code{"saliva"}. All by default.
#' @return A \code{data.table} with \code{cpg_id}, \code{target_gene},
#'   \code{context}, \code{direction} (+1/-1: the sign of the observed
#'   association, or of the SMR estimate for saliva), \code{tss_dist},
#'   \code{n_cohorts}, \code{cohort_agreement}, \code{cohorts}, \code{source},
#'   \code{evidence} (\code{measured_eqtm} or \code{smr}) and \code{flags}
#'   (saliva: \code{mhc}, \code{inv17} for the 17q21.31 inversion). The
#'   attribute \code{"contexts"} describes each context and its agreement
#'   with the measured table.
#' @seealso \code{\link{cpgd_annotate}} to attach these signs to a result.
#' @examples
#' x <- cpgd_reported_contexts("liver")
#' table(x$n_cohorts)
#' attr(x, "contexts")
#' @export
cpgd_reported_contexts <- function(context = NULL) {
  if (is.null(.cpgd_env$reported_contexts)) {
    X <- .cpgd_data("reported_contexts")
    X[, "target_gene" := toupper(trimws(as.character(get("target_gene"))))]
    X[, "flags" := as.character(get("flags"))]
    X[which(get("flags") == ""), "flags" := NA_character_]
    .cpgd_env$reported_contexts <- X
  }
  X <- data.table::copy(.cpgd_env$reported_contexts)
  if (!is.null(context)) {
    want <- match.arg(context, CPGD_REPORTED_CONTEXTS$context, several.ok = TRUE)
    # computed outside `[`, where `context` would resolve to the table's own column
    keep <- X[["context"]] %in% want
    X <- X[keep]
  }
  data.table::setattr(X, "contexts", CPGD_REPORTED_CONTEXTS)
  X[]
}

#' Region flags: CpGs in CoRSIVs, imprinted regions and allele-specific methylation
#'
#' Per-CpG flags for three kinds of region where methylation is known to be
#' set by something other than the tissue's regulatory state:
#' \describe{
#'   \item{in_corsiv}{the CpG lies in a correlated region of systemic
#'     interindividual variation (CoRSIVs 2.0): methylation there is
#'     established before gastrulation and is shared across tissues}
#'   \item{in_imprinted_or_parental_asm}{imprinted, or allele-specific
#'     methylation by parent of origin}
#'   \item{in_asm}{allele-specific methylation (sequence dependent or parental;
#'     Rosenski et al. 2025)}
#' }
#' The flags are information. In four checks (measured cohorts against each
#' other, SMR against measured blood, BLUEPRINT against measured blood, SMR
#' against BLUEPRINT) sign agreement inside these regions differed from
#' outside by 0.06 at most, in both directions and with no consistent loss
#' (largest drop 0.03), so no flag changes a direction or an accuracy
#' (\code{tools/build_region_flags.py}).
#'
#' @return A \code{data.table} with \code{cpg_id}, \code{in_corsiv},
#'   \code{in_imprinted_or_parental_asm} and \code{in_asm} (logical), one row
#'   per flagged CpG (66,636); a CpG that is absent carries no flag.
#' @seealso \code{\link{cpgd_annotate}} to attach the flags to a result.
#' @examples
#' f <- cpgd_region_flags()
#' colSums(f[, -1])
#' @export
cpgd_region_flags <- function() {
  if (is.null(.cpgd_env$region_flags)) {
    X <- .cpgd_data("cpg_region_flags")
    for (cc in setdiff(names(X), "cpg_id")) data.table::set(X, j = cc, value = as.logical(X[[cc]]))
    data.table::setkeyv(X, "cpg_id")
    .cpgd_env$region_flags <- X
  }
  data.table::copy(.cpgd_env$region_flags)
}

#' Attach reported tissue contexts and region flags to a result
#'
#' Adds information columns to any table with a \code{cpg_id} column, such as
#' the output of \code{\link{cpg_gene_pairs}} or
#' \code{\link{cpg_expression_direction}}. Nothing else in the table is
#' touched: \code{best_direction} and its accuracy stay as they are.
#'
#' @param x A data frame with a \code{cpg_id} column. Context signs are
#'   matched on the pair when \code{x} also has \code{target_gene} (or
#'   \code{requested_gene}); otherwise every context sign of the CpG is listed
#'   with its gene.
#' @param contexts Add \code{context_directions}: the signs of
#'   \code{\link{cpgd_reported_contexts}} for the row, as
#'   \code{"liver:-1|retina:+1"} (with the gene, \code{"liver:GENE:-1"}, when
#'   \code{x} names no gene), and \code{context_n}, their number.
#' @param region_flags Add \code{in_corsiv},
#'   \code{in_imprinted_or_parental_asm} and \code{in_asm}
#'   (\code{\link{cpgd_region_flags}}; \code{FALSE} for a CpG without a flag).
#' @return \code{x} as a \code{data.table}, with the added columns, in the
#'   same row order.
#' @examples
#' x <- data.frame(cpg_id = c("cg00000165", "cg00000029"), target_gene = c("X", "RBL2"))
#' cpgd_annotate(x)
#' @export
cpgd_annotate <- function(x, contexts = TRUE, region_flags = TRUE) {
  if (!is.data.frame(x) || !"cpg_id" %in% names(x)) {
    stop("`x` must be a data frame with a `cpg_id` column.", call. = FALSE)
  }
  out <- data.table::as.data.table(x)
  out <- data.table::copy(out)
  if (isTRUE(contexts)) {
    C <- cpgd_reported_contexts()
    gcol <- intersect(c("target_gene", "requested_gene"), names(out))
    sgn <- ifelse(C$direction > 0, "+1", "-1")
    if (length(gcol)) {
      C[, ".lab" := paste0(get("context"), ":", sgn)]
      A <- C[, list(context_directions = paste(get(".lab"), collapse = "|"), context_n = .N),
             by = c("cpg_id", "target_gene")]
      key <- data.table::data.table(cpg_id = out$cpg_id,
                                    target_gene = toupper(trimws(as.character(out[[gcol[1]]]))))
      hit <- A[key, on = c("cpg_id", "target_gene")]
    } else {
      C[, ".lab" := paste0(get("context"), ":", get("target_gene"), ":", sgn)]
      A <- C[, list(context_directions = paste(get(".lab"), collapse = "|"), context_n = .N),
             by = "cpg_id"]
      hit <- A[data.table::data.table(cpg_id = out$cpg_id), on = "cpg_id"]
    }
    data.table::set(out, j = "context_directions", value = hit$context_directions)
    data.table::set(out, j = "context_n", value = ifelse(is.na(hit$context_n), 0L, hit$context_n))
  }
  if (isTRUE(region_flags)) {
    Fl <- cpgd_region_flags()
    hit <- Fl[data.table::data.table(cpg_id = out$cpg_id), on = "cpg_id"]
    for (cc in c("in_corsiv", "in_imprinted_or_parental_asm", "in_asm")) {
      data.table::set(out, j = cc, value = hit[[cc]] %in% TRUE)
    }
  }
  out[]
}
