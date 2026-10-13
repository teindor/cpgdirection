#' Verify the installation against a shipped reference
#'
#' Scores 40 CpGs whose expected output was computed when the package was built
#' and compares the results field by field. Run this once after installing; if
#' it fails, do not use the results.
#'
#' The check is on \code{probability_plus1}, the model's actual output, together
#' with the assigned gene and evidence tier. Direction is reported but does not
#' gate the result, because it is the sign of the probability and would only
#' restate it.
#'
#' @param regenerate Optional path. Writes a fresh reference file from the
#'   installed package instead of testing against the shipped one. Use after a
#'   deliberate change of behaviour, never to make a failure go away.
#' @return Invisibly, \code{TRUE} on success. The report is emitted with
#'   \code{message()}, so \code{suppressMessages()} silences it.
#' @examplesIf cpgd_has_data("lookup_blood_hg19")
#' cpgd_selftest()
#' @export
cpgd_selftest <- function(regenerate = NULL) {

  ref_path <- system.file("extdata", "selftest_reference.csv", package = "cpgdirection")
  if (!nzchar(ref_path) || !file.exists(ref_path)) {
    message("reference file not found; cannot verify"); return(invisible(FALSE))
  }
  # An empty CSV field in a character column comes back as "" from fread, not
  # NA, while the package returns NA_character_. Without this the two are
  # compared as different values and every abstaining row is scored a mismatch.
  ref <- data.table::fread(ref_path, showProgress = FALSE,
                           na.strings = c("NA", ""))

  message("cpgdirection self-test")
  message(sprintf("  lookup table: %s", basename(.cpgd_lookup_path("blood"))))
  info <- tryCatch(cpgd_lookup_info("blood"), error = function(e) NULL)
  if (!is.null(info)) {
    message(sprintf("  %s pairs, %s CpGs, build %s",
                format(info$pairs, big.mark = ","),
                format(info$cpgs, big.mark = ","), info$genome_build))
  }

  got <- .cpgd_one_tissue(ref$input, tissue = "blood", verbose = FALSE)

  if (!is.null(regenerate)) {
    cols <- c("input", "cpg_id", "target_gene", "tss_dist", "status",
              "direction", "probability_plus1", "confidence", "evidence_tier")
    data.table::fwrite(got[, intersect(cols, names(got)), with = FALSE], regenerate)
    message(sprintf("  wrote a new reference to %s", regenerate))
    return(invisible(TRUE))
  }

  cols <- c("cpg_id", "target_gene", "direction", "probability_plus1",
            "confidence", "evidence_tier")
  m <- merge(ref[, intersect(cols, names(ref)), with = FALSE],
             got[, intersect(cols, names(got)), with = FALSE],
             by = "cpg_id", suffixes = c("_ref", "_got"))

  # Missing on one side only is a MISMATCH, not an unknown. The previous version
  # let NA propagate through `==`, so a single such row made the whole score NA
  # and the test reported failure without a number to act on.
  blank_to_na <- function(v) {
    if (is.character(v)) v[!is.na(v) & !nzchar(trimws(v))] <- NA_character_
    v
  }
  agree <- function(a, b) {
    a <- blank_to_na(a); b <- blank_to_na(b)
    both_na <- is.na(a) & is.na(b)
    eq <- !is.na(a) & !is.na(b) & (a == b)
    mean(both_na | eq)
  }
  near <- function(a, b, tol = 1e-4) {
    a <- suppressWarnings(as.numeric(a)); b <- suppressWarnings(as.numeric(b))
    both_na <- is.na(a) & is.na(b)
    eq <- !is.na(a) & !is.na(b) & (abs(a - b) <= tol)
    mean(both_na | eq)
  }

  ok_g <- agree(m$target_gene_ref, m$target_gene_got)
  ok_p <- near(m$probability_plus1_ref, m$probability_plus1_got)
  ok_c <- near(m$confidence_ref, m$confidence_got)
  ok_t <- agree(m$evidence_tier_ref, m$evidence_tier_got)
  ok_d <- agree(m$direction_ref, m$direction_got)

  message(sprintf("  compared %d CpGs", nrow(m)))
  message(sprintf("    gene assignment     %6.1f%%", 100 * ok_g))
  message(sprintf("    probability         %6.1f%%   <- the model's output", 100 * ok_p))
  message(sprintf("    confidence          %6.1f%%", 100 * ok_c))
  message(sprintf("    evidence tier       %6.1f%%", 100 * ok_t))
  message(sprintf("    direction           %6.1f%%   (informational)", 100 * ok_d))

  pass <- isTRUE(ok_g > 0.99) && isTRUE(ok_p > 0.99) && isTRUE(ok_t > 0.99)
  message(if (pass) "  SELF-TEST PASSED" else
          "  SELF-TEST FAILED - do not use these results.")

  if (pass && ok_d < 0.99) {
    message(sprintf(
      "  Note: direction agrees on %.0f%% while the probability agrees on 100%%.\n",
      100 * ok_d),
      "  The reference records direction as missing for rows it abstained on;\n",
      "  this version keeps the value and signals the abstention in `call`.\n",
      "  The model is unchanged. Regenerate with cpgd_selftest(regenerate=...)\n",
      "  if you want the reference to match current reporting.")
  }

  message("  tissues available:")
  for (t in CPGD_TISSUES) {
    i <- tryCatch(cpgd_lookup_info(t), error = function(e) NULL)
    if (!is.null(i)) message(sprintf("    %-18s %s pairs, %s CpGs", t,
        format(i$pairs, big.mark = ","), format(i$cpgs, big.mark = ",")))
  }
  invisible(pass)
}
