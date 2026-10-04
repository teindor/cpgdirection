#' Write a results table to HTML
#'
#' Renders a \code{cpg_expression_direction()} result as a formatted HTML table
#' and writes it to disk. Uses \pkg{gt} when it is installed, and falls back to
#' a plain styled table when it is not, so a long run never fails at the last
#' step for want of an optional package.
#'
#' @param x A result from \code{\link{cpg_expression_direction}}.
#' @param file File name. A relative name is written inside \code{dir}.
#' @param dir Directory to write into. Defaults to the working directory, so
#'   \code{setwd()} controls where results land.
#' @param title Table title.
#' @param subtitle Table subtitle. \code{NULL} generates one from the data.
#' @param compact Show the interpretive columns only. \code{FALSE} writes every
#'   column, which is wide but complete. Default \code{TRUE}.
#' @param open Open the file in a browser when done. Default \code{FALSE}.
#'
#' @return Invisibly, the path written.
#' @examplesIf requireNamespace("gt", quietly = TRUE) && cpgd_has_data("lookup_blood_hg19")
#' res <- cpg_expression_direction(c("cg02079741_TC21_POMC", "cg00000029"),
#'                                 universal = FALSE, verbose = FALSE)
#' cpgd_report(res, "hpa_panel.html", dir = tempdir())
#' @export
cpgd_report <- function(x,
                        file = "cpgdirection_results.html",
                        dir = getwd(),
                        title = "CpG to expression direction",
                        subtitle = NULL,
                        compact = TRUE,
                        open = FALSE) {

  if (!is.data.frame(x)) stop("`x` must be a cpgdirection result.", call. = FALSE)
  path <- if (grepl("^(/|[A-Za-z]:)", file)) file else file.path(dir, file)
  if (!grepl("\\.html?$", path, ignore.case = TRUE)) path <- paste0(path, ".html")
  d <- data.table::as.data.table(x)

  full <- "best_direction" %in% names(d)
  if (isTRUE(compact)) {
    want <- if (full)
      c("input", "cpg_id", "requested_gene", "best_direction", "best_evidence",
        "best_confidence", "best_tier", "best_expected_accuracy",
        "direction_uncertain", "best_direction_filled", "best_direction_flipped",
        # smr_gene and smr_gene_match are not optional detail. A causal
        # direction whose subject is unstated is ambiguous, and on real panels
        # roughly three quarters of SMR rows concern a different gene from the
        # rest of the row. Omitting them from the compact report - which is what
        # most users actually read - would present those directions as though
        # they were about the requested gene.
        "smr_direction", "smr_gene", "smr_gene_match", "smr_tier",
        "smr_agreement", "smr_heidi_status",
        "onco_direction", "onco_gene", "onco_gene_match", "onco_tier",
        "onco_agreement",
        "n_tissues_calling", "dist_universal",
        "dist_dir_blood", "dist_dir_nasal", "dist_dir_solid",
        "measured_genes", "note")
    else
      c("input", "cpg_id", "target_gene", "gene_source", "tss_dist", "status",
        "direction", "confidence", "evidence_tier", "expected_accuracy",
        "call", "note")
    d <- d[, intersect(want, names(d)), with = FALSE]
  }

  # Deliberately not a name beginning with a dot: gt passes column names through
  # as argument names, where a leading dot risks partial-matching gt's own
  # dot-prefixed parameters.
  dircol <- if (full) "best_direction" else "direction"
  if (dircol %in% names(d)) {
    data.table::set(d, j = "direction_call", value = data.table::fcase(
      d[[dircol]] ==  1, "+1  higher meth -> HIGHER expr",
      d[[dircol]] == -1, "-1  higher meth -> LOWER expr",
      default = "no call"))
    nm <- setdiff(names(d), "direction_call")
    at <- which(nm == dircol)
    data.table::setcolorder(d, c(nm[seq_len(at)], "direction_call",
                                 nm[-seq_len(at)]))
  }

  if (is.null(subtitle)) {
    n <- nrow(d)
    k <- if (dircol %in% names(d)) sum(!is.na(d[[dircol]])) else NA_integer_
    subtitle <- sprintf(
      "%d CpGs; %s directions assigned. Every direction is conditional on the CpG being an eQTM.",
      n, if (is.na(k)) "?" else as.character(k))
  }

  ok <- requireNamespace("gt", quietly = TRUE)
  if (ok) {
    .cpgd_report_gt(d, path, title, subtitle, dircol)
  } else {
    message("Package 'gt' is not installed; writing a plain HTML table instead.\n",
            "  install.packages(\"gt\") for the formatted version.")
    .cpgd_report_plain(d, path, title, subtitle)
  }
  message("wrote ", path)
  if (isTRUE(open)) utils::browseURL(path)
  invisible(path)
}


.cpgd_report_gt <- function(d, path, title, subtitle, dircol) {
  df <- as.data.frame(d)
  g <- gt::gt(df)
  g <- gt::tab_header(g, title = gt::md(paste0("**", title, "**")),
                      subtitle = subtitle)

  if ("direction_call" %in% names(df)) {
    pos <- which(df[[dircol]] == 1)
    neg <- which(df[[dircol]] == -1)
    non <- which(is.na(df[[dircol]]))
    if (length(pos))
      g <- gt::tab_style(g, gt::cell_fill(color = "#d8f0dc"),
                         gt::cells_body(columns = "direction_call", rows = pos))
    if (length(neg))
      g <- gt::tab_style(g, gt::cell_fill(color = "#fbdcdc"),
                         gt::cells_body(columns = "direction_call", rows = neg))
    if (length(non))
      g <- gt::tab_style(g, list(gt::cell_fill(color = "#f0f0f0"),
                                 gt::cell_text(color = "#888888")),
                         gt::cells_body(columns = "direction_call", rows = non))
  }

  # Evidence strength is the column readers skip, so make it impossible to skip.
  if ("best_evidence" %in% names(df)) {
    shade <- c(measured = "#bfe3c6", catalogue_consensus = "#d9ead8",
               catalogue_single = "#f2efd6", distance_only = "#fae6cf",
               distance_tissue_conflict = "#e6dcf0",
               distance_uninformative = "#f0f0f0", tissue_conflict = "#f6ddd0",
               no_evidence = "#f0f0f0")
    for (lv in names(shade)) {
      rw <- which(df$best_evidence == lv)
      if (length(rw))
        g <- gt::tab_style(g, gt::cell_fill(color = shade[[lv]]),
                           gt::cells_body(columns = "best_evidence", rows = rw))
    }
  }

  for (nm in intersect(c("best_confidence", "confidence", "p_universal"), names(df)))
    g <- gt::fmt_number(g, columns = nm, decimals = 3)
  for (nm in intersect(c("dist_universal", "tss_dist"), names(df)))
    g <- gt::fmt_number(g, columns = nm, decimals = 0, use_seps = TRUE)

  g <- gt::sub_missing(g, missing_text = "\u2014")
  g <- gt::opt_row_striping(g)
  g <- gt::tab_options(g, table.font.size = gt::px(12),
                       heading.title.font.size = gt::px(18),
                       data_row.padding = gt::px(3))
  if ("direction_uncertain" %in% names(df)) {
    rw <- which(df$direction_uncertain %in% TRUE)
    if (length(rw))
      g <- gt::tab_style(g, gt::cell_text(style = "italic", color = "#777777"),
                         gt::cells_body(columns = "direction_uncertain", rows = rw))
    g <- gt::tab_source_note(g, gt::md(
      "*Uncertain rows.* `best_direction_filled` and `best_direction_flipped` bracket the rows whose sign is nearly free to flip. They are a **sensitivity analysis, not a menu** \u2014 run the analysis under both and report whether the conclusion survives. Picking the column that fits better is selecting on the outcome."))
  }
  g <- gt::tab_source_note(g, gt::md(
    "*Direction, not existence.* Every accuracy figure is conditional on the CpG being an eQTM; the catalogues behind the models contain only significant associations. Read `best_evidence` alongside `best_direction`: `measured` is worth ~1.00, `distance_only` ~0.62."))
  g <- gt::tab_source_note(g, gt::md(
    paste0("Generated by cpgdirection ", utils::packageVersion("cpgdirection"),
           " on ", format(Sys.Date()), ". Lookup tables are hg19.")))
  gt::gtsave(g, filename = path)
}


.cpgd_report_plain <- function(d, path, title, subtitle) {
  esc <- function(v) {
    v <- ifelse(is.na(v), "\u2014", as.character(v))
    v <- gsub("&", "&amp;", v, fixed = TRUE)
    v <- gsub("<", "&lt;",  v, fixed = TRUE)
    gsub(">", "&gt;", v, fixed = TRUE)
  }
  df <- as.data.frame(d)
  hdr <- paste0("<th>", esc(names(df)), "</th>", collapse = "")
  body <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td>", esc(unlist(lapply(df[i, , drop = FALSE], as.character))),
                          "</td>", collapse = ""), "</tr>"), character(1))
  html <- c(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    paste0("<title>", esc(title), "</title>"),
    "<style>body{font-family:system-ui,sans-serif;margin:2rem;color:#222}",
    "table{border-collapse:collapse;font-size:12px}",
    "th,td{border:1px solid #ddd;padding:4px 7px;text-align:left}",
    "th{background:#f4f4f4}tr:nth-child(even){background:#fafafa}",
    "p.sub{color:#555}</style></head><body>",
    paste0("<h2>", esc(title), "</h2>"),
    paste0("<p class='sub'>", esc(subtitle), "</p>"),
    "<table>", paste0("<tr>", hdr, "</tr>"), body, "</table>",
    "<p class='sub'><em>Direction, not existence.</em> Every accuracy figure is",
    "conditional on the CpG being an eQTM.</p>",
    paste0("<p class='sub'>cpgdirection ", utils::packageVersion("cpgdirection"),
           ", ", format(Sys.Date()), "</p>"),
    "</body></html>")
  writeLines(html, path)
}


#' Is the balance of +1 and -1 calls what it should be?
#'
#' A sanity check for a result that looks suspiciously one-sided. Inverse
#' associations are the majority everywhere, so a run returning mostly \code{-1}
#' is expected; a run returning \emph{only} \code{-1} may still be correct, and
#' this says which.
#'
#' Called with no argument it reports the balance in the shipped tables. Called
#' on a result it compares that result's catalogue-derived calls against the
#' relevant base rate with an exact binomial test.
#'
#' @param x Optional result from \code{\link{cpg_expression_direction}}.
#' @param tissue Which tissue's base rate to test against. Default
#'   \code{"blood"}.
#' @return Invisibly, a list of the figures printed.
#'
#' @section Why distance-only calls are excluded:
#' The pooled distance curve falls below 0.5 out to roughly 19 kb and rises
#' above it only beyond roughly 300 kb, with an uninformative band between. A
#' panel whose CpGs all sit within a few tens of kilobases of their target
#' therefore \strong{cannot} receive a \code{+1} from the distance layer, no
#' matter how many CpGs it contains. That is a property of the panel's geometry,
#' not evidence about its CpGs, so those calls are not tested here.
#'
#' @examplesIf cpgd_has_data("lookup_blood_hg19")
#' cpgd_direction_balance()
#' @export
cpgd_direction_balance <- function(x = NULL, tissue = "blood") {

  tissue <- match.arg(tissue, CPGD_TISSUES)
  L <- .cpgd_lookup(tissue)
  dd <- suppressWarnings(as.numeric(L$direction))
  base_pos <- mean(dd > 0, na.rm = TRUE)

  M <- tryCatch(cpgd_measured_eqtms(), error = function(e) NULL)
  meas_pos <- if (is.null(M)) NA_real_ else
    mean(suppressWarnings(as.numeric(M$direction[M$tissue == tissue])) > 0, na.rm = TRUE)

  cat(sprintf("\ncpgdirection: balance of directions in %s\n", tissue))
  cat(sprintf("  predicted in the shipped table : %.1f%% +1\n", 100 * base_pos))
  if (!is.na(meas_pos))
    cat(sprintf("  measured (ground truth)        : %.1f%% +1\n", 100 * meas_pos))
  if (!is.na(meas_pos) && meas_pos - base_pos > 0.05) {
    cat(sprintf(
      "\n  The model predicts +1 less often than it truly occurs (%.1f%% vs %.1f%%).\n",
      100 * base_pos, 100 * meas_pos),
      "  It shrinks towards the majority class, so a +1 is conservative and a\n",
      "  -1 carries less weight than its count suggests.\n", sep = "")
  }

  res <- list(tissue = tissue, predicted_frac_pos = base_pos,
              measured_frac_pos = meas_pos)

  if (!is.null(x) && is.data.frame(x)) {
    d <- data.table::as.data.table(x)
    dircol <- if ("best_direction" %in% names(d)) "best_direction" else "direction"
    if (!dircol %in% names(d)) { cat("\n  (no direction column in `x`)\n"); return(invisible(res)) }
    if ("best_evidence" %in% names(d)) {
      keep <- d$best_evidence %in% c("measured", "catalogue_consensus", "catalogue_single")
      excl <- sum(d$best_evidence == "distance_only", na.rm = TRUE)
    } else {
      keep <- rep(TRUE, nrow(d)); excl <- 0L
    }
    v <- suppressWarnings(as.numeric(d[[dircol]][keep]))
    v <- v[!is.na(v)]
    cat(sprintf("\n  your run: %d catalogue-derived calls", length(v)))
    if (excl > 0) cat(sprintf(" (%d distance-only calls excluded)", excl))
    cat("\n")
    if (length(v)) {
      npos <- sum(v > 0)
      cat(sprintf("    +1 %d    -1 %d    (%.1f%% +1)\n", npos, length(v) - npos,
                  100 * npos / length(v)))
      bt <- stats::binom.test(npos, length(v), p = base_pos)
      cat(sprintf("    exact binomial vs the %.1f%% base rate: p = %.3g\n",
                  100 * base_pos, bt$p.value))
      cat(if (bt$p.value < 0.05)
            "    -> more one-sided than the base rate explains; look at the panel.\n"
          else
            "    -> consistent with the base rate; no anomaly to explain.\n")
      res$n_calls <- length(v); res$n_pos <- npos; res$p_value <- bt$p.value
    }
    if (excl > 0) {
      cat(sprintf(
        "\n  The %d distance-only calls are excluded on purpose, and for a stronger\n", excl),
        "  reason than sampling. The blood distance curve peaks at 0.449 and never\n",
        "  reaches 0.5 at any distance, while a distance call requires all three\n",
        "  tissue curves to agree. distance_only therefore CANNOT return +1 for any\n",
        "  CpG, ever. Testing those calls against a base rate would be circular.\n", sep = "")
    }
  }
  invisible(res)
}
