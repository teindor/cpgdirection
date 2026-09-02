# Internal helpers. Not exported.

.cpgd_env <- new.env(parent = emptyenv())

CPGD_TISSUES <- c("blood", "nasal_epithelium", "solid_tissue")

# The evidence levels, in ladder order, defined ONCE. They were previously
# written out inline in the print method and edited by hand as layers were
# added, which duplicated two of them and made factor() refuse the vector.
# A constant cannot drift from itself.
CPGD_EVIDENCE <- c(
  "measured",
  "smr_high",
  "catalogue_consensus",
  "smr_moderate",
  "catalogue_single",
  "smr_weak",
  "distance_only",
  "distance_targeted",
  "catalogue_targeted",
  "targeted_last_resort",
  "distance_tissue_conflict",
  "distance_uninformative",
  "tissue_conflict",
  "no_evidence")
stopifnot(!anyDuplicated(CPGD_EVIDENCE))

# Path of the BUNDLED lookup file, or NA when this (thin) installation holds
# none and the table comes from the cpgdirectionData ExperimentHub package
# instead. See R/data_backend.R for the resolution order.
.cpgd_lookup_path <- function(tissue = "blood") {
  p <- system.file("extdata", sprintf("lookup_%s_hg19.csv.gz", tissue),
                   package = "cpgdirection")
  if (!nzchar(p) || !file.exists(p)) return(NA_character_)
  p
}

# Read once per session and cache. The table is ~1.6M rows; fread takes a few
# seconds, and repeating that on every call would make interactive use painful.
.cpgd_lookup <- function(tissue = "blood") {
  key <- paste0("lookup_", tissue)
  if (!is.null(.cpgd_env[[key]])) return(.cpgd_env[[key]])
  L <- .cpgd_data(sprintf("lookup_%s_hg19", tissue))
  req <- c("cpg_id", "target_gene", "tss_dist", "status", "direction",
           "probability_plus1", "confidence", "evidence_tier")
  miss <- setdiff(req, names(L))
  if (length(miss)) {
    stop("The lookup table is missing expected columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  data.table::setkeyv(L, "cpg_id")
  .cpgd_env[[key]] <- L
  L
}

# Accepts a character vector, a data.frame, or a path to a text/CSV file.
# Strips identifiers to their canonical cg form and, for Project Alpha-style
# column names (cgXXXXXXX_TC21_GENE), parses the gene out of the name.
#
# `dedupe` controls the cpg_id-only collapse at the end. The legacy one-row-per
# CpG interfaces keep the historical behaviour (dedupe = TRUE). The pair-
# preserving workflow (cpg_gene_pairs) MUST pass dedupe = FALSE: collapsing by
# cpg_id alone folds "cg123_GENE_A" and "cg123_GENE_B" into one row and loses a
# biological co-effect. The parser parses; it does not decide which biological
# records are redundant.
.cpgd_parse_input <- function(cpgs, genes = NULL, dedupe = TRUE) {

  given <- NA_character_
  origin <- "none"

  if (is.data.frame(cpgs)) {
    nm <- names(cpgs)
    cc <- nm[grepl("^(cpg|probe|cg|id)", nm, ignore.case = TRUE)][1]
    if (is.na(cc)) cc <- nm[1]
    gc_ <- nm[grepl("^(gene|symbol|target)", nm, ignore.case = TRUE)][1]
    raw <- as.character(cpgs[[cc]])
    if (!is.na(gc_)) { given <- toupper(trimws(as.character(cpgs[[gc_]]))); origin <- "user" }
  } else if (is.character(cpgs) && length(cpgs) == 1L && file.exists(cpgs)) {
    if (grepl("\\.(csv|tsv)(\\.gz)?$", cpgs, ignore.case = TRUE)) {
      return(.cpgd_parse_input(.cpgd_fread(cpgs), genes,
                               dedupe = dedupe))
    }
    raw <- trimws(readLines(cpgs, warn = FALSE))
    raw <- raw[nzchar(raw) & !startsWith(raw, "#")]
  } else {
    raw <- trimws(as.character(cpgs))
  }

  if (!is.null(genes)) {
    if (length(genes) != length(raw)) {
      stop("`genes` must be the same length as `cpgs` (", length(raw), ").", call. = FALSE)
    }
    given <- toupper(trimws(as.character(genes))); origin <- "user"
  }

  # Extract the cg identifier from ANYWHERE in the string. Callers prefix names
  # for all sorts of reasons ("Zcg00335286_TC21_MC2R" from a z-scored column),
  # and taking the first underscore field would silently drop every such row.
  # Fully vectorised: regexpr gives position and length, substr extracts in place.
  # regmatches() would return only the matching elements, a shorter vector, and
  # realigning it to the input is easy to get wrong.
  pos <- regexpr("cg[0-9]{6,}", raw, ignore.case = TRUE)
  len <- attr(pos, "match.length")
  cpg <- ifelse(pos > 0L, tolower(substr(raw, pos, pos + len - 1L)), NA_character_)

  parts <- strsplit(raw, "_", fixed = TRUE)
  # Gene candidates: everything AFTER the token that carries the cg identifier,
  # as a whole and as its last token. Panel names carry things like
  # "Zcg00335286_TC21_MC2R" and "Zcg01961214_TC21_LINC02210_CRHR1", where the
  # informative symbol sits at the end.
  #
  # Locating the cg token is not the same as dropping positions 1 and 2. An
  # earlier version did the latter, which silently returned NA for the plainest
  # form a user would ever type -- "cg00039463_TRAP1" has only two tokens, so
  # removing the first two left nothing and the gene was lost. Probe-type codes
  # (TC21, BC11) are dropped by name rather than by position, and anything that
  # survives is validated against the TSS gene list by the caller, so a spurious
  # candidate costs nothing while a missing one costs a whole evidence layer.
  .drop_code <- function(p) p[!grepl("^([TB]C[0-9]+|[0-9]+)$", p)]
  tail_toks <- lapply(parts, function(p) {
    i <- grep("cg[0-9]{6,}", p, ignore.case = TRUE)[1]
    if (is.na(i)) i <- 1L
    tl <- if (i < length(p)) p[(i + 1L):length(p)] else character(0)
    .drop_code(tl)
  })
  tail_full <- vapply(tail_toks, function(tl)
    if (length(tl)) toupper(paste(tl, collapse = "_")) else NA_character_,
    character(1))
  tail_last <- vapply(tail_toks, function(tl)
    if (length(tl)) toupper(tl[length(tl)]) else NA_character_, character(1))

  if (all(is.na(given))) {
    gene  <- tail_full
    gene2 <- tail_last
    origin <- "parsed"
  } else {
    gene  <- given
    gene2 <- given
  }

  d <- data.table::data.table(input = raw, cpg_id = cpg,
                              given_gene = gene, given_gene2 = gene2,
                              gene_origin = origin)
  # Position of the row in the caller's input, assigned BEFORE any row is
  # removed or collapsed, so every surviving row can be traced back to the
  # physical input column it came from.
  d[, "input_id" := .I]
  d <- d[!is.na(d$cpg_id), ]
  for (col in c("given_gene", "given_gene2")) {
    bad <- which(d[[col]] %in% c("", "NA", "NAN", "NONE"))
    if (length(bad)) data.table::set(d, i = bad, j = col, value = NA_character_)
  }
  # A panel can carry the same probe several times under different suffixes
  # (cg09527270_TC11 / _TC12 / _TC13). The legacy interfaces deduplicate, but
  # record the raw count so callers can say how many inputs collapsed rather
  # than losing them silently. The pair workflow passes dedupe = FALSE and
  # receives every parsed row.
  n_raw <- nrow(d)
  if (isTRUE(dedupe)) d <- unique(d, by = "cpg_id")
  data.table::setattr(d, "n_raw", n_raw)
  d
}

#' Show what the identifier parser extracted
#'
#' Diagnostic. Given the identifiers you intend to pass to
#' \code{\link{cpg_expression_direction}}, shows the CpG and gene the parser
#' pulls out of each, so a naming convention that does not parse can be seen
#' rather than guessed at.
#'
#' @param cpgs Character vector of identifiers, as you would pass them.
#' @return A \code{data.table} with \code{input}, \code{cpg_id},
#'   \code{given_gene} and \code{given_gene2}.
#' @examples
#' cpgd_parse_check(c("Zcg00335286_TC21_MC2R", "cg00000029", "not_a_probe"))
#' @export
cpgd_parse_check <- function(cpgs) {
  d <- .cpgd_parse_input(cpgs)
  if (!nrow(d)) {
    message("Nothing parsed. Identifiers need 'cg' followed by 6+ digits.")
    return(invisible(d))
  }
  d[, c("input", "cpg_id", "given_gene", "given_gene2"), with = FALSE]
}
