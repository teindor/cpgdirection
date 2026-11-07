# Internal helpers. Not exported.

.cpgd_env <- new.env(parent = emptyenv())

CPGD_TISSUES <- c("blood", "nasal_epithelium", "solid_tissue")

# The evidence levels, in ladder order, defined ONCE. They were previously
# written out inline in the print method and edited by hand as layers were
# added, which duplicated two of them and made factor() refuse the vector.
# A constant cannot drift from itself.
CPGD_EVIDENCE <- c(
  "measured",
  "smr_replicated",
  "smr_high",
  "catalogue_consensus",
  "smr_moderate",
  "catalogue_single",
  "onco_consensus",
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

# Share of HELIX-tested pairs with nominal p < 0.05, by rung (5% expected under
# the null). The complement is "not detected at n = 823", not "absent".
CPGD_SMR_DETECTION_CHILD <- c(replicated = 0.40, high = 0.40, moderate = 0.20,
                              weak = 0.28, reported = 0.21)

# child_blood_helix values in smr_directions (Hub version with the column):
# sig = HELIX sigPair; nominal = p < 0.05, not sig; null = tested, p >= 0.05;
# untested. Informational only: never a rung condition.
CPGD_CHILD_BLOOD_LEVELS <- c("sig", "nominal", "null", "untested")
CPGD_CHILD_DETECTED <- c("sig", "nominal")

# Who the quoted blood figures describe, by age_group (shown by the print
# methods). Nasal and solid figures do not depend on age.
CPGD_AGE_POPULATION <- c(
  adult = "adult blood (MESA monocytes and T cells as the yardstick)",
  child = "child blood (HELIX, 6-11 y, n = 823, where it scores the source on 100+ pairs; adult otherwise)",
  any   = "adult and child blood side by side")

# The accuracy key of an SMR row that supplies a direction: which validated
# stratum its sign belongs to. `src` is the ladder rung (smr_replicated,
# smr_high, smr_moderate, smr_weak); anything else gives NA. Within moderate,
# eQTLGen pairs whose cell-type FinnGen arm disagrees (discordant_celltype S2)
# are their own stratum: they validate well below the rest of the rung.
.cpgd_smr_key <- function(src, tier, outcome, replication) {
  eq <- outcome %in% "eQTLGen"
  data.table::fcase(
    src %in% "smr_replicated" & tier %in% "S1", "replicated_S1",
    src %in% "smr_replicated" & !is.na(tier), "replicated_S2",
    src %in% "smr_high", "high",
    src %in% "smr_moderate" & !eq, "moderate_finngen",
    src %in% "smr_moderate" & tier %in% "S3", "moderate_s3_replicated",
    src %in% "smr_moderate" & replication %in% "discordant_celltype",
      "moderate_discordant_celltype",
    src %in% "smr_moderate", "moderate",
    src %in% "smr_weak" & replication %in% "discordant", "weak_discordant",
    src %in% "smr_weak", "weak",
    default = NA_character_)
}

# require_detected acts only with age_group = "child"; anywhere else it is
# reported and ignored rather than silently honoured or silently dropped.
.cpgd_require_detected <- function(age_group, require_detected) {
  if (!is.logical(require_detected) || length(require_detected) != 1L ||
      is.na(require_detected)) {
    stop("`require_detected` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(require_detected) && !identical(age_group, "child")) {
    warning("require_detected = TRUE applies only with age_group = \"child\"; ",
            "ignored.", call. = FALSE)
    return(FALSE)
  }
  require_detected
}

# An SMR resource built before the HELIX column cannot answer the child
# questions. Say so once, return the column as NA, and refuse to apply
# require_detected (it would otherwise remove the whole SMR layer, silently).
# Returns the require_detected value to use.
.cpgd_check_child_column <- function(S, age_group, require_detected) {
  if (identical(age_group, "adult") || is.null(S) || !nrow(S) ||
      "child_blood_helix" %in% names(S)) return(require_detected)
  warning("this smr_directions resource predates the child_blood_helix ",
          "column: it is returned as NA",
          if (isTRUE(require_detected)) ", and require_detected cannot be applied (ignored)" else "",
          ". Update cpgdirectionData.", call. = FALSE)
  FALSE
}

# The child-blood column of an SMR table, NA where the table predates it.
.cpgd_smr_child <- function(S) {
  if ("child_blood_helix" %in% names(S)) as.character(S$child_blood_helix)
  else rep(NA_character_, nrow(S))
}

# Effective SMR rung for a table row. Tables built before the replication
# rebuild carry no smr_rung column; they map tier -> rung as the ladder always
# did, so a stale cached table degrades to the old behaviour, never to an error.
.cpgd_smr_rung <- function(S) {
  if ("smr_rung" %in% names(S)) return(as.character(S$smr_rung))
  tier <- as.character(S$smr_tier)
  data.table::fcase(tier == "S1", "high", tier == "S2", "moderate",
                    tier == "S3", "weak", default = NA_character_)
}
.cpgd_smr_col <- function(S, col, default) {
  if (col %in% names(S)) as.character(S[[col]]) else rep(default, nrow(S))
}

# Verbose summaries go through message() so that they can be silenced with
# suppressMessages() and never mix with a function's return value on stdout;
# the print methods (which cat()) are reused by capturing their output.
.cpgd_message_print <- function(x) {
  message(paste(utils::capture.output(print(x)), collapse = "\n"))
  invisible(x)
}

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
