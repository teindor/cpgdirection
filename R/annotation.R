# One place to ask where a CpG is and what the package knows about it (1.1.0).

#' Where a CpG is, in hg19, hg38 and T2T-CHM13, and what the package knows about it
#'
#' One row per probe with its coordinates in the genome build you work in, the
#' arrays that carry it, its gene annotation and its probe quality flags.
#' It covers every probe ID of the four Illumina methylation arrays (HM27,
#' 450K, EPIC v1, EPIC v2; 1,084,297 IDs), not only the EPIC v2 probes the
#' direction layers are keyed on.
#'
#' @section Where the coordinates come from:
#' \strong{hg19.} The package's own table (\code{\link{cpgd_cpg_positions}},
#' the one the direction engine uses) wherever it holds the probe, so the two
#' always agree. For probes it does not hold (450K and EPIC v1 probes that
#' EPIC v2 dropped) the Illumina 450K manifest, else the Illumina EPIC v1
#' coordinates (build 37 in both).
#'
#' \strong{hg38.} The Zhou-lab (SeSAMe) hg38 manifests, taking the manifest
#' of the first array generation that carries the probe (HM450, then EPIC,
#' then EPIC v2, then HM27): the 1-based position of the first base of the
#' CpG on the forward strand (\code{CpG_beg} + 1). A 450K probe therefore
#' sits exactly where the SeSAMe HM450 manifest puts it. Where no manifest
#' places a probe, the position of the package's EPIC v2 annotation layer is
#' used. Probes on alternate or unplaced contigs keep that contig name.
#'
#' \strong{T2T-CHM13.} From \code{\link{cpgd_probe_qc}}; EPIC v2 probes only.
#'
#' \code{source} (or \code{source_hg19}, \code{source_hg38}) says which table
#' placed each probe. A coordinate the package does not hold is \code{NA},
#' never lifted over or guessed.
#'
#' @section Coverage and ambiguous probes:
#' hg38: every probe ID but 19 (8 of the 450K array, 14 of EPIC v1, 3 of EPIC
#' v2), which no manifest maps. hg19: 98.9\% (EPIC v2 99.95\%, 450K 99.96\%,
#' EPIC v1 98.8\%); about 12,200 probe IDs have no hg19 coordinate here, most
#' of them EPIC v1 probes that EPIC v2 dropped and HM27 probes not on the
#' 450K array.
#'
#' \code{hg38_ambiguous} is \code{TRUE} for 1,094 probes that two sources
#' place differently: probes whose sequence maps to more than one place, so
#' that the manifests of different arrays disagree, and EPIC v2 replicate
#' probes that map to different places. They are placed by the priority
#' above; treat their coordinate, in any build, with caution.
#'
#' @section What this does not change:
#' The direction engine works in hg19: its models and distance curves were
#' fitted on hg19 CpG-to-TSS distances, and \code{tss_dist} everywhere in
#' the package is an hg19 distance. The hg38 and CHM13 coordinates returned
#' here are for locating and joining your data; no direction is recomputed
#' in another build.
#'
#' @param cpgs Probe identifiers, as your data name them: bare 450K or EPIC
#'   v1 IDs, or EPIC v2 IDs with their replicate suffix
#'   (\code{cg00000029_TC21}); a gene suffix (\code{cg00000029_RBL2}) is
#'   ignored too, and case does not matter. \code{NULL} returns every probe.
#' @param build One or more of \code{"hg19"}, \code{"hg38"}, \code{"chm13"}.
#'   With one build the coordinate columns are named \code{chr} and
#'   \code{pos}; with several, \code{chr_<build>} and \code{pos_<build>}.
#'   All three by default.
#' @param genes Add \code{genes}: the genes the EPIC v2 annotation layer
#'   assigns to the CpG (any of its four sources), separated by \code{";"}.
#' @param qc Add the probe quality flags of \code{\link{cpgd_probe_qc}}
#'   (\code{qc_exclude}, \code{masked_general}, \code{cross_hybridizing}).
#' @return A \code{data.table}, one row per element of \code{cpgs} in the
#'   order given (or one per probe): \code{input} (the ID exactly as given,
#'   when \code{cpgs} is given), \code{cpg_id} (without suffix),
#'   \code{known} (the ID is a probe of one of the four arrays), the
#'   coordinate columns and their \code{source}, for hg38 also
#'   \code{strand_hg38} (where the EPIC v2 annotation gives it) and
#'   \code{hg38_ambiguous}, \code{array} (the arrays that carry the probe,
#'   for example \code{"450K;EPICv1;EPICv2"}) with \code{on_hm27},
#'   \code{on_450k}, \code{on_epic_v1}, \code{on_epic_v2}, and the optional
#'   gene and quality columns. The region flags of \code{cpgd_region_flags()}
#'   are added when the package holds them. A probe that cannot be placed has
#'   \code{NA} coordinates.
#' @examples
#' cpgd_cpg_annotation(c("cg00000029", "cg00000029_TC21", "cg_unknown"), build = "hg38")
#' @export
cpgd_cpg_annotation <- function(cpgs = NULL, build = c("hg19", "hg38", "chm13"),
                                genes = TRUE, qc = TRUE) {
  build <- match.arg(build, c("hg19", "hg38", "chm13"), several.ok = TRUE)
  if (is.null(.cpgd_env$cpg_coordinates)) {
    X <- .cpgd_data("cpg_coordinates")
    for (cc in c("chr_hg19", "source_hg19", "chr_hg38", "source_hg38")) {
      v <- as.character(X[[cc]]); v[which(v == "")] <- NA_character_
      data.table::set(X, j = cc, value = v)
    }
    data.table::set(X, j = "hg38_ambiguous", value = as.logical(X[["hg38_ambiguous"]]) %in% TRUE)
    # probe IDs are matched without regard to case
    data.table::set(X, j = ".key", value = tolower(X[["cpg_id"]]))
    data.table::setkeyv(X, ".key")
    .cpgd_env$cpg_coordinates <- X
  }
  X <- .cpgd_env$cpg_coordinates
  if (is.null(cpgs)) {
    out <- data.table::copy(X)
    out[, "known" := TRUE]
  } else {
    input <- as.character(cpgs)
    id <- sub("_.*$", "", trimws(input))
    out <- X[data.table::data.table(.key = tolower(id)), on = ".key"]
    out[, "known" := !is.na(get("array"))]
    i <- which(is.na(out$cpg_id))
    data.table::set(out, i = i, j = "cpg_id", value = trimws(input)[i])   # unknown: as given
    data.table::set(out, i = i, j = "hg38_ambiguous", value = FALSE)
    out[, "input" := input]
  }
  keys <- data.table::data.table(.key = out$.key)
  lower_key <- function(D) {
    D <- data.table::as.data.table(D)
    data.table::set(D, j = ".key", value = tolower(D[["cpg_id"]]))
    D
  }
  # hg19: the package's own table (the one the direction engine uses) comes
  # first; the bundled table only fills probes it does not hold
  if ("hg19" %in% build) {
    P <- tryCatch(cpgd_cpg_positions(), error = function(e) NULL)
    if (!is.null(P)) {
      h <- lower_key(P)[keys, on = ".key"]
      i <- which(!is.na(h$pos))
      data.table::set(out, i = i, j = "chr_hg19", value = as.character(h$chr[i]))
      data.table::set(out, i = i, j = "pos_hg19", value = as.integer(h$pos[i]))
      data.table::set(out, i = i, j = "source_hg19", value = "cpgdirection")
    }
  }
  # hg38: the bundled table holds the position wherever it is not simply the
  # value of the EPIC v2 annotation layer; that layer supplies the rest
  A <- NULL
  if ("hg38" %in% build || isTRUE(genes)) {
    A <- .cpgd_epicv2_summary()
    if (!is.null(A)) A <- A[keys, on = ".key"]
  }
  out[, "strand_hg38" := NA_character_]
  if ("hg38" %in% build && !is.null(A)) {
    i <- which(is.na(out$pos_hg38) & !is.na(A$position))
    data.table::set(out, i = i, j = "chr_hg38", value = as.character(A$chr[i]))
    data.table::set(out, i = i, j = "pos_hg38", value = as.integer(A$position[i]))
    data.table::set(out, i = i, j = "source_hg38", value = "cpgdirection")
    data.table::set(out, j = "strand_hg38", value = A$strand)
  }
  Q <- NULL
  if ("chm13" %in% build || isTRUE(qc)) {
    Q <- tryCatch(cpgd_probe_qc(), error = function(e) NULL)
    if (!is.null(Q)) Q <- lower_key(Q)[keys, on = ".key"]
  }
  if ("chm13" %in% build) {
    data.table::set(out, j = "chr_chm13",
                    value = if (is.null(Q)) NA_character_ else as.character(Q$chm13_chr))
    data.table::set(out, j = "pos_chm13",
                    value = if (is.null(Q)) NA_integer_ else as.integer(Q$chm13_pos))
  }
  ar <- out$array
  has <- function(a) !is.na(ar) & grepl(a, ar, fixed = TRUE)
  data.table::set(out, j = "on_hm27", value = has("HM27"))
  data.table::set(out, j = "on_450k", value = has("450K"))
  data.table::set(out, j = "on_epic_v1", value = has("EPICv1"))
  data.table::set(out, j = "on_epic_v2", value = has("EPICv2"))
  one <- length(build) == 1L
  if (one) {
    # one build asked for: plain `chr`, `pos` (and `source`), ready to join on
    data.table::setnames(out, paste0(c("chr_", "pos_"), build), c("chr", "pos"))
    if (build %in% c("hg19", "hg38")) data.table::setnames(out, paste0("source_", build), "source")
  }
  cols <- c(if (!is.null(cpgs)) "input", "cpg_id", "known",
            if (one) c("chr", "pos", if (build %in% c("hg19", "hg38")) "source")
            else c(unlist(lapply(build, function(b) paste0(c("chr_", "pos_"), b))),
                   intersect(paste0("source_", build), c("source_hg19", "source_hg38"))),
            if ("hg38" %in% build) c("strand_hg38", "hg38_ambiguous"),
            "array", "on_hm27", "on_450k", "on_epic_v1", "on_epic_v2")
  if (isTRUE(genes)) {
    data.table::set(out, j = "genes",
                    value = if (is.null(A)) NA_character_ else A$genes)
    cols <- c(cols, "genes")
  }
  if (isTRUE(qc)) {
    for (cc in c("qc_exclude", "masked_general", "cross_hybridizing")) {
      data.table::set(out, j = cc, value = if (is.null(Q) || !cc %in% names(Q)) NA else as.logical(Q[[cc]]))
    }
    cols <- c(cols, "qc_exclude", "masked_general", "cross_hybridizing")
  }
  Fl <- tryCatch(.cpgd_data("cpg_region_flags", required = FALSE), error = function(e) NULL)
  if (!is.null(Fl)) {
    h <- lower_key(Fl)[keys, on = ".key"]
    for (cc in setdiff(names(Fl), c("cpg_id", ".key"))) {
      data.table::set(out, j = cc, value = as.logical(h[[cc]]) %in% TRUE)
      cols <- c(cols, cc)
    }
  }
  out <- out[, cols, with = FALSE]
  data.table::setattr(out, "build", build)
  out[]
}

# The EPIC v2 probe-gene annotation, one row per CpG: its hg38 position and
# strand and the genes any annotation source assigns to it. Built once per
# session from the (CpG x gene) layer.
.cpgd_epicv2_summary <- function() {
  if (!is.null(.cpgd_env$epicv2_summary)) return(.cpgd_env$epicv2_summary)
  A <- tryCatch(.cpgd_data("epicv2_probe_gene_annotation", required = FALSE),
                error = function(e) NULL)
  if (is.null(A) || !nrow(A)) return(NULL)
  A <- data.table::as.data.table(A)
  ok <- !is.na(A$position) & grepl("^chr", A$chr)
  U <- unique(A[ok, c("cpg_id", "chr", "position", "strand"), with = FALSE], by = "cpg_id")
  st <- as.character(U$strand); st[which(st == "")] <- NA_character_
  data.table::set(U, j = "strand", value = st)
  hasg <- !is.na(A$target_gene) & nzchar(A$target_gene)
  G <- A[hasg, list(genes = paste(unique(get("target_gene")), collapse = ";")), by = "cpg_id"]
  S <- merge(U, G, by = "cpg_id", all = TRUE)
  data.table::set(S, j = ".key", value = tolower(S[["cpg_id"]]))
  data.table::setkeyv(S, ".key")
  .cpgd_env$epicv2_summary <- S
  S
}
