#' cpgdirection: the direction of a CpG's effect on gene expression
#'
#' Assigns a direction to the association between methylation at a CpG site and
#' the expression of a nearby gene, for sites with no measured record, and
#' reports how much each answer is worth.
#'
#' Almost everything is called through the \pkg{data.table}, \pkg{stats},
#' \pkg{utils} and \pkg{tools} namespaces explicitly. What must be imported is
#' the non-standard evaluation symbols, plus the three data.table functions used
#' bare inside \code{[.data.table} expressions, where a \code{::} prefix would
#' be evaluated in the wrong frame. The gene annotation packages
#' (\pkg{TxDb.Hsapiens.UCSC.hg19.knownGene}, \pkg{org.Hs.eg.db},
#' \pkg{GenomicFeatures}, \pkg{GenomicRanges}, \pkg{AnnotationDbi}) are
#' Suggests: only the distance-only layer needs them, it checks for them with
#' \code{requireNamespace()} and says what to install when they are missing,
#' and their data objects are reached through \code{asNamespace()}.
#'
#' Code provenance: the design, data assembly, validation and all
#' methodological decisions are the author's. The R code was drafted and
#' edited with the assistance of Claude (Anthropic) and was run, tested and
#' checked by the author; no code was copied from other software projects.
#' Assisted-by: Claude (Anthropic).
#'
#' @keywords internal
#' @importFrom data.table := .SD .I .N fifelse fcase setorderv %chin%
"_PACKAGE"
