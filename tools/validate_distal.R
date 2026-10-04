#!/usr/bin/env Rscript
# Validate the distal_links layer (ENCODE-rE2G) against the package's measured
# eQTMs and SMR pairs. A link carries no sign, so the questions are about
# TARGET assignment, plus one check that the links do not smuggle in a sign:
#
#   Q1 coverage   How many reference CpGs / pairs have any blood-class link?
#   Q2 targets    For reference CpGs with links: precision (linked genes that are
#                 a reference gene for that CpG) and recall (reference pairs that
#                 are linked), by rE2G score and by number of biosamples, against
#                 the manifest annotation on the SAME CpGs. Precision is bounded
#                 above by reference incompleteness (an eQTM study tests a window
#                 of genes, not all genes), so compare it with the manifest's
#                 precision, not with 1.
#   Q3 distance   Do links reach pairs the manifest misses, and at what distance?
#   Q4 sign       Share of negative reference signs among linked pairs, by element
#                 class and score bin, against the reference base rate. The layer
#                 is documented as sign-free; this shows whether that is true in
#                 the data (a promoter-class link is expected to look like a
#                 near-TSS pair, which is distance information the ladder already
#                 has, not new direction evidence).
#
# Reference sets, evaluated separately: measured blood, nasal, solid; SMR S1+S2.
#
# Usage
#   Rscript tools/validate_distal.R --data-dir <hub_upload dir> \
#       --out tools/reports/validate_distal.md

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
data_dir <- get_arg("--data-dir")
out_md   <- get_arg("--out", "tools/reports/validate_distal.md")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)
# package code: an installed cpgdirection, else load_all() from the package root,
# else (no Bioconductor deps on this machine) source R/ directly -- the data
# accessors used here need only data.table.
# Prefer the package root (the code under validation), not whatever older
# cpgdirection happens to be installed in the library.
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else if (requireNamespace("cpgdirection", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cpgdirection))
} else stop("run from the cpgdirection package root or install cpgdirection")
if (!exists("cpgd_distal_links")) stop("cpgdirection < 2.99.6 loaded; run from the package root")

wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(c - h, c + h)
}
fmt_ci <- function(k, n) { if (n == 0) return("-"); ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2]) }
md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  h <- paste0("| ", paste(names(d), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  r <- apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(h, s, r), collapse = "\n")
}

# ---- data --------------------------------------------------------------------------
D <- as.data.table(cpgd_distal_links())
if (is.null(D) || !nrow(D)) stop("distal_links layer not available; set --data-dir")
D[, target_gene := toupper(trimws(target_gene))]
# per CpG x gene within the blood class (what a blood question consults)
L <- D[tissue_class == "blood", list(
  n_biosamples = uniqueN(biosample),
  score_max = max(score, na.rm = TRUE),
  promoter_class = any(grepl("promoter", element_class, ignore.case = TRUE))),
  by = c("cpg_id", "target_gene")]
L[, score_bin := cut(score_max, c(0, 0.5, 0.9, 1), labels = c("0.2-0.5", "0.5-0.9", ">=0.9"), include.lowest = TRUE)]
L[, bs_bin := cut(n_biosamples, c(0, 2, 9, Inf), labels = c("1-2", "3-9", ">=10"))]
message(sprintf("blood-class links: %s pairs, %s CpGs", format(nrow(L), big.mark = ","), format(uniqueN(L$cpg_id), big.mark = ",")))

A <- tryCatch(as.data.table(cpgd_manifest_genes()), error = function(e) NULL)
if (!is.null(A)) {
  A <- unique(A[, list(cpg_id, target_gene = toupper(trimws(target_gene)))])
  A <- A[nzchar(target_gene)]
}

M <- as.data.table(cpgd_measured_eqtms())
M <- M[, list(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
       by = c("cpg_id", "tissue", "direction", "tss_dist")]
M <- unique(M[nzchar(target_gene)])
S <- tryCatch(as.data.table(cpgd_smr_directions()), error = function(e) NULL)
refs <- list(
  blood = M[tissue == "blood",            list(cpg_id, target_gene, ref_dir = direction, dist = abs(tss_dist))],
  nasal = M[tissue == "nasal_epithelium", list(cpg_id, target_gene, ref_dir = direction, dist = abs(tss_dist))],
  solid = M[tissue == "solid_tissue",     list(cpg_id, target_gene, ref_dir = direction, dist = abs(tss_dist))])
if (!is.null(S) && "smr_tier" %in% names(S)) {
  dcol <- intersect(c("cpg_gene_dist", "tss_dist"), names(S))[1]
  refs$smr_S1S2 <- S[smr_tier %in% c("S1", "S2"),
                     list(cpg_id, target_gene = toupper(target_gene), ref_dir = direction,
                          dist = if (!is.na(dcol)) abs(get(dcol)) else NA_real_)]
}
refs <- lapply(refs, function(r) unique(r[!is.na(target_gene) & nzchar(target_gene)]))

# ---- evaluation ----------------------------------------------------------------------
pr <- function(R, X) {
  # precision/recall of a proposal table X (cpg_id, target_gene) against reference R,
  # restricted to CpGs present in BOTH (so the two proposal sources are compared on
  # the same CpGs when called with the same R)
  cp <- intersect(unique(R$cpg_id), unique(X$cpg_id))
  Rc <- R[cpg_id %in% cp]; Xc <- X[cpg_id %in% cp]
  hit_x <- Xc[R, on = c("cpg_id", "target_gene"), nomatch = NULL]
  list(n_cpg = length(cp),
       proposals = nrow(Xc), precision = fmt_ci(nrow(unique(hit_x)), nrow(Xc)),
       ref_pairs = nrow(Rc), recall = fmt_ci(nrow(unique(hit_x)), nrow(Rc)))
}

report <- c("# distal_links (ENCODE-rE2G, blood class) vs peripheral references", "",
            sprintf("Layer: %s blood-class CpG-gene links over %s CpGs (26 biosamples; brain class not used here).",
                    format(nrow(L), big.mark = ","), format(uniqueN(L$cpg_id), big.mark = ",")),
            "Precision = linked genes that are a reference gene for that CpG; recall = reference pairs that are linked;",
            "both on the CpGs the two sources share. Wilson 95% CI. The manifest annotation is the comparator.", "")
for (nm in names(refs)) {
  R <- refs[[nm]]
  message(sprintf("[%s] %d reference pairs", nm, nrow(R)))
  report <- c(report, sprintf("## Reference: %s (%s pairs, %s CpGs)", nm,
                              format(nrow(R), big.mark = ","), format(uniqueN(R$cpg_id), big.mark = ",")), "")
  # Q1
  cov_cpg  <- mean(unique(R$cpg_id) %in% L$cpg_id)
  J <- merge(R, L, by = c("cpg_id", "target_gene"))
  cov_pair <- nrow(J) / nrow(R)
  report <- c(report, "### Q1 coverage", "",
              sprintf("- reference CpGs with >= 1 blood-class link: %.1f%%", 100 * cov_cpg),
              sprintf("- reference pairs that are linked: %.1f%% (%s)", 100 * cov_pair, format(nrow(J), big.mark = ",")), "")
  # Q2
  q2 <- list(rE2G_all = pr(R, L[, list(cpg_id, target_gene)]))
  for (b in levels(L$score_bin)) q2[[paste0("rE2G_score_", b)]] <- pr(R, L[score_bin == b, list(cpg_id, target_gene)])
  for (b in levels(L$bs_bin))    q2[[paste0("rE2G_biosamples_", b)]] <- pr(R, L[bs_bin == b, list(cpg_id, target_gene)])
  if (!is.null(A)) {
    q2$manifest_all <- pr(R, A)
    # the fair comparison: manifest restricted to the CpGs that have links
    q2$manifest_on_linked_cpgs <- pr(R[cpg_id %in% L$cpg_id], A)
    q2$rE2G_on_manifest_cpgs   <- pr(R[cpg_id %in% A$cpg_id], L[, list(cpg_id, target_gene)])
  }
  t2 <- rbindlist(lapply(names(q2), function(k) c(list(proposal = k), q2[[k]])))
  report <- c(report, "### Q2 target proposal quality", "", md_table(t2), "")
  print(t2)
  # Q3
  if (!is.null(A) && any(!is.na(R$dist))) {
    inA <- R[A, on = c("cpg_id", "target_gene"), nomatch = NULL][, list(cpg_id, target_gene)]
    R[, in_manifest := paste(cpg_id, target_gene) %in% paste(inA$cpg_id, inA$target_gene)]
    R[, linked := paste(cpg_id, target_gene) %in% paste(L$cpg_id, L$target_gene)]
    t3 <- R[, list(n = .N, median_dist = sprintf("%.0f", median(dist, na.rm = TRUE)),
                   p90_dist = sprintf("%.0f", quantile(dist, 0.9, na.rm = TRUE))),
            by = c("in_manifest", "linked")][order(in_manifest, linked)]
    report <- c(report, "### Q3 distance profile of reference pairs by manifest / link status", "", md_table(t3), "",
                "Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.", "")
    print(t3)
  }
  # Q4
  if (nrow(J)) {
    base_neg <- mean(R$ref_dir == -1)
    t4a <- J[, list(n = .N, neg_share = fmt_ci(sum(ref_dir == -1), .N)), by = "promoter_class"][order(promoter_class)]
    t4b <- J[, list(n = .N, neg_share = fmt_ci(sum(ref_dir == -1), .N)), by = "score_bin"][order(score_bin)]
    t4c <- J[, list(n = .N, neg_share = fmt_ci(sum(ref_dir == -1), .N)), by = "bs_bin"][order(bs_bin)]
    report <- c(report, sprintf("### Q4 sign among linked pairs (reference base rate of -1: %.3f)", base_neg), "",
                "By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):", "",
                md_table(t4a), "", "By rE2G score:", "", md_table(t4b), "", "By number of biosamples:", "", md_table(t4c), "")
    print(t4a); print(t4b)
  }
}
report <- c(report, "## Reading the result", "",
            "- The layer earns its keep as a TARGET source if rE2G precision on the shared CpGs is at or above the manifest's",
            "  and its recall adds pairs the manifest lacks (Q2, Q3).",
            "- It must NOT be promoted to a direction rung on the strength of Q4: a negative-share above base rate in",
            "  promoter-class links is the TSS-distance effect, already on the ladder as distance_only.",
            "- Score bins should show a precision gradient; if they do not, `distal_score_max` is provenance, not a confidence.", "")
dir.create(dirname(out_md), showWarnings = FALSE, recursive = TRUE)
writeLines(report, out_md)
message("wrote ", out_md)
