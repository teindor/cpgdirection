#!/usr/bin/env Rscript
# Fold the osteoarthritis joint-tissue eQTMs (Kreitmaier et al. 2026; cartilage and synovium, full statistics
# from the MSK Knowledge Portal; tools/ingest_kreitmaier2026_oa.R) into the bundled measured layer as a
# solid-tissue source (PI ruling 11 Nov 2026: pairs significant at the study's own threshold only; the
# source is not widened to sub-threshold pairs; fat pad is not included).
#
#   inst/extdata/measured_eqtms.csv.gz  + cartilage and synovium pairs, tissue = solid_tissue,
#                                         cohorts Kreitmaier2026_cartilage / Kreitmaier2026_synovium
#                                         (one donor group, .cpgd_donor_group()), source = kreitmaier2026_msk
#
# Cohort rule as in tools/build_nov2026_extdata.R and inst/scripts/build_measured_eqtms.R: one vote per
# cohort per pair, sign by majority, ties dropped; n_cohorts / cohort_agreement / cohorts recomputed.
#
# Usage (package root):
#   Rscript tools/build_kreitmaier2026_extdata.R [--add tools/out/measured_eqtms_add_kreitmaier2026.csv.gz]
# Writes inst/extdata/measured_eqtms.csv.gz and tools/reports/build_kreitmaier2026_extdata.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
add_f <- get_arg("--add", "tools/out/measured_eqtms_add_kreitmaier2026.csv.gz")
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)), file.exists(add_f))
big <- function(x) format(x, big.mark = ",")
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }

E0 <- fread("inst/extdata/measured_eqtms.csv.gz")
cols <- names(E0)
if (any(E0$source == "kreitmaier2026_msk")) stop("the measured table already holds kreitmaier2026_msk rows; restore the previous table before rebuilding")
E <- copy(E0); E[, target_gene := toupper(trimws(target_gene))]
stopifnot(!anyDuplicated(E[, .(cpg_id, target_gene, tissue)]))
K <- fread(add_f); K[, target_gene := toupper(trimws(target_gene))]
stopifnot(all(K$tissue == "solid_tissue"), all(K$direction %in% c(-1L, 1L)), all(K$source == "kreitmaier2026_msk"),
          all(abs(K$tss_dist) <= 250000), !anyDuplicated(K[, .(cpg_id, target_gene)]))
votes <- function(D, existing) {
  D <- copy(D)
  D[, n := as.integer(n_cohorts)]
  D[, k := as.integer(round(n * cohort_agreement))]
  D[, `:=`(plus = fifelse(direction > 0, k, n - k), minus = fifelse(direction > 0, n - k, k))]
  D[, .(cpg_id, target_gene, tissue, plus, minus, tss_dist, cohorts, source, existing = existing)]
}
V <- rbindlist(list(votes(E, TRUE), votes(K, FALSE)), use.names = TRUE)
setorder(V, cpg_id, target_gene, tissue, -existing)
old <- E[, .(cpg_id, target_gene, tissue, old = direction)]
M <- V[, .(plus = sum(plus), minus = sum(minus), tss_dist = tss_dist[!is.na(tss_dist)][1],
           cohorts = paste(unique(unlist(strsplit(cohorts, "|", fixed = TRUE))), collapse = "|"),
           source = source[1], n_src = .N),
       by = .(cpg_id, target_gene, tissue)]
M[, n_cohorts := plus + minus]
kk <- paste(K$cpg_id, K$target_gene)
in_k <- M$tissue == "solid_tissue" & paste(M$cpg_id, M$target_gene) %in% kk
n_new <- sum(in_k & M$n_src == 1L)
n_shared <- sum(in_k & M$n_src == 2L)
n_tie <- sum(M$plus == M$minus)
M <- M[plus != minus]
M[, direction := as.integer(sign(plus - minus))]
M[, cohort_agreement := round(pmax(plus, minus) / n_cohorts, 3)]
chk <- merge(M[, .(cpg_id, target_gene, tissue, direction)], old, by = c("cpg_id", "target_gene", "tissue"))
n_flip <- sum(chk$direction != chk$old)
ks <- merge(K[, .(cpg_id, target_gene, kd = direction)], old[tissue == "solid_tissue"], by = c("cpg_id", "target_gene"))
grp <- function(x) vapply(strsplit(x, "|", fixed = TRUE), function(v) {
  g <- ifelse(startsWith(v, "GTEx"), "GTEx", ifelse(startsWith(v, "Kreitmaier2026"), "Kreitmaier2026", v)); length(unique(g)) }, integer(1))
rep2_before <- sum(grp(E[tissue == "solid_tissue"]$cohorts) >= 2L & E[tissue == "solid_tissue"]$cohort_agreement >= 0.99)
M <- M[, cols, with = FALSE]
setorder(M, tissue, cpg_id, target_gene)
stopifnot(!anyDuplicated(M[, .(cpg_id, target_gene, tissue)]), all(M$direction %in% c(-1L, 1L)))
rep2_after <- sum(grp(M[tissue == "solid_tissue"]$cohorts) >= 2L & M[tissue == "solid_tissue"]$cohort_agreement >= 0.99)
fwrite(M, "inst/extdata/measured_eqtms.csv.gz", compress = "gzip")
tab <- merge(E0[, .(before = .N), by = tissue], M[, .(after = .N, multi_cohort = sum(n_cohorts > 1)), by = tissue], by = "tissue", all = TRUE)
message("measured table:"); print(tab)
rep <- c("# Osteoarthritis joint-tissue eQTMs folded into the measured layer", "",
         sprintf("Run %s. Add file: `%s` (%s pairs, %s CpGs, %s genes).", format(Sys.time(), "%Y-%m-%d %H:%M"), add_f,
                 big(nrow(K)), big(uniqueN(K$cpg_id)), big(uniqueN(K$target_gene))), "", md(tab), "",
         sprintf("- new solid-tissue pairs: %s; pairs already measured in solid tissue: %s (same sign %s, opposite %s)",
                 big(n_new), big(n_shared), big(sum(ks$kd == ks$old)), big(sum(ks$kd != ks$old))),
         sprintf("- pairs dropped as cohort ties: %s; existing pairs whose sign changed: %s", big(n_tie), big(n_flip)),
         sprintf("- solid-tissue pairs seen in 2+ donor groups with every cohort agreeing: %s before, %s after", big(rep2_before), big(rep2_after)),
         sprintf("- sources after: %s", paste(sprintf("%s %s", names(table(M$source)), big(as.integer(table(M$source)))), collapse = "; ")), "")
dir.create("tools/reports", showWarnings = FALSE)
writeLines(rep, "tools/reports/build_kreitmaier2026_extdata.md")
cat(paste(rep, collapse = "\n"), "\n")
