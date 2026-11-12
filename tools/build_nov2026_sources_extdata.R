#!/usr/bin/env Rscript
# Fold the four measured sources admitted on 12 Nov 2026 (tools/ingest_nov2026_measured_sources.R) into the
# bundled measured layer. Each cleared the evidence bar (sign agreement >= 0.85 with the pairs the package
# already measures) and each is one cohort and one donor group.
#
#   inst/extdata/measured_eqtms.csv.gz  + kidney tubule (Liu 2022)       solid_tissue  Liu2022_kidney         liu2022_kidney
#                                       + placenta (Delahaye 2018)       solid_tissue  Delahaye2018_placenta  delahaye2018_placenta
#                                       + adipose (Voisin 2016)          solid_tissue  Voisin2016_adipose     voisin2016_adipose
#                                       + whole blood, GTP (Kennedy 2018) blood        Kennedy2018_GTP_blood  kennedy2018_gtp
#
# Cohort rule as in tools/build_kreitmaier2026_extdata.R and inst/scripts/build_measured_eqtms.R: one vote per
# cohort per pair, sign by majority, ties dropped; n_cohorts / cohort_agreement / cohorts recomputed. A pair
# already in the table keeps its `source` and `tss_dist`.
#
# Usage (package root):
#   Rscript tools/build_nov2026_sources_extdata.R [--add tools/out/measured_eqtms_add_nov2026_sources.csv.gz]
# Writes inst/extdata/measured_eqtms.csv.gz and tools/reports/build_nov2026_sources_extdata.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
add_f <- get_arg("--add", "tools/out/measured_eqtms_add_nov2026_sources.csv.gz")
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)), file.exists(add_f))
big <- function(x) format(x, big.mark = ",")
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(trimws(x), collapse = " | "), " |"))) }
new_src <- c(liu2022_kidney = "solid_tissue", delahaye2018_placenta = "solid_tissue", voisin2016_adipose = "solid_tissue", kennedy2018_gtp = "blood")

E0 <- fread("inst/extdata/measured_eqtms.csv.gz")
cols <- names(E0)
if (any(E0$source %in% names(new_src)) || any(grepl("Liu2022_kidney|Delahaye2018_placenta|Voisin2016_adipose|Kennedy2018_GTP_blood", E0$cohorts)))
  stop("the measured table already holds these sources; restore the previous table before rebuilding")
E <- copy(E0); E[, target_gene := toupper(trimws(target_gene))]
stopifnot(!anyDuplicated(E[, .(cpg_id, target_gene, tissue)]))
K <- fread(add_f); K[, target_gene := toupper(trimws(target_gene))]
setorder(K, source, cpg_id, target_gene)
stopifnot(all(K$source %in% names(new_src)), all(K$tissue == new_src[K$source]), all(K$direction %in% c(-1L, 1L)), all(K$n_cohorts == 1L),
          all(abs(K$tss_dist) <= 250000), !anyDuplicated(K[, .(source, cpg_id, target_gene)]))
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
           source = source[1], had = any(existing)),
       by = .(cpg_id, target_gene, tissue)]
M[, n_cohorts := plus + minus]
n_tie <- M[plus == minus, .N, by = tissue]
M <- M[plus != minus]
M[, direction := as.integer(sign(plus - minus))]
M[, cohort_agreement := round(pmax(plus, minus) / n_cohorts, 3)]
chk <- merge(M[, .(cpg_id, target_gene, tissue, direction)], old, by = c("cpg_id", "target_gene", "tissue"))
flips <- chk[direction != old, .N, by = tissue]
per <- rbindlist(lapply(names(new_src), function(s) {
  k <- K[source == s]; x <- merge(k[, .(cpg_id, target_gene, tissue, kd = direction)], old, by = c("cpg_id", "target_gene", "tissue"))
  data.table(source = s, tissue = unname(new_src[s]), pairs = big(nrow(k)), already_measured_in_tissue = big(nrow(x)),
             same_sign = big(sum(x$kd == x$old)), opposite = big(sum(x$kd != x$old)), new_to_tissue = big(nrow(k) - nrow(x))) }))
grp <- function(x) vapply(strsplit(x, "|", fixed = TRUE), function(v) {
  g <- ifelse(startsWith(v, "GTEx"), "GTEx", ifelse(startsWith(v, "Kreitmaier2026"), "Kreitmaier2026", ifelse(startsWith(v, "EVA-PR"), "EVAPR",
       ifelse(startsWith(v, "MESA"), "MESA", ifelse(startsWith(v, "HELIX"), "HELIX", v))))); length(unique(g)) }, integer(1))
rep2 <- function(D) D[, .(n = sum(grp(cohorts) >= 2L & cohort_agreement >= 0.99)), by = tissue]
r0 <- rep2(E); M <- M[, cols, with = FALSE]
setorder(M, tissue, cpg_id, target_gene)
stopifnot(!anyDuplicated(M[, .(cpg_id, target_gene, tissue)]), all(M$direction %in% c(-1L, 1L)), nrow(M) >= nrow(E0) - sum(n_tie$N))
r1 <- rep2(M)
fwrite(M, "inst/extdata/measured_eqtms.csv.gz", compress = "gzip")
tab <- merge(E0[, .(before = .N), by = tissue], M[, .(after = .N, multi_cohort = sum(n_cohorts > 1)), by = tissue], by = "tissue", all = TRUE)
tab <- merge(tab, merge(r0[, .(tissue, two_groups_agree_before = n)], r1[, .(tissue, two_groups_agree_after = n)], by = "tissue"), by = "tissue", all.x = TRUE)
message("measured table:"); print(tab); print(per)
one <- function(d) if (nrow(d)) paste(sprintf("%s %s", d$tissue, big(d$N)), collapse = "; ") else "none"
rep <- c("# Kidney, placenta, adipose and GTP whole-blood eQTMs folded into the measured layer", "",
         sprintf("Run %s. Add file: `%s` (%s rows, %s CpGs, %s genes).", format(Sys.time(), "%Y-%m-%d %H:%M"), add_f,
                 big(nrow(K)), big(uniqueN(K$cpg_id)), big(uniqueN(K$target_gene))), "", md(tab), "",
         "`two_groups_agree`: pairs seen in two or more donor groups with every cohort agreeing.", "", md(per), "",
         sprintf("- pairs dropped as cohort ties: %s", one(n_tie)),
         sprintf("- existing pairs whose sign changed: %s", one(flips)),
         sprintf("- sources after: %s", paste(sprintf("%s %s", names(table(M$source)), trimws(big(as.integer(table(M$source))))), collapse = "; ")), "")
dir.create("tools/reports", showWarnings = FALSE)
writeLines(rep, "tools/reports/build_nov2026_sources_extdata.md")
cat(paste(rep, collapse = "\n"), "\n")
