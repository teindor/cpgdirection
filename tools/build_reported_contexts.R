#!/usr/bin/env Rscript
# Reported tissue contexts (PI rulings of 11-12 Nov 2026): CpG-gene signs measured or estimated in a tissue the
# package has no voting class for. They are REPORTED: never part of the measured table, never a vote, never a
# source of best_direction. Each failed the admission bar for a voting source (sign agreement >= 0.85 with the
# pairs the package already measures), or is a single SMR context.
#
#   retina             Advani et al. 2024, eQTM_significant_pairs_retina.txt (top CpG per gene, adjusted p <= 0.05)
#   liver              Bonder et al. 2014 (BMC Genomics, Supplementary Table 5) and the Horvath liver pairs of
#                      Wang et al. 2021 (Table S2): two independent cohorts; a pair both hold with the same sign is
#                      replicated within liver (n_cohorts = 2)
#   pancreatic_islets  Olsson et al. 2014 (PLoS Genet, S25 Table a)
#   fibroblasts        Wagner et al. 2014 (Genome Biol, Additional file 13; hg18 positions mapped to array CpGs)
#   saliva             SMR: GoDMC cis-mQTL instruments against saliva cis-eQTLs (lead 15); every pair kept, with
#                      flags for the MHC and the 17q21.31 inversion
#
# The published tables are parsed by tools/parse_reported_sources.py into one long table (--long: source, tissue,
# cpg_id, target_gene, direction, p). Here: genes named as the catalogue names them (cpgdirectionData gene table,
# else the unique Ensembl GRCh37 symbol), CpGs placed by the package position table (else the Illumina 450K
# manifest, build 37), pairs on one chromosome within the window, one vote per cohort, majority sign, ties dropped.
#
# Usage (package root):
#   Rscript tools/build_reported_contexts.R --long <sources_long.csv.gz> --wang <measured_eqtms_add_wang2021.csv.gz> \
#     --saliva <saliva_outcome_godmc_pairs_250kb.csv.gz> --gene-table <gene_table.csv.gz> --ens <ens_grch37.csv.gz> \
#     [--hm450 <GPL13534 manifest csv.gz>] [--data-dir <dir>] [--window 250000]
# Writes inst/extdata/reported_contexts.csv.gz and tools/reports/build_reported_contexts.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
long_f <- get_arg("--long"); wang_f <- get_arg("--wang"); sal_f <- get_arg("--saliva")
gt_f <- get_arg("--gene-table"); ens_f <- get_arg("--ens"); hm_f <- get_arg("--hm450"); data_dir <- get_arg("--data-dir")
window <- as.numeric(get_arg("--window", "250000"))
for (f in list(long_f, wang_f, sal_f, gt_f, ens_f)) stopifnot(!is.null(f), file.exists(path.expand(f)))
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)))
suppressMessages(devtools::load_all(".", quiet = TRUE))
wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f] (n=%s)", k / n, ci[1], ci[2], format(n, big.mark = ",")) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(trimws(x), collapse = " | "), " |"))) }
big <- function(x) format(x, big.mark = ",")

# ---- genes and positions (as tools/ingest_nov2026_measured_sources.R) --------------------------
GT <- fread(gt_f)
EN <- unique(fread(ens_f), by = "ensgene")
EN[, `:=`(chr = paste0("chr", sub("^chr", "", as.character(chr))), sym = toupper(trimws(as.character(symbol))))]
EN[, tss := fifelse(strand == 1, as.numeric(start), as.numeric(end))]
EN <- EN[nzchar(sym) & sym != "NA" & chr %in% paste0("chr", c(1:22, "X", "Y"))]
EN[, n_sym := .N, by = sym]
E2 <- EN[!ensgene %in% GT$gene_id & n_sym == 1L & !sym %in% toupper(GT$display), .(gene_id = ensgene, target_gene = sym, gchr = chr, tss)]
GM <- rbind(GT[, .(gene_id, target_gene = toupper(display), gchr = chr, tss = as.numeric(tss))], E2)
stopifnot(!anyDuplicated(GM$target_gene))
P <- as.data.table(cpgd_cpg_positions())[, .(cpg_id, chr, pos = as.numeric(pos))]
if (!is.null(hm_f)) {
  H <- suppressWarnings(fread(path.expand(hm_f), skip = "IlmnID,", fill = TRUE, select = c("Name", "Genome_Build", "CHR", "MAPINFO"), colClasses = "character"))
  H <- H[grepl("^c[gh]\\d", Name) & Genome_Build == "37" & CHR %in% c(1:22, "X", "Y") & grepl("^\\d+$", MAPINFO),
         .(cpg_id = Name, chr = paste0("chr", CHR), pos = as.numeric(MAPINFO))]
  P <- rbind(P, unique(H, by = "cpg_id")[!cpg_id %in% P$cpg_id])
}

# ---- measured contexts ------------------------------------------------------------------------
ctx_of <- c(advani2024_retina = "retina", bonder2014_liver = "liver", olsson2014_islets = "pancreatic_islets", wagner2014_fibroblast = "fibroblasts")
coh_of <- c(advani2024_retina = "Advani2024_retina", bonder2014_liver = "Bonder2014_liver", olsson2014_islets = "Olsson2014_islets", wagner2014_fibroblast = "Wagner2014_fibroblasts")
L <- fread(long_f)[source %in% names(ctx_of) & direction %in% c(-1, 1)]
n_read <- L[, .(rows_read = .N), by = source]
L <- L[, .(sym = toupper(trimws(unlist(strsplit(as.character(target_gene), "\\s*///\\s*|\\s*;\\s*")))), cpg_id, direction), by = .(source, row = seq_len(nrow(L)))]
L <- merge(L[, row := NULL], GM, by.x = "sym", by.y = "target_gene")[, target_gene := sym][, sym := NULL]
W <- merge(L, P, by = "cpg_id")
W <- W[chr == gchr][, tss_dist := as.integer(round(tss - pos))][abs(tss_dist) <= window]
C <- W[, .(direction = { s <- sum(direction); if (s == 0) NA_real_ else sign(s) }, tss_dist = tss_dist[1]), by = .(source, cpg_id, target_gene)][!is.na(direction)]
C[, `:=`(context = unname(ctx_of[source]), cohort = unname(coh_of[source]))]
# the Horvath liver pairs of Wang 2021 (already in the package's gene names and coordinates)
HL <- fread(wang_f)[cohorts == "Horvath_liver" & direction %in% c(-1, 1) & abs(tss_dist) <= window,
                    .(source = "wang2021_S2", cpg_id, target_gene = toupper(trimws(target_gene)), direction = as.numeric(direction),
                      tss_dist = as.integer(tss_dist), context = "liver", cohort = "Horvath_liver")]
C <- rbind(C, HL, use.names = TRUE)
liver_x <- merge(C[cohort == "Bonder2014_liver", .(cpg_id, target_gene, a = direction)], C[cohort == "Horvath_liver", .(cpg_id, target_gene, b = direction)], by = c("cpg_id", "target_gene"))
R <- C[, .(plus = sum(direction > 0), minus = sum(direction < 0), tss_dist = tss_dist[1], cohorts = paste(sort(unique(cohort)), collapse = "|"),
           source = paste(sort(unique(source)), collapse = "|")), by = .(context, cpg_id, target_gene)]
n_tie <- R[plus == minus, .N]; R <- R[plus != minus]
R[, `:=`(direction = as.integer(sign(plus - minus)), n_cohorts = plus + minus, cohort_agreement = round(pmax(plus, minus) / (plus + minus), 3),
         evidence = "measured_eqtm", flags = NA_character_)]

# ---- saliva (SMR) -------------------------------------------------------------------------------
S <- fread(sal_f)
stopifnot(all(c("cpg_id", "target_gene", "direction", "tss_dist", "mhc", "inv17") %in% names(S)), all(abs(S$tss_dist) <= window))
SV <- S[, .(context = "saliva", cpg_id, target_gene = toupper(trimws(target_gene)), tss_dist = as.integer(round(tss_dist)), cohorts = "GoDMC_x_saliva_eQTL",
            source = "saliva_smr_lead15", direction = as.integer(direction), n_cohorts = 1L, cohort_agreement = 1, evidence = "smr",
            flags = fifelse(as.logical(mhc) & as.logical(inv17), "mhc|inv17", fifelse(as.logical(mhc), "mhc", fifelse(as.logical(inv17), "inv17", NA_character_))))]
SV <- unique(SV, by = c("cpg_id", "target_gene"))
cols <- c("cpg_id", "target_gene", "context", "direction", "tss_dist", "n_cohorts", "cohort_agreement", "cohorts", "source", "evidence", "flags")
OUT <- rbind(R[, cols, with = FALSE], SV[, cols, with = FALSE])
setorder(OUT, context, cpg_id, target_gene)
stopifnot(!anyDuplicated(OUT[, .(context, cpg_id, target_gene)]), all(OUT$direction %in% c(-1L, 1L)), all(abs(OUT$tss_dist) <= window))

# ---- agreement with the measured table (the admission test; why these are reported, not voting) --
M <- as.data.table(cpgd_measured_eqtms()); M[, target_gene := toupper(trimws(target_gene))]
anyM <- M[, .(s = sum(direction)), by = .(cpg_id, target_gene)][s != 0][, r := sign(s)]
tab <- OUT[, { x <- merge(.SD, anyM, by = c("cpg_id", "target_gene"))
               .(pairs = big(.N), cpgs = big(uniqueN(cpg_id)), genes = big(uniqueN(target_gene)), negative = sprintf("%.3f", mean(direction == -1)),
                 two_cohorts_agree = big(sum(n_cohorts >= 2L & cohort_agreement == 1)), agreement_with_measured = fmt(sum(x$direction == x$r), nrow(x))) },
           by = .(context, evidence), .SDcols = c("cpg_id", "target_gene", "direction", "n_cohorts", "cohort_agreement")]
print(tab)
fwrite(OUT, "inst/extdata/reported_contexts.csv.gz", compress = "gzip")
rep <- c("# Reported tissue contexts", "", sprintf("Run %s. Window %s bp. Rows %s. Never a vote and never a source of best_direction.", format(Sys.time(), "%Y-%m-%d %H:%M"), big(window), big(nrow(OUT))), "",
         md(tab), "", "`agreement_with_measured`: sign agreement with every pair of the measured table (any tissue, majority sign); the bar for a voting source is 0.85.", "",
         sprintf("- rows read from the parsed tables: %s", paste(sprintf("%s %s", n_read$source, big(n_read$rows_read)), collapse = "; ")),
         sprintf("- liver: Bonder 2014 and Horvath liver share %s pairs, same sign in %s; pairs dropped as cohort ties: %s", big(nrow(liver_x)), big(sum(liver_x$a == liver_x$b)), big(n_tie)),
         sprintf("- saliva flags: MHC %s, 17q21.31 inversion %s", big(sum(grepl("mhc", OUT$flags))), big(sum(grepl("inv17", OUT$flags)))), "")
dir.create("tools/reports", showWarnings = FALSE)
writeLines(rep, "tools/reports/build_reported_contexts.md")
cat(paste(rep, collapse = "\n"), "\n")
