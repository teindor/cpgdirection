#!/usr/bin/env Rscript
# Lead 3 (Nov 2026 scan, b2): T2T-CHM13 / pangenome probe uniqueness for EPIC v2
# (and 450K), from Dong ... Korthauer 2025, Cell Rep 44:115755, MIT-licensed sets at
# github.com/functionalepigenomics/Illumina_Infinium_HumanMethylation_BeadChips_Annotation
# (files dated 6 Jan 2025). Cross-checked against the Zhou-lab EPICv2 general mask
# (InfiniumAnnotationV1, hg38 manifest + mask) and validated against the package's
# own direction sources.
#
# Sandbox result (3 Nov 2026, EPICv2, 933,252 cg probes in the Zhou hg38 manifest):
#   T2T-ambiguous probes            74,719  of which Zhou-masked 25,293 (34%), Zhou-clean 49,426
#   Pan+T2T-ambiguous probes       111,971  of which Zhou-masked 29,092 (26%), Zhou-clean 82,879
#   Zhou-masked but T2T-unambiguous  6,533  (mostly SNP-in-probe masks, 1,080 M_nonuniq)
#   per CpG: 73,943 CpGs have EVERY probe T2T-ambiguous; 48,914 of those are not Zhou-masked.
# So the flag is large and mostly new. It enters cpgd_probe_qc() as INFORMATIONAL
# columns; whether it should feed qc_exclude is what section 4 below tests.
#
# Usage (package root):
#   Rscript tools/build_probe_uniqueness.R [--dir ~/probe_uniqueness] [--data-dir <dir>]
# Writes <dir>/epicv2_probe_uniqueness.csv.gz (per cpg_id; candidate Hub resource),
#        <dir>/dong2025_probe_uniqueness.csv.gz (per probe, EPICv2 + 450K),
#        <dir>/probe_uniqueness_t2t.md (cross-tabs and validation).

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
dir <- path.expand(get_arg("--dir", "~/probe_uniqueness"))
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)
dir.create(file.path(dir, "dl"), showWarnings = FALSE, recursive = TRUE)
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else suppressPackageStartupMessages(library(cpgdirection))

# ---- 1. downloads --------------------------------------------------------------------
dong <- "https://raw.githubusercontent.com/functionalepigenomics/Illumina_Infinium_HumanMethylation_BeadChips_Annotation/main"
zhou <- "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/EPICv2"
sets <- data.table(
  array = c(rep("EPICv2", 5), "HM450K", "HM450K"),
  set   = c("t2t", "pangenome_xpop", "pangenome_eur", "pan_t2t", "pan_t2t_eur", "t2t", "pan_t2t"),
  path  = c("T2T/T2T_EPICv2_annotate.bed.tar.xz", "Pangenome/Pangenome_EPICv2_annotate.bed.tar.xz",
            "Pangenome/Pangenome_EUR_EPICv2_annotate.bed.tar.xz", "Pangenome+T2T/Pan+T2T_EPICv2_annotate.bed.tar.xz",
            "Pangenome+T2T/Pan+T2T_EUR_EPICv2_annotate.bed.tar.xz", "T2T/T2T_HM450K_annotate.bed.tar.xz",
            "Pangenome+T2T/Pan+T2T_HM450K_annotate.bed.tar.xz"))
sets[, file := file.path(dir, "dl", basename(path))]
sets[, bed := sub("\\.tar\\.xz$", "", file)]
fetch <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 1000) return(invisible())
  message("downloading ", basename(dest))
  st <- system2("curl", c("-sS", "-L", "--fail", "--retry", "3", "-o", shQuote(dest), shQuote(url)))
  if (st != 0) stop("download failed: ", url)
}
for (k in seq_len(nrow(sets))) {
  fetch(file.path(dong, sets$path[k]), sets$file[k])
  if (!file.exists(sets$bed[k])) {
    system2("tar", c("-xJf", shQuote(sets$file[k]), "-C", shQuote(file.path(dir, "dl"))))
    unlink(list.files(file.path(dir, "dl"), "^\\._", full.names = TRUE))
  }
}
fetch(file.path(zhou, "EPICv2.hg38.manifest.tsv.gz"), file.path(dir, "dl", "EPICv2.hg38.manifest.tsv.gz"))
fetch(file.path(zhou, "EPICv2.hg38.mask.tsv.gz"),     file.path(dir, "dl", "EPICv2.hg38.mask.tsv.gz"))

# ---- 2. per-probe table -------------------------------------------------------------------
read_bed <- function(f) fread(f, sep = "\t", select = c("CpG", "Chr", "Position", "All_freq", "EUR_freq", "Repeat"),
                              na.strings = "NA", showProgress = FALSE)
P <- rbindlist(lapply(unique(sets$array), function(arr) {
  s <- sets[array == arr]
  ids <- lapply(seq_len(nrow(s)), function(k) read_bed(s$bed[k]))
  names(ids) <- s$set
  base <- ids[["t2t"]][, .(probe_id = CpG, chm13_chr = Chr, chm13_pos = Position,
                           snp_maf05_all = !is.na(All_freq), snp_maf05_eur = !is.na(EUR_freq),
                           repeat_class = sub("^.*#", "", Repeat))]
  allp <- unique(unlist(lapply(ids, `[[`, "CpG")))
  base <- rbindlist(list(base, data.table(probe_id = setdiff(allp, base$probe_id))), fill = TRUE)
  for (st in c("t2t", "pangenome_xpop", "pangenome_eur", "pan_t2t", "pan_t2t_eur"))
    base[, (st) := if (st %in% names(ids)) probe_id %in% ids[[st]]$CpG else NA]
  base[, array := arr]
  base[, cpg_id := sub("_.*$", "", probe_id)]
  setcolorder(base, c("array", "probe_id", "cpg_id"))
  message(sprintf("%s: %s probes in any set; T2T %s, cross-population pangenome %s, Pan+T2T %s", arr,
                  format(nrow(base), big.mark = ","), format(sum(base$t2t), big.mark = ","),
                  format(sum(base$pangenome_xpop, na.rm = TRUE), big.mark = ","), format(sum(base$pan_t2t), big.mark = ",")))
  base
}), use.names = TRUE, fill = TRUE)
fwrite(P, file.path(dir, "dong2025_probe_uniqueness.csv.gz"), compress = "gzip")

# ---- 3. EPICv2 universe from the Zhou manifest; cross-tab with the Zhou mask -----------------
man  <- fread(file.path(dir, "dl", "EPICv2.hg38.manifest.tsv.gz"), select = c("Probe_ID", "CpG_chrm"), showProgress = FALSE)
mask <- fread(file.path(dir, "dl", "EPICv2.hg38.mask.tsv.gz"), showProgress = FALSE)
E <- merge(man[grepl("^cg", Probe_ID), .(probe_id = Probe_ID)],
           mask[, .(probe_id = Probe_ID, zhou_general = M_general == "TRUE", maskUniq)], by = "probe_id", all.x = TRUE)
E[is.na(zhou_general), zhou_general := FALSE]
E <- merge(E, P[array == "EPICv2", .(probe_id, t2t, pangenome_xpop, pan_t2t, pan_t2t_eur, snp_maf05_all, chm13_chr, chm13_pos)],
           by = "probe_id", all.x = TRUE)
for (st in c("t2t", "pangenome_xpop", "pan_t2t", "pan_t2t_eur")) E[is.na(get(st)), (st) := FALSE]
E[, cpg_id := sub("_.*$", "", probe_id)]
ct <- function(flag) { x <- E[, .N, by = .(unambiguous = get(flag), zhou_general)][order(-unambiguous, zhou_general)]
  x[, set := flag]; x }
xt <- rbindlist(lapply(c("t2t", "pangenome_xpop", "pan_t2t"), ct))
message("\n== EPICv2 probes: uniqueness vs Zhou general mask"); print(xt)

# per CpG: every-replicate rule for the exclusion-grade reading, any-replicate for the informational one
C <- E[, .(n_probes = .N,
           n_t2t_unambiguous = sum(t2t), n_pan_t2t_unambiguous = sum(pan_t2t),
           n_pangenome_xpop_unambiguous = sum(pangenome_xpop), n_pan_t2t_eur_unambiguous = sum(pan_t2t_eur),
           t2t_ambiguous = all(!t2t), pan_t2t_ambiguous = all(!pan_t2t),
           pangenome_xpop_ambiguous = all(!pangenome_xpop), pan_t2t_eur_ambiguous = all(!pan_t2t_eur),
           t2t_snp_maf05 = any(snp_maf05_all %in% TRUE),
           chm13_chr = chm13_chr[which(t2t)[1]], chm13_pos = chm13_pos[which(t2t)[1]],
           zhou_general_all = all(zhou_general)), by = cpg_id]
C[, uniqueness_source := "Dong2025_CellRep_T2T_pangenome"]
Q <- as.data.table(cpgd_probe_qc())
C <- merge(C, Q[, .(cpg_id, qc_exclude, masked_general, cross_hybridizing, ch_blat_predicted, mapping_flagged)],
           by = "cpg_id", all.x = TRUE)
message(sprintf("\nEPICv2 CpGs %s; T2T-ambiguous (every probe) %s; Pan+T2T-ambiguous %s; of the T2T-ambiguous, qc_exclude already %s, Garvan cross-hyb %s, BLAT-predicted %s",
                format(nrow(C), big.mark = ","), format(sum(C$t2t_ambiguous), big.mark = ","), format(sum(C$pan_t2t_ambiguous), big.mark = ","),
                format(sum(C$t2t_ambiguous & C$qc_exclude %in% TRUE), big.mark = ","),
                format(sum(C$t2t_ambiguous & C$cross_hybridizing %in% TRUE), big.mark = ","),
                format(sum(C$t2t_ambiguous & C$ch_blat_predicted %in% TRUE), big.mark = ",")))
C[, class := fcase(qc_exclude %in% TRUE, "qc_exclude",
                   t2t_ambiguous & pan_t2t_ambiguous, "ambiguous_T2T_and_pangenome",
                   t2t_ambiguous, "ambiguous_T2T_only",
                   pan_t2t_ambiguous, "ambiguous_pangenome_only",
                   default = "unambiguous")]
cls <- C[, .N, by = class][order(-N)]; print(cls)

# ---- 4. does ambiguity predict worse direction evidence? -------------------------------------
wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2]) }
M <- as.data.table(cpgd_measured_eqtms()); M[, gene := toupper(trimws(target_gene))]
S <- as.data.table(cpgd_smr_directions());  S[, gene := toupper(trimws(target_gene))]
if (!"smr_rung" %in% names(S)) S[, smr_rung := fcase(smr_tier == "S1", "high", smr_tier == "S2", "moderate", default = "weak")]
# (a) measured blood vs SMR (blood), by probe class -- two independent sign sources
V <- merge(M[tissue == "blood", .(cpg_id, gene, m_dir = direction)],
           S[smr_rung %in% c("replicated", "high", "moderate"), .(cpg_id, gene, s_dir = direction, smr_rung)], by = c("cpg_id", "gene"))
V <- merge(V, C[, .(cpg_id, class)], by = "cpg_id", all.x = TRUE); V[is.na(class), class := "not_on_EPICv2"]
va <- V[, .(pairs = .N, cpgs = uniqueN(cpg_id), agree = fmt(sum(m_dir == s_dir), .N)), by = class][order(class)]
message("\n== measured blood vs SMR sign agreement, by probe class"); print(va)
# (b) measured blood vs measured nasal / solid (cross-tissue replication of a measured sign)
W <- merge(M[tissue == "blood", .(cpg_id, gene, b_dir = direction)],
           M[tissue != "blood", .(cpg_id, gene, o_dir = direction, tissue)], by = c("cpg_id", "gene"))
W <- merge(W, C[, .(cpg_id, class)], by = "cpg_id", all.x = TRUE); W[is.na(class), class := "not_on_EPICv2"]
wa <- W[, .(pairs = .N, agree = fmt(sum(b_dir == o_dir), .N)), by = class][order(class)]
message("\n== measured blood vs measured other-tissue sign agreement, by probe class"); print(wa)
# (c) does ambiguity change how often a CpG has any measured eQTM / any SMR pair (coverage, not accuracy)
cov <- C[, .(cpgs = .N, with_measured = mean(cpg_id %in% M$cpg_id), with_smr = mean(cpg_id %in% S$cpg_id)), by = class][order(class)]
message("\n== coverage by class"); print(cov)
# (d) multi-cohort measured pairs: does cohort agreement drop for ambiguous probes?
if ("cohort_agreement" %in% names(M)) {
  mc <- merge(M[n_cohorts >= 2, .(cpg_id, gene, cohort_agreement)], C[, .(cpg_id, class)], by = "cpg_id")
  ca <- mc[, .(pairs = .N, mean_cohort_agreement = sprintf("%.3f", mean(cohort_agreement, na.rm = TRUE))), by = class][order(class)]
  message("\n== cohort agreement of multi-cohort measured pairs, by class"); print(ca)
} else ca <- NULL

# ---- 5. write ------------------------------------------------------------------------------------
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
out <- C[, .(cpg_id, n_probes, n_t2t_unambiguous, n_pan_t2t_unambiguous, n_pangenome_xpop_unambiguous, n_pan_t2t_eur_unambiguous,
             t2t_ambiguous, pan_t2t_ambiguous, pangenome_xpop_ambiguous, pan_t2t_eur_ambiguous, t2t_snp_maf05,
             chm13_chr, chm13_pos, uniqueness_source)]
fwrite(out, file.path(dir, "epicv2_probe_uniqueness.csv.gz"), compress = "gzip")
rep <- c("# EPICv2 probe uniqueness under T2T-CHM13 and the pangenome (Dong et al. 2025) -- build and validation", "",
         sprintf("Run %s. Sources: functionalepigenomics annotation repo (MIT), Zhou InfiniumAnnotationV1 EPICv2 hg38 manifest + mask.", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
         "## Probes: uniqueness vs Zhou general mask", "", md(xt), "",
         "## CpGs by class (every-replicate rule)", "", md(cls), "",
         "## Measured blood vs SMR sign agreement, by class", "", md(va), "",
         "## Measured blood vs measured other-tissue sign agreement, by class", "", md(wa), "",
         "## Coverage by class", "", md(cov), "",
         if (!is.null(ca)) c("## Cohort agreement of multi-cohort measured pairs, by class", "", md(ca), "") else character(),
         "## Reading", "",
         "If the ambiguous classes agree with the unambiguous class within their Wilson intervals, the flags stay informational",
         "(columns on cpgd_probe_qc(); no change to qc_exclude). If an ambiguous class is clearly worse on (a) AND (b), that class",
         "is a candidate for qc_exclude in a later version, after the same check on the Project Alpha arrays (probe-level",
         "detection p / bead counts) -- the Dong result is in silico mapping, not measured failure.", "")
writeLines(rep, file.path(dir, "probe_uniqueness_t2t.md"))
message("\nwrote ", file.path(dir, "epicv2_probe_uniqueness.csv.gz"), ", dong2025_probe_uniqueness.csv.gz and probe_uniqueness_t2t.md")
