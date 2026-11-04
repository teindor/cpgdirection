#!/usr/bin/env Rscript
# Lead 6 (Nov 2026 scan, c3): does the package's direction evidence behave differently
# for stimulus-responsive genes (cytokines, chemokines, interferon-stimulated, immediate-
# early) than for the rest? Dwaraka ... Burd (bioRxiv 10.64898/2026.07.14.738581) found,
# in 333 T-cell-enriched blood samples, that methylation predicts lineage/identity genes
# well and stimulus-responsive genes poorly, from distal enhancer CpGs rather than
# promoters. If that holds in our own sources, a gene-class covariate should scale
# confidence down for those genes (reported first; it earns a place in the ladder only
# if the agreement gap is clear and replicates across the checks below).
#
# Everything here comes from tables the package already ships plus GO.db / org.Hs.eg.db;
# nothing is downloaded. Usage (package root):
#   Rscript tools/gene_class_covariate.R [--out tools/reports/gene_class_covariate.md] [--data-dir <dir>]

suppressPackageStartupMessages({ library(data.table); library(AnnotationDbi); library(org.Hs.eg.db); library(GO.db) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
out_file <- get_arg("--out", "tools/reports/gene_class_covariate.md")
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)
if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else suppressPackageStartupMessages(library(cpgdirection))

# ---- gene classes ---------------------------------------------------------------------------
go_genes <- function(go_id) {
  ids <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = go_id, keytype = "GOALL", columns = "SYMBOL"),
                  error = function(e) NULL)
  if (is.null(ids)) character() else unique(toupper(ids$SYMBOL))
}
classes <- list(
  cytokine      = go_genes("GO:0005125"),   # cytokine activity
  chemokine     = go_genes("GO:0008009"),   # chemokine activity
  ifn_response  = go_genes("GO:0034340"),   # response to type I interferon
  immediate_early = c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "EGR4",
                      "IER2", "IER3", "IER5", "NR4A1", "NR4A2", "NR4A3", "ATF3", "DUSP1", "DUSP2", "DUSP5",
                      "KLF2", "KLF4", "KLF6", "ZFP36", "BTG2", "CYR61", "CCN1", "CTGF", "CCN2", "MYC", "ARC",
                      "NPAS4", "GADD45B", "PPP1R15A", "SOCS3", "TNFAIP3", "NFKBIA", "NFKBIZ", "CXCL8", "IL6",
                      "TNF", "CCL2", "CCL3", "CCL4", "PTGS2", "SRF", "ELK1", "MAFF", "HES1", "ID1", "ID2", "ID3"),
  checkpoint_ligand = c("CD274", "PDCD1LG2", "CD80", "CD86", "LGALS9", "TNFRSF14", "PVR", "NECTIN2", "VSIR", "CD276", "VTCN1", "HHLA2"))
stim <- unique(unlist(classes))
message(sprintf("stimulus-responsive set: %s genes (cytokine %d, chemokine %d, IFN %d, IEG %d, checkpoint %d)",
                length(stim), length(classes$cytokine), length(classes$chemokine), length(classes$ifn_response),
                length(classes$immediate_early), length(classes$checkpoint_ligand)))
# lineage / identity reference set: haematopoietic lineage markers via GO "immune system development" would be
# too broad; use the complement ("other") as the reference and report the canonical T/B/myeloid identity genes
# separately so the Dwaraka contrast is visible.
identity <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD2", "CD5", "CD7", "LCK", "ZAP70", "TCF7", "LEF1", "BCL11B",
              "GATA3", "TBX21", "RORC", "FOXP3", "IL7R", "CCR7", "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "EBF1",
              "CD14", "LYZ", "CSF1R", "CEBPA", "CEBPB", "SPI1", "ITGAM", "FCGR3A", "NCAM1", "KLRD1", "NKG7", "GNLY",
              "PTPRC", "RUNX1", "RUNX3", "IKZF1", "IKZF3", "ETS1", "MPO", "ELANE", "CSF3R", "FCGR3B")
gene_class <- function(g) fcase(g %in% stim, "stimulus_responsive", g %in% identity, "lineage_identity", default = "other")

wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
fmt <- function(k, n) if (n == 0) "-" else { ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2]) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }

M <- as.data.table(cpgd_measured_eqtms()); M[, gene := toupper(trimws(target_gene))]; M[, class := gene_class(gene)]
S <- as.data.table(cpgd_smr_directions());  S[, gene := toupper(trimws(target_gene))]; S[, class := gene_class(gene)]
if (!"smr_rung" %in% names(S)) S[, smr_rung := fcase(smr_tier == "S1", "high", smr_tier == "S2", "moderate", default = "weak")]

# (1) independent-source sign agreement: measured blood vs blood SMR, by gene class
V <- merge(M[tissue == "blood", .(cpg_id, gene, class, m_dir = direction, tss_dist)],
           S[smr_rung %in% c("replicated", "high", "moderate"), .(cpg_id, gene, s_dir = direction)], by = c("cpg_id", "gene"))
t1 <- V[, .(pairs = .N, genes = uniqueN(gene), agree = fmt(sum(m_dir == s_dir), .N)), by = class][order(class)]
# (2) cross-tissue agreement of measured signs (blood vs nasal / solid), by class
W <- merge(M[tissue == "blood", .(cpg_id, gene, class, b_dir = direction)],
           M[tissue != "blood", .(cpg_id, gene, o_dir = direction)], by = c("cpg_id", "gene"))
t2 <- W[, .(pairs = .N, genes = uniqueN(gene), agree = fmt(sum(b_dir == o_dir), .N)), by = class][order(class)]
# (3) cohort agreement within tissue for multi-cohort pairs, by class
t3 <- if ("cohort_agreement" %in% names(M)) M[n_cohorts >= 2, .(pairs = .N, mean_cohort_agreement = sprintf("%.3f", mean(cohort_agreement))), by = .(tissue, class)][order(tissue, class)] else NULL
# (4) where the measured CpGs sit: share of measured pairs beyond 10 kb of the TSS (Dwaraka: distal for the well-predicted genes)
t4 <- M[, .(pairs = .N, share_beyond_10kb = sprintf("%.3f", mean(abs(tss_dist) > 10000)), share_negative = sprintf("%.3f", mean(direction == -1))), by = .(tissue, class)][order(tissue, class)]
# (5) SMR instrument strength by class: share of pairs at the replicated/high rungs
t5 <- S[, .(pairs = .N, genes = uniqueN(gene), share_replicated_or_high = sprintf("%.3f", mean(smr_rung %in% c("replicated", "high")))), by = class][order(class)]
# (6) coverage: of the genes in each class that are expressed-and-annotated anywhere in the package, how many have ANY measured blood pair
ann <- tryCatch(as.data.table(cpgd_probe_gene_annotation()), error = function(e) NULL)
t6 <- if (!is.null(ann)) { ann[, gene := toupper(trimws(target_gene))]; G <- unique(ann[, .(gene)]); G[, class := gene_class(gene)]
  G[, .(genes = .N, with_measured_blood = sprintf("%.3f", mean(gene %in% M[tissue == "blood", gene])),
        with_smr = sprintf("%.3f", mean(gene %in% S$gene))), by = class][order(class)] } else NULL

for (x in list(t1, t2, t3, t4, t5, t6)) if (!is.null(x)) print(x)
rep <- c("# Gene-class covariate: stimulus-responsive vs lineage vs other genes, in the package's own evidence", "",
         sprintf("Run %s. Classes: GO cytokine activity, chemokine activity, response to type I interferon, an immediate-early list and checkpoint ligands (%d genes); %d lineage/identity genes; everything else 'other'.",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), length(stim), length(identity)), "",
         "## 1. Measured blood vs blood SMR sign agreement", "", md(t1), "",
         "## 2. Measured blood vs measured other-tissue sign agreement", "", md(t2), "",
         if (!is.null(t3)) c("## 3. Cohort agreement of multi-cohort measured pairs", "", md(t3), "") else character(),
         "## 4. Position and sign of measured pairs", "", md(t4), "",
         "## 5. SMR rung mix", "", md(t5), "",
         if (!is.null(t6)) c("## 6. Coverage of annotated genes", "", md(t6), "") else character(),
         "## Reading", "",
         "A covariate is justified if stimulus-responsive genes are clearly lower on 1 AND 2 (non-overlapping Wilson intervals",
         "with 'other'), in which case it enters as a reported column (gene_class) and a confidence note, not a rung change;",
         "the stress-response pathways DMSA tests are exactly where it would bite. Lower coverage alone (6) is not evidence of",
         "wrong signs and changes nothing.", "")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
writeLines(rep, out_file); message("wrote ", out_file)
