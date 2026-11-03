#!/usr/bin/env Rscript
# tools/finngen_smr_celltype.R -- cell-type-specific blood SMR directions:
# GoDMC clumped cis-mQTL instruments (exposure) x FinnGen single-nucleus
# multiome cis-eQTLs (outcome), one context per cell type, in the schema of
# the package's bulk `smr_directions` layer plus a `context` column.
#
# Exposure (mQTL). godmc_instruments.csv.gz, produced by godmc_extract.sh from
#   GoDMC assoc_meta_all.csv.gz (Min et al. 2021, Nat Genet; hg19): rows with
#   clumped == TRUE and p < 1e-5. Columns: cpg, snp (chr:pos:SNP), beta_a1, se,
#   pval, samplesize, allele1, allele2, freq_a1, cistrans, num_studies, hetpval.
#   Effects are oriented to allele1. Only cis rows (cistrans == TRUE) are used.
#   SNP rsIDs come from the 1000G EUR PLINK .bim (hg19) by chr:pos + allele set;
#   palindromic strand-ambiguous pairs and allele mismatches are refused.
#   CpG hg19 positions come from cpgdirection::cpgd_cpg_positions().
# Outcome (eQTL). FinnGen BESDs built by tools/finngen_to_smr.R (rsID, hg19).
# LD reference. The same 1000G EUR panel as the bulk layer.
#
# Settings. SMR v1.03 defaults plus --peqtl-smr 5e-8 (the package's stated
#   threshold); HEIDI defaults (--heidi-min-m 3, --heidi-max-m 20) match the
#   "fewer than three instruments" / "11-20 instruments" wording of the bulk
#   layer's documentation. HEIDI is reported, never applied (package policy).
#
# Tiers (as in ?cpgd_smr_directions). For each CpG-gene pair, every clumped
#   GoDMC instrument of the CpG that has an eQTL record for the gene gives a
#   Wald-ratio direction sign(b_eQTL / b_mQTL) after allele alignment; an
#   instrument counts when its Wald p < 5e-8.
#   S1: > 1 counting instrument, all agreeing.  S2: exactly one.  S3: counting
#   instruments disagree (direction = majority; ties dropped).
#   p_SMR, p_HEIDI, nsnp_HEIDI and top_instrument come from SMR's top-SNP test.
#
# Usage
#   Rscript tools/finngen_smr_celltype.R besd --instruments <godmc_instruments.csv.gz> \
#       --bim <1000G_EUR.bim> --out <prefix>
#   Rscript tools/finngen_smr_celltype.R run --mqtl <prefix> --bfile <1000G_EUR> \
#       --eqtl <finngen BESD prefix> --context blood_<ct> --out <prefix>
#   Output of "run": <out>.smr, <out>.instruments.txt, <out>.directions.csv.gz:
#   cpg_id target_gene ensembl_id direction smr_tier p_SMR n_instruments
#   instrument_agreement cpg_gene_dist top_instrument p_HEIDI nsnp_HEIDI
#   heidi_status context

suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
stopifnot(mode %in% c("besd", "run"))
comp <- c(A = "T", C = "G", G = "C", T = "A")

if (mode == "besd") {
  inst <- get_arg("--instruments"); bim <- get_arg("--bim"); out <- get_arg("--out")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

  g <- fread(inst)
  message(nrow(g), " GoDMC clumped instruments read")
  g <- g[cistrans == TRUE | cistrans == "TRUE"]
  message(nrow(g), " cis rows kept")
  g[, c("chr", "pos", "type") := tstrsplit(snp, ":", fixed = TRUE)]
  g <- g[type == "SNP"]
  g[, chr := sub("^chr", "", chr)]; g[, pos := as.integer(pos)]
  g[, allele1 := toupper(allele1)]; g[, allele2 := toupper(allele2)]

  b <- fread(bim, header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  b[, chr := as.character(chr)]
  m <- merge(g, b[, .(chr, pos, rsid, bim_a1 = a1, bim_a2 = a2)], by = c("chr", "pos"),
             all = FALSE, allow.cartesian = TRUE)
  same  <- (m$allele1 == m$bim_a1 & m$allele2 == m$bim_a2) |
           (m$allele1 == m$bim_a2 & m$allele2 == m$bim_a1)
  snv   <- nchar(m$allele1) == 1L & nchar(m$allele2) == 1L
  palin <- snv & m$allele1 == comp[m$allele2]
  message(sum(same), " allele-identical; ", sum(!same & palin), " palindromic refused; ",
          sum(!same & !palin), " allele mismatches refused; ",
          nrow(g) - length(unique(m$snp)), " instrument SNPs absent from the .bim")
  m <- m[same]
  m <- m[!(rsid %in% rsid[duplicated(unique(m[, .(snp, rsid)])$rsid)])]

  suppressPackageStartupMessages(library(cpgdirection))
  P <- cpgdirection::cpgd_cpg_positions()
  m <- merge(m, P[, .(cpg = cpg_id, cpg_chr = sub("^chr", "", as.character(chr)), cpg_pos = pos)],
             by = "cpg", all = FALSE)
  message(nrow(m), " instrument-CpG rows with rsID and CpG position")

  q <- m[, .(SNP = rsid, Chr = chr, BP = pos, A1 = allele1, A2 = allele2, Freq = freq_a1,
             Probe = cpg, Probe_Chr = cpg_chr, Probe_bp = cpg_pos, Gene = cpg,
             Orientation = "+", b = beta_a1, SE = se, p = pval)]
  qfile <- paste0(out, ".query.txt")
  fwrite(q, qfile, sep = "\t")
  fwrite(m[, .(cpg, rsid, chr, pos, allele1, allele2, beta_a1, se, pval, cpg_chr, cpg_pos)],
         paste0(out, ".instruments.csv.gz"), compress = "gzip")
  cmd <- sprintf("smr --qfile %s --make-besd --out %s", shQuote(qfile), shQuote(out))
  message("running: ", cmd)
  if (system(cmd) != 0) stop("smr --make-besd failed")
  message("mQTL BESD written: ", out, ".{besd,esi,epi}")
  quit(status = 0)
}

if (mode == "run") {
  mqtl <- get_arg("--mqtl"); bfile <- get_arg("--bfile"); eqtl <- get_arg("--eqtl")
  ctx  <- get_arg("--context"); out <- get_arg("--out")
  p_thr <- as.numeric(get_arg("--peqtl-smr", "5e-8"))
  extract <- get_arg("--extract-snp", "")
  extract_opt <- if (nzchar(extract)) paste(" --extract-snp", shQuote(extract)) else ""

  cmd <- sprintf("smr --bfile %s --beqtl-summary %s --beqtl-summary %s%s --peqtl-smr %s --thread-num 4 --out %s",
                 shQuote(bfile), shQuote(mqtl), shQuote(eqtl), extract_opt, format(p_thr, scientific = TRUE), shQuote(out))
  message("running: ", cmd)
  if (system(cmd) != 0) stop("smr failed")
  S <- fread(paste0(out, ".smr"))
  message(nrow(S), " CpG-gene pairs tested by SMR")

  I <- fread(paste0(mqtl, ".instruments.csv.gz"))
  snpfile <- paste0(out, ".instrument_snps.txt")
  writeLines(unique(I$rsid), snpfile)
  cmd2 <- sprintf("smr --beqtl-summary %s --extract-snp %s --query 1 --out %s",
                  shQuote(eqtl), shQuote(snpfile), shQuote(paste0(out, ".instruments")))
  message("running: ", cmd2)
  if (system(cmd2) != 0) stop("smr --query failed")
  E <- fread(paste0(out, ".instruments.txt"))
  setnames(E, c("SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp",
                "Gene", "Orientation", "b", "SE", "p"))

  W <- merge(I[, .(cpg, rsid, allele1, allele2, b_m = beta_a1, se_m = se, cpg_pos)],
             E[, .(rsid = SNP, A1, A2, Probe, Gene, Probe_bp, b_e = b, se_e = SE)],
             by = "rsid", allow.cartesian = TRUE)
  flip <- W$allele1 == W$A2 & W$allele2 == W$A1
  keep <- (W$allele1 == W$A1 & W$allele2 == W$A2) | flip
  n_refused <- sum(!keep)
  W <- W[keep]; W[flip[keep], b_e := -b_e]
  W[, wald := b_e / b_m]
  W[, wald_se := sqrt((se_e / b_m)^2 + (b_e * se_m / b_m^2)^2)]
  W[, wald_p := 2 * pnorm(-abs(wald / wald_se))]
  W <- W[wald_p < p_thr]
  message(n_refused, " instrument-gene comparisons refused on allele grounds; ",
          nrow(W), " instrument-gene Wald ratios at p < ", p_thr)

  A <- W[, .(n_instruments = .N,
             n_pos = sum(wald > 0), n_neg = sum(wald < 0),
             cpg_gene_dist = Probe_bp[1] - cpg_pos[1]), by = .(cpg, Probe, Gene)]
  A[, direction := fifelse(n_pos > n_neg, 1L, fifelse(n_neg > n_pos, -1L, NA_integer_))]
  A[, instrument_agreement := pmax(n_pos, n_neg) / n_instruments]
  A[, smr_tier := fifelse(n_instruments > 1 & instrument_agreement == 1, "S1",
                   fifelse(n_instruments == 1, "S2", "S3"))]
  A <- A[!is.na(direction)]

  S2 <- S[, .(cpg = Expo_ID, Probe = Outco_ID, top_instrument = topSNP, p_SMR, p_HEIDI, nsnp_HEIDI)]
  R <- merge(A, S2, by = c("cpg", "Probe"), all.x = TRUE)
  R[, heidi_status := fifelse(is.na(p_HEIDI), "not_tested", fifelse(p_HEIDI >= 0.01, "pass", "fail"))]
  R[, context := ctx]
  setnames(R, c("cpg", "Gene", "Probe"), c("cpg_id", "target_gene", "ensembl_id"))
  R <- R[, .(cpg_id, target_gene, ensembl_id, direction, smr_tier, p_SMR, n_instruments,
             instrument_agreement, cpg_gene_dist, top_instrument, p_HEIDI, nsnp_HEIDI,
             heidi_status, context)]
  fwrite(R, paste0(out, ".directions.csv.gz"), compress = "gzip")
  message(ctx, ": ", nrow(R), " pairs -- S1 ", sum(R$smr_tier == "S1"), ", S2 ",
          sum(R$smr_tier == "S2"), ", S3 ", sum(R$smr_tier == "S3"),
          "; HEIDI pass/fail/not_tested ", sum(R$heidi_status == "pass"), "/",
          sum(R$heidi_status == "fail"), "/", sum(R$heidi_status == "not_tested"))
}
