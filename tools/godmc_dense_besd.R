#!/usr/bin/env Rscript
# tools/godmc_dense_besd.R -- dense GoDMC cis-mQTL exposure BESD for SMR + HEIDI.
# Reads every cis row of GoDMC assoc_meta_all.csv.gz (Min et al. 2021; hg19;
# all meta-analysed associations at p < 1e-5), so that HEIDI has the mQTL
# effects of the SNPs in LD around each top instrument. Rows with clumped ==
# TRUE are the instruments and are written to <out>.instruments.csv.gz in the
# layout tools/finngen_smr_celltype.R expects. rsIDs from the 1000G EUR .bim
# (hg19) by chr:pos + allele set; CpG positions from cpgd_cpg_positions().
# Usage: Rscript tools/godmc_dense_besd.R --assoc <assoc_meta_all.csv.gz> \
#          --bim <1000G_EUR.bim> --out <prefix>
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) { i <- match(flag, args); if (is.na(i)) stop("missing ", flag); args[i + 1L] }
assoc <- get_arg("--assoc"); bim <- get_arg("--bim"); out <- get_arg("--out")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

message("streaming cis rows from ", basename(assoc), " (10-20 min)")
awk <- "NR>1 { gsub(/\\r/,\"\"); if ($11==\"TRUE\" || $11==\"\\\"TRUE\\\"\") print $1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$7\",\"$8\",\"$9\",\"$NF }"
g <- fread(cmd = sprintf("gzip -dc %s | awk -F, '%s'", shQuote(assoc), awk), header = FALSE,
           col.names = c("cpg", "snp", "beta_a1", "se", "pval", "allele1", "allele2", "freq_a1", "clumped"),
           showProgress = FALSE)
for (col in c("cpg", "snp", "allele1", "allele2", "clumped")) set(g, j = col, value = gsub('"', "", g[[col]]))
message(nrow(g), " cis rows; ", sum(g$clumped == "TRUE"), " clumped instruments")
g[, c("chr", "pos", "type") := tstrsplit(snp, ":", fixed = TRUE)]
g <- g[type == "SNP"]; g[, chr := sub("^chr", "", chr)]; g[, pos := as.integer(pos)]
g[, allele1 := toupper(allele1)]; g[, allele2 := toupper(allele2)]

b <- fread(bim, header = FALSE, col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
b[, chr := as.character(chr)]
m <- merge(g, b[, .(chr, pos, rsid, bim_a1 = a1, bim_a2 = a2)], by = c("chr", "pos"), all = FALSE, allow.cartesian = TRUE)
same <- (m$allele1 == m$bim_a1 & m$allele2 == m$bim_a2) | (m$allele1 == m$bim_a2 & m$allele2 == m$bim_a1)
message(sum(same), " allele-identical of ", nrow(m), " position matches; ",
        length(unique(g$snp)) - length(unique(m$snp)), " SNPs absent from the .bim")
m <- m[same]
dup <- unique(m[, .(snp, rsid)])[, .N, by = rsid][N > 1, rsid]
m <- m[!(rsid %in% dup)]
m <- m[!duplicated(m[, .(cpg, rsid)])]

suppressPackageStartupMessages(library(cpgdirection))
P <- cpgdirection::cpgd_cpg_positions()
m <- merge(m, P[, .(cpg = cpg_id, cpg_chr = sub("^chr", "", as.character(chr)), cpg_pos = pos)], by = "cpg", all = FALSE)
message(nrow(m), " SNP-CpG rows with rsID and CpG position; ",
        sum(m$clumped == "TRUE"), " of them instruments; ", length(unique(m$cpg)), " CpGs")

q <- m[, .(SNP = rsid, Chr = chr, BP = pos, A1 = allele1, A2 = allele2, Freq = freq_a1,
           Probe = cpg, Probe_Chr = cpg_chr, Probe_bp = cpg_pos, Gene = cpg, Orientation = "+",
           b = beta_a1, SE = se, p = pval)]
qfile <- paste0(out, ".query.txt")
fwrite(q, qfile, sep = "\t")
fwrite(m[clumped == "TRUE", .(cpg, rsid, chr, pos, allele1, allele2, beta_a1, se, pval, cpg_chr, cpg_pos)],
       paste0(out, ".instruments.csv.gz"), compress = "gzip")
cmd <- sprintf("smr --qfile %s --make-besd --out %s", shQuote(qfile), shQuote(out))
message("running: ", cmd)
if (system(cmd) != 0) stop("smr --make-besd failed")
message("dense mQTL BESD written: ", out, ".{besd,esi,epi}")
