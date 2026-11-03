#!/usr/bin/env Rscript
# tools/finngen_to_smr.R -- turn FinnGen single-nucleus multiome cis-eQTL
# nominal statistics (v1, GRCh38, no rsIDs) into SMR BESD inputs keyed on rsID
# and hg19 positions, so they can be run against the package's GoDMC blood
# mQTL instruments with the same SMR v1.03 settings as the bulk layer.
#
# Source
#   FinnGen multiome summary statistics v1 (Kanai et al. 2025, medRxiv
#   doi:10.1101/2025.11.25.25340489; Nature 2026 doi:10.1038/s41586-026-11078-2).
#   Public after the access form at finngen.fi/en/access_results;
#   gs://finngen-public-data-multiome/ ; README.md shipped with the release.
#   eQTL/cis_nominal/finngen_multiome_v1.eQTL.{l1|l2}.{CELL_TYPE}.cis_nominal.tsv.gz
#   columns: #CHR POS cell_type phenotype_id MarkerID AF_Allele2 BETA SE p.value
#   MarkerID = {CHR}_{POS}_{REF}_{ALT}; BETA and AF are oriented to ALT (allele 2).
#   Variants within 1 Mb of the TSS, minor allele count > 20; SAIGE-QTL.
#   All coordinates GRCh38. Genes are Ensembl IDs.
#
# mode "map"   build the variant dictionary once (FinnGen MarkerID -> rsID,
#              hg19 position, A1/A2) by lifting every FinnGen variant from hg38
#              to hg19 and matching chr + pos + allele set to the 1000G EUR
#              PLINK .bim used as the SMR LD reference. Strand flips of
#              palindromic (A/T, C/G) pairs are refused, not assumed.
# mode "query" per cell type: stream the nominal file chromosome by chromosome
#              through tabix, join the dictionary and the hg19 gene table
#              (EnsDb.Hsapiens.v75), write an SMR "query" file, run
#              `smr --make-besd`.
#
# Dependencies: data.table, rtracklayer, GenomicRanges, EnsDb.Hsapiens.v75,
#               tabix (htslib) and smr (v1.03) on PATH; hg38ToHg19.over.chain.
#
# Usage
#   Rscript tools/finngen_to_smr.R map --nominal <PBMC cis_nominal.tsv.gz> \
#       --bim <1000G_EUR.bim> --chain <hg38ToHg19.over.chain> --out <map.tsv.gz>
#   Rscript tools/finngen_to_smr.R query --nominal <cell type cis_nominal.tsv.gz> \
#       --map <map.tsv.gz> --out <prefix>      # writes <prefix>.query.txt + BESD
#
# Query columns: SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b SE p
#   SNP = rsID; Chr/BP = hg19 (from the .bim); A1 = FinnGen ALT (effect allele),
#   A2 = REF; Freq = AF_Allele2; Probe = Ensembl ID; Probe_bp = hg19 TSS;
#   Gene = HGNC symbol; Orientation = strand.

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
stopifnot(mode %in% c("map", "query"))

comp <- c(A = "T", C = "G", G = "C", T = "A")

if (mode == "map") {
  suppressPackageStartupMessages({ library(rtracklayer); library(GenomicRanges) })
  nominal <- get_arg("--nominal"); bim <- get_arg("--bim")
  chain   <- get_arg("--chain");   out <- get_arg("--out")

  message("reading FinnGen variant universe from ", basename(nominal))
  v <- fread(cmd = sprintf("gzip -dc %s | cut -f5", shQuote(nominal)),
             header = TRUE, col.names = "MarkerID", showProgress = FALSE)
  v <- unique(v)
  v[, c("chr38", "pos38", "ref", "alt") := tstrsplit(MarkerID, "_", fixed = TRUE)]
  v[, pos38 := as.integer(pos38)]
  v[, chr38 := ifelse(grepl("^chr", chr38), chr38, paste0("chr", chr38))]
  message(nrow(v), " unique variants")

  message("liftOver hg38 -> hg19")
  gr <- GRanges(v$chr38, IRanges(v$pos38, v$pos38))
  lo <- liftOver(gr, import.chain(chain))
  n  <- lengths(lo)
  keep <- n == 1L
  v <- v[keep]; lo <- unlist(lo[keep])
  v[, chr19 := sub("^chr", "", as.character(seqnames(lo)))]
  v[, pos19 := start(lo)]
  message(sum(!keep), " variants dropped at liftOver; ", nrow(v), " kept")

  message("reading .bim ", basename(bim))
  b <- fread(bim, header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  b[, chr := as.character(chr)]
  setkey(b, chr, pos)

  m <- merge(v, b[, .(chr19 = chr, pos19 = pos, rsid, bim_a1 = a1, bim_a2 = a2)],
             by = c("chr19", "pos19"), all = FALSE, allow.cartesian = TRUE)
  same  <- (m$ref == m$bim_a1 & m$alt == m$bim_a2) | (m$ref == m$bim_a2 & m$alt == m$bim_a1)
  flipc <- (comp[m$ref] == m$bim_a1 & comp[m$alt] == m$bim_a2) |
           (comp[m$ref] == m$bim_a2 & comp[m$alt] == m$bim_a1)
  palin <- m$ref == comp[m$alt]
  snv <- nchar(m$ref) == 1L & nchar(m$alt) == 1L
  flipc[!snv] <- FALSE; palin[!snv] <- FALSE
  ok <- same
  message(sum(same), " allele-identical; ", sum(flipc & !palin & !same),
          " matched on the other strand (non-palindromic); ",
          sum(palin & !same), " palindromic strand-ambiguous refused; ",
          sum(!ok & !palin), " allele mismatches refused")
  m <- m[ok]
  m <- m[!duplicated(MarkerID)]
  m <- m[!(rsid %in% rsid[duplicated(rsid)])]
  out_tab <- m[, .(MarkerID, rsid, chr19, pos19, A1 = alt, A2 = ref)]
  fwrite(out_tab, out, sep = "\t", compress = "gzip")
  message("wrote ", nrow(out_tab), " mapped variants to ", out)
  quit(status = 0)
}

if (mode == "query") {
  suppressPackageStartupMessages({ library(EnsDb.Hsapiens.v75) })
  nominal <- get_arg("--nominal"); mapf <- get_arg("--map"); out <- get_arg("--out")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

  map <- fread(mapf); setkey(map, MarkerID)

  g <- as.data.table(genes(EnsDb.Hsapiens.v75, return.type = "data.frame",
                           columns = c("gene_id", "gene_name", "seq_name",
                                       "gene_seq_start", "gene_seq_end", "seq_strand")))
  g <- g[seq_name %in% c(1:22, "X")]
  g[, tss := ifelse(seq_strand == 1L, gene_seq_start, gene_seq_end)]
  g <- g[, .(phenotype_id = gene_id, Gene = gene_name, Probe_Chr = seq_name,
             Probe_bp = tss, Orientation = ifelse(seq_strand == 1L, "+", "-"))]
  setkey(g, phenotype_id)

  chrs <- c(1:22, "X")
  qfile <- paste0(out, ".query.txt")
  if (file.exists(qfile)) file.remove(qfile)
  first <- TRUE; n_tot <- 0L; n_kept <- 0L
  for (ch in chrs) {
    rows <- NULL
    for (region in c(paste0("chr", ch), as.character(ch))) {
      rows <- tryCatch(fread(cmd = sprintf("tabix %s %s", shQuote(nominal), region),
                             header = FALSE, showProgress = FALSE), error = function(e) NULL)
      if (!is.null(rows) && nrow(rows)) break
    }
    if (is.null(rows) || !nrow(rows)) next
    setnames(rows, c("CHR", "POS", "cell_type", "phenotype_id", "MarkerID",
                     "AF_Allele2", "BETA", "SE", "p.value"))
    n_tot <- n_tot + nrow(rows)
    rows[, phenotype_id := sub("\\..*$", "", phenotype_id)]
    x <- map[rows, on = "MarkerID", nomatch = NULL]
    x <- g[x, on = "phenotype_id", nomatch = NULL]
    x <- x[is.finite(BETA) & is.finite(SE) & SE > 0]
    q <- x[, .(SNP = rsid, Chr = chr19, BP = pos19, A1, A2, Freq = AF_Allele2,
               Probe = phenotype_id, Probe_Chr, Probe_bp, Gene, Orientation,
               b = BETA, SE, p = p.value)]
    n_kept <- n_kept + nrow(q)
    fwrite(q, qfile, sep = "\t", append = !first, col.names = first)
    first <- FALSE
    message("chr", ch, ": ", nrow(rows), " pairs read, ", nrow(q), " written")
  }
  message("total ", n_tot, " pairs read; ", n_kept, " with rsID + hg19 gene written to ", qfile)

  cmd <- sprintf("smr --qfile %s --make-besd --out %s", shQuote(qfile), shQuote(out))
  message("running: ", cmd)
  status <- system(cmd)
  if (status != 0) stop("smr --make-besd failed")
  message("BESD written: ", out, ".{besd,esi,epi}")
}
