#!/usr/bin/env Rscript
# CpG coordinates and array membership for every probe of the four Illumina methylation arrays (HM27, 450K,
# EPIC v1, EPIC v2), so that cpgd_cpg_annotation() can place any array probe in hg19 and hg38.
#
# hg19. The package's own table (Hub resource cpg_positions_hg19, EPIC v2 probes) is what the direction engine
#   uses; it is never changed and not repeated here. This table adds hg19 for probes that table does not hold
#   (450K and EPIC v1 probes that EPIC v2 dropped): the Illumina 450K manifest v1.1 (GEO GPL13534, build 37),
#   else the Illumina EPIC v1 coordinates (CHR, MAPINFO, build 37).
# hg38. One position per probe, from the Zhou-lab (SeSAMe / InfiniumAnnotation) hg38 manifests, taking the
#   manifest of the first array generation that carries the probe: HM450, then EPIC, then EPICv2, then HM27;
#   position = CpG_beg + 1 (1-based first base of the CpG on the forward strand). Where no manifest places the
#   probe, the position of the package's EPIC v2 annotation layer is used. The table stores an hg38 position
#   only where it is not simply the annotation layer's value (cpgd_cpg_annotation() falls back to that layer).
#   `hg38_ambiguous` marks probes that two sources place differently (multi-mapping probes; EPIC v2 replicate
#   probes mapping to different places): they are placed by the same priority, and flagged.
# arrays. Which arrays carry the probe, for every probe.
# Probe IDs are matched without regard to case (the package stores ch. and nv- IDs in lower case) and written
# in the manifests' spelling.
#
# Usage (package root):
#   Rscript tools/build_cpg_coordinates.R --hm450-hg19 <GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz> \
#     --epic-hg19 <tsv.gz: cgid, CHR, MAPINFO> --hm450-hg38 <tsv.gz: probe_id, chr_hg38, beg_hg38> \
#     --epic-hg38 <same columns> --epicv2-hg38 <same columns, replicate suffix removed, distinct rows> \
#     --hm27-hg38 <same columns> [--check <tsv.gz: cpg, chr, pos; hg38 positions to reproduce>] [--data-dir <dir>]
# Writes inst/extdata/cpg_coordinates.csv.gz and tools/reports/build_cpg_coordinates.md.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
f450_19 <- get_arg("--hm450-hg19"); fe_19 <- get_arg("--epic-hg19"); fe_38 <- get_arg("--epic-hg38")
f450_38 <- get_arg("--hm450-hg38"); f27_38 <- get_arg("--hm27-hg38"); fv2_38 <- get_arg("--epicv2-hg38")
chk_f <- get_arg("--check"); data_dir <- get_arg("--data-dir")
for (f in list(f450_19, fe_19, fe_38, f450_38, f27_38, fv2_38)) stopifnot(!is.null(f), file.exists(path.expand(f)))
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir)
stopifnot(file.exists("DESCRIPTION"), grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1)))
suppressMessages(devtools::load_all(".", quiet = TRUE))
big <- function(x) format(x, big.mark = ",")
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(trimws(x), collapse = " | "), " |"))) }
ok_id <- function(x) grepl("^(cg[0-9]|ch\\.|rs[0-9]|nv-)", x, ignore.case = TRUE)
std_chr <- paste0("chr", c(1:22, "X", "Y", "M"))

# ---- what the package holds (EPIC v2), keyed in lower case -------------------------------------
P19 <- as.data.table(cpgd_cpg_positions())[, .(key = tolower(cpg_id), chr19 = chr, pos19 = as.numeric(pos))]
A <- cpgdirection:::.cpgd_data("epicv2_probe_gene_annotation")
P38 <- unique(A[!is.na(position) & grepl("^chr", chr), .(key = tolower(cpg_id), chr38 = chr, pos38 = as.numeric(position))], by = "key")
v2_pkg <- unique(c(P19$key, tolower(A$cpg_id)))

# ---- outside sources ---------------------------------------------------------------------------
H19 <- suppressWarnings(fread(path.expand(f450_19), skip = "IlmnID,", fill = TRUE, select = c("Name", "Genome_Build", "CHR", "MAPINFO"), colClasses = "character"))
H19 <- H19[ok_id(Name)]; ids450_ill <- unique(H19$Name)
H19 <- unique(H19[Genome_Build == "37" & CHR %in% c(1:22, "X", "Y", "MT", "M") & grepl("^\\d+$", MAPINFO),
                  .(key = tolower(Name), chr19 = paste0("chr", sub("^MT$", "M", CHR)), pos19 = as.numeric(MAPINFO))], by = "key")
E19 <- fread(path.expand(fe_19), colClasses = "character"); setnames(E19, c("id", "CHR", "MAPINFO"))
E19 <- unique(E19[ok_id(id) & grepl("^\\d+$", MAPINFO), .(key = tolower(id), chr19 = paste0("chr", sub("^chr", "", CHR)), pos19 = as.numeric(MAPINFO))], by = "key")
E19 <- E19[chr19 %in% std_chr]
z38 <- function(f) { z <- fread(path.expand(f), colClasses = "character", na.strings = c("NA", ""))[ok_id(probe_id)]
  ids <- unique(z$probe_id)
  p <- unique(z[grepl("^chr", chr_hg38) & grepl("^\\d+$", beg_hg38), .(key = tolower(probe_id), chr38 = chr_hg38, pos38 = as.numeric(beg_hg38) + 1)])
  multi <- unique(p$key[duplicated(p$key)])            # one ID, several places (EPIC v2 replicates)
  list(ids = ids, pos = p[!key %chin% multi], multi = multi) }
Z450 <- z38(f450_38); ZE <- z38(fe_38); ZV2 <- z38(fv2_38); Z27 <- z38(f27_38)

# ---- each outside source against the package, on shared probes ----------------------------------
chk <- function(lab, X, Pk, a, b) { m <- merge(X, Pk, by = "key"); same <- m[[paste0(a, ".x")]] == m[[paste0(a, ".y")]] & m[[paste0(b, ".x")]] == m[[paste0(b, ".y")]]
  off1 <- m[[paste0(a, ".x")]] == m[[paste0(a, ".y")]] & abs(m[[paste0(b, ".x")]] - m[[paste0(b, ".y")]]) == 1
  data.table(source = lab, shared_with_package = big(nrow(m)), identical = sprintf("%.4f", mean(same)), one_base_apart = big(sum(off1)), elsewhere = big(sum(!same & !off1)), .ok = mean(same)) }
cal <- rbind(chk("Illumina 450K manifest, hg19", H19, P19, "chr19", "pos19"), chk("Illumina EPIC v1 coordinates, hg19", E19, P19, "chr19", "pos19"),
             chk("Zhou HM450 manifest, hg38", Z450$pos, P38, "chr38", "pos38"), chk("Zhou EPIC manifest, hg38", ZE$pos, P38, "chr38", "pos38"),
             chk("Zhou EPICv2 manifest, hg38", ZV2$pos, P38, "chr38", "pos38"), chk("Zhou HM27 manifest, hg38", Z27$pos, P38, "chr38", "pos38"))
print(cal[, !".ok"]); stopifnot(all(cal$.ok > 0.99))

# ---- probes, spelling, membership ------------------------------------------------------------------
sp <- unique(data.table(id = c(ZV2$ids, ZE$ids, Z450$ids, ids450_ill, Z27$ids)))[, key := tolower(id)]
sp <- unique(sp, by = "key")
pk <- setdiff(v2_pkg, sp$key); pk <- pk[ok_id(pk)]
X <- rbind(sp, setnames(data.table(id = pk, kk = pk), "kk", "key")); setnames(X, "id", "cpg_id")   # data.table(key =) would set a key
on450 <- unique(tolower(c(Z450$ids, ids450_ill))); on27 <- tolower(Z27$ids); onE <- tolower(ZE$ids); onV2 <- unique(c(v2_pkg, tolower(ZV2$ids)))
X[, array := { a <- cbind(ifelse(key %chin% on27, "HM27", NA), ifelse(key %chin% on450, "450K", NA), ifelse(key %chin% onE, "EPICv1", NA), ifelse(key %chin% onV2, "EPICv2", NA))
              apply(a, 1, function(r) paste(r[!is.na(r)], collapse = ";")) }]
stopifnot(!anyDuplicated(X$key), all(nzchar(X$array)))

# ---- hg19: supplement only -------------------------------------------------------------------------
S19 <- rbind(cbind(H19, source_hg19 = "illumina_450k_manifest"), cbind(E19[!key %chin% H19$key], source_hg19 = "illumina_epic_v1_manifest"))
h19 <- S19[X[, .(key)], on = "key"]; need19 <- !(X$key %chin% P19$key)
X[, `:=`(chr_hg19 = fifelse(need19, h19$chr19, NA_character_), pos_hg19 = as.integer(fifelse(need19, h19$pos19, NA_real_)),
         source_hg19 = fifelse(need19, h19$source_hg19, NA_character_))]

# ---- hg38: one position per probe by priority; ambiguity ---------------------------------------------
src <- list(zhou_hm450_hg38 = Z450$pos, zhou_epic_hg38 = ZE$pos, zhou_epicv2_hg38 = ZV2$pos, zhou_hm27_hg38 = Z27$pos, cpgdirection = P38)
L <- rbindlist(lapply(names(src), function(s) cbind(src[[s]], source = s, rank = match(s, names(src)))))
setorder(L, key, rank)
amb <- L[, .(n = uniqueN(paste(chr38, pos38))), by = key][n > 1L, key]
amb <- union(amb, ZV2$multi)
F38 <- L[, .SD[1L], by = key]
f <- F38[X[, .(key)], on = "key"]; p <- P38[X[, .(key)], on = "key"]
same_as_pkg <- !is.na(f$pos38) & !is.na(p$pos38) & f$chr38 == p$chr38 & f$pos38 == p$pos38
X[, `:=`(chr_hg38 = fifelse(same_as_pkg, NA_character_, f$chr38), pos_hg38 = as.integer(fifelse(same_as_pkg, NA_real_, f$pos38)),
         source_hg38 = fifelse(same_as_pkg, NA_character_, f$source), hg38_ambiguous = fifelse(key %chin% amb, TRUE, NA))]
X[, `:=`(.has19 = key %chin% P19$key | !is.na(pos_hg19), .has38 = !is.na(f$pos38), .f38chr = f$chr38, .f38pos = f$pos38)]

# ---- reports --------------------------------------------------------------------------------------
per <- rbindlist(lapply(c("HM27", "450K", "EPICv1", "EPICv2"), function(a) X[grepl(a, array, fixed = TRUE),
  .(array = a, probe_ids = big(.N), with_hg19 = big(sum(.has19)), with_hg38 = big(sum(.has38)), no_hg38 = big(sum(!.has38)), hg38_ambiguous = big(sum(hg38_ambiguous %in% TRUE)))]))
tot <- X[, .(array = "all", probe_ids = big(.N), with_hg19 = big(sum(.has19)), with_hg38 = big(sum(.has38)), no_hg38 = big(sum(!.has38)), hg38_ambiguous = big(sum(hg38_ambiguous %in% TRUE)))]
print(rbind(per, tot))
kind <- X[, .(probe_ids = .N, no_hg38 = sum(!.has38), on_alt_contig = sum(.has38 & !.f38chr %chin% std_chr)), by = .(type = tolower(substr(cpg_id, 1, 2)))]
src_tab <- X[, .N, by = .(source_hg38 = fifelse(is.na(source_hg38), fifelse(.has38, "cpgdirection (annotation layer; not stored)", "none"), source_hg38))][order(-N)]
chk_lines <- NULL
if (!is.null(chk_f)) {
  T0 <- fread(path.expand(chk_f)); setnames(T0, 1:3, c("cpg", "chr", "pos")); T0[, key := tolower(cpg)]
  m <- merge(T0, X[, .(key, .f38chr, .f38pos)], by = "key", all.x = TRUE)
  n_ok <- sum(m$chr == m$.f38chr & m$pos == m$.f38pos, na.rm = TRUE)
  chk_lines <- c("## Reference positions reproduced", "", sprintf("`%s`: %s probes; identical chromosome and position for %s; different or missing for %s.",
                 basename(chk_f), big(nrow(T0)), big(n_ok), big(nrow(T0) - n_ok)), "")
  message(chk_lines[3]); stopifnot(n_ok == nrow(T0))
}
out <- X[order(key), .(cpg_id, array, chr_hg19, pos_hg19, source_hg19, chr_hg38, pos_hg38, source_hg38, hg38_ambiguous)]
fwrite(out, "inst/extdata/cpg_coordinates.csv.gz", compress = "gzip")
kind[, `:=`(probe_ids = big(probe_ids), no_hg38 = big(no_hg38), on_alt_contig = big(on_alt_contig))]; src_tab[, N := big(N)]
rep <- c("# CpG coordinates and array membership (HM27, 450K, EPIC v1, EPIC v2)", "", sprintf("Run %s. %s probe IDs.", format(Sys.time(), "%Y-%m-%d %H:%M"), big(nrow(out))), "",
         "## Outside sources against the package, on the probes both hold", "", md(cal[, !".ok"]), "",
         "## Coverage by array", "", md(rbind(per, tot)), "",
         "`no_hg38`: no source places the probe on hg38 (the manifests give no mapping). `hg38_ambiguous`: two sources place the probe differently.", "",
         "## By type of probe ID", "", md(kind), "", "## Source of the hg38 position", "", md(src_tab), "",
         sprintf("hg19 positions supplied here (the rest are the package's): %s. hg38 positions stored here: %s.", big(sum(!is.na(out$pos_hg19))), big(sum(!is.na(out$pos_hg38)))), "", chk_lines)
dir.create("tools/reports", showWarnings = FALSE)
writeLines(rep, "tools/reports/build_cpg_coordinates.md")
cat(sprintf("wrote inst/extdata/cpg_coordinates.csv.gz (%s rows)\n", big(nrow(out))))
