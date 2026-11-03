#!/usr/bin/env Rscript
# Rebuild the SMR layer with outcome-study replication -- v2 (2026-11-03).
#
# v2 supersedes the 03-Nov v1 builder after lead 1 (tools/explore_smr_pairs.R):
# the FinnGen side now takes its sign from SMR's own top-cis-mQTL test where that
# test passes p < 5e-8 (sign(b_SMR)), which (a) adds the pairs SMR tested that had
# no clumped GoDMC instrument (PBMC: 36,518 new pairs, 0.899 vs measured blood,
# n = 782) and (b) is the admission rule for every FinnGen-only pair (0.889 on the
# instrument-table pairs, n = 675). HEIDI is shipped and conditions nothing; the
# v1 "S1 or S2 + HEIDI pass, provisional" rung is withdrawn.
#
# Inputs
#   * the shipped smr_directions table (GoDMC x eQTLGen; 413,982 pairs) via the
#     package accessor;
#   * smr/bychr/smrd_<ctx>_chr<c>.smr -- SMR v1.03 two-BESD results per context and
#     chromosome (exposure GoDMC dense cis mQTL, outcome FinnGen L1 cis-eQTL);
#   * smr/smrd_blood_<ctx>.directions.csv.gz -- the clumped-instrument Wald tables
#     (tools/finngen_smr_celltype.R), used for pairs SMR did not pass at 5e-8.
#
# Output: ONE table, smr_directions.csv.gz, one row per CpG-gene pair:
#   original columns (cpg_id, target_gene, direction, smr_tier, p_SMR, n_instruments,
#   instrument_agreement, cpg_gene_dist, top_instrument, p_HEIDI, nsnp_HEIDI,
#   heidi_status) unchanged in meaning, plus
#   outcome            eQTLGen | FinnGen_<context>
#   replication        concordant | discordant | eqtlgen_only | finngen_only |
#                      concordant_celltype | discordant_celltype
#   finngen_direction, finngen_tier, finngen_p_SMR, finngen_heidi_status, finngen_contexts
#   finngen_evidence   smr_top_snp | instrument_vote | both | NA   (what the FinnGen sign rests on)
#   smr_rung           replicated | high | moderate | weak | reported
#
# Rung rules (all validated against the measured blood eQTMs, 2-3 Nov 2026):
#   eQTLGen S1/S2, FinnGen PBMC sign agrees        -> replicated  (S1 0.993 n=299; S2 0.958 n=3,060)
#   eQTLGen S1/S2, FinnGen PBMC sign disagrees     -> weak, eQTLGen sign kept (0.767 n=223)
#   eQTLGen S3                                     -> weak
#   eQTLGen row without a FinnGen PBMC pair        -> high / moderate / weak by tier
#   FinnGen PBMC-only with top-mQTL SMR p < 5e-8   -> moderate (0.899 n=782 SMR-tested; 0.889 n=675 instrument pairs)
#   FinnGen PBMC-only without SMR p < 5e-8         -> reported
#   FinnGen cell-type-only (any)                   -> reported (SMR-tested 0.768 n=228; instrument 0.70 n=177)
#
# Usage (from the cpgdirection package root)
#   Rscript tools/build_smr_replication.R \
#       --smr-dir  "$HOME/finngen_multiome/smr/bychr" \
#       --inst-glob "$HOME/finngen_multiome/smr/smrd_blood_*.directions.csv.gz" \
#       --out-dir tools/reports/smr_replication_build [--data-dir <hub_upload dir>]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L]
}
smr_dir   <- path.expand(get_arg("--smr-dir", "~/finngen_multiome/smr/bychr"))
inst_glob <- get_arg("--inst-glob", "~/finngen_multiome/smr/smrd_blood_*.directions.csv.gz")
out_dir   <- get_arg("--out-dir", "tools/reports/smr_replication_build")
data_dir  <- get_arg("--data-dir")
P_ADMIT   <- 5e-8
if (!is.null(data_dir)) options(cpgdirection.data_dir = data_dir, cpgdirectionData.local = data_dir)

if (file.exists("DESCRIPTION") && grepl("^Package: cpgdirection", readLines("DESCRIPTION", 1))) {
  ok <- tryCatch({ suppressMessages(devtools::load_all(".", quiet = TRUE)); TRUE }, error = function(e) FALSE)
  if (!ok) for (f in list.files("R", full.names = TRUE)) sys.source(f, envir = globalenv())
} else if (requireNamespace("cpgdirection", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cpgdirection))
} else stop("run from the cpgdirection package root or install cpgdirection")

wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(c - h, c + h)
}
fmt_ci <- function(k, n) { if (n == 0) return("-"); ci <- wilson(k, n); sprintf("%.3f [%.3f, %.3f]", k / n, ci[1], ci[2]) }
md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  paste(c(paste0("| ", paste(names(d), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
          apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))), collapse = "\n")
}
read_gz <- function(f) fread(cmd = sprintf("gzip -dc %s", shQuote(f)))
gkey <- function(g) gsub("_", "-", toupper(trimws(g)))
big <- function(x) format(x, big.mark = ",")
heidi_of <- function(p) fifelse(is.na(p), "not_tested", fifelse(p >= 0.05, "pass", "fail"))
ctx_order <- c("PBMC", "CD4_T", "CD8_T", "Mono", "NK", "B", "DC", "other_T", "other")

# ---- 1. the shipped eQTLGen layer ---------------------------------------------------
E <- as.data.table(cpgd_smr_directions())
stopifnot(all(c("cpg_id", "target_gene", "direction", "smr_tier", "p_SMR", "n_instruments",
                "instrument_agreement", "cpg_gene_dist", "top_instrument",
                "p_HEIDI", "nsnp_HEIDI", "heidi_status") %in% names(E)))
if ("outcome" %in% names(E)) {
  message("input table already carries 'outcome': rebuilding from its eQTLGen rows only")
  E <- E[outcome == "eQTLGen"]
  drop <- intersect(c("outcome", "replication", "finngen_direction", "finngen_tier", "finngen_p_SMR",
                      "finngen_heidi_status", "finngen_contexts", "finngen_evidence", "smr_rung"), names(E))
  if (length(drop)) E[, (drop) := NULL]
}
E[, target_gene := toupper(trimws(target_gene))]
E[, .gkey := gkey(target_gene)]
stopifnot(!anyDuplicated(E[, .(cpg_id, .gkey)]))
message(sprintf("eQTLGen layer: %s pairs", big(nrow(E))))

# ---- 2a. FinnGen SMR-tested pairs at p < 5e-8, per context --------------------------
T <- rbindlist(lapply(ctx_order, function(ct) {
  files <- Sys.glob(file.path(smr_dir, sprintf("smrd_%s_chr*.smr", ct)))
  if (!length(files)) return(NULL)
  S <- rbindlist(lapply(files, function(f) {
    d <- fread(f, showProgress = FALSE)
    d[!is.na(b_SMR) & b_SMR != 0 & !is.na(p_SMR) & p_SMR < P_ADMIT,
      .(cpg_id = Expo_ID, target_gene = toupper(trimws(Outco_Gene)), ensembl_id = Outco_ID,
        direction = sign(b_SMR), p_SMR, p_HEIDI, nsnp_HEIDI, top_instrument = topSNP,
        cpg_gene_dist = Outco_bp - Expo_bp)]
  }), use.names = TRUE)
  S <- S[nzchar(target_gene)]
  S[, .gkey := gkey(target_gene)]
  setorder(S, cpg_id, .gkey, p_SMR)
  S <- S[!duplicated(S[, .(cpg_id, .gkey)])]
  S[, ctx := ct]
  message(sprintf("  %-8s SMR-tested pairs at p<%g: %s", ct, P_ADMIT, big(nrow(S))))
  S
}), use.names = TRUE)
if (!nrow(T)) stop("no .smr files found under ", smr_dir)
T[, heidi_status := heidi_of(p_HEIDI)]
T[, evidence := "smr_top_snp"]
T[, `:=`(smr_tier = "S2", n_instruments = 1L, instrument_agreement = 1)]

# ---- 2b. FinnGen instrument-vote pairs (the directions tables) ---------------------
I <- rbindlist(lapply(Sys.glob(path.expand(inst_glob)), read_gz), fill = TRUE)
if (!nrow(I)) stop("no instrument tables match ", inst_glob)
I[, target_gene := toupper(trimws(target_gene))]
I <- I[!is.na(direction) & direction != 0 & nzchar(target_gene)]
I[, .gkey := gkey(target_gene)]
I[, ctx := sub("^blood_", "", context)]
I[, heidi_status := heidi_of(p_HEIDI)]
I[, evidence := "instrument_vote"]
I <- I[, .(cpg_id, target_gene, ensembl_id, .gkey, ctx, direction = as.numeric(direction), smr_tier,
           p_SMR, n_instruments, instrument_agreement, cpg_gene_dist, top_instrument,
           p_HEIDI, nsnp_HEIDI, heidi_status, evidence)]

# ---- 2c. one FinnGen record per pair and context: SMR top-SNP where it passed, else vote
TI <- merge(T[, .(cpg_id, .gkey, ctx, smr_dir_ = direction)],
            I[, .(cpg_id, .gkey, ctx, inst_dir_ = direction)], by = c("cpg_id", ".gkey", "ctx"), all = TRUE)
both_key <- TI[!is.na(smr_dir_) & !is.na(inst_dir_), paste(cpg_id, .gkey, ctx)]
F <- rbindlist(list(T[, names(I), with = FALSE], I[!paste(cpg_id, .gkey, ctx) %in% paste(T$cpg_id, T$.gkey, T$ctx)]),
               use.names = TRUE)
F[paste(cpg_id, .gkey, ctx) %in% both_key, evidence := "both"]
sanity <- TI[!is.na(smr_dir_) & !is.na(inst_dir_), .(n = .N, agree = fmt_ci(sum(smr_dir_ == inst_dir_), .N)), by = ctx][order(match(ctx, ctx_order))]
message("\n== sign(b_SMR) vs instrument vote on pairs both cover (p_SMR < 5e-8):"); print(sanity)

# union across contexts: PBMC first, then the cell types in a fixed order
F[, ctx_rank := match(ctx, ctx_order)]; F[is.na(ctx_rank), ctx_rank := length(ctx_order) + 1L]
setorder(F, cpg_id, .gkey, ctx_rank, p_SMR, na.last = TRUE)
ctxs <- F[, .(finngen_contexts = paste(unique(ctx[order(ctx_rank)]), collapse = ";")), by = .(cpg_id, .gkey)]
U <- F[!duplicated(F[, .(cpg_id, .gkey)])]
U <- merge(U, ctxs, by = c("cpg_id", ".gkey"))
U[, finngen_source := fifelse(ctx == "PBMC", "PBMC", "celltype")]
U[, admissible := finngen_source == "PBMC" & !is.na(p_SMR) & p_SMR < P_ADMIT & evidence %in% c("smr_top_snp", "both")]
message(sprintf("\nFinnGen union: %s pairs (%s PBMC, %s cell-type only); %s PBMC pairs admissible (top-mQTL SMR p < %g)",
                big(nrow(U)), big(sum(U$finngen_source == "PBMC")), big(sum(U$finngen_source != "PBMC")),
                big(sum(U$admissible)), P_ADMIT))

# ---- 3. annotate the eQTLGen rows ---------------------------------------------------
Uf <- U[, .(cpg_id, .gkey, finngen_direction = direction, finngen_tier = smr_tier, finngen_p_SMR = p_SMR,
            finngen_heidi_status = heidi_status, finngen_contexts, finngen_evidence = evidence, finngen_source)]
E <- merge(E, Uf, by = c("cpg_id", ".gkey"), all.x = TRUE)
E[, outcome := "eQTLGen"]
E[, replication := fcase(
  is.na(finngen_direction), "eqtlgen_only",
  finngen_source == "PBMC" & direction == finngen_direction, "concordant",
  finngen_source == "PBMC" & direction != finngen_direction, "discordant",
  direction == finngen_direction, "concordant_celltype",
  default = "discordant_celltype")]
E[, smr_rung := fcase(
  replication == "concordant" & smr_tier %in% c("S1", "S2"), "replicated",
  replication == "discordant" & smr_tier %in% c("S1", "S2"), "weak",
  smr_tier == "S1", "high",
  smr_tier == "S2", "moderate",
  default = "weak")]

# ---- 4. FinnGen-only rows (pairs the eQTLGen layer lacks) ---------------------------
ek <- paste(E$cpg_id, E$.gkey)
N <- U[!paste(cpg_id, .gkey) %in% ek]
N <- N[, .(cpg_id, target_gene, .gkey, direction, smr_tier, p_SMR, n_instruments, instrument_agreement,
           cpg_gene_dist, top_instrument, p_HEIDI, nsnp_HEIDI, heidi_status,
           outcome = paste0("FinnGen_", ctx), replication = "finngen_only",
           finngen_direction = direction, finngen_tier = smr_tier, finngen_p_SMR = p_SMR,
           finngen_heidi_status = heidi_status, finngen_contexts, finngen_evidence = evidence,
           smr_rung = fifelse(admissible, "moderate", "reported"))]

# ---- 5. assemble ----------------------------------------------------------------------
keep <- c("cpg_id", "target_gene", "direction", "smr_tier", "p_SMR", "n_instruments",
          "instrument_agreement", "cpg_gene_dist", "top_instrument", "p_HEIDI", "nsnp_HEIDI",
          "heidi_status", "outcome", "replication", "finngen_direction", "finngen_tier",
          "finngen_p_SMR", "finngen_heidi_status", "finngen_contexts", "finngen_evidence", "smr_rung")
E[, finngen_source := NULL]
OUT <- rbindlist(list(E[, c(keep, ".gkey"), with = FALSE], N[, c(keep, ".gkey"), with = FALSE]), use.names = TRUE)
stopifnot(!anyDuplicated(OUT[, .(cpg_id, .gkey)]))
OUT[, .gkey := NULL]
setorder(OUT, cpg_id, target_gene)

comp <- OUT[, .N, by = .(outcome = fifelse(outcome == "eQTLGen", "eQTLGen", sub("_.*$", "", outcome)),
                         replication, finngen_evidence, smr_rung)][order(outcome, replication, finngen_evidence, smr_rung)]
rung_tot <- OUT[, .(rows = .N, cpgs = uniqueN(cpg_id)), by = smr_rung][order(-rows)]
message("\n== rows by outcome / replication / evidence / rung"); print(comp)
message("\n== rows by rung"); print(rung_tot)

# ---- 6. validation of every rung against the measured blood eQTMs -------------------
M <- as.data.table(cpgd_measured_eqtms())[tissue == "blood"]
M <- unique(M[, list(target_gene = toupper(trimws(unlist(strsplit(target_gene, ";", fixed = TRUE))))),
              by = .(cpg_id, ref_dir = direction)])
M[, .gkey := gkey(target_gene)]
V <- merge(OUT[, .(cpg_id, .gkey = gkey(target_gene), direction, smr_tier, heidi_status, outcome,
                   replication, finngen_evidence, smr_rung)],
           unique(M[, .(cpg_id, .gkey, ref_dir)]), by = c("cpg_id", ".gkey"))
val_rung <- V[, .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N)), by = smr_rung][order(smr_rung)]
val <- V[, .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N),
             baseline = sprintf("%.3f", max(mean(ref_dir == 1), mean(ref_dir == -1)))),
         by = .(smr_rung, replication, smr_tier)][order(smr_rung, replication, smr_tier)]
val_fo <- V[replication == "finngen_only", .(n = .N, agree = fmt_ci(sum(direction == ref_dir), .N)),
            by = .(outcome = sub("_.*$", "", outcome), finngen_evidence, heidi_status, smr_rung)][order(outcome, finngen_evidence, heidi_status)]
message("\n== validation vs measured blood, by rung"); print(val_rung)
message("\n== validation vs measured blood, by rung x replication x tier"); print(val)
message("\n== FinnGen-only rows by evidence x HEIDI (HEIDI is reported, never applied)"); print(val_fo)

# ---- 7. write -------------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "smr_directions.csv.gz")
fwrite(OUT, out_file, compress = "gzip")
report <- c("# SMR layer with outcome-study replication -- build report (v2)", "",
            sprintf("Built %s from %s eQTLGen pairs, %s FinnGen SMR-tested pairs (p < %g) and %s FinnGen instrument-vote pairs; union %s FinnGen pairs.",
                    format(Sys.time(), "%Y-%m-%d %H:%M"), big(nrow(E)), big(nrow(T)), P_ADMIT, big(nrow(I)), big(nrow(U))),
            sprintf("Output: `%s`, %s rows.", out_file, big(nrow(OUT))), "",
            "## sign(b_SMR) vs instrument vote, pairs both cover", "", md_table(sanity), "",
            "## Rows by outcome, replication, FinnGen evidence and rung", "", md_table(comp), "",
            "## Rows by rung", "", md_table(rung_tot), "",
            "## Validation against measured blood eQTMs (Wilson 95% CI)", "",
            "### By rung", "", md_table(val_rung), "",
            "### By rung x replication x tier", "", md_table(val), "",
            "### FinnGen-only rows by evidence x HEIDI status", "", md_table(val_fo), "",
            "Rung rules are in the header of tools/build_smr_replication.R; the accuracy strings in",
            "R/utils.R (CPGD_SMR_ACCURACY) must match the 'By rung' table above.", "")
writeLines(report, file.path(out_dir, "build_smr_replication.md"))
message("\nwrote ", out_file, " and ", file.path(out_dir, "build_smr_replication.md"))
