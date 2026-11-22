#!/usr/bin/env Rscript
## The catalogue-model figures of cpgd_accuracy_table() and of the
## `expected_accuracy` / `best_expected_accuracy` columns: the agreement of the
## shipped models' calls with measured cis-eQTMs from cohorts outside the
## models' training, with Wilson 95% intervals.
##
## Yardsticks (one sign per CpG-gene pair: the majority of the yardstick's
## cohorts, ties dropped):
##   adult blood  MESA monocytes and T cells (Wang et al. 2021, Table S2), not a
##                training cohort of the blood model.
##   child blood  HELIX cell-adjusted eQTMs (Ruiz-Arenas et al. 2022, Dryad),
##                restricted to pairs the blood lookup does not hold as
##                DIRECT_eQTM.
##   nasal        EVA-PR nasal epithelium (eQTM Atlas exports), restricted to
##                pairs the nasal lookup does not hold as DIRECT_eQTM. The nasal
##                model was trained on EVA-PR, so this figure is same-study.
##   solid        FUSION skeletal muscle and GTEx solid tissues (eQTM Atlas
##                exports), each scored in turn and pooled; a pair measured in
##                both counts once (FUSION).
##
## Sources scored (the package's call rule: tier A confidence >= 0.30, tier B
## >= 0.50; PREDICTED rows only, as the ladders form them):
##   catalogue_model_<tissue>_<A|B>   one tissue model's call
##   catalogue_consensus_<A|B>        2+ tissue models calling and agreeing
##                                    (worst tier)
##   catalogue_single_<tissue>_<A|B>  exactly one tissue model calling
##
## Inputs (all are open data):
##   --atlas-long  the per-cohort Atlas table written by build_measured_eqtms.R
##                 (columns cohort, tissue, cpg_id, target_gene, direction,
##                 tss_dist)
##   --helix       HELIX pairs (columns cpg_id, target_gene, direction)
##   --mesa        MESA pairs (columns tissue, cohorts, cpg_id, target_gene,
##                 direction); rows with tissue == "blood" are used
##   --data-dir    optional; the folder named by options(cpgdirection.data_dir)
##   --out         CSV to write (default: a file in tempdir())
##
## Usage:
##   Rscript validate_catalogue_accuracy.R --atlas-long <f> --helix <f>
##     --mesa <f> [--data-dir <dir>] [--out <csv>]

suppressPackageStartupMessages({
  library(data.table)
  library(cpgdirection)
})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1L]
}
data_dir <- get_arg("--data-dir")
if (!is.null(data_dir)) options(cpgdirection.data_dir = path.expand(data_dir))
atlas_long <- get_arg("--atlas-long")
helix_file <- get_arg("--helix")
mesa_file <- get_arg("--mesa")
out_csv <- get_arg("--out", tempfile("catalogue_accuracy_", fileext = ".csv"))
stopifnot(!is.null(atlas_long), !is.null(helix_file), !is.null(mesa_file),
          file.exists(atlas_long), file.exists(helix_file), file.exists(mesa_file))

gk <- function(x) gsub("_", "-", toupper(trimws(x)))
donor_group <- function(cohort) {
  cohort <- trimws(as.character(cohort))
  fcase(startsWith(cohort, "GTEx"), "GTEx",
        startsWith(cohort, "EVA-PR"), "EVAPR",
        grepl("Framingham", cohort, fixed = TRUE), "FHS",
        startsWith(cohort, "HELIX"), "HELIX",
        startsWith(cohort, "MESA"), "MESA",
        default = cohort)
}
wilson <- function(k, n, z = stats::qnorm(0.975)) {
  p <- k / n
  d <- 1 + z^2 / n
  c0 <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(p = p, lo = pmax(0, c0 - h), hi = pmin(1, c0 + h))
}

## measured signs per cohort, with donor groups
A <- fread(atlas_long)
A <- A[!is.na(direction) & direction != 0,
       .(cohort, tissue, cpg_id, g = gk(target_gene),
         d = sign(as.numeric(direction)), grp = donor_group(cohort))]
H <- fread(helix_file)[, .(cohort = "HELIX_2022_Dryad_adj", tissue = "blood",
                           cpg_id, g = gk(target_gene),
                           d = sign(as.numeric(direction)), grp = "HELIX")]
W <- fread(mesa_file)[tissue == "blood"]
W <- W[, .(cohort = cohorts, tissue, cpg_id, g = gk(target_gene),
           d = sign(as.numeric(direction)),
           grp = donor_group(sub("[|].*$", "", cohorts)))]
C <- unique(rbindlist(list(A, H, W), use.names = TRUE)[!is.na(d) & d != 0],
            by = c("cohort", "cpg_id", "g"))
C[, pk := paste(cpg_id, g)]
yard <- function(X) X[, .(y = {
  s <- sum(d)
  if (s > 0) 1 else if (s < 0) -1 else NA_real_
}), by = pk][!is.na(y)]

## the shipped models' calls
tissues <- c("blood", "nasal_epithelium", "solid_tissue")
short <- c(blood = "blood", nasal_epithelium = "nasal", solid_tissue = "solid")
L <- lapply(setNames(tissues, tissues), function(t)
  as.data.table(cpgdirection:::.cpgd_lookup(t)))
direct <- lapply(L, function(x)
  unique(paste(x$cpg_id[x$status == "DIRECT_eQTM"],
               gk(x$target_gene[x$status == "DIRECT_eQTM"]))))
calls <- rbindlist(lapply(tissues, function(t) {
  x <- L[[t]]
  x <- data.table(pk = paste(x$cpg_id, gk(x$target_gene)),
                  d = as.numeric(x$direction), cf = as.numeric(x$confidence),
                  ti = as.character(x$evidence_tier), st = as.character(x$status))
  setorderv(x, c("pk", "cf"), order = c(1L, -1L))
  x <- x[, .SD[1L], by = pk]
  x[, call := fcase(is.na(d), NA_real_,
                    ti == "A" & !is.na(cf) & cf >= 0.30, d,
                    ti == "B" & !is.na(cf) & cf >= 0.50, d,
                    default = NA_real_)]
  x[!is.na(call) & st == "PREDICTED", .(pk, call, ti, tissue = t)]
}))
MOD <- calls[, .(pk, d = call,
                 source = paste0("catalogue_model_", short[tissue], "_", ti))]
CAT <- calls[, .(n_call = .N, n_pos = sum(call == 1), n_neg = sum(call == -1),
                 worst = max(ti), t1 = tissue[1], ti1 = ti[1], d1 = call[1]),
             by = pk]
CAT[, source := fcase(
  n_call >= 2L & (n_pos == n_call | n_neg == n_call),
  paste0("catalogue_consensus_", worst),
  n_call == 1L, paste0("catalogue_single_", short[t1], "_", ti1),
  default = NA_character_)]
CAT[, d := fifelse(n_call == 1L, d1, fifelse(n_pos > n_neg, 1, -1))]
SRC <- rbindlist(list(MOD, CAT[!is.na(source), .(pk, d, source)]))

## score on each yardstick
folds <- list(
  list(name = "MESA adult blood", y = yard(C[grp == "MESA" & tissue == "blood"]),
       drop = NULL),
  list(name = "HELIX child blood (non-catalogue pairs)", y = yard(C[grp == "HELIX"]),
       drop = "blood"),
  list(name = "EVA-PR nasal (non-catalogue pairs)",
       y = yard(C[grp == "EVAPR" & tissue == "nasal_epithelium"]),
       drop = "nasal_epithelium"),
  list(name = "FUSION + GTEx solid", y = {
    fu <- yard(C[grp == "FUSION"])
    gt <- yard(C[grp == "GTEx" & tissue == "solid_tissue"])
    rbind(fu, gt[!(pk %chin% fu$pk)])
  }, drop = NULL))
res <- rbindlist(lapply(folds, function(f) {
  Y <- f$y
  if (!is.null(f$drop)) Y <- Y[!(pk %chin% direct[[f$drop]])]
  m <- merge(SRC, Y, by = "pk")
  m[, .(n = .N, k = sum(d == y)), by = source][, yardstick := f$name][]
}))
w <- wilson(res$k, res$n)
res[, `:=`(agree = round(w$p, 3), lo = round(w$lo, 3), hi = round(w$hi, 3))]
setcolorder(res, c("yardstick", "source", "n", "k", "agree", "lo", "hi"))
setorderv(res, c("yardstick", "source"))
fwrite(res, out_csv)
print(res)
message("written: ", out_csv)
