#!/usr/bin/env Rscript
# Ruling 3 follow-up (4 Nov 2026): how each evidence source transports across tissue and age.
# For every source rung (blood SMR by rung x replication x tier; brain SMR by tier; measured blood by cohort)
# it reports sign agreement, with a Wilson CI and the majority-class baseline, against every independent
# measured reference the package now has:
#   child blood   measured rows whose cohorts include HELIX      (age shift, same tissue)
#   adult blood   measured blood rows without HELIX              (reference tissue)
#   nasal         measured nasal epithelium (EVA-PR)              (tissue shift, airway)
#   solid         measured solid_tissue                           (tissue shift, mixed organs)
#   adult cortex  brain_measured_eqtms.csv.gz (ROSMAP aged DLPFC) (tissue shift, brain)
# A source is never scored against itself: HELIX rows are dropped from the measured-blood source when the
# reference is child blood, and so on. Read-only; changes nothing in the package.
#
# Usage (package root, after the November data build):
#   Rscript tools/validate_ladder_transport.R [--min-n 200]
# Writes tools/reports/validate_ladder_transport.md and tools/out/ladder_transport.csv.

suppressPackageStartupMessages({ library(data.table); devtools::load_all(".", quiet = TRUE) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1L] }
min_n <- as.integer(get_arg("--min-n", "200"))
wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d; c(c - h, c + h) }
md <- function(d) { d <- as.data.frame(d); c(paste0("| ", paste(names(d), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"), apply(d, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))) }
key <- function(D) D[, target_gene := toupper(trimws(target_gene))][]

# ---- references ---------------------------------------------------------------------------------
M <- key(fread("inst/extdata/measured_eqtms.csv.gz"))
M[, helix := grepl("HELIX", cohorts)]
refs <- list(
  child_blood  = M[tissue == "blood" & helix, .(cpg_id, target_gene, r = direction)],
  adult_blood  = M[tissue == "blood" & !helix, .(cpg_id, target_gene, r = direction)],
  nasal        = M[grepl("nasal|airway", tissue), .(cpg_id, target_gene, r = direction)],
  solid        = M[grepl("solid", tissue), .(cpg_id, target_gene, r = direction)])
bf <- "inst/extdata/brain_measured_eqtms.csv.gz"
if (file.exists(bf)) refs$adult_cortex <- key(fread(bf))[, .(cpg_id, target_gene, r = direction)]
message("references: ", paste(sprintf("%s %s", names(refs), vapply(refs, nrow, 1L)), collapse = "; "))

# ---- sources ------------------------------------------------------------------------------------
S <- key(fread(path.expand("~/cpgd_local/smr_directions.csv.gz")))
if (!"smr_rung" %in% names(S)) stop("cpgd_smr_directions() has no smr_rung column: run after the rulings v2c build")
if (!"replication" %in% names(S)) S[, replication := NA_character_]
src <- list()
src$blood_smr_rung <- S[, .(cpg_id, target_gene, s = direction, level = smr_rung)]
src$blood_smr_rung_repl_tier <- S[, .(cpg_id, target_gene, s = direction, level = paste(smr_rung, replication, smr_tier, sep = "/"))]
BD <- tryCatch(key(as.data.table(cpgd_brain_directions())), error = function(e) NULL)
if (!is.null(BD)) src$brain_smr_tier <- BD[, .(cpg_id, target_gene, s = direction, level = paste0("brain_", tier))]
mb <- M[tissue == "blood", .(cpg_id, target_gene, s = direction, helix, cohorts)]

score <- function(x, source, reference) x[, { b <- max(mean(r == 1), mean(r == -1)); k <- sum(s == r); ci <- wilson(k, .N)
  .(source = source, reference = reference, n = .N, agree = round(k / .N, 3), lo = round(ci[1], 3), hi = round(ci[2], 3), baseline = round(b, 3)) }, by = level]

R <- list()
for (rn in names(refs)) {
  for (sn in names(src)) R[[length(R) + 1]] <- score(merge(src[[sn]], refs[[rn]], by = c("cpg_id", "target_gene")), sn, rn)
  if (rn != "adult_blood") {
    m <- if (rn == "child_blood") mb[helix == FALSE] else mb
    R[[length(R) + 1]] <- score(merge(m[, .(cpg_id, target_gene, s, level = "measured_blood")], refs[[rn]], by = c("cpg_id", "target_gene")), "measured_blood", rn)
  }
}
R <- rbindlist(R)
setcolorder(R, c("source", "level", "reference"))
dir.create("tools/out", showWarnings = FALSE, recursive = TRUE); dir.create("tools/reports", showWarnings = FALSE)
fwrite(R, "tools/out/ladder_transport.csv")

# ---- wide summary and the ruling-3 checks ------------------------------------------------------
cell <- function(a, lo, hi, n) ifelse(n < 30, "-", sprintf("%.3f [%.3f, %.3f] (n=%s)", a, lo, hi, format(n, big.mark = ",")))
R[, txt := cell(agree, lo, hi, n)]
W1 <- dcast(R[source %in% c("blood_smr_rung", "measured_blood", "brain_smr_tier")], source + level ~ reference, value.var = "txt", fill = "-")
lad <- c("replicated", "high", "moderate", "weak", "reported")
ord <- R[source == "blood_smr_rung" & level %in% lad & n >= min_n][, level := factor(level, lad)][order(reference, level)]
order_ok <- ord[, .(rungs_scored = paste(level, collapse = " > "), monotone = all(diff(agree) <= 0.02)), by = reference]
mod <- R[source == "blood_smr_rung" & level == "moderate", .(reference, moderate = agree, n_mod = n)]
mea <- R[source == "measured_blood", .(reference, measured = agree, n_meas = n)]
gap <- merge(mod, mea, by = "reference", all = TRUE)[, gap_measured_minus_moderate := round(measured - moderate, 3)][]
print(W1); print(order_ok); print(gap)

rep <- c("# Ladder transport across tissue and age", "", sprintf("Run %s. Rungs scored only where n >= %d for the order check.", format(Sys.time(), "%Y-%m-%d %H:%M"), min_n), "",
  "## Agreement by source and reference (Wilson 95% CI)", "", md(W1), "",
  "## Does the blood SMR rung order hold in each reference? (monotone allows 0.02 slack)", "", md(order_ok), "",
  "## Ruling 3: measured blood vs SMR moderate, by reference", "", md(gap), "",
  "Reading: ruling 3 (keep smr_moderate below catalogue_consensus / measured evidence) stands if the measured-minus-moderate gap",
  "stays positive and widens away from blood. Full long table, including rung x replication x tier, in tools/out/ladder_transport.csv.", "")
writeLines(rep, "tools/reports/validate_ladder_transport.md")
message("wrote tools/reports/validate_ladder_transport.md and tools/out/ladder_transport.csv")
