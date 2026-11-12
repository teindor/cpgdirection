# The SMR layer with outcome-study replication (1.1.0: eQTLGen, FinnGen and
# INTERVAL; rung rule D; per-pair accuracy). Every test skips on a table built
# before the INTERVAL rebuild, so an older cached cpgdirectionData resource does
# not fail the suite.

.smr_rep_table <- function() {
  S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
  if (is.null(S) || !nrow(S) ||
      !all(c("smr_rung", "replication", "outcome", "interval_status", "smr_accuracy") %in% names(S))) return(NULL)
  S
}

test_that("the replication columns are well-formed and follow the rung rule", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  expect_setequal(unique(S$smr_rung), c("replicated", "high", "moderate", "weak", "reported"))
  expect_true(all(S$replication %in% c("concordant", "discordant", "eqtlgen_only", "finngen_only",
                                       "concordant_celltype", "discordant_celltype", "interval_only")))
  expect_true(all(S$interval_status %in% c("agrees", "disagrees", "absent", "interval_only")))
  expect_true(all(S$heidi_status %in% c("pass", "fail", "not_tested")))
  # one row per pair, whatever the outcome study
  expect_false(anyDuplicated(S[, c("cpg_id", "target_gene")]) > 0)
  # every pair carries a held-out accuracy
  expect_true(all(!is.na(S$smr_accuracy) & S$smr_accuracy >= 0 & S$smr_accuracy <= 1))
  adm <- S$interval_evidence %in% c("both", "smr_top_snp")
  rp <- S$smr_rung == "replicated"
  # replicated: tier S1 or S2, two bulk outcomes with the same sign, none dissenting
  expect_true(all(S$smr_tier[rp] %in% c("S1", "S2")))
  expect_false(any(S$interval_status[rp] == "disagrees"))
  expect_true(all(S$replication[rp] == "concordant" | (S$interval_status[rp] == "agrees" & adm[rp])))
  expect_true(all(S$outcome[rp] %in% c("eQTLGen", "FinnGen_PBMC")))
  # only an admissible INTERVAL record promotes: a vote-only or conflicting record never makes a pair replicated
  expect_true(all(S$replication[rp & !adm] == "concordant"))
  # any INTERVAL dissent demotes to weak (two against one) or reported
  dis <- S$interval_status == "disagrees"
  expect_true(all(S$smr_rung[dis] %in% c("weak", "reported")))
  expect_true(all(S$smr_rung[dis & S$replication == "concordant"] == "weak"))
  # a pair seen in INTERVAL alone: by its tier when the record is admissible, weak otherwise
  io <- S$interval_status == "interval_only"
  expect_true(all(S$outcome[io] == "INTERVAL" & S$replication[io] == "interval_only"))
  expect_true(all(S$smr_rung[io & adm & S$smr_tier == "S1"] == "high"))
  expect_true(all(S$smr_rung[io & adm & S$smr_tier == "S2"] == "moderate"))
  expect_true(all(S$smr_rung[io & !adm] == "weak"))
  # where INTERVAL is absent the eQTLGen-FinnGen rules stand as before
  ab <- S$interval_status == "absent"
  expect_true(all(S$smr_rung[ab & S$replication == "discordant" & S$smr_tier == "S1"] == "weak"))
  expect_true(all(S$smr_rung[ab & S$replication == "discordant" & S$smr_tier %in% c("S2", "S3")] == "reported"))
  expect_true(all(S$smr_rung[ab & S$replication == "concordant" & S$smr_tier == "S3"] == "moderate"))
  fo <- S[S$replication == "finngen_only", ]
  expect_true(all(startsWith(fo$outcome, "FinnGen")))
  # cell-type-only rows are reported unless an admissible INTERVAL record agrees
  ct <- fo$outcome != "FinnGen_PBMC"
  expect_true(all(fo$smr_rung[ct & fo$interval_status == "absent"] == "reported"))
  expect_true(all(fo$smr_rung[ct] %in% c("reported", "moderate")))
  # the eQTLGen sign is kept on discordant rows
  d <- S[S$replication == "discordant", ]
  expect_true(all(d$direction != d$finngen_direction))
  expect_true(all(d$outcome == "eQTLGen"))
})

test_that("FinnGen PBMC-only admission is the top-cis-mQTL SMR threshold where INTERVAL is silent", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  pb <- S[S$replication == "finngen_only" & S$outcome == "FinnGen_PBMC" & S$interval_status == "absent", ]
  skip_if(!nrow(pb), "no FinnGen PBMC-only rows")
  adm <- !is.na(pb$p_SMR) & pb$p_SMR < 5e-8 & pb$finngen_evidence %in% c("smr_top_snp", "both")
  expect_true(all(pb$smr_rung[adm] == "moderate"))
  expect_true(all(pb$smr_rung[!adm] == "reported"))
  expect_true(any(adm))
  # HEIDI conditions nothing: admitted rows span pass, fail and not_tested, and
  # a HEIDI-failing pair is admitted whenever its SMR p clears the threshold
  expect_true(all(c("pass", "fail") %in% pb$heidi_status[adm]))
  expect_true(all(pb$smr_rung[adm & pb$heidi_status == "fail"] == "moderate"))
  # the SMR-tested sign is the sign of b_SMR at the top cis-mQTL (S2, one instrument)
  expect_true(all(pb$smr_tier[pb$finngen_evidence == "smr_top_snp"] == "S2"))
  expect_true(all(pb$n_instruments[pb$finngen_evidence == "smr_top_snp"] == 1L))
})

test_that("replicated rows lose only to sources ranked above them; reported rows never answer", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  rep1 <- S[S$smr_rung == "replicated", ][1:50, ]
  rep1 <- rep1[!is.na(rep1$cpg_id), ]
  skip_if(!nrow(rep1), "no replicated rows")
  p <- cpg_gene_pairs(rep1$cpg_id, genes = rep1$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  # the row's accuracy key is the band of its per-pair accuracy; whatever beat
  # it is ranked above that band in blood
  o <- cpgdirection:::.cpgd_rung_order("blood")
  band <- cpgdirection:::.cpgd_smr_key("smr_replicated", p$smr_accuracy)
  expect_identical(band, cpgdirection:::.cpgd_smr_band(p$smr_accuracy))
  sk <- paste0("smr_", band)
  beat <- !is.na(band) & p$best_rung != sk & !is.na(p$best_direction)
  expect_true(all(match(p$best_rung[beat], o) < match(sk[beat], o)))
  expect_true(any(p$best_evidence == "smr_replicated"))
  w <- p$best_evidence == "smr_replicated"
  expect_true(all(p$smr_rung[w] == "replicated"))
  expect_true(all(p$best_rung[w] == sk[w]))
  expect_true(all(grepl("per-pair accuracy", p$best_expected_accuracy[w], fixed = TRUE)))
  # a pair below the accuracy floor supplies nothing, whatever its rung
  expect_false(any(w & p$smr_accuracy < cpgdirection:::CPGD_SMR_FLOOR))

  rp <- S[S$smr_rung == "reported", ][1:50, ]
  rp <- rp[!is.na(rp$cpg_id), ]
  skip_if(!nrow(rp), "no reported rows")
  q <- cpg_gene_pairs(rp$cpg_id, genes = rp$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  # the direction is visible in the smr_* columns but never promoted
  expect_true(all(!is.na(q$smr_direction)))
  expect_false(any(q$best_evidence %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak")))
})

test_that("FinnGen-only admitted rows reach smr_moderate and carry their band's figure", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  m <- S[S$replication == "finngen_only" & S$smr_rung == "moderate", ][1:50, ]
  m <- m[!is.na(m$cpg_id), ]
  skip_if(!nrow(m), "no admitted FinnGen-only rows")
  p <- cpg_gene_pairs(m$cpg_id, genes = m$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  w <- p[p$best_evidence == "smr_moderate", ]
  skip_if(!nrow(w), "every admitted row was outranked by a stronger layer")
  expect_true(all(startsWith(w$smr_outcome, "FinnGen")))
  expect_true(all(startsWith(w$best_rung, "smr_acc_")))
  expect_true(all(w$smr_accuracy >= cpgdirection:::CPGD_SMR_FLOOR))
  expect_true(all(grepl("per-pair accuracy", w$best_expected_accuracy, fixed = TRUE)))
  expect_false(any(grepl("PROVISIONAL|HEIDI", w$best_expected_accuracy)))
})

test_that("discordant rows carry the eQTLGen sign at the weak rung", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  d <- S[S$replication == "discordant" & S$smr_tier == "S1", ][1:50, ]
  d <- d[!is.na(d$cpg_id), ]
  skip_if(!nrow(d), "no discordant rows")
  p <- cpg_gene_pairs(d$cpg_id, genes = d$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  w <- p[p$best_evidence == "smr_weak", ]
  skip_if(!nrow(w), "every discordant row was outranked by a stronger layer")
  expect_true(all(w$smr_replication == "discordant"))
  expect_true(all(startsWith(w$best_rung, "smr_acc_")))
  m <- merge(w, d, by.x = c("cpg_id", "target_gene"), by.y = c("cpg_id", "target_gene"))
  expect_true(all(m$best_direction == m$direction))
})

test_that("the CpG-level interface surfaces the rung columns", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without the INTERVAL columns")
  r <- cpg_expression_direction(S$cpg_id[S$smr_rung == "replicated"][1], verbose = FALSE)
  expect_true(all(c("smr_rung", "smr_replication", "smr_outcome", "smr_accuracy") %in% names(r)))
  expect_true("smr_replicated" %in% CPGD_EVIDENCE)
})
