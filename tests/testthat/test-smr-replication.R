# The SMR layer with outcome-study replication (1.1.0, v2 build). Every test
# skips on a table built before the rebuild, so an older cached cpgdirectionData
# resource degrades to the old behaviour rather than failing the suite.

.smr_rep_table <- function() {
  S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
  if (is.null(S) || !nrow(S) || !all(c("smr_rung", "replication", "outcome") %in% names(S))) return(NULL)
  S
}

test_that("the replication columns are well-formed", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  expect_setequal(setdiff(unique(S$smr_rung), NA),
                  intersect(c("replicated", "high", "moderate", "weak", "reported"), unique(S$smr_rung)))
  expect_true(all(S$replication %in% c("concordant", "discordant", "eqtlgen_only", "finngen_only",
                                       "concordant_celltype", "discordant_celltype")))
  # one row per pair, whatever the outcome study
  expect_false(anyDuplicated(S[, c("cpg_id", "target_gene")]) > 0)
  # the rung rules, as the builder states them
  expect_true(all(S$smr_tier[S$smr_rung == "replicated"] %in% c("S1", "S2")))
  expect_true(all(S$replication[S$smr_rung == "replicated"] == "concordant"))
  expect_true(all(S$smr_rung[S$replication == "discordant" & S$smr_tier == "S1"] == "weak"))
  expect_true(all(S$smr_rung[S$replication == "discordant" & S$smr_tier %in% c("S2", "S3")] == "reported"))
  expect_true(all(S$smr_rung[S$replication == "concordant" & S$smr_tier == "S3"] == "moderate"))
  fo <- S[S$replication == "finngen_only", ]
  expect_true(all(fo$outcome != "eQTLGen"))
  # cell-type-only rows are reported, never admitted
  expect_true(all(fo$smr_rung[fo$outcome != "FinnGen_PBMC"] == "reported"))
  # the eQTLGen sign is kept on discordant rows
  d <- S[S$replication == "discordant", ]
  expect_true(all(d$direction != d$finngen_direction))
  expect_true(all(d$outcome == "eQTLGen"))
})

test_that("FinnGen PBMC-only admission is the top-cis-mQTL SMR threshold, nothing else", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  skip_if(!"finngen_evidence" %in% names(S), "table built before the v2 rule")
  pb <- S[S$replication == "finngen_only" & S$outcome == "FinnGen_PBMC", ]
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
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  rep1 <- S[S$smr_rung == "replicated", ][1:50, ]
  rep1 <- rep1[!is.na(rep1$cpg_id), ]
  skip_if(!nrow(rep1), "no replicated rows")
  p <- cpg_gene_pairs(rep1$cpg_id, genes = rep1$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  # whatever beat a replicated row is ranked above its stratum in blood
  o <- cpgdirection:::.cpgd_rung_order("blood")
  sk <- paste0("smr_", cpgdirection:::.cpgd_smr_key("smr_replicated", p$smr_tier,
                                                     p$smr_outcome, p$smr_replication))
  beat <- p$best_rung != sk & !is.na(p$best_direction)
  expect_true(all(match(p$best_rung[beat], o) < match(sk[beat], o)))
  expect_true(any(p$best_evidence == "smr_replicated"))
  expect_true(all(p$smr_rung[p$best_evidence == "smr_replicated"] == "replicated"))

  rp <- S[S$smr_rung == "reported", ][1:50, ]
  rp <- rp[!is.na(rp$cpg_id), ]
  skip_if(!nrow(rp), "no reported rows")
  q <- cpg_gene_pairs(rp$cpg_id, genes = rp$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  # the direction is visible in the smr_* columns but never promoted
  expect_true(all(!is.na(q$smr_direction)))
  expect_false(any(q$best_evidence %in% c("smr_replicated", "smr_high", "smr_moderate", "smr_weak")))
})

test_that("FinnGen PBMC-only admitted rows reach smr_moderate with the FinnGen accuracy string", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  m <- S[S$replication == "finngen_only" & S$smr_rung == "moderate", ][1:50, ]
  m <- m[!is.na(m$cpg_id), ]
  skip_if(!nrow(m), "no admitted FinnGen-only rows")
  p <- cpg_gene_pairs(m$cpg_id, genes = m$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  w <- p[p$best_evidence == "smr_moderate", ]
  skip_if(!nrow(w), "every admitted row was outranked by a stronger layer")
  expect_true(all(w$smr_outcome == "FinnGen_PBMC"))
  expect_true(all(grepl("FinnGen PBMC outcome only", w$best_expected_accuracy, fixed = TRUE)))
  expect_false(any(grepl("PROVISIONAL|HEIDI", w$best_expected_accuracy)))
})

test_that("discordant rows carry the eQTLGen sign at the weak rung with the discordant accuracy string", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  d <- S[S$replication == "discordant" & S$smr_tier == "S1", ][1:50, ]
  d <- d[!is.na(d$cpg_id), ]
  skip_if(!nrow(d), "no discordant rows")
  p <- cpg_gene_pairs(d$cpg_id, genes = d$target_gene, gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, probe_qc = "ignore")
  w <- p[p$best_evidence == "smr_weak", ]
  skip_if(!nrow(w), "every discordant row was outranked by a stronger layer")
  expect_true(all(w$smr_replication == "discordant"))
  expect_true(all(grepl("disagree", w$best_expected_accuracy, fixed = TRUE)))
  m <- merge(w, d, by.x = c("cpg_id", "target_gene"), by.y = c("cpg_id", "target_gene"))
  expect_true(all(m$best_direction == m$direction))
})

test_that("the CpG-level interface surfaces the rung columns", {
  S <- .smr_rep_table(); skip_if(is.null(S), "SMR table without replication columns")
  r <- cpg_expression_direction(S$cpg_id[S$smr_rung == "replicated"][1], verbose = FALSE)
  expect_true(all(c("smr_rung", "smr_replication", "smr_outcome") %in% names(r)))
  expect_true("smr_replicated" %in% CPGD_EVIDENCE)
})
