## The catalogue lookups once paired a third of their blood and solid rows
## with a gene on another chromosome (fixed in cpgdirectionData 0.99.5). Every
## row's gene must sit on its CpG's chromosome, by the EPIC v2 annotation.

test_that("every lookup row pairs a CpG with a gene on its own chromosome", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19") &&
              cpgd_has_data("epicv2_probe_gene_annotation"), "data layers unavailable")
  A <- data.table::as.data.table(cpgdirection:::.cpgd_data("epicv2_probe_gene_annotation"))
  A <- A[!is.na(A$chr) & nzchar(A$chr), ]
  A$g <- toupper(A$target_gene)
  N <- A[, list(n = .N), by = c("g", "chr")]
  N <- N[order(N$g, -N$n), ]
  top <- N[!duplicated(N$g), ]
  tot <- N[, list(tot = sum(get("n"))), by = "g"]
  top <- merge(top, tot, by = "g")
  top <- top[top$n / top$tot >= 0.9, ]
  P <- data.table::as.data.table(cpgdirection:::.cpgd_data("cpg_positions_hg19"))
  for (t in CPGD_TISSUES) {
    L <- data.table::as.data.table(cpgdirection:::.cpgd_lookup(t))
    cc <- P$chr[match(L$cpg_id, P$cpg_id)]
    gg <- top$chr[match(toupper(L$target_gene), top$g)]
    ok <- !is.na(cc) & !is.na(gg)
    expect_gt(mean(ok), 0.8, label = paste(t, "share of rows checkable"))
    expect_identical(sum(cc[ok] != gg[ok]), 0L, label = paste(t, "rows on another chromosome"))
    expect_true(all(abs(L$tss_dist) <= 250000), label = paste(t, "|tss_dist| <= 250 kb"))
  }
})
