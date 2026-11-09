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

## Since cpgdirectionData 0.99.7 every name resolves to one Ensembl gene and a
## gene's names are one candidate. Before, the cohorts' different names of one
## gene (renamed HGNC symbols, clone names, ENSG IDs) were separate genes.
test_that("one gene is one candidate, under one name", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19"), "data layers unavailable")
  same <- list(c("AFDN", "MLLT4"), c("CRHR1-IT1", "LINC02210"), c("DGCR14", "ESS2"),
               c("ADA2", "CECR1"))
  for (t in CPGD_TISSUES) {
    g <- unique(cpgdirection:::.cpgd_lookup(t)$target_gene)
    for (s in same) expect_lte(sum(s %in% g), 1L, label = paste(t, paste(s, collapse = "/")))
  }
  g <- unique(cpgdirection:::.cpgd_lookup("blood")$target_gene)
  expect_true(all(c("AFDN", "CRHR1-IT1", "DGCR14", "CECR1") %in% g))
})

## Genes first seen in the GTEx and Framingham Atlas exports (not hg19) were
## once placed at the wrong end on the minus strand: GTEx gives a CpG as TSS +
## a strand-signed distance, and Framingham's TSS is the gene's leftmost
## coordinate. Two minus-strand genes, one from each export, must sit at their
## Ensembl GRCh37 TSS: the highest coordinate, which is the 5' end on the minus
## strand.
test_that("minus-strand genes from the non-hg19 Atlas exports sit at their 5' end", {
  skip_if_not(cpgd_has_data("lookup_blood_hg19") && cpgd_has_data("cpg_positions_hg19"),
              "data layers unavailable")
  L <- data.table::as.data.table(cpgdirection:::.cpgd_lookup("blood"))
  P <- data.table::as.data.table(cpgdirection:::.cpgd_data("cpg_positions_hg19"))
  tss <- c(ADCY5 = 123168605, FTCD = 47575481)   # chr3 (GTEx), chr21 (Framingham)
  for (g in names(tss)) {
    x <- L[L$target_gene == g, ]
    expect_gt(nrow(x), 100)
    expect_equal(as.numeric(unique(P$pos[match(x$cpg_id, P$cpg_id)] - x$tss_dist)), tss[[g]],
                 label = paste(g, "TSS"))
  }
  x <- L[L$cpg_id == "cg12793116" & L$target_gene == "ADCY5", ]
  expect_equal(as.integer(x$tss_dist), -63243L)
})
