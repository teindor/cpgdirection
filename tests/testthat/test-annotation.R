# cpgd_cpg_annotation(): coordinates in hg19, hg38 and CHM13, array membership,
# genes and probe QC for every probe of the four arrays (1.1.0).

test_that("the bundled coordinate table covers the four arrays and never repeats a probe", {
  X <- cpgdirection:::.cpgd_data("cpg_coordinates")
  expect_identical(names(X)[1:9], c("cpg_id", "array", "chr_hg19", "pos_hg19", "source_hg19", "chr_hg38",
                                    "pos_hg38", "source_hg38", "hg38_ambiguous"))
  expect_false(anyDuplicated(tolower(X$cpg_id)) > 0)
  expect_gt(nrow(X), 1080000L)
  expect_true(all(unlist(strsplit(unique(X$array), ";", fixed = TRUE)) %in% c("HM27", "450K", "EPICv1", "EPICv2")))
  expect_gt(sum(grepl("450K", X$array, fixed = TRUE)), 485000L)
  expect_gt(sum(grepl("EPICv1", X$array, fixed = TRUE)), 865000L)
  expect_gt(sum(grepl("EPICv2", X$array, fixed = TRUE)), 930000L)
})

test_that("hg19 is the package's own table wherever it holds the probe", {
  skip_if_not(cpgd_has_data("cpg_positions_hg19"), "position table unavailable")
  P <- cpgd_cpg_positions()
  set.seed(1); ids <- c("cg00000029", "cg00000109", sample(P$cpg_id, 2000))
  a <- cpgd_cpg_annotation(ids, build = "hg19", genes = FALSE, qc = FALSE)
  expect_identical(a$input, ids)
  expect_true(all(a$known & a$on_epic_v2))
  expect_identical(a$pos, as.integer(P$pos[match(ids, P$cpg_id)]))
  expect_identical(a$chr, as.character(P$chr[match(ids, P$cpg_id)]))
  expect_true(all(a$source == "cpgdirection"))
})

test_that("one build gives chr and pos; several give one pair of columns per build", {
  ids <- c("cg00000029", "cg00000109", "cg00000155")
  a <- cpgd_cpg_annotation(ids)
  expect_true(all(c("chr_hg19", "pos_hg19", "chr_hg38", "pos_hg38", "chr_chm13", "pos_chm13", "source_hg19",
                    "source_hg38", "array", "genes", "qc_exclude") %in% names(a)))
  expect_false(any(c("chr", "pos", "source") %in% names(a)))
  expect_true(all(!is.na(a$pos_hg38)))
  # the builds differ: hg38 is not a copy of hg19
  expect_false(any(a$pos_hg38 == a$pos_hg19))
  expect_true(all(a$chr_hg38 == a$chr_hg19))
  expect_true(grepl("RBL2", a$genes[1], fixed = TRUE))
  expect_identical(a$array[1], "450K;EPICv1;EPICv2")
  h <- cpgd_cpg_annotation(ids, build = "hg38", genes = FALSE, qc = FALSE)
  # (the region flags follow when the package holds them)
  expect_identical(setdiff(names(h), c("in_corsiv", "in_imprinted_or_parental_asm", "in_asm")),
                   c("input", "cpg_id", "known", "chr", "pos", "source", "strand_hg38", "hg38_ambiguous",
                     "array", "on_hm27", "on_450k", "on_epic_v1", "on_epic_v2"))
  expect_identical(h$pos, a$pos_hg38); expect_identical(h$chr, a$chr_hg38)
  expect_identical(attr(h, "build"), "hg38")
  expect_true(is.integer(h$pos))
  k <- cpgd_cpg_annotation(ids, build = "chm13", genes = FALSE, qc = FALSE)
  expect_true(all(c("chr", "pos") %in% names(k)) && !"source" %in% names(k))
})

test_that("hg38 follows the SeSAMe manifest of the first array that carries the probe", {
  X <- cpgdirection:::.cpgd_data("cpg_coordinates")
  # a 450K probe EPIC v2 dropped: placed by the manifests, absent from the engine's hg19 table
  old <- X[!grepl("EPICv2", X$array, fixed = TRUE) & grepl("450K", X$array, fixed = TRUE) &
             !is.na(X$pos_hg19) & !is.na(X$pos_hg38)][1:5]
  a <- cpgd_cpg_annotation(old$cpg_id, build = c("hg19", "hg38", "chm13"))
  expect_true(all(a$known & a$on_450k & !a$on_epic_v2))
  expect_identical(a$pos_hg19, as.integer(old$pos_hg19))
  expect_identical(a$pos_hg38, as.integer(old$pos_hg38))
  expect_true(all(a$source_hg19 %in% c("illumina_450k_manifest", "illumina_epic_v1_manifest")))
  expect_true(all(a$source_hg38 == "zhou_hm450_hg38"))
  expect_true(all(is.na(a$pos_chm13)))
  skip_if_not(cpgd_has_data("cpg_positions_hg19"), "position table unavailable")
  expect_false(any(old$cpg_id %chin% cpgd_cpg_positions()$cpg_id))
  # probes two sources place differently are flagged, and still placed by the priority
  amb <- X[X$hg38_ambiguous %in% TRUE & !is.na(X$pos_hg38)][1:5]
  b <- cpgd_cpg_annotation(amb$cpg_id, build = "hg38", genes = FALSE, qc = FALSE)
  expect_true(all(b$hg38_ambiguous))
  expect_identical(b$pos, as.integer(amb$pos_hg38))
  expect_false(any(cpgd_cpg_annotation(c("cg00000029", "nope"), build = "hg38")$hg38_ambiguous))
})

test_that("input handling: suffixes, case, order, unknown IDs, column switches", {
  inp <- c("cg00000109_TC21", "not_a_probe", "cg00000029_RBL2", "CG00000029", "cg00000029")
  a <- cpgd_cpg_annotation(inp, build = "hg38", genes = FALSE, qc = FALSE)
  expect_identical(a$input, inp)
  expect_identical(a$cpg_id, c("cg00000109", "not_a_probe", "cg00000029", "cg00000029", "cg00000029"))
  expect_identical(a$known, c(TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_true(is.na(a$pos[2]) && is.na(a$chr[2]) && !a$on_epic_v2[2] && is.na(a$array[2]))
  expect_identical(a$pos[3], a$pos[5]); expect_identical(a$pos[4], a$pos[5])
  expect_false(any(c("chr_hg19", "pos_chm13", "genes", "qc_exclude") %in% names(a)))
  expect_error(cpgd_cpg_annotation("cg00000029", build = "hg18"), "should be one of")
  expect_identical(nrow(cpgd_cpg_annotation(character(0))), 0L)
  # the non-CpG probes of EPIC v2, in the manifest's spelling or the package's lower case
  X <- cpgdirection:::.cpgd_data("cpg_coordinates")
  ch <- X$cpg_id[grepl("^ch\\.", X$cpg_id) & grepl("EPICv2", X$array, fixed = TRUE)][1]
  z <- cpgd_cpg_annotation(c(ch, tolower(ch), paste0(ch, "_BC21")), build = "hg38", genes = FALSE, qc = FALSE)
  expect_true(all(z$known) && all(!is.na(z$pos)) && length(unique(z$pos)) == 1L)
  expect_identical(z$cpg_id, rep(ch, 3))
})

test_that("with no input every probe is returned once, almost all placed on hg38, and the cache is not handed out", {
  all1 <- cpgd_cpg_annotation(build = "hg38", genes = FALSE, qc = FALSE)
  expect_false(anyDuplicated(tolower(all1$cpg_id)) > 0)
  expect_gt(nrow(all1), 1080000L)
  skip_if_not(cpgd_has_data("epicv2_probe_gene_annotation"), "annotation layer unavailable")
  for (a in c("on_450k", "on_epic_v1", "on_epic_v2")) expect_lt(sum(is.na(all1$pos[all1[[a]]])), 20L)
  expect_lt(sum(is.na(all1$pos)), 20L)
  all1[, pos := 0L]
  expect_true(cpgd_cpg_annotation("cg00000029", build = "hg38")$pos > 0L)
})
