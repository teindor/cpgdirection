# The pair-preserving API: automatic target discovery, many-to-many candidate
# union, pair-level direction resolution, and pair invariance.
#
# Fixture tests inject small source tables through the `sources` argument so
# they are hermetic: they exercise the discovery/union/resolution logic without
# depending on the packaged resources or the Bioconductor annotation stack.
# Real-resource tests are guarded and run wherever the installed package is
# complete.

library(data.table)

# ---- fixtures --------------------------------------------------------------

fx_manifest <- function(...) {
  d <- data.table(...)
  d[, array := "EPICv2"]
  d[, annotation_source := "Illumina_EPICv2"]
  if (!"refgene_group" %in% names(d)) d[, refgene_group := "TSS200"]
  d
}

fx_lookup_row <- function(cpg, gene, direction = -1, tier = "A",
                          conf = 0.8, prob = 0.1, status = "PREDICTED") {
  data.table(cpg_id = cpg, target_gene = gene, tss_dist = 500,
             status = status, direction = direction,
             probability_plus1 = prob, confidence = conf,
             evidence_tier = tier)
}

fx_measured_row <- function(cpg, gene, direction = 1, tissue = "blood") {
  data.table(cpg_id = cpg, target_gene = gene, tissue = tissue,
             direction = direction, tss_dist = 500)
}

fx_smr_row <- function(cpg, gene, direction = -1, tier = "S1", p = 1e-10) {
  data.table(cpg_id = cpg, target_gene = gene, direction = direction,
             smr_tier = tier, p_SMR = p, n_instruments = 3L,
             instrument_agreement = 1, cpg_gene_dist = 1000,
             top_instrument = "rs1", p_HEIDI = 0.5, nsnp_HEIDI = 3L,
             heidi_status = "pass")
}

CG1 <- "cg10000001"
CG2 <- "cg10000002"

fx_qc_row <- function(cpg, masked_general = FALSE, masked_partial = FALSE,
                      reasons = "", pos_ok = TRUE) {
  data.table(cpg_id = cpg, n_probes = 1L,
             n_masked_general = as.integer(masked_general),
             mask_reasons = reasons, masked_general = masked_general,
             masked_partial = masked_partial, pos_hg19_verified = pos_ok)
}

# empty-but-typed stand-ins: "this source exists and holds nothing"
EMPTY_SRC <- list(
  manifest = fx_manifest(cpg_id = "cgX", target_gene = "X")[0],
  lookup   = list(blood = fx_lookup_row("cgX", "X")[0]),
  measured = fx_measured_row("cgX", "X")[0],
  smr      = fx_smr_row("cgX", "X")[0],
  probe_qc = fx_qc_row("cgX")[0])

src_with <- function(...) {
  s <- EMPTY_SRC
  mod <- list(...)
  for (nm in names(mod)) s[[nm]] <- mod[[nm]]
  s
}

# ---- Test 1: bare manifest CpG discovers its gene, no genes= required ------

test_that("a bare CpG auto-discovers its manifest gene", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = CG1,
                                               target_gene = "GENEA")))
  expect_equal(p$target_gene, "GENEA")
  expect_true(p$has_manifest)
  expect_match(p$mapping_sources, "EPICv2_manifest")
})

# ---- Test 2: manifest multi-gene annotation stays many-to-many -------------

test_that("a multi-gene manifest annotation yields one row per gene", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = c(CG1, CG1),
                                               target_gene = c("GENEA", "GENEB"))))
  expect_equal(nrow(p), 2L)
  expect_setequal(p$target_gene, c("GENEA", "GENEB"))
  expect_true(all(p$is_coeffect))
  expect_true(all(p$n_targets_for_cpg == 2L))
})

# ---- Tests 3-4: measured and SMR augment the manifest targets --------------

test_that("measured eQTM targets augment manifest targets", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
                        measured = fx_measured_row(CG1, "GENEB", direction = 1)))
  expect_setequal(p$target_gene, c("GENEA", "GENEB"))
  b <- p[p$target_gene == "GENEB", ]
  expect_equal(b$mapping_primary, "measured_eQTM")
  expect_equal(b$best_evidence, "measured")
  expect_equal(b$best_direction, 1)
})

test_that("SMR targets augment manifest and measured targets", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
                        measured = fx_measured_row(CG1, "GENEB"),
                        smr      = fx_smr_row(CG1, "GENEC", direction = -1,
                                              tier = "S1")))
  expect_setequal(p$target_gene, c("GENEA", "GENEB", "GENEC"))
  cc <- p[p$target_gene == "GENEC", ]
  expect_equal(cc$best_evidence, "smr_high")
  expect_equal(cc$best_direction, -1)
})

# ---- Test 5: an input suffix augments rather than overrides ----------------

test_that("a suffix gene augments the auto-discovered targets by default", {
  p <- cpg_gene_pairs(paste0(CG1, "_GENEC"), universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = c(CG1, CG1),
                                               target_gene = c("GENEA", "GENEB"))))
  expect_setequal(p$target_gene, c("GENEA", "GENEB", "GENEC"))
  cc <- p[p$target_gene == "GENEC", ]
  expect_equal(cc$mapping_primary, "input_annotation")
})

test_that("annotation_mode = 'strict' restricts an annotated input to its gene", {
  p <- cpg_gene_pairs(c(paste0(CG1, "_GENEA"), CG2),
                      annotation_mode = "strict",
                      universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(
                          cpg_id = c(CG1, CG1, CG2),
                          target_gene = c("GENEA", "GENEB", "GENEX"))))
  expect_setequal(p[p$cpg_id == CG1, ]$target_gene, "GENEA")
  expect_setequal(p[p$cpg_id == CG2, ]$target_gene, "GENEX")
})

# ---- pair invariance -------------------------------------------------------

test_that("a fixed pair resolves identically however the gene was discovered", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    lookup   = list(blood = fx_lookup_row(CG1, "GENEA", direction = -1,
                                          tier = "A", conf = 0.8)))
  bare     <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                             sources = srcs)
  suffixed <- cpg_gene_pairs(paste0(CG1, "_GENEA"), universal = FALSE,
                             verbose = FALSE, sources = srcs)
  pairwise <- cpg_gene_pairs(CG1, genes = "GENEA", gene_mode = "pairwise",
                             universal = FALSE, verbose = FALSE, sources = srcs)
  rows <- lapply(list(bare, suffixed, pairwise),
                 function(p) p[p$target_gene == "GENEA", ])
  for (r in rows) {
    expect_equal(r$best_direction, rows[[1]]$best_direction)
    expect_equal(r$best_evidence,  rows[[1]]$best_evidence)
    expect_equal(r$best_confidence, rows[[1]]$best_confidence)
    expect_equal(r$direction_tier, rows[[1]]$direction_tier)
  }
  # only the mapping provenance may differ
  expect_match(suffixed[suffixed$target_gene == "GENEA", ]$mapping_sources,
               "input_annotation")
  expect_match(pairwise[pairwise$target_gene == "GENEA", ]$mapping_sources,
               "requested_pair")
})

# ---- parser: no cpg_id-only collapse in the pair workflow ------------------

test_that("two suffixed inputs on one CpG both survive parsing", {
  q <- cpgdirection:::.cpgd_parse_input(
    c(paste0(CG1, "_GENEA"), paste0(CG1, "_GENEB")), dedupe = FALSE)
  expect_equal(nrow(q), 2L)
  expect_equal(q$input_id, 1:2)
  expect_setequal(q$given_gene, c("GENEA", "GENEB"))

  p <- cpg_gene_pairs(c(paste0(CG1, "_GENEA"), paste0(CG1, "_GENEB")),
                      universal = FALSE, verbose = FALSE,
                      sources = EMPTY_SRC)
  expect_setequal(p$target_gene, c("GENEA", "GENEB"))
})

test_that("the legacy parser default still collapses technical replicates", {
  q <- cpgdirection:::.cpgd_parse_input(
    c(paste0(CG1, "_TC11"), paste0(CG1, "_TC12")))
  expect_equal(nrow(q), 1L)
  expect_equal(attr(q, "n_raw"), 2L)
})

# ---- filter and pairwise semantics -----------------------------------------

test_that("gene_mode = 'filter' restricts targets discovered first", {
  srcs <- src_with(manifest = fx_manifest(
    cpg_id = rep(CG1, 3), target_gene = c("GENEA", "GENEB", "GENEC")))
  p <- cpg_gene_pairs(CG1, genes = c("GENEA", "GENEC"), gene_mode = "filter",
                      universal = FALSE, verbose = FALSE, sources = srcs)
  expect_setequal(p$target_gene, c("GENEA", "GENEC"))
})

test_that("gene_mode = 'pairwise' evaluates exactly the supplied pairs", {
  srcs <- src_with(manifest = fx_manifest(
    cpg_id = rep(CG1, 3), target_gene = c("GENEA", "GENEB", "GENEC")))
  p <- cpg_gene_pairs(cpgs = c(CG1, CG1), genes = c("GENEA", "GENEB"),
                      gene_mode = "pairwise",
                      universal = FALSE, verbose = FALSE, sources = srcs)
  expect_equal(nrow(p), 2L)
  expect_setequal(p$target_gene, c("GENEA", "GENEB"))
})

test_that("gene_mode is never inferred from matching vector lengths", {
  srcs <- src_with(manifest = fx_manifest(
    cpg_id = c(CG1, CG2), target_gene = c("GENEA", "GENEB")))
  # two cpgs, two genes, default filter mode: genes act as a SET filter,
  # so GENEB is kept for CG2 even though genes[1] pairs positionally with CG1
  p <- cpg_gene_pairs(c(CG1, CG2), genes = c("GENEB", "GENEA"),
                      universal = FALSE, verbose = FALSE, sources = srcs)
  expect_setequal(p$pair_id, c(paste0(CG1, "|GENEA"), paste0(CG2, "|GENEB")))
})

test_that("pairwise mode requires genes", {
  expect_error(cpg_gene_pairs(CG1, gene_mode = "pairwise", verbose = FALSE),
               "pairwise")
})

# ---- opposite directions and pair-specific abstention ----------------------

test_that("opposite directions across target genes are both retained", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        measured = rbind(
                          fx_measured_row(CG1, "GENEA", direction = -1),
                          fx_measured_row(CG1, "GENEB", direction = 1))))
  expect_equal(nrow(p), 2L)
  expect_equal(p[p$target_gene == "GENEA", ]$best_direction, -1)
  expect_equal(p[p$target_gene == "GENEB", ]$best_direction, 1)
})

test_that("abstention is pair-specific, never CpG-global", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = c(CG1, CG1),
                                               target_gene = c("GENEA", "GENEB")),
                        measured = fx_measured_row(CG1, "GENEA", direction = 1)))
  a <- p[p$target_gene == "GENEA", ]
  b <- p[p$target_gene == "GENEB", ]
  expect_true(a$usable)
  expect_false(b$usable)
  expect_equal(b$best_evidence, "no_evidence")
  expect_true(nzchar(b$abstain_reason))
})

# ---- evidence ladder inside the pair ---------------------------------------

test_that("the pair ladder keeps the package's evidence precedence", {
  # measured beats SMR S1 beats catalogue consensus
  srcs <- src_with(
    measured = fx_measured_row(CG1, "GENEA", direction = 1),
    smr      = fx_smr_row(CG1, "GENEA", direction = -1, tier = "S1"),
    lookup   = list(
      blood            = fx_lookup_row(CG1, "GENEA", direction = -1),
      nasal_epithelium = fx_lookup_row(CG1, "GENEA", direction = -1)))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_equal(p$best_evidence, "measured")
  expect_equal(p$best_direction, 1)

  srcs$measured <- EMPTY_SRC$measured
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_equal(p$best_evidence, "smr_high")
  expect_equal(p$best_direction, -1)

  srcs$smr <- EMPTY_SRC$smr
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_equal(p$best_evidence, "catalogue_consensus")
  expect_equal(p$best_direction, -1)
})

test_that("evidence for a different gene is never promoted into a pair", {
  # SMR knows GENEB; the requested pair is GENEA: its row must not borrow it
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    smr      = fx_smr_row(CG1, "GENEB", direction = -1, tier = "S1"))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  a <- p[p$target_gene == "GENEA", ]
  expect_true(is.na(a$smr_direction))
  expect_false(a$best_evidence %in% c("smr_high", "smr_moderate", "smr_weak"))
  b <- p[p$target_gene == "GENEB", ]
  expect_equal(b$best_evidence, "smr_high")
})

test_that("catalogue tissues disagreeing on a pair is a conflict, not a coin flip", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    lookup = list(
      blood            = fx_lookup_row(CG1, "GENEA", direction = -1),
      nasal_epithelium = fx_lookup_row(CG1, "GENEA", direction = 1)))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_equal(p$best_evidence, "tissue_conflict")
  expect_true(is.na(p$best_direction))
  expect_false(p$usable)
})

# ---- distance layer inside the pair ----------------------------------------

test_that("the distance layer resolves a pair no catalogue holds", {
  skip_if_not(file.exists(system.file("extdata", "distance_curves.csv",
                                      package = "cpgdirection")))
  srcs <- src_with(manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"))
  srcs$cpg_pos  <- data.table(cpg_id = CG1, chr = "chr1", pos = 1000000L)
  srcs$gene_tss <- data.table(gene = "GENEA", chr = "chr1", tss = 1000500L,
                              strand = "+")
  p <- cpg_gene_pairs(CG1, universal = TRUE, verbose = FALSE, sources = srcs)
  expect_true(p$best_evidence %in%
                c("distance_only", "distance_uninformative",
                  "distance_tissue_conflict"))
  expect_equal(p$abs_dist, 500)
})

# ---- warnings --------------------------------------------------------------

test_that("a CpG with no mapping warns about the resources, not 'supply genes'", {
  expect_warning(
    cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                   sources = EMPTY_SRC),
    "EPIC-v2 manifest, packaged catalogue, measured eQTM, SMR")
})

# ---- provenance bookkeeping ------------------------------------------------

test_that("all mapping sources supporting one pair are preserved", {
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE,
                      sources = src_with(
                        manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
                        measured = fx_measured_row(CG1, "GENEA"),
                        smr      = fx_smr_row(CG1, "GENEA")))
  expect_equal(nrow(p), 1L)
  expect_true(p$has_manifest & p$has_measured_eqtm & p$has_smr)
  expect_match(p$mapping_sources, "measured_eQTM")
  expect_match(p$mapping_sources, "SMR")
  expect_match(p$mapping_sources, "EPICv2_manifest")
  expect_equal(p$mapping_strength, "measured")
})

test_that("input rows are traceable from pairs", {
  p <- cpg_gene_pairs(c(paste0("Z", CG1, "_TC21_GENEA"), CG2),
                      universal = FALSE, verbose = FALSE,
                      sources = src_with(manifest = fx_manifest(
                        cpg_id = c(CG1, CG2), target_gene = c("GENEA", "GENEB"))))
  expect_equal(p[p$cpg_id == CG1, ]$input_id, 1L)
  expect_equal(p[p$cpg_id == CG2, ]$input_id, 2L)
  expect_match(p[p$cpg_id == CG1, ]$input, "TC21_GENEA")
})

# ---- real packaged resources ----------------------------------------------

test_that("bare EPIC v2 CpGs discover manifest genes from the packaged table", {
  skip_on_cran()
  skip_if_not(cpgdirection:::.cpgd_has_data("epicv2_probe_gene_annotation"),
              "annotation layer unavailable")
  m <- cpgd_manifest_genes()
  expect_true(all(c("cpg_id", "target_gene", "array", "annotation_source",
                    "refgene_group") %in% names(m)))
  expect_true("CRHBP" %in% m[m$cpg_id == "cg26261055", ]$target_gene)

  p <- suppressWarnings(cpg_gene_pairs("cg26261055", universal = FALSE,
                                       verbose = FALSE))
  expect_true("CRHBP" %in% p$target_gene)
  expect_true(p[p$target_gene == "CRHBP", ]$has_manifest)
})

test_that("pair results are invariant on real resources too", {
  skip_on_cran()
  skip_if_not(cpgdirection:::.cpgd_has_data("epicv2_probe_gene_annotation"),
              "annotation layer unavailable")
  a <- suppressWarnings(cpg_gene_pairs("cg26261055", universal = FALSE,
                                       verbose = FALSE))
  b <- suppressWarnings(cpg_gene_pairs("cg26261055_CRHBP", universal = FALSE,
                                       verbose = FALSE))
  ra <- a[a$target_gene == "CRHBP", ]
  rb <- b[b$target_gene == "CRHBP", ]
  expect_equal(ra$best_direction, rb$best_direction)
  expect_equal(ra$best_evidence,  rb$best_evidence)
  expect_equal(ra$best_confidence, rb$best_confidence)
})

test_that("cpg_expression_direction() is unchanged by the pair machinery", {
  skip_on_cran()
  r <- cpg_expression_direction("cg00000029", tissue = "blood", verbose = FALSE)
  expect_equal(nrow(r), 1L)
  expect_true(all(c("target_gene", "status", "call") %in% names(r)))
})

# ---- the manifest builder --------------------------------------------------

test_that("cpgd_build_manifest_genes() parses an Illumina-format CSV", {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "[Heading]",
    "Descriptor File Name,EPIC-8v2-0_TEST.csv",
    "[Assay]",
    paste("IlmnID,Name,CHR,MAPINFO,Strand_FR,UCSC_RefGene_Name,",
          "UCSC_RefGene_Group,GencodeV41_Name,GencodeV41_Group", sep = ""),
    "cg10000001_TC21,cg10000001,chr1,1000,F,GENEA;GENEA;GENEB,TSS200;TSS1500;Body,,",
    "cg10000002_BC11,cg10000002,chr2,2000,R,,,GENEC;ENSG00000000001.1,exon_1;TSS200",
    "cg10000003_TC11,cg10000003,chr3,3000,F,,,,",
    "[Controls]",
    "ctl1,negative,0,0"), f)
  out <- tempfile(fileext = ".csv.gz")
  suppressMessages(cpgd_build_manifest_genes(f, out = out))
  a <- read.csv(out, stringsAsFactors = FALSE)  # base R: no R.utils on the builders

  # every semicolon-delimited gene preserved, never just the first
  expect_setequal(a[a$cpg_id == "cg10000001", ]$target_gene,
                  c("GENEA", "GENEB"))
  # transcript groups pooled per gene
  expect_equal(a[a$cpg_id == "cg10000001" & a$target_gene == "GENEA",
                 ]$refgene_group, "TSS1500;TSS200")
  # Gencode supplements where UCSC is empty; Ensembl-ID entries are dropped
  expect_equal(a[a$cpg_id == "cg10000002", ]$target_gene, "GENEC")
  expect_equal(a[a$cpg_id == "cg10000002", ]$refgene_source, "GencodeV41")
  # unannotated probes and control rows contribute nothing
  expect_false("cg10000003" %in% a$cpg_id)
  expect_true(all(grepl("^cg", a$cpg_id)))
  expect_true(all(a$annotation_source == "Illumina_EPICv2"))
})

# ---- probe QC --------------------------------------------------------------

test_that("a fully masked probe is excluded by default", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = c(CG1, CG2),
                           target_gene = c("GENEA", "GENEB")),
    probe_qc = fx_qc_row(CG1, masked_general = TRUE,
                         reasons = "M_nonuniq"))
  p <- suppressMessages(
    cpg_gene_pairs(c(CG1, CG2), universal = FALSE, verbose = FALSE,
                   sources = srcs))
  expect_false(CG1 %in% p$cpg_id)                 # wrong-locus probe is out
  expect_true(CG2 %in% p$cpg_id)                  # clean CpG unaffected
  expect_equal(attr(p, "qc_excluded_cpgs"), CG1)
})

test_that("probe_qc = 'flag' keeps masked probes, marked", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = fx_qc_row(CG1, masked_general = TRUE, reasons = "M_nonuniq"))
  p <- cpg_gene_pairs(CG1, probe_qc = "flag", universal = FALSE,
                      verbose = FALSE, sources = srcs)
  expect_true(CG1 %in% p$cpg_id)
  expect_true(p$probe_masked)
  expect_equal(p$probe_mask_reasons, "M_nonuniq")
})

test_that("a partially masked CpG is kept and flagged, never excluded", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = fx_qc_row(CG1, masked_partial = TRUE))
  p <- suppressMessages(
    cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs))
  expect_true(CG1 %in% p$cpg_id)
  expect_true(p$probe_masked_partial)
  expect_false(p$probe_masked)
})

test_that("probe_qc = 'ignore' skips QC entirely", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = fx_qc_row(CG1, masked_general = TRUE))
  p <- cpg_gene_pairs(CG1, probe_qc = "ignore", universal = FALSE,
                      verbose = FALSE, sources = srcs)
  expect_true(CG1 %in% p$cpg_id)
  expect_false("probe_masked" %in% names(p))
})

test_that("excluding every submitted CpG returns an empty flagged result", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = fx_qc_row(CG1, masked_general = TRUE))
  expect_warning(
    p <- suppressMessages(
      cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)),
    "probe QC")
  expect_equal(nrow(p), 0L)
  expect_equal(attr(p, "qc_excluded_cpgs"), CG1)
})

test_that("an unverified hg19 position is visible on the pair row", {
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = fx_qc_row(CG1, pos_ok = FALSE))
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_false(p$pos_hg19_verified)
})

test_that("the packaged annotation carries multi-source agreement", {
  skip_on_cran()
  skip_if_not(cpgdirection:::.cpgd_has_data("epicv2_probe_gene_annotation"),
              "annotation layer unavailable")
  m <- cpgd_manifest_genes()
  skip_if_not("n_annotation_sources" %in% names(m))
  expect_true(all(c("annotation_source", "n_annotation_sources") %in% names(m)))
  expect_true(any(m$n_annotation_sources >= 2))
  crhbp <- m[m$cpg_id == "cg26261055" & m$target_gene == "CRHBP", ]
  expect_true(nrow(crhbp) == 1L && crhbp$n_annotation_sources >= 2)
})

test_that("the packaged probe QC table loads and flags a known fraction", {
  skip_on_cran()
  skip_if_not(cpgdirection:::.cpgd_has_data("epicv2_probe_qc"),
              "probe QC layer unavailable")
  qc <- cpgd_probe_qc()
  expect_true(all(c("masked_general", "masked_partial",
                    "pos_hg19_verified") %in% names(qc)))
  frac <- mean(qc$masked_general)
  expect_gt(frac, 0.005)   # the mask is real
  expect_lt(frac, 0.15)    # and not runaway
})

test_that("qc_exclude, when present, is the exclusion criterion", {
  qcrow <- fx_qc_row(CG1)
  qcrow[, qc_exclude := TRUE]   # e.g. Garvan WGBS cross-hybridization,
  # not the Zhou mask
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = c(CG1, CG2),
                           target_gene = c("GENEA", "GENEB")),
    probe_qc = qcrow)
  p <- suppressMessages(
    cpg_gene_pairs(c(CG1, CG2), universal = FALSE, verbose = FALSE,
                   sources = srcs))
  expect_false(CG1 %in% p$cpg_id)
  expect_true(CG2 %in% p$cpg_id)
})

test_that("Garvan flags travel onto pair rows when present in the QC table", {
  qcrow <- fx_qc_row(CG1)
  qcrow[, c("cross_hybridizing", "mapping_flagged", "qc_exclude") :=
          list(FALSE, TRUE, FALSE)]
  srcs <- src_with(
    manifest = fx_manifest(cpg_id = CG1, target_gene = "GENEA"),
    probe_qc = qcrow)
  p <- cpg_gene_pairs(CG1, universal = FALSE, verbose = FALSE, sources = srcs)
  expect_false(p$cross_hybridizing)
  expect_true(p$mapping_flagged)
})
