## cpgd_build_gene_tss(source = "manifest", genome_build =)

.tss_manifest <- function(cols, rows) {
    f <- tempfile(fileext = ".csv")
    writeLines(c(paste(cols, collapse = ","), rows), f)
    f
}

test_that("an hg38 manifest gives an hg38 table and is refused as hg19", {
    f <- .tss_manifest(
        c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "Genome_Build"),
        c(
            "cg1_TC21,1,1000,GENEA,GRCh38", "cg2_TC21,1,3000,GENEA,GRCh38",
            "cg3_TC21,2,500,GENEB,GRCh38"
        )
    )
    d <- suppressMessages(cpgd_build_gene_tss(
        source = "manifest", manifest_file = f,
        genome_build = "hg38"
    ))
    expect_identical(attr(d, "genome_build"), "hg38")
    expect_identical(d[d$gene == "GENEA", ]$tss, 2000L)
    expect_identical(sort(d$chr), c("chr1", "chr2"))
    expect_error(
        suppressMessages(cpgd_build_gene_tss(
            source = "manifest",
            manifest_file = f
        )),
        "This manifest is hg38"
    )
})

test_that("build-labelled columns are chosen by genome_build", {
    f <- .tss_manifest(
        c(
            "IlmnID", "CHR_hg19", "MAPINFO_hg19", "CHR_hg38", "Start_hg38",
            "UCSC_RefGene_Name"
        ),
        c("cg1_TC21,1,1000,1,1500,GENEA", "cg2_TC21,3,700,3,900,GENEB")
    )
    d19 <- suppressMessages(cpgd_build_gene_tss(source = "manifest", manifest_file = f))
    d38 <- suppressMessages(cpgd_build_gene_tss(
        source = "manifest", manifest_file = f,
        genome_build = "hg38"
    ))
    expect_identical(attr(d19, "genome_build"), "hg19")
    expect_identical(d19[d19$gene == "GENEA", ]$tss, 1000L)
    expect_identical(d38[d38$gene == "GENEA", ]$tss, 1500L)
    expect_identical(d38[d38$gene == "GENEB", ]$tss, 900L)
})

test_that("hg38 is refused when the manifest carries no hg38 coordinates", {
    f <- .tss_manifest(
        c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name"),
        c("cg1_TC21,1,1000,GENEA")
    )
    expect_error(
        cpgd_build_gene_tss(
            source = "manifest", manifest_file = f,
            genome_build = "hg38"
        ),
        "carries no hg38 coordinates"
    )
    d <- suppressMessages(cpgd_build_gene_tss(source = "manifest", manifest_file = f))
    expect_identical(attr(d, "genome_build"), "hg19")
})

test_that("the TxDb source is hg19 only", {
    expect_error(
        cpgd_build_gene_tss(source = "txdb", genome_build = "hg38"),
        "hg19 only"
    )
})
