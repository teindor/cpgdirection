# Regression tests for two bugs found in 2.2.3, both of which returned a
# plausible answer rather than an error, and so would not have been caught by
# anything already in the suite.
#
# 1. The parser required more than two underscore-separated tokens, so the
#    plainest identifier a user would type -- "cg00039463_TRAP1" -- lost its
#    gene entirely and the whole distance layer went unavailable.
#
# 2. The SMR layer keyed its join on requested_gene only. With no gene supplied
#    the key was NA, so the layer could never fire in the mode most people run.
#    Every S1/S2/S3 figure was correct, and none of it reached the user.

test_that("a two-token cg_GENE identifier yields the gene", {
    d <- cpgdirection:::.cpgd_parse_input("cg00039463_TRAP1")
    expect_equal(d$cpg_id[1], "cg00039463")
    expect_true("TRAP1" %in% toupper(c(d$given_gene[1], d$given_gene2[1])))
})

test_that("probe-type codes are dropped, not counted as genes", {
    d <- cpgdirection:::.cpgd_parse_input(
        c("Zcg00335286_TC21_MC2R", "cg09527270_TC11", "cg00039463_BC11_CREBBP")
    )
    expect_equal(d$cpg_id, c("cg00335286", "cg09527270", "cg00039463"))
    expect_equal(toupper(d$given_gene2[1]), "MC2R")
    expect_true(is.na(d$given_gene[2])) # TC11 alone is not a gene
    expect_equal(toupper(d$given_gene2[3]), "CREBBP")
})

test_that("the SMR layer is reported when no gene is supplied", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    # Take a CpG the shipped table definitely covers rather than hard-coding one:
    # the extdata file can be regenerated, and a literal would rot silently.
    cpg <- S$cpg_id[1]
    r <- cpg_expression_direction(cpg, verbose = FALSE)

    # Reported unconditionally. Whether it ALSO supplies best_direction depends on
    # smr_gene_match, which is the next test; the layer being visible at all is
    # what regressed in 2.2.3.
    expect_true(all(c("smr_tier", "smr_gene", "smr_gene_match") %in% names(r)))
    expect_false(is.na(r$smr_tier[1]))
    expect_false(is.na(r$smr_direction[1]))
    expect_false(is.na(r$smr_gene[1]))
})

test_that("only gene-matched SMR rows may supply best_direction", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    r <- cpg_expression_direction(utils::head(unique(S$cpg_id), 300), verbose = FALSE)
    smr_won <- r$best_evidence %in% c("smr_high", "smr_moderate", "smr_weak")
    # Every row where SMR won must be a gene match. The converse need not hold:
    # a matched row can still be outranked by a measured eQTM.
    expect_true(all(r$smr_gene_match[smr_won] %in% TRUE))
    # And agreement is only computed where the two concern the same gene.
    expect_true(all(is.na(r$smr_agreement[!(r$smr_gene_match %in% TRUE)])))
})

test_that("a supplied gene is never silently replaced by another", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    cpg <- S$cpg_id[1]
    # A symbol that cannot be this CpG's SMR target. The direction may still be
    # REPORTED for whatever gene the instrument reaches, but it must be flagged as
    # a non-match and must not become best_direction.
    r <- cpg_expression_direction(cpg, genes = "ZZZ_NOT_A_GENE", verbose = FALSE)
    expect_false(r$smr_gene_match[1])
    expect_false(r$best_evidence[1] %in% c("smr_high", "smr_moderate", "smr_weak"))
    expect_true(is.na(r$smr_agreement[1]))
})

test_that("HEIDI columns are shipped and surfaced", {
    # The manuscript states that p_HEIDI, nsnp_HEIDI and a three-level status are
    # shipped with the package. This test is what makes that sentence true.
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    expect_true(all(c("p_HEIDI", "nsnp_HEIDI", "heidi_status") %in% names(S)))
    expect_setequal(
        setdiff(unique(S$heidi_status), NA),
        c("pass", "fail", "not_tested")
    )
    # not_tested must be exactly the rows with too few instruments, never a
    # relabelled pass.
    expect_true(all(is.na(S$p_HEIDI[S$heidi_status == "not_tested"])))
    expect_true(all(S$p_HEIDI[S$heidi_status == "pass"] >= 0.05))
    expect_true(all(S$p_HEIDI[S$heidi_status == "fail"] < 0.05))

    r <- cpg_expression_direction(S$cpg_id[1], verbose = FALSE)
    expect_true(all(c("smr_heidi_status", "smr_p_heidi") %in% names(r)))
})

test_that("nothing in the package conditions on HEIDI", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    # Tiers must be independent of HEIDI status: if a future edit quietly starts
    # filtering, S1 would stop appearing among HEIDI-failing pairs.
    expect_true(nrow(S[S$smr_tier == "S1" & S$heidi_status == "fail", ]) > 0)
    expect_true(nrow(S[S$smr_tier == "S1" & S$heidi_status == "pass", ]) > 0)
})

test_that("the compact report carries the columns needed to read an SMR row", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")
    skip_if_not_installed("gt")

    r <- cpg_expression_direction(utils::head(unique(S$cpg_id), 20), verbose = FALSE)
    f <- tempfile(fileext = ".html")
    on.exit(unlink(f), add = TRUE)
    suppressMessages(cpgd_report(r, file = f))
    txt <- paste(readLines(f, warn = FALSE), collapse = " ")
    # A direction without its gene is ambiguous; the report must not drop it.
    for (col in c("smr_gene", "smr_gene_match", "smr_heidi_status")) {
        expect_true(grepl(col, txt, fixed = TRUE), info = col)
    }
})

test_that("smr_in_table separates absent CpGs from unmatched genes", {
    S <- tryCatch(cpgd_smr_directions(), error = function(e) NULL)
    skip_if(is.null(S) || !nrow(S), "SMR table unavailable")

    cpg <- S$cpg_id[1]
    r <- cpg_expression_direction(cpg, genes = "ZZZ_NOT_A_GENE", verbose = FALSE)
    expect_true("smr_in_table" %in% names(r))
    expect_true(r$smr_in_table[1]) # present in the table, just not for that gene
})
