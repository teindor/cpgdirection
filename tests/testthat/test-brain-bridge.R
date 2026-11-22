# The brain bridge (2.3.0). These tests pin the claims the Paper 2 manuscript
# makes about the shipped module, so a drifting data file fails loudly.

test_that("brain directions table is shipped and shaped", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    B <- cpgd_brain_directions()
    expect_true(all(c(
        "cpg_id", "target_gene", "direction", "tier",
        "heidi_status"
    ) %in% names(B)))
    expect_gt(nrow(B), 200000)
    expect_setequal(setdiff(unique(B$tier), NA), c("T1", "T2", "T3"))
    # autosomal pipeline: no sex-chromosome artefact route into the bridge
    if ("cpg_chr" %in% names(B)) {
        expect_false(any(grepl("^(X|Y)$", as.character(B$cpg_chr))))
    }
})

test_that("the deliverable matches the manuscript counts", {
    D <- cpgd_bridge_deliverable()
    expect_equal(nrow(D), 728L)
    expect_true(all(c("grade_v4", "strict_qc", "saliva_ok") %in% names(D)))
    app <- D[D$grade_v4 == "A++" & D$strict_qc == TRUE]
    expect_equal(nrow(app), 205L)
    expect_equal(sum(app$saliva_ok == TRUE), 184L)
})

test_that("cg06846259 returns the POMC story", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    r <- cpg_brain_bridge("cg06846259", tissue = "saliva")
    expect_equal(nrow(r), 1L)
    expect_true(grepl("POMC", r$t1_genes[1]))
    expect_true(grepl("-1", r$t1_directions[1]))
    expect_true(r$bridge_usable[1] %in% c(TRUE, "TRUE"))
})

test_that("tissue argument switches the usable flag", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    D <- cpgd_bridge_deliverable()
    # a CpG usable from blood but not saliva
    cand <- D[D$blood_ok == TRUE & D$saliva_ok == FALSE]$cpg_id
    skip_if(length(cand) == 0)
    a <- cpg_brain_bridge(cand[1], tissue = "blood")
    b <- cpg_brain_bridge(cand[1], tissue = "saliva")
    expect_true(a$bridge_usable[1] %in% c(TRUE, "TRUE"))
    expect_false(b$bridge_usable[1] %in% c(TRUE, "TRUE"))
})

test_that("unknown CpGs come back scored or NA, never dropped", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    r <- cpg_brain_bridge(c("cg06846259", "cg00000029"))
    expect_equal(nrow(r), 2L)
})

test_that("every input format the grand function takes works here too", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    forms <- list(
        bare   = "cg06846259",
        panel  = "Zcg06846259_TC21_POMC",
        twotok = "cg06846259_POMC",
        frame  = data.frame(probe = "cg06846259", gene = "POMC")
    )
    for (nm in names(forms)) {
        r <- cpg_brain_bridge(forms[[nm]])
        expect_equal(r$cpg_id[1], "cg06846259", info = nm)
    }
})

test_that("co-regulation columns exist and keep their evidence classes apart", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    r <- cpg_brain_bridge("cg06846259")
    expect_true(all(c("co_up", "co_down", "hic_genes") %in% names(r)))
    # POMC is a -1 partner: it must sit in co_down, never in co_up
    expect_true(grepl("POMC", r$co_down[1]))
    expect_true(is.na(r$co_up[1]) || !grepl("POMC", r$co_up[1]))

    # a multi-partner CpG splits its partners by sign, consistently with t1_genes
    B <- cpgd_brain_directions()
    t1 <- B[B$tier == "T1"]
    multi <- names(which(table(t1$cpg_id) >= 3))[1]
    skip_if(is.na(multi))
    m <- cpg_brain_bridge(multi)
    listed <- unlist(strsplit(na.omit(c(m$co_up[1], m$co_down[1])), ";"))
    expect_setequal(listed, unlist(strsplit(m$t1_genes[1], ";")))
})

test_that("a gene named in the input is answered, not discarded", {
    skip_if_not(cpgd_has_data("brain_directions"), "brain directions unavailable")
    r <- cpg_brain_bridge("cg06846259_POMC")
    expect_true("requested_gene_direction" %in% names(r))
    expect_equal(as.integer(r$requested_gene_direction[1]), -1L)
    expect_equal(r$requested_gene_tier[1], "T1")
    # and a gene that is NOT a partner of this CpG comes back NA, never swapped
    r2 <- cpg_brain_bridge("cg06846259_CREBBP")
    expect_true(is.na(r2$requested_gene_direction[1]))
})
