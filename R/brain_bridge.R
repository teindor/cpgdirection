# The brain bridge: what a peripheral measurement licenses about the brain.
# New in 2.3.0. Everything here is built from open resources (IMAGE-CpG, AMAZE,
# GSE59685, Sommerer et al., Gibbs GSE15745, Brain-mMeta x BrainMeta v2);
# components derived from CC-BY-NC sources are deliberately NOT shipped.

#' Brain causal directions (tier T1/T2)
#'
#' The genome-wide brain layer: 220,383 CpG-gene pairs with an SMR causal
#' direction (Brain-mMeta mQTLs as exposure, BrainMeta v2 cis-eQTLs as outcome),
#' graded by a score validated out-of-sample in blood and against the
#' independent Gibbs brain catalogue.
#'
#' @section What the tiers are worth in brain:
#' Validated against measured brain eQTMs (GSE15745; a cohort verified absent
#' from both SMR arms), tier T1 holds 92-96\% across every significance rung of
#' the validation catalogue; tier T2 falls to 67-77\%. Accuracies are best
#' supported within ~100 kb (the validating array is promoter-weighted) and are
#' conditional on an association existing. Agreement between independently
#' derived blood and brain directions is only ~65\%, so do NOT substitute a
#' peripheral eQTM direction for these.
#'
#' @return A \code{data.table}: \code{cpg_id}, \code{target_gene},
#'   \code{direction}, \code{tier}, \code{score}, \code{n_instruments},
#'   \code{instrument_agreement}, \code{cpg_gene_dist}, \code{hla},
#'   \code{inv17q21}, \code{p_SMR}, \code{p_HEIDI}, \code{heidi_status}.
#' @examplesIf cpgd_has_data("brain_directions")
#' b <- cpgd_brain_directions()
#' head(b)
#' @export
cpgd_brain_directions <- function() {
    if (!is.null(.cpgd_env$brain)) {
        return(.cpgd_env$brain)
    }
    B <- .cpgd_data("brain_directions")
    data.table::setkeyv(B, "cpg_id")
    .cpgd_env$brain <- B
    B
}

#' The graded periphery-to-brain deliverable
#'
#' 728 CpGs where a tier-T1 brain causal direction coexists with robust
#' peripheral-brain methylation concordance (rho >= 0.5 on BOTH IMAGE-CpG
#' arrays, in at least two peripheral tissues, outside the extended MHC), with
#' per-cohort support flags and grades. Grade A++ (n = 205; 184 saliva-usable)
#' is supported by four independent cohorts.
#' @return A \code{data.table}, one row per CpG; see the paper for column
#'   definitions. \code{strict_qc == FALSE} rows sit on probes failing the
#'   stringent probe filter and should be treated as provisional.
#' @examples
#' d <- cpgd_bridge_deliverable()
#' head(d)
#' @export
cpgd_bridge_deliverable <- function() {
    if (!is.null(.cpgd_env$bridge)) {
        return(.cpgd_env$bridge)
    }
    D <- .cpgd_data("bridge_deliverable")
    data.table::setkeyv(D, "cpg_id")
    .cpgd_env$bridge <- D
    D
}

#' What does a peripheral CpG measurement license about the brain?
#'
#' For each input CpG: the brain causal direction and target gene(s) where a
#' tier-T1/T2 estimate exists, the peripheral-brain concordance evidence for
#' the tissue you measured, the deliverable grade where the CpG reaches the
#' graded set, and otherwise a genome-wide synthesised saliva-bridge score.
#'
#' @section The three caveats this function will keep repeating:
#' \enumerate{
#'   \item Directions are BRAIN-derived. Peripheral eQTM directions agree with
#'     brain only ~65\% of the time and are never used here.
#'   \item Peripheral-brain concordance couplings rest on n = 27 living-brain
#'     donors; per-CpG estimates have a cross-resource reliability of r ~ 0.2-0.4.
#'     Grades exist precisely because single estimates cannot be trusted alone.
#'   \item Sex-chromosome and discontinued probes produce replicable artefactual
#'     concordance; \code{strict_qc} is the guard. This function never returns
#'     an X/Y probe as bridge-usable.
#' }
#'
#' @param cpgs CpG identifiers in any of the formats
#'   \code{cpg_expression_direction()} accepts.
#' @param tissue Which peripheral tissue you measured: "saliva" (default),
#'   "blood" or "buccal".
#' @return A \code{data.table} of class \code{cpgd_bridge}: per CpG,
#'   \code{bridge_grade}, \code{bridge_usable} (concordance evidence in YOUR
#'   tissue), \code{rho_epic}/\code{rho_450k} (IMAGE-CpG, n = 27),
#'   \code{amaze_replicated}, \code{hannon_cortical_rho} (blood, n = 75),
#'   \code{sommerer_buccal_R} (n = 120), \code{ensemble_pctile} (genome-wide
#'   synthesised saliva score), \code{strict_qc}, and the brain layer's
#'   \code{t1_genes}/\code{t1_directions}. Where a gene was named in the input,
#'   \code{requested_gene}, \code{requested_gene_direction} and
#'   \code{requested_gene_tier} answer for that gene specifically.
#'
#' @section Co-regulation columns:
#'   \code{co_up} and \code{co_down} list T1 partner genes this CpG's
#'   methylation moves up and down respectively - genes within one column rise
#'   and fall together (causal co-direction, from the brain layer itself); a
#'   CpG with both columns filled is bivalent and its partners move in opposite
#'   directions. \code{hic_genes} lists genes contacted by the same
#'   adult-cortex enhancer the CpG sits in (PsychENCODE): candidate
#'   co-regulation by physical contact, with no direction attached - the
#'   evidence class that corroborated causal targets 2.5x over nearest genes.
#' @examplesIf cpgd_has_data("saliva_bridge_scores")
#' cpg_brain_bridge("cg06846259", tissue = "saliva") # POMC, grade A++
#' @export
cpg_brain_bridge <- function(cpgs, tissue = c("saliva", "blood", "buccal")) {
    tissue <- match.arg(tissue)
    q <- .cpgd_parse_input(cpgs)
    if (!nrow(q)) stop("No CpG identifiers recognised.", call. = FALSE)

    out <- data.table::data.table(cpg_id = unique(q$cpg_id))
    # Genes carried in the input (panel names like cg..._TC21_MC2R, or a gene
    # column in a data.frame) are honoured, not silently discarded: the
    # requested gene is echoed back with its own T1 answer, while t1_genes
    # still lists every partner so nothing is hidden.
    gq <- q[!duplicated(q$cpg_id), c("cpg_id", "given_gene", "given_gene2")]
    out <- merge(out, gq, by = "cpg_id", all.x = TRUE)

    ## brain layer: T1/T2 genes and directions per CpG
    B <- cpgd_brain_directions()
    bt <- B[B$tier == "T1"]
    agg <- bt[, list(
        t1_genes = paste(get("target_gene"), collapse = ";"),
        t1_directions = paste(get("direction"), collapse = ";")
    ),
    by = "cpg_id"
    ]
    out <- merge(out, agg, by = "cpg_id", all.x = TRUE)
    t2 <- unique(B[B$tier == "T2"]$cpg_id)
    out[, "has_t2_direction" := out$cpg_id %in% t2]

    ## deliverable
    D <- cpgd_bridge_deliverable()
    ok_col <- paste0(tissue, "_ok")
    am_raw <- paste0(tissue, "_amaze_raw")
    am_adj <- paste0(tissue, "_amaze_adj")
    rho_e <- paste0("rho_", tissue, "_epic")
    rho_4 <- paste0("rho_", tissue, "_450k") # buccal has no 450K column
    # NOT the deliverable's `genes`/`directions`: those duplicate t1_genes /
    # t1_directions and duplicated columns invite people to trust whichever they
    # noticed first.
    dcols <- intersect(
        c(
            "cpg_id", "grade_v4", "strict_qc", ok_col, am_raw, am_adj,
            rho_e, rho_4, "hannon_cortical_raw", "sommerer_buccal_R"
        ),
        names(D)
    )
    out <- merge(out, D[, dcols, with = FALSE], by = "cpg_id", all.x = TRUE)
    data.table::setnames(out,
        old = intersect(c(
            "grade_v4", ok_col, am_raw, am_adj, rho_e, rho_4,
            "hannon_cortical_raw"
        ), names(out)),
        new = intersect(
            c(
                "bridge_grade", "bridge_usable", "amaze_replicated_raw",
                "amaze_replicated_adj", "rho_epic", "rho_450k",
                "hannon_cortical_rho"
            ),
            c(
                "bridge_grade", "bridge_usable", "amaze_replicated_raw",
                "amaze_replicated_adj", "rho_epic", "rho_450k",
                "hannon_cortical_rho"
            )[
                c(
                    "grade_v4", ok_col, am_raw, am_adj, rho_e, rho_4,
                    "hannon_cortical_raw"
                ) %in% names(out)
            ]
        ),
        skip_absent = TRUE
    )

    ## genome-wide synthesised score (saliva-bridge; reported for every tissue,
    ## labelled for what it is)
    S <- .cpgd_data("saliva_bridge_scores", required = FALSE)
    if (!is.null(S)) {
        S <- S[, c("cpg_id", "ensemble_pctile"), with = FALSE]
        out <- merge(out, S, by = "cpg_id", all.x = TRUE)
    }

    ## co-regulation, where knowledge exists - two evidence classes, kept apart:
    ##
    ## (1) causal co-direction, from our own brain layer: genes moved by THIS
    ##     CpG's methylation in the same direction rise and fall together.
    ##     co_up = T1 partners at +1, co_down = T1 partners at -1. A CpG with
    ##     entries in both is bivalent - its partners move in OPPOSITE directions,
    ##     which is co-regulation but not co-activation, and listing them in one
    ##     column would assert the wrong thing.
    cu <- bt[bt$direction == 1L,
        list(co_up = paste(get("target_gene"), collapse = ";")),
        by = "cpg_id"
    ]
    cd <- bt[bt$direction == -1L,
        list(co_down = paste(get("target_gene"), collapse = ";")),
        by = "cpg_id"
    ]
    out <- merge(out, cu, by = "cpg_id", all.x = TRUE)
    out <- merge(out, cd, by = "cpg_id", all.x = TRUE)

    ## (2) physical co-targets, from adult-cortex Hi-C: genes contacted by the
    ##     same enhancer the CpG sits in. Candidate co-regulation only - contact
    ##     is not direction - but this evidence class corroborated the causal
    ##     target 2.5x over the nearest gene, so it earns a labelled column.
    HC <- .cpgd_data("hic_cotargets", required = FALSE)
    if (!is.null(HC)) {
        out <- merge(out, HC, by = "cpg_id", all.x = TRUE)
    }

    ## if the input named a gene, answer for that gene specifically
    has_req <- !all(is.na(out$given_gene)) || !all(is.na(out$given_gene2))
    if (has_req) {
        B1 <- cpgd_brain_directions()
        B1 <- B1[B1$tier %in% c("T1", "T2")]
        key <- data.table::data.table(
            cpg_id = out$cpg_id,
            g1 = toupper(out$given_gene), g2 = toupper(out$given_gene2)
        )
        hit <- function(gcol) {
            m <- merge(key[, c("cpg_id", gcol), with = FALSE],
                B1[, list(
                    cpg_id = get("cpg_id"),
                    g = toupper(get("target_gene")),
                    d = get("direction"), tr = get("tier")
                )],
                by.x = c("cpg_id", gcol), by.y = c("cpg_id", "g")
            )
            m[!duplicated(m$cpg_id)]
        }
        h <- hit("g1")
        h2 <- hit("g2")
        out[, "requested_gene" := data.table::fifelse(
            !is.na(out$given_gene),
            out$given_gene, out$given_gene2
        )]
        idx1 <- match(out$cpg_id, h$cpg_id)
        idx2 <- match(out$cpg_id, h2$cpg_id)
        out[, "requested_gene_direction" := data.table::fifelse(
            !is.na(idx1), h$d[idx1],
            data.table::fifelse(!is.na(idx2), h2$d[idx2], NA_integer_)
        )]
        out[, "requested_gene_tier" := data.table::fifelse(
            !is.na(idx1), h$tr[idx1],
            data.table::fifelse(!is.na(idx2), h2$tr[idx2], NA_character_)
        )]
    }
    out[, c("given_gene", "given_gene2") := NULL]

    data.table::setattr(out, "tissue", tissue)
    data.table::setattr(out, "class", c("cpgd_bridge", class(out)))
    out[]
}

#' @export
print.cpgd_bridge <- function(x, ...) {
    tis <- attr(x, "tissue")
    n <- nrow(x)
    arrow <- function(d) {
        if (is.na(d)) {
            ""
        } else if (as.integer(d) > 0) {
            "raises"
        } else {
            "lowers"
        }
    }

    if (n <= 10) {
        ## card layout: the answer first, evidence beneath it, one CpG at a time
        for (i in seq_len(n)) {
            r <- x[i]
            head_bits <- c(
                r$cpg_id,
                if (!is.na(r$bridge_grade)) sprintf("grade %s", r$bridge_grade) else "not in the graded set",
                if (isTRUE(r$strict_qc %in% c(TRUE, "TRUE"))) {
                    "strict QC"
                } else if (!is.na(r$bridge_grade)) "PROVISIONAL: discontinued probe" else NULL,
                if (isTRUE(r$bridge_usable %in% c(TRUE, "TRUE"))) {
                    sprintf("usable from %s", tis)
                } else if (!is.na(r$bridge_grade)) sprintf("NOT robust from %s", tis) else NULL
            )
            cat(paste(head_bits, collapse = "  |  "), "\n")

            if ("requested_gene" %in% names(x) && !is.na(r$requested_gene)) {
                if (!is.na(r$requested_gene_direction)) {
                    cat(sprintf(
                        "  %s: in brain, higher methylation here %s %s expression (%s)\n",
                        r$requested_gene, arrow(r$requested_gene_direction),
                        r$requested_gene, r$requested_gene_tier
                    ))
                } else {
                    cat(sprintf(
                        "  %s: not a brain partner of this CpG (no substitution made)\n",
                        r$requested_gene
                    ))
                }
            }
            if (!is.na(r$t1_genes)) {
                # One plain statement per direction, never a sign convention the reader
                # has to decode: "-POMC" cost a question; a sentence would not have.
                g <- strsplit(r$t1_genes, ";", fixed = TRUE)[[1]]
                d <- strsplit(r$t1_directions, ";", fixed = TRUE)[[1]]
                lo <- g[d != "1"]
                hi <- g[d == "1"]
                parts <- c(
                    if (length(lo)) sprintf("lowers %s expression", paste(lo, collapse = ", ")),
                    if (length(hi)) sprintf("raises %s expression", paste(hi, collapse = ", "))
                )
                cat(sprintf(
                    "  brain (T1): higher methylation here %s\n",
                    paste(parts, collapse = " and ")
                ))
            } else if (isTRUE(r$has_t2_direction)) {
                cat("  brain: T2 evidence only (see cpgd_brain_directions())\n")
            } else {
                cat("  brain: no tiered causal direction\n")
            }

            cw <- character(0)
            if (!is.na(r$rho_epic)) {
                cw <- c(cw, sprintf(
                    "%s rho %.2f/%.2f (two arrays, n=27)", tis,
                    r$rho_epic,
                    if ("rho_450k" %in% names(x) && !is.na(r$rho_450k)) r$rho_450k else NA
                ))
            }
            if ("amaze_replicated_raw" %in% names(x) && !is.na(r$bridge_grade)) {
                cw <- c(cw, sprintf(
                    "AMAZE %s",
                    if (isTRUE(r$amaze_replicated_raw %in% c(TRUE, "TRUE")) ||
                        isTRUE(r$amaze_replicated_adj %in% c(TRUE, "TRUE"))) {
                        "replicated"
                    } else {
                        "not replicated"
                    }
                ))
            }
            if (!is.na(r$hannon_cortical_rho)) {
                cw <- c(cw, sprintf("blood-brain %.2f (n=75)", r$hannon_cortical_rho))
            }
            if ("sommerer_buccal_R" %in% names(x) && !is.na(r$sommerer_buccal_R)) {
                cw <- c(cw, sprintf("buccal-brain %.2f (n=120)", r$sommerer_buccal_R))
            }
            if (length(cw)) cat("  concordance: ", paste(cw, collapse = " | "), "\n", sep = "")
            if ("ensemble_pctile" %in% names(x) && !is.na(r$ensemble_pctile)) {
                cat(sprintf("  synthesised bridge score: %.1f percentile\n", r$ensemble_pctile))
            }
            co <- character(0)
            if (!is.na(r$co_up)) co <- c(co, paste0("up-together: ", r$co_up))
            if (!is.na(r$co_down)) co <- c(co, paste0("down-together: ", r$co_down))
            if (length(co) == 2) co <- c(co, "(bivalent CpG: the two sets move oppositely)")
            if (length(co)) cat("  co-regulated: ", paste(co, collapse = "  "), "\n", sep = "")
            if ("hic_genes" %in% names(x) && !is.na(r$hic_genes)) {
                cat("  same-enhancer contacts (Hi-C, no direction): ", r$hic_genes, "\n", sep = "")
            }
            cat("\n")
        }
        cat("Caveats: directions are brain-derived (peripheral eQTM directions agree\n")
        cat("with brain only ~65%); couplings rest on n = 27 donors - trust grades\n")
        cat("over any single rho. Full columns: as.data.frame(x).\n")
        return(invisible(x))
    }

    ## many CpGs: summary + compact table, no 20-column dump
    cat(sprintf("cpg_brain_bridge: %d CpGs, measured tissue = %s\n", n, tis))
    ng <- sum(!is.na(x$bridge_grade))
    if (ng) {
        tb <- table(x$bridge_grade)
        cat(sprintf(
            "  graded deliverable: %d  (%s)\n", ng,
            paste(sprintf("%s %d", names(tb), tb), collapse = ", ")
        ))
    } else {
        cat("  graded deliverable: 0\n")
    }
    cat(sprintf(
        "  usable from %s: %d   T1 brain direction: %d   T2 only: %d\n",
        tis, sum(x$bridge_usable %in% c(TRUE, "TRUE")),
        sum(!is.na(x$t1_genes)),
        sum(x$has_t2_direction & is.na(x$t1_genes))
    ))
    if ("ensemble_pctile" %in% names(x)) {
        k <- sum(is.na(x$bridge_grade) & x$ensemble_pctile >= 95, na.rm = TRUE)
        cat(sprintf("  ungraded but top-5%% synthesised bridge score: %d\n", k))
    }
    show <- intersect(c(
        "cpg_id", "bridge_grade", "bridge_usable", "t1_genes",
        "t1_directions", "ensemble_pctile", "strict_qc"
    ), names(x))
    cat("\n")
    print(data.table::as.data.table(x)[, show, with = FALSE], nrows = 20)
    cat("\nFull columns: as.data.frame(x). Caveats: see ?cpg_brain_bridge.\n")
    invisible(x)
}
