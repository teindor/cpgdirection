# cpgdirection 1.1.0 (in development -- lands after Bioconductor accepts 0.99.x)

* **`cpgd_cpg_annotation()`: where a CpG is, in hg19, hg38 and T2T-CHM13,
  for every probe of the four arrays.**
  - **What.** One row per probe ID, as your data name it (bare 450K and
    EPIC v1 IDs, EPIC v2 IDs with their replicate suffix): `chr` and `pos`
    for the build asked for (`build = "hg38"`; with several builds,
    `chr_<build>` and `pos_<build>`), the arrays that carry the probe
    (`array`, for example `"450K;EPICv1;EPICv2"`, with `on_hm27`, `on_450k`,
    `on_epic_v1`, `on_epic_v2`), the genes of the EPIC v2 annotation layer,
    and the probe quality flags. `source` says which table placed the probe.
  - **Why.** The coordinates were already in the package but in three
    tables, and only for EPIC v2 probes. Many groups work in hg38, and
    about 154,000 probe IDs of the 450K and EPIC v1 arrays that EPIC v2
    dropped had no position at all.
  - **hg38.** From the Zhou-lab (SeSAMe) hg38 manifests, taking the
    manifest of the first array generation that carries the probe (HM450,
    EPIC, EPIC v2, HM27; `CpG_beg` + 1), so a 450K probe sits exactly where
    the SeSAMe HM450 manifest puts it. Every probe ID is placed but 19 that
    no manifest maps. `hg38_ambiguous` flags 1,094 probes that two sources
    place differently (multi-mapping probes).
  - **hg19.** The package's own table wherever it holds the probe (always
    equal to `cpgd_cpg_positions()`), else the Illumina 450K and EPIC v1
    coordinates; 98.9% of all probe IDs. T2T-CHM13: EPIC v2 probes only.
  - **New table.** `cpg_coordinates` (1,084,297 probe IDs): array
    membership for every probe and the positions the package's other
    resources do not hold. A coordinate the package does not hold is `NA`,
    never lifted over.
  - **Unchanged.** The direction engine works in hg19 (its models and
    distance curves were fitted on hg19 distances). `cpgd_cpg_positions()`
    and every direction are as before.

* **The catalogue's training chain is rebuilt (cpgdirectionData 0.99.7),
  and the ladder's figures with it.**
  - **Why.** 0.99.5 corrected the gene placement but kept the original
    pools and saved models, whose build was never committed. Those pools
    had:
    - wrong-end distances for Framingham's minus-strand genes (2,654 blood
      rows, measured from a point a median 15.8 kb from the Ensembl TSS);
    - non-hg19 TSSs in 156 TCGA rows;
    - one gene under several names, each treated as its own gene;
    - a pseudo-gene, "NAN", lumping 3,935 HELIX rows from 1,461 unnamed
      transcript clusters.

    The calibrators were fitted to models unlike the shipped ones
    (scikit-learn's automatic early stopping), and neither they nor the
    solid model could be regenerated.
  - **What 0.99.7 does.**
    - Each identified name (11,005 of 11,478) resolves to one Ensembl
      GRCh37 gene, at one strand-aware hg19 TSS that the training distances
      and the lookups share.
    - The three models are retrained from committed, deterministic scripts,
      with early stopping off so each calibrator matches its model.
    - The solid model's tier-B calls abstain: they agree 0.41 with the solid
      yardsticks, below chance.
    - The candidate rule, tiers, thresholds and `DIRECT_eQTM` rule are
      unchanged.
  - **Rows.** Blood 1,970,053, nasal 1,500,942, solid 1,184,630.
  - **Validation.** On independent yardsticks with genes matched by Ensembl
    ID, blood and nasal hold, with differences the size refitting alone
    produces. Solid is lower: GTEx solid -0.021 at ladder level. The solid
    model is unstable; over 11 retrainings its tier-A calls agreed 0.60 to
    0.69 with FUSION.

  This is a correction, not a gain in accuracy.
* **`R/rung_accuracy_table.R` is regenerated.** The procedure first
  reproduced the committed table exactly on the 0.99.5 lookups.
  - Adult blood: catalogue_consensus_A 0.779 [0.745] to 0.793 [0.759];
    catalogue_model_blood_A 0.713 to 0.715; catalogue_single_nasal_A 0.629
    to 0.615.
  - Child blood: catalogue_model_nasal_A 0.597 to 0.550. Its yardstick
    leaves out the blood catalogue's `DIRECT_eQTM` pairs, which now include
    more HELIX pairs.
  - Nasal epithelium: catalogue_single_blood_A now has 47 pairs of its own
    (29 before), so it carries its own figure (0.468) instead of the pooled
    one (0.703). It leaves the order: a point below 0.5 never answers.
  - Solid tissue: catalogue_model_solid_A 0.675 [0.652] to 0.639 [0.619];
    catalogue_single_solid_A 0.668 to 0.627.
  - Elsewhere the orders change by local swaps.
* **`cpgd_selftest()`'s reference is regenerated.** 26 of 40 rows change:
  - probabilities move, because the models are retrained, and two calls
    change: cg10073842-MAGEL2 is now -1 and cg10556198-CRHR1-IT1 +1;
  - ten `tss_dist` values move by 1 bp, because EVA-PR's plus-strand TSSs
    were 0-based;
  - one moves by 28 bp: RPL32P1 no longer takes the Zhou TSS, which is
    rounded to 100 bp;
  - cg12793116-ADCY5 is now at -63,243 bp (the strand fix);
  - four CpGs requesting CRHR1/LINC02210 now answer for CRHR1-IT1, LINC02210
    being the same gene;
  - cg22759407-PDE4D is now `UNMAPPED`: PDE4D's only blood training row
    lies 1.26 Mb from its TSS, so PDE4D is no longer a blood candidate.
* **README and vignette.**
  - The model table shows the new AUCs. Blood's saved 0.793 is now 0.848
    because the pseudo-gene is gone; on identified genes it is 0.848 before
    and 0.851 after.
  - The `best_evidence` ranges follow the new table (catalogue_consensus
    0.79-0.83, catalogue_single 0.55-0.90).
  - The tier table gives the current values. Blood tier A had been quoted as
    0.75, solid tier A as 0.71 and blood tier B as 0.89; they are 0.72, 0.64
    and 0.88. Solid tier B now abstains.
  - The solid model's instability is stated.
* **Counts.** `DESCRIPTION`, `CITATION.cff` and `.zenodo.json` now give
  3,322,411 distinct pairs over 759,540 CpGs.
* **New tests.**
  - A gene's names are one candidate (AFDN/MLLT4, CRHR1-IT1/LINC02210,
    DGCR14/ESS2, CECR1/ADA2).
  - Two minus-strand genes, one from GTEx and one from Framingham, sit at
    their 5' end.
* `test-accuracy-values.R`'s comments follow the new figures.
* Requires cpgdirectionData (>= 0.99.7).
* **The catalogue lookups are rebuilt (cpgdirectionData 0.99.5), and the
  ladder's figures with them.** The shipped lookups paired CpGs with genes on
  another chromosome in 450,553 of the 1,408,361 blood rows whose gene
  chromosome could be checked, 288,290 of 829,224 solid and 52,327 of
  1,395,873 nasal. Their build placed each gene at its correct TSS position
  but on a wrong chromosome for 46% of the blood lookup's Atlas genes (54% of
  solid, 7% of nasal), and took the TSS of 934 genes from the Framingham and
  GTEx Atlas exports, which are not hg19 (a median 326 kb off). The rebuild
  reuses the saved models, pools, features, candidate rule and call
  thresholds and changes only the gene placement; with the shipped placement
  it regenerates the shipped tables in every column. Every gene it uses has
  its chromosome confirmed by the EPIC v2 annotation and HGNC or, where
  neither has the gene, by HELIX, the measured eQTMs and its own cis
  training partners; no rebuilt row pairs a CpG with a gene on another
  chromosome. Rows: blood 2,051,017, nasal 1,523,103, solid 1,188,710.
  `R/rung_accuracy_table.R` is regenerated from them: catalogue calls that
  did not change keep their accuracy (0.754 against MESA in blood); the new
  calls, mostly for genes that were misplaced, agree 0.643, so the blood
  figures fall (single blood A 0.751 [0.705] to 0.692 [0.661]; consensus A
  0.800 [0.764] to 0.779 [0.745]) and catalogue_single_blood_A now ranks
  below smr_moderate_finngen in blood. In solid tissue the solid model's
  single tier-A calls fall from 0.778 [0.730] (325 pairs) to 0.668 [0.634]
  (792) and now rank below distance_only; in adult blood they rise from 0.39
  to 0.577 [0.512] and enter the order, below distance_only. The child-blood
  and nasal yardsticks leave out each catalogue's DIRECT_eQTM pairs, which
  changed, so the figures of other sources move there, and with them the
  cross-tissue figures some sources borrow. The rule that orders the ladder
  is unchanged. `cpgd_selftest()`'s reference is regenerated from the
  rebuilt blood lookup, and the accuracy ranges in the README and the
  vignette are those of the new table. `DESCRIPTION`, `CITATION.cff` and
  `.zenodo.json` gave "3,105,697 pairs over 761,445 CpGs", a count no
  current data reproduce; they now give the three rebuilt catalogues'
  3,455,039 distinct pairs over 761,385 CpGs. New test: no lookup row pairs
  a CpG with a gene on another chromosome.
* **The evidence ladder is ranked by accuracy, per tissue.** Every source
  is now scored in each target tissue against a measured yardstick that is
  independent of it (`tools/validate_rung_accuracy.R`, report in
  `tools/reports/validate_rung_accuracy.md`): adult blood against MESA
  monocytes and T cells (Wang et al. 2021; not a catalogue training cohort,
  not an SMR arm), child blood against HELIX pairs outside the blood
  catalogue's training set, nasal epithelium against EVA-PR, solid tissue
  against FUSION and GTEx, each held out in turn. The table
  (`CPGD_RUNG_ACCURACY`, generated into `R/rung_accuracy_table.R`) orders the
  ladder in the call's tissue (`target_tissue`, `cpg_gene_pairs(tissue = )`,
  blood otherwise) and gives every call its figure. Sources are ranked by
  the lower bound of their 95% interval in the tissue -- the number a call
  returns by default -- and by nothing else.
  - Sources are finer than before: measured signs by the tissue they were
    measured in, and replicated only where 2+ independent donor groups agree
    (the strata of EVA-PR, MESA's two cell types and GTEx's tissues are one
    study each; nasal epithelium has one, so no nasal sign is replicated);
    SMR by stratum (a new stratum for eQTLGen pairs whose FinnGen cell-type
    arm disagrees, 0.71 in adult blood against 0.80 for the rest of
    moderate); catalogue single calls by which tissue's model made them. The
    new column `best_rung` names the source that answered; `best_evidence`
    keeps the coarse label.
  - Catalogue sources are scored on their predictions only: a DIRECT_eQTM
    row is a training cohort's measured sign, which the measured rungs
    already carry. Distance is scored only against yardsticks its curves were
    not fitted on (MESA and the solid folds). Measured data are the shipped
    ones (no liver), and a solid pair measured in both FUSION and GTEx counts
    once.
  - The order is the lower bound alone: no tier, stratum or parent
    placement, so a small sample sits where its interval puts it. A source
    whose point estimate in a tissue is below 0.5, or that no yardstick
    scores at all, never answers there. Catalogue conflicts still abstain
    rather than be settled by a source ranked at or below the distance
    curves. Age never changes the order.
  - What changes in adult blood (lower bounds): replicated SMR S1 (0.89)
    now outranks a single measured blood study (0.87), which edges SMR high
    (0.86); the tumour consensus (0.82) and SMR replicated S2 (0.82) move
    above SMR moderate (0.77) and the catalogue consensus (0.76, now formed
    from predicted calls); single calls of the nasal model (0.59) fall below
    the distance curves (0.66) and of the solid-tissue model (point 0.39)
    leave the blood ladder. Directions change where these sources disagree
    on a pair.
  - Every direction has a weight, in every output; only an abstention has
    none. A declared target tissue's last-resort sign is its distance curve's
    sign and carries that curve's validated figure (`best_rung` =
    `distance_curve_<tissue>`); a source no yardstick scores (the
    solid-tissue model's tier-B predictions) carries the neutral 0.5.
  - One catalogue call rule in every result: the pair resolver now applies
    the single-tissue rule (tier A at confidence 0.30, tier B at 0.50,
    DIRECT_eQTM rows alike; before, DIRECT rows always called in pairs), and
    a lone model's call is a single call even with `min_tissues = 1`. In
    both ladders a calling DIRECT_eQTM row is a measured sign of its tissue
    (it answers through the measured rung where the measured layer has no
    sign), a consensus or conflict is between predicted calls on the same
    gene, and models naming different genes are single calls about their own
    genes.
  - The accuracy text of every call is rebuilt from the same table, e.g.
    `0.93 [0.89, 0.95] in adult blood: SMR tier S1, same sign in eQTLGen and
    FinnGen PBMC (vs MESA, n=244)`. The per-rung string constants of 1.0
    (`CPGD_SMR_ACCURACY`, `CPGD_MEASURED_ACCURACY`, `CPGD_CATALOGUE_ACCURACY`,
    ...) are gone. Sources with fewer than 30 pairs in a tissue quote their
    cross-tissue figure; child blood quotes HELIX where HELIX scores the
    source on 100+ pairs.
  - The brain bridge's blood-SMR fallback is ordered by its accuracy in
    cortex: replicated > high > weak > moderate.
  - The per-tissue measured signs behind the ladder are read for the
    queried CpGs only, and grouped without `get()` inside `j`. Grouping the
    whole measured table (93,818 records) that way exhausted 8 GB in
    `cpg_gene_pairs()` when the Matrix package was loaded, as it is in any
    session that loads dmsa. Results are unchanged.
* **One accuracy number per call, for engines that weight by it.** Every
  accuracy text has a numeric twin: `best_accuracy` in
  `cpg_expression_direction()` and `cpg_gene_pairs()`, `accuracy` in the
  single-tissue result and in `cpg_direction_universal()`, and
  `requested_gene_accuracy` in `cpg_brain_bridge()`: the probability that the
  called sign is right, given coupling. New argument
  `accuracy = c("lower", "median")` on the five direction functions:
  - `"lower"` (default): the lower bound of the Wilson 95% interval of the
    source's agreement with the tissue's yardstick (the cluster-bootstrap
    interval for brain T1/T2). The number a call returns is never below 0.5,
    so a small sample cannot turn into a weight that reverses the sign;
    `cpgd_accuracy_values()$lower` shows the bound as computed.
  - `"median"`: the point estimate.
  - The default comes from `getOption("cpgdirection.accuracy", "lower")`, so
    `options(cpgdirection.accuracy = "median")` switches every call in a
    session, including those DMSA makes. With `age_group = "any"` the number
    is the smaller of the adult and child figures. `attr(result, "accuracy")`
    records the choice.
  - `cpgd_accuracy_values()` (new, exported) lists every quoted text with its
    number, n and source; the numbers are looked up by the exact text.
  - For a weighting engine: P(d = +1) = `best_accuracy` where `best_direction`
    is +1 and `1 - best_accuracy` where it is -1; the expected-sign multiplier
    is `best_direction * (2 * best_accuracy - 1)`.
* **Brain T1/T2 accuracy with the same-person overlap removed**
  (`tools/validate_brain_overlap.R`, 4 Nov 2026). Expected agreement against
  an adult cortex cohort sharing no donor with either SMR arm: T1 0.964
  [0.926, 0.998] (in-sample ROSMAP 0.995), T2 0.943 [0.899, 0.984] (in-sample
  0.939: no inflation), estimated without comparing any two donor-sharing
  datasets. `CPGD_BRAIN_ACCURACY` and `requested_gene_best_accuracy` quote
  these figures; the ladder order is unchanged.
* **New argument `age_group = c("adult", "child", "any")`** on
  `cpgd_smr_directions()`, `cpg_expression_direction()`,
  `cpg_direction_all_tissues()`, `cpg_gene_pairs()` and `cpg_brain_bridge()`
  (spec approved 4 Nov 2026). Age moves no direction and no rung -- child
  blood (HELIX, 6-11 y) scores the SMR rungs higher than adult blood (S1
  0.96 vs 0.93, S2 0.96 vs 0.85, high 0.97 vs 0.91, moderate 0.90 vs 0.80,
  weak 0.64 vs 0.58) and 1-25-year cortex
  matches adult cortex at 0.975 -- so the argument changes only which
  validation figures are quoted, which age-specific columns are returned,
  and an optional filter. `"adult"` gives the same output as the default.
  - `"child"`: blood figures in `best_expected_accuracy` quote child blood
    for every source HELIX scores on 100+ pairs outside the catalogue's
    training set -- measured, SMR, tumour and the nasal model's single
    calls (SMR S1 0.963, n = 887; high 0.966, n = 770; S2 0.960, n = 2,188;
    moderate 0.902, n = 2,304; weak 0.635, n = 192) -- and the adult (MESA)
    figure for the rest; the outputs gain `child_blood_helix` /
    `smr_child_blood_helix` (`sig`, `nominal`, `null`, `untested`), and the
    print methods name the population the figures describe. `"any"` quotes
    adult and child side by side.
  - `require_detected = FALSE` (new, child only): when `TRUE`, SMR pairs whose
    `child_blood_helix` is not `sig`/`nominal` are removed before discovery
    and resolution; measured, catalogue, manifest and distance evidence are
    never filtered, and a message reports how many SMR pairs were dropped.
    With any other `age_group` it warns and is ignored.
  - The documentation now says what the HELIX tested-but-null run showed:
    an SMR accuracy is conditional on the pair being coupled in the tissue,
    and most SMR pairs show no detectable coupling in child blood at n = 823
    (detected at p < 0.05: replicated 0.40, high 0.40, moderate 0.20, weak
    0.28, reported 0.21). `null` means not detected, not absent, and never
    changes a rung.
  - Coverage limits are stated: blood in children 6-11 and adults; cortex
    1-25 y and adults; no adolescent blood, no fetal or infant cortex, no
    age data for nasal epithelium or solid tissue. Project Alpha users
    should set `age_group = "child"`.
  - Requires the cpgdirectionData version of `smr_directions` that carries
    `child_blood_helix` (241,044 of 549,316 pairs tested); against an older
    resource the column comes back `NA` with a warning and `require_detected`
    is ignored rather than emptying the SMR layer.

* **Measured brain cortex as the top rung of the brain bridge.** New exported
  accessor `cpgd_brain_measured()` reads `brain_measured_eqtms` (one row per
  pair per cohort: ROSMAP aged DLPFC, 8,815 pairs, Wang et al. 2021 S2;
  Jaffe DLPFC age-adjusted, 383 pairs, lead 18b) and collapses it with one
  vote per cohort, majority sign, ties dropped (9,012 pairs; 174 measured in
  both cohorts, which agree on 174 of 180 shared pairs). `cpg_brain_bridge()`
  gains `measured_brain_genes`, `measured_brain_directions`,
  `measured_brain_n_cohorts` and `measured_brain_cohorts`, a print line for
  them, and -- where the input names a gene -- `requested_gene_best_direction`,
  `requested_gene_best_evidence` and `requested_gene_best_accuracy` from the
  brain-target ladder: measured brain > brain T1 > brain T2 > measured blood
  > blood SMR (ruling 3 transport result, 4 Nov 2026). Blood evidence is a
  labelled fallback (0.84 measured, 0.67-0.77 SMR against adult cortex) and
  never ranks above brain SMR. `requested_gene_direction` and
  `requested_gene_tier` keep their meaning (brain SMR layer only).
* **Brain T1/T2 accuracy, now with the 450K references and their overlap
  stated.** T1 0.995 (n = 1,521) and T2 0.939 (n = 2,904) against ROSMAP
  cortex are in-sample -- ROSMAP is in both SMR arms (Brain-mMeta mQTLs and
  BrainMeta v2 eQTLs) -- and the Jaffe DLPFC figures (T1 0.969, T2 0.918) are
  ROSMAP-free but share their 450K donors with the Brain-mMeta exposure arm.
  GSE15745 (0.92-0.96 / 0.67-0.77) stays the one reference free of both arms
  (`CPGD_BRAIN_ACCURACY`). The distance caveat is dropped (figures unchanged
  within 100 kb).
* **Age caveat of the brain bridge revised** (lead 18b): with age regressed
  out, 1-25-year cortex eQTM signs match adult ROSMAP cortex at 0.975 and Jaffe
  adults at 0.980; the earlier 0.73 developmental figure was age confounding.
  Brain directions are validated in postnatal cortex from 1 year to old age;
  fetal and infant (< 1 y) cortex is untested. The wording follows
  `age_group` and is returned in `attr(result, "brain_caveat")`.
* **SMR documentation.** `replicated` and `high` were described as
  equivalent (they differ by 0.00-0.03 in all five references; since the
  ranked ladder above, the per-tissue validation orders them); the
  child-blood column and the cross-tissue transport table
  (`tools/validate_ladder_transport.R`) are in `?cpgd_smr_directions`.
* **Documentation of the November data build.** `?cpgd_measured_eqtms` names
  the HELIX Dryad (`HELIX_2022_Dryad_adj`, signed, cell-adjusted) and MESA
  (`MESA_monocytes`, `MESA_Tcells`) blood cohorts and the new `source`
  values; `?cpgd_probe_qc` documents the informational T2T-CHM13 / pangenome
  columns (Dong et al. 2025) and why `qc_exclude` is unchanged;
  `inst/scripts/make-extdata.md` covers the build and the new
  `brain_measured_eqtms.csv.gz`.
* Internal: SMR accuracy strings are filled by one helper for both ladders
  (`.cpgd_fill_smr_accuracy()`), and internal callers read the full SMR table
  through `.cpgd_smr_table()`; the exported `cpgd_smr_directions()` is the
  age-aware view and returns a copy, so editing its result never reaches the
  cached table (likewise `cpgd_brain_measured()`). With `require_detected`,
  `smr_in_table` is still judged on the unfiltered table and the print method
  says the pairs were removed rather than that no instrument exists. A column
  subset of a `cpgd_bridge` result now prints as a plain table instead of
  failing.
* New tests: `test-age-group.R`, `test-brain-measured.R`.

* **The SMR rung is rebuilt with outcome-study replication.** The GoDMC
  instruments were run against a second, independent outcome study -- the
  FinnGen immune-multiome single-nucleus cis-eQTLs (PBMC and eight L1 cell
  types; `tools/finngen_smr_celltype.R`, `tools/godmc_dense_besd.R`), SMR
  1.03 per chromosome with the dense GoDMC cis mQTL set as exposure -- and
  `smr_directions` now carries, per pair, `outcome`, `replication`, the
  FinnGen side (`finngen_direction`, `finngen_tier`, `finngen_p_SMR`,
  `finngen_heidi_status`, `finngen_contexts`, `finngen_evidence`) and
  `smr_rung`, which both ladders act on (`tools/build_smr_replication.R`).
  Validation against the measured blood eQTMs
  (`tools/validate_celltype_smr.R`, `tools/reports/validate_celltype_smr.md`,
  `tools/reports/build_smr_replication.md`, 2-3 Nov 2026):

  - New evidence level **`smr_replicated`**, between `measured` and
    `smr_high`: eQTLGen S1/S2 pairs whose sign FinnGen PBMC reproduces --
    S1 0.986 (n = 1,795), S2 0.971 (n = 3,387). 68,536 pairs over 41,759 CpGs
    move from the 0.85 rung to a 0.97 rung.
  - eQTLGen pairs whose sign FinnGen PBMC **contradicts**: S1 drop to `smr_weak`
    with the eQTLGen sign kept (0.92, n = 111 -- the larger study is usually
    right, but n is below the 200-pair bar for a rung); S2 and S3 become
    `reported` with no sign called (0.59, n = 302 and 0.38, n = 68, neither
    beating its majority-class baseline). Conflicting evidence abstains rather
    than guesses, as ACMG/AMP treats conflicting criteria.
  - eQTLGen S3 pairs (instruments disagree) whose majority sign FinnGen PBMC
    **reproduces** rise from `smr_weak` to `smr_moderate` (0.855, n = 214;
    Wilson 0.80-0.90), with their own accuracy string. Unreplicated S3 stays
    `smr_weak` (0.70, n = 227).
  - Pairs the eQTLGen layer lacks enter from FinnGen PBMC at
    `smr_moderate` when the SMR test at the pair's top cis-mQTL reaches
    p < 5e-8 (sign = sign(b_SMR); 0.904, n = 1,798, Wilson 0.89-0.92 -- above
    the rung's historical 0.84-0.86). That is 78,228 PBMC-only pairs admitted
    (39,442 with both an SMR test and an instrument vote, 38,786 SMR-tested
    only). `finngen_evidence` records
    whether a row's sign is the SMR estimate at the top cis-mQTL
    (`smr_top_snp`), the vote of the pair's instruments
    (`instrument_vote`), or both.
  - FinnGen PBMC pairs below that threshold and every cell-type-only pair
    (0.744, n = 519) are `reported`: visible in the `smr_*` columns, never
    `best_direction`. The table has 549,316 rows: 413,982 eQTLGen + 135,334
    FinnGen-only (rungs: moderate 375,385, replicated 68,536, reported 57,106,
    weak 28,201, high 20,088).
  - HEIDI conditions nothing. Its v1 role as a gate on the FinnGen-only
    admission is withdrawn: that separation rested on 110 pairs, and the
    top-cis-mQTL SMR threshold reaches 0.90 on 1,798 without it (HEIDI-passing
    admitted pairs: 0.95-0.97 on n = 150, noted, below the rung bar). HEIDI status
    on FinnGen rows is still shipped, recomputed at the table's contract
    (pass = p >= 0.05).
  - Cell types agree with PBMC on 99% of shared pairs and PBMC is right
    81-94% of the time when they differ; a cell type therefore never
    overrides PBMC and contributes coverage only (`finngen_contexts`).
    `blood_other` is pooled into the union.
  - `CPGD_SMR_ACCURACY` (R/utils.R) holds the rung strings;
    `.cpgd_smr_rung()` maps tier -> rung for tables built before this
    version, so a stale cached resource keeps the old behaviour.
  - New Hub version of `smr_directions` in cpgdirectionData (row count from
    the v2 build report).

* **Undermining record** (Undermind workspace "cpgdirection monthly horizon
  scan", 2 Nov 2026): no MR/SMR guideline treats cross-outcome concordance as
  a formal tier, but false-sign work on TWAS (Gerlach et al. 2025) and the
  MR+colocalisation review (Zuber et al. 2022) support demoting
  non-replicating signs rather than discarding them; Wallace et al. 2025
  explain why HEIDI-type heterogeneity tests reject true colocalisations
  more often as sample size grows, which is why the package never conditions
  on HEIDI; ClinGen practice (Pejaver et al. 2022, Brnich et al. 2019)
  supports stating each rung's validated n and lower confidence bound, which
  the accuracy strings do.
# cpgdirection 0.99.0

**Bioconductor submission numbering.** This version continues 2.99.8 (the last
GitHub-numbered release) and carries the changes requested in the Bioconductor
review (BiocContributions #174). The number drops to 0.99.0 because that is
what Bioconductor requires of a new package; it becomes 1.0.0 at the next
Bioconductor release. The code and data are those of 2.99.8 plus the items
below.

* **Dependencies.** `TxDb.Hsapiens.UCSC.hg19.knownGene`, `org.Hs.eg.db`,
  `GenomicFeatures`, `GenomicRanges` and `AnnotationDbi` move from Imports to
  Suggests: only the distance-only layer uses them, it checks for them with
  `requireNamespace()` and says what to install. `cpgdirectionData` stays in
  Imports and is now called directly, without a `requireNamespace()` guard.
  `tools` is declared for `tools::R_user_dir()`.

* **Nothing is written to the working directory by default.** `cpgd_report()`
  defaults to `dir = tempdir()`; `html_out = TRUE` writes the report into
  `tempdir()` and reports the path; `cpgd_build_gene_tss()` defaults to
  `out = NULL`, returning the table and writing only when a path is given
  (`cpgd_gene_tss()` keeps its cache under `tools::R_user_dir()`).

* **`message()` instead of `cat()`** outside print methods: `cpgd_selftest()`,
  `cpgd_direction_balance()`, `cpg_direction_universal()` and the `verbose`
  summaries of `cpg_expression_direction()`, `cpg_gene_pairs()` and the full
  ladder now report through `message()`, so `suppressMessages()` silences them
  and stdout carries only return values.

* **hg38 manifests are refused.** `cpgd_build_gene_tss(source = "manifest")`
  uses hg19-labelled columns when a re-annotated manifest carries them,
  refuses a manifest whose columns or `Genome_Build` say hg38 (EPIC v2's own
  manifest included) with a message naming the alternatives, and otherwise says
  that the coordinates are taken as hg19. Every layer in the package is hg19.

* **Vignette** rebuilt in the Bioconductor shape (`BiocStyle::html_document`,
  table of contents): an Introduction placing the package in the ecosystem
  (versus ELMER and MethReg, which infer methylation-expression relationships
  from a user's own matched data), an Installation section, the worked
  examples, and session information.

* **README** describes the current installation path (Bioconductor, or GitHub
  with the data package) and the ExperimentHub data layers; the local-tarball
  instructions and the bundled-data description of the old fat build are gone.

* **`inst/scripts/`** documents how each file in `inst/extdata` was made and
  from which sources (`make-extdata.md`); the measured-table builder moves
  there from `tools/`.

* **Code provenance** is stated in `DESCRIPTION`, the README and the package
  help page: the design, data assembly, validation and all methodological
  decisions are the author's; the R code was drafted and edited with the
  assistance of Claude (Anthropic) and was run, tested and checked by the
  author (Assisted-by: Claude (Anthropic)). `PUSH_INSTRUCTIONS.md`, a note
  from the preparation of the first version, is removed from the repository.

# cpgdirection 2.99.8

* **The measured layer is rebuilt from the eQTM Atlas, with cohort provenance.**
  `inst/scripts/build_measured_eqtms.R` (new) reads the per-cohort exports of the
  eQTM Atlas (Sriram et al. 2026; hg19) and keeps cis pairs (|CpG - TSS| <=
  250 kb, the window of the earlier build; decided from coordinates because
  the Atlas labels every FUSION row "trans") with a usable sign that pass each
  cohort's own threshold. Eighteen cohort files were read: EVA-PR nasal
  epithelium (four subsets) and white blood cells, GTEx whole blood and seven
  organs, FUSION muscle, Colonomics colon, Framingham whole blood. HELIX and
  EvoImmunoPop carry no sign in the export and the Framingham export is a
  10,000-row cap, so the earlier build's Framingham + HELIX blood pairs are
  kept as `cohorts = "packaged_v1"`. The measured table grows from 26,633 to
  ~75,000 pairs (blood ~26,500, nasal ~11,000, solid ~37,500) and gains
  `n_cohorts`, `cohort_agreement`, `cohorts` and `source`.

  - **Cross-cohort sign agreement** on shared pairs is 0.89-1.00 for every
    cohort pair with >= 20 shared pairs, within and across tissues (EVA-PR WBC
    vs Framingham 0.992 on 2,869; FUSION vs GTEx muscle 0.993 on 299; GTEx lung
    vs ovary 0.985 on 2,227). The measured rung's expected accuracy is now
    stated from this replication rather than assumed: `~1.00` where 2+ cohorts
    agree, `0.96-0.99` for a single-cohort pair (`CPGD_MEASURED_ACCURACY`).
    Both ladders report `measured_n_cohorts` and `measured_cohort_agreement`.
  - **Tumour pairs leave the measured rung.** The 1,704 solid-tissue pairs of
    the earlier build are tumour-derived (the solid catalogue is documented as
    tumour tissue) and agree with the normal-tissue cohorts on only 0.58 of
    shared pairs; tumour evidence is the Onco-eQTM layer. `solid_tissue` now
    means normal muscle, lung, ovary, colon, breast, prostate and testis,
    pooled for the direction call with the organ kept in the long table.
  - The far-cis pairs the Atlas calls "cis" up to 1 Mb are not ingested: beyond
    100 kb they are majority-positive and single-cohort, a different population
    from the near-TSS pairs the rung was built on.

* **distal_links gains scE2G cell-type links.** `tools/build_distal_links.R`
  reads the scE2G_multiome "pairs" files served by the Broad
  (`data.broadinstitute.org/scE2Gpredictions/PBMC_9_Cell_Types/scE2G/`),
  streaming each 650-770 MB file through `gzip | awk` so only links at or above
  the model's calling threshold (`E2G.Score.qnorm >= 0.177`,
  `--sce2g-min-score`) reach R, and names the biosample from the file's
  `CellType` column (PBMC/BMMC clusters are `tissue_class = "blood"`). The nine
  PBMC cell types (cDC, CD14 monocytes x2, CD16 monocytes, B, naive and memory
  CD4, central-memory CD8, NK) join the 26 bulk ENCODE-rE2G biosamples, so a
  blood link can be read per cell type.
  The rebuilt layer holds 10,205,122 rows, 1,289,742 CpG-gene pairs and
  525,222 CpGs; scE2G contributes 2,168,868 rows (225-263 thousand per cell
  type; scores 0.177-1.00, median 0.78).

* **The Onco-eQTM rung re-validated against the rebuilt measured layer.** With
  tumour pairs gone, the solid reference is normal tissue from non-TCGA donors
  and becomes a fourth independent check: O1 agrees with it on 0.874 of 1,976
  pairs (baseline 0.672), with blood 0.830 (n = 1,955) and nasal 0.937
  (n = 1,500) unchanged. `CPGD_ONCO_ACCURACY` now quotes all four references.

# cpgdirection 2.99.7

* **The tumour consensus earns a rung.** `tools/validate_onco.R` (new, with
  `tools/validate_distal.R`) scores every Onco-eQTM tier against the package's
  independent peripheral references. Tier O1 (>= 5 cancer types, one sign)
  agrees with measured blood eQTMs 0.82 (n = 1,455; majority baseline 0.72),
  with nasal epithelium 0.93 (n = 1,221; baseline 0.71) and with SMR tier S1
  0.70 (n = 1,033; baseline 0.60). Near the TSS the agreement is mostly the
  shared -1 base rate (0.88 vs 0.87 at <= 1.5 kb); the information is in
  gene-body and 1.5-100 kb pairs (blood 0.76-0.80 vs baselines 0.56-0.58,
  nasal 0.94 vs 0.56-0.60), where the layer also calls +1 at the reference
  rate and is right 72-87% of the time -- the one thing no current rung does.
  Against SMR S2 the agreement is 0.61 (n = 8,499; baseline 0.54) and against S3 0.47 (n = 429; baseline 0.53), consistent with the placement above `smr_weak`. O2 is a few points above baseline, O3/O4 are at chance. The solid-tissue
  reference (0.97) is excluded from the decision because it may share TCGA
  material with Onco-eQTM.

  - New evidence level `onco_consensus`, between `catalogue_single` and
    `smr_weak`, supplied by tier O1 only, in both `cpg_gene_pairs()` and
    `cpg_expression_direction()`; `direction_tier` "O1"; expected accuracy
    "0.70-0.93 (Onco-eQTM O1 ...)" with the three validated figures.
  - The tumour layer is now reported on every run (`onco_direction`,
    `onco_tier`, `onco_n_cancers`, `onco_sign_agreement`, `onco_median_r`,
    `onco_agreement`); `include_onco = TRUE` in `cpg_gene_pairs()` only
    governs whether tumour-only TARGETS enter the candidate set.
    `cpg_expression_direction()` keys and gates the layer exactly like SMR
    (`onco_gene`, `onco_gene_match`; the strongest tumour pair is named even
    when no gene matched, and barred from `best_direction`).
  - `onco_agreement` is `NA` where the tumour sign IS the best direction:
    agreement is a statement about two independent signs.
  - `--append-solid` in `tools/build_onco_eqtm.R` is answered: no. The rung is
    the right home for this evidence; the measured layer is not.

* **distal_links validated as a target source, not a direction source.** On
  the blood CpGs both sources cover, rE2G recall exceeds the manifest's (0.46
  vs 0.43) and the high-confidence subset (score >= 0.9, or >= 10 biosamples)
  is more precise than the manifest (0.52 vs 0.43); the 1,567 blood pairs the
  links find that the manifest lacks sit at a median 30 kb (p90 115 kb). The
  sign check came out as it should: non-promoter links at the base rate,
  promoter-class links negative -- the TSS-distance effect already on the
  ladder. No direction rung; `include = "distal"` stays opt-in.

# cpgdirection 2.99.6

* **Two new opt-in evidence layers, both kept OUT of `best_direction`.** The
  October 2026 horizon scan turned up resources that propose targets or carry
  signs from contexts the peripheral ladder does not answer for. They enter as
  reported columns, with their own provenance, so that agreement or
  disagreement with the peripheral call is visible instead of silently
  resolved by precedence -- the same rule brain SMR has followed since 2.3.

  - `cpgd_distal_links()` and `cpg_gene_pairs(include = "distal")`: CpG-to-gene
    links from ENCODE-rE2G enhancer-gene predictions (Gschwind et al. 2026,
    Nature), intersected with each probe's hg38 coordinate and keyed on the
    cg ID. The shipped layer (Hub resource `distal_links`) holds the
    thresholded rE2G tables of 47 released ENCODE annotation sets covering 26
    resting peripheral-blood cell types (monocytes, B, CD4, CD8, naive and
    memory T, regulatory T, NK, dendritic cells, GM12878, Th1/Th17) and brain
    regions and cell types (dorsolateral prefrontal and frontal cortex,
    caudate, putamen, temporal and frontal gyri, bulk brain, astrocytes,
    glutamatergic neurons): 8,036,254 rows, 1,149,157 CpG-gene pairs,
    514,549 CpGs, 19,551 genes; rE2G scores 0.20-1.00. Activated and
    stimulated cell types are excluded because the eQTM reference data are
    unstimulated blood. ENCODE ships an ABC-model table beside every rE2G
    table; the builder recognises its header and skips it (`--keep-abc`
    keeps it as `source = "ABC"`) so one score scale applies. A distal link
    is a TARGET proposal with no sign: `mapping_source = "distal_link"`
    (priority below every direction-bearing source and below the manifest),
    `has_distal`, `distal_sources`, `distal_n_biosamples`,
    `distal_score_max`, `distal_biosamples`. It is not counted in
    `n_annotation_sources`, which tallies positional tracks. A blood
    question consults blood-class biosamples only. Built by
    `tools/build_distal_links.R` from files fetched with
    `tools/fetch_encode_re2g.py` and mapped by `tools/distal_biosamples.csv`.
    The builder also accepts scE2G single-cell maps (Sheth et al. 2026,
    Nature Genetics) and Human Cell Epigenome Atlas loops (Zhou et al. 2026,
    Science); neither is in this release's table.
  - `cpgd_onco_eqtm()` and `cpg_gene_pairs(include_onco = TRUE)`: the
    Onco-eQTM cross-tumour consensus (Korra, Nishana and Kumar 2026, NAR
    Genomics and Bioinformatics; 27 TCGA cancer types, 6,880 tumours, 450K,
    hg19; CC BY-NC 4.0), pooled per CpG-gene pair from the per-cancer cis
    tables of the Onco-eQTM download centre (5,233,328 CpG-gene-cancer rows,
    r recovered from the t statistic and the per-cancer n). 2,464,990 pairs
    with tiers O1 (>= 5 cancers, one sign; 113,560), O2 (2-4 cancers, one
    sign; 566,351), O3 (cancers disagree; 549,150) and O4 (one cancer;
    1,235,929). Sign sanity check: promoter pairs (|TSS distance| < 1.5 kb)
    are 86% negative, gene-body pairs 39% negative -- the textbook pattern.
    Reported as `onco_direction`, `onco_tier`, `onco_n_cancers`,
    `onco_sign_agreement` and `onco_agreement` (does the tumour sign match
    the peripheral `best_direction` of the same pair?); `mapping_source =
    "tumour_eQTM"` ranks below every peripheral source. Tumour evidence
    never enters `best_direction`. Built by `tools/build_onco_eqtm.R`; its
    `--append-solid` option is deliberately opt-in and labelled: a measured
    record is the top rung and tissue-blind, and 2.5 million tumour pairs
    would otherwise outrank blood SMR for a saliva study.

* **Exeter re-annotation v3.0 (GENCODE 49).** `tools/build_exeter_v3_annotation.R`
  swaps the Exeter gene track of the packaged EPIC v2 annotation union from
  v2.0 (GENCODE 47, Zenodo 15181885) to v3.0 (GENCODE 49, Zenodo 20704849,
  June 2026), leaves the other tracks byte-identical, recomputes
  `annotation_source`/`n_annotation_sources`, and writes a diff report
  (`tools/reports/exeter_v3_diff.md`). The rebuilt table has 1,006,212
  pairs over 793,233 CpGs (was 996,697 / 791,550): 9,647 pairs gained, 132
  lost, 14,214 TSS200/TSS1500 relabels. v3.0 drops the Promoter_2000bp /
  Enhancer_5000bp columns, so those labels are lifted out of the v2.0
  feature strings into a dedicated `exeter_regulatory` column (416,200
  pairs) and the regulatory track is kept. The gene track is named
  `Exeter_GENCODEv49`; no package code keys on the old name.

* New Hub resources registered in cpgdirectionData 0.99.2: `distal_links`,
  `onco_eqtm_consensus`. Installations without them behave exactly as before
  (the accessors return `NULL` with a message; the defaults of
  `cpg_gene_pairs()` do not touch them).

# cpgdirection 2.99.5

* **Build fix (Bioconductor #174).** Every platform failed `R CMD check` on
  2.99.4 because `data.table::fread()` on a `.csv.gz` path silently requires
  the `R.utils` package, which the Bioconductor builders do not have (and
  this package never declared). All compressed reads now go through one
  internal helper (`.cpgd_fread()`) that decompresses with base R
  (`gzfile()`) before parsing; the tables it returns are identical to the
  previous ones. No new dependency. A test that read its own `.gz` output
  with `fread()` now uses `read.csv()`.

# cpgdirection 2.6.0

* **Bioconductor-ready data backend.** All data layers now resolve through
  one internal backend (`R/data_backend.R`) with three tiers, tried in
  order: `options(cpgdirection.data_dir=)` (an explicit local directory, for
  air-gapped machines and pinned analyses), the package's own bundled
  `inst/extdata` (the full GitHub build -- unchanged behaviour, no network
  ever needed), and the new `cpgdirectionData` ExperimentHub companion
  package (the thin Bioconductor build). The fat and thin builds run
  identical code and return identical tables; only the transport differs.
  `tools/make_thin_build.R` produces the submission tarball (extdata drops
  from ~100 MB to 4.4 MB; the eight large layers move to the Hub).
* `cpgd_lookup_info()$file` now reports the actual data source
  (`"cpgdirectionData (ExperimentHub)"` on a thin install).
* Reporting polish: printing a column subset of a `cpgd_pairs` table shows a
  plain table instead of a summary with misleading zeros, and in
  `gene_mode = "filter"` a CpG whose discovered targets were all removed by
  the gene filter is no longer miscounted as "no target found" -- discovery
  is judged before the filter.

# cpgdirection 2.5.0

* **The mapping database is now multi-source, cross-validated, and
  probe-QC-aware.** A CpG assigned to the wrong gene poisons every downstream
  direction, so the annotation attacks mis-assignment from several independent
  angles at once.

  - The packaged EPIC v2 probe-to-gene table is rebuilt as the UNION of five
    annotation tracks -- Illumina UCSC_RefGene, Illumina GencodeV41, the Zhou
    lab (SeSAMe) transcript-aware GENCODE v41 annotation, the Exeter GENCODE 47
    re-annotation (Mallabar-Rimmer et al. 2025, Zenodo 15181885), and the
    Exeter distance-based regulatory assignments (Promoter_2000bp /
    Enhancer_5000bp, with an In_GeneHancer flag). Coverage grows from 349,835
    to 791,550 CpGs with at least one gene (996,697 CpG-gene pairs, 38,950
    genes). Every pair carries per-track support flags and
    `n_annotation_sources`; 84.4% of pairs are corroborated by >= 2
    independent tracks and 54.8% by >= 3. `zhou_dist_tss` records the signed
    distance to the nearest annotated TSS of the assigned gene.
  - NEW probe QC resource and accessor `cpgd_probe_qc()`, from three
    independent angles: the Zhou general mask (degenerate/non-unique mapping,
    SNP-contaminated extension base) collapsed to CpG level over replicate
    probes; the Garvan experimental evidence (Peters et al. 2024, BMC
    Genomics) -- WGBS-confirmed cross-hybridization, BLAT-predicted
    off-targets, mapping mismatches; and an independent verification of the
    packaged hg19 coordinates against the Zhou hg19 manifest (99.5% exact;
    3,632 CpGs flagged unverified). The two cross-hybridization angles
    corroborate each other: 7,535 of the 11,754 WGBS-confirmed off-target
    CpGs are also sequence-masked by Zhou, and each catches probes the other
    misses.
  - `cpg_gene_pairs()` gains `probe_qc = c("exclude", "flag", "ignore")`.
    Default `"exclude"`: the 35,657 CpGs (3.85%) whose EVERY replicate probe
    is Zhou-masked or WGBS-confirmed cross-hybridizing (`qc_exclude`) are
    dropped -- a probe reading the wrong locus is the one failure no
    direction evidence can repair -- with the exclusions reported, listed in
    `attr(result, "qc_excluded_cpgs")`, and recoverable via
    `probe_qc = "flag"`. A CpG with at least one clean replicate is kept and
    flagged `probe_masked_partial` instead. Pair rows carry `probe_masked`,
    `probe_mask_reasons`, `cross_hybridizing`, `mapping_flagged` and
    `pos_hg19_verified`.

* Annotation provenance and probe reliability remain SEPARATE facts from
  direction evidence: `annotation_source`/`n_annotation_sources` say who
  assigns the gene, `cpgd_probe_qc()` says whether the probe can be trusted,
  and `best_evidence` says what the direction rests on.

# cpgdirection 2.4.0

* **Automatic CpG-to-gene discovery and a pair-preserving API.** A bare EPIC
  v2 CpG identifier is now sufficient: `cpg_gene_pairs("cg26261055")`
  discovers every supported target gene automatically and returns a separate
  direction/evidence record for every CpG x gene pair.

  - `cpg_gene_pairs()` -- the new many-to-many interface. Candidate targets
    are the UNION of the packaged EPIC v2 manifest annotation, the tissue
    catalogue, measured eQTMs and SMR (brain targets on request via
    `include_brain = TRUE`; genes parsed from input suffixes augment the set
    by default). Direction is resolved WITHIN each `cpg_id + target_gene`
    pair by the same evidence ladder as `cpg_expression_direction()`
    (`.cpgd_resolve_pair_direction()`), so a fixed pair returns the same
    answer whether the gene came from the manifest, an input suffix, or
    `genes=`. Opposite directions across a CpG's genes are retained;
    abstention is pair-specific.
  - `genes = NULL` means auto-discover, never "cannot calculate direction".
    `genes` with `gene_mode = "filter"` (default) restricts discovered
    targets; `gene_mode = "pairwise"` evaluates explicit `cpgs[i] x genes[i]`
    hypotheses. The mode is never inferred from vector lengths.
  - `cpgd_manifest_genes()` -- new packaged long-form EPIC v2 probe-to-gene
    annotation (402,098 CpG-gene rows, 349,835 CpGs, 29,175 genes) built from
    the official Illumina EPIC-8v2-0 manifest RefGene fields; every
    semicolon-delimited gene becomes its own row, never just the first.
    `cpgd_build_manifest_genes()` rebuilds it from a manifest CSV. Manifest
    annotation is a mapping source, not proof of regulation, and its
    provenance stays distinct from measured/SMR evidence.
  - Mapping provenance (`mapping_sources`, `mapping_primary`,
    `mapping_strength`, `has_*`) is kept separate from direction evidence
    (`best_evidence`, `direction_tier`) on every row.

* The input parser gained `input_id` and a `dedupe` argument. The legacy
  one-row-per-CpG interfaces still collapse technical replicates by `cpg_id`;
  the pair workflow parses every row and never deduplicates by `cpg_id`
  alone, so `cg123_GENE_A` and `cg123_GENE_B` both survive.

* `cpg_expression_direction()` is unchanged and remains the convenience
  one-result-per-CpG interface; its documentation now points many-to-many
  consumers (DMSA included) at `cpg_gene_pairs()`.

# cpgdirection 2.3.0

* **The brain bridge.** Three new data layers and one user-facing function
  answer the question most users actually have: what does a saliva or blood
  measurement license about the brain?

  - `cpgd_brain_directions()` — 220,383 brain CpG-gene pairs with SMR causal
    directions (Brain-mMeta x BrainMeta v2), tiered by a score validated
    out-of-sample in blood (T1 96.3%, T2 83.0%) and against the independent
    Gibbs brain catalogue (T1 92-96% across every rung).
  - `cpgd_bridge_deliverable()` — the graded 728-CpG intersection of a T1 brain
    direction with robust peripheral-brain concordance; grade A++ (205 CpGs,
    184 saliva-usable) is supported by four independent cohorts.
  - `cpg_brain_bridge(my_cpgs, tissue = "saliva")` — per CpG: brain direction
    and target(s), concordance evidence for YOUR tissue, deliverable grade,
    genome-wide synthesised bridge score, and strict-QC flags.

* The print method repeats three caveats on purpose: directions are
  brain-derived (peripheral eQTM directions agree with brain only ~65%);
  concordance couplings rest on n = 27 living-brain donors; sex-chromosome and
  discontinued probes produce replicable artefactual concordance and are
  guarded by `strict_qc`.

* Licence note: shipped bridge data derive from open resources (IMAGE-CpG,
  AMAZE, GSE59685, Sommerer, Gibbs). Components derived from CC-BY-NC sources
  (the pediatric coupling table) are deliberately excluded from the package and
  live in the OSF deposit.

* **Co-regulation columns, two evidence classes kept apart.** `co_up`/`co_down`
  list the T1 partner genes a CpG's methylation moves in each direction —
  causal co-direction from the brain layer itself, with bivalent CpGs visibly
  split rather than merged. `hic_genes` lists genes contacted by the same
  adult-cortex enhancer (PsychENCODE; 23,194 CpGs): candidate co-regulation by
  physical contact, no direction attached. TF-regulon co-activation is
  deliberately NOT shipped: the motif-to-polarity chain performed at its
  majority baseline in brain and does not meet the package's evidence bar.

* Gene-aware input: a gene named in the identifier or a data.frame column is
  echoed back with its own answer (`requested_gene`, `requested_gene_direction`,
  `requested_gene_tier`) — `NA` when the gene is not a partner of that CpG,
  never a silently substituted neighbour.

* +12.6 MB of extdata; tests pin the manuscript's headline counts (728 / 205 /
  184), the POMC example, all input formats, and the co_up/co_down sign split,
  so silent data drift fails the build.

# cpgdirection 2.2.5

* **HEIDI results are now shipped.** `smr_directions.csv.gz` gains `p_HEIDI`,
  `nsnp_HEIDI` and `heidi_status`; `cpg_expression_direction()` surfaces the
  latter two as `smr_heidi_status` and `smr_p_heidi`. Row count is unchanged at
  413,982.

* **Nothing conditions on them, deliberately.** Across all 413,982 pairs, 5.4%
  passed HEIDI, 78.1% failed, and 16.4% had fewer than three instruments after
  LD pruning and could not be tested. Against the directly measured eQTMs,
  passing pairs were concordant 80.8% of the time and failing pairs 87.3%
  (p = 0.002) — the opposite of what a useful filter looks like.

  The cause is power, not biology. HEIDI's ability to reject rises with the
  number of instruments, and so does the precision of the Wald ratio.
  Concordance by instruments surviving pruning runs 69.9% (<3), 74.8% (3-5),
  74.0% (6-10), 88.0% (11-20), and mean instrument count is 17.3 among failing
  pairs against 12.0 among passing ones. With both QTL resources above 30,000
  individuals, HEIDI detects heterogeneity far below the level that would flip a
  sign, so filtering on it would discard the best-instrumented and most accurate
  pairs first.

  `not_tested` is kept strictly apart from `pass` throughout. Recording an
  untested pair as passing would convert absence of evidence into evidence,
  which is the error this whole exercise exists to avoid.

* The tier structure is unchanged and S1 > S2 > S3 holds within every HEIDI
  stratum.

* `cpgd_smr_directions()` reads the new columns defensively, so a stale cached
  table yields `NA` rather than an error.

* **`cpgd_report()` was dropping the columns that make an SMR row readable.**
  Its compact column list is hand-written, so `smr_gene`, `smr_gene_match` and
  `smr_heidi_status` never reached the HTML — the output most users actually
  look at. On a real panel about three quarters of SMR rows concern a different
  gene from the rest of the row, so a report showing `smr_direction` without
  `smr_gene` presents those directions as though they were about the requested
  gene. All three are now included, and a test renders a report and greps for
  them.

# cpgdirection 2.2.4

* **The SMR layer could never fire unless a gene was named.** Its join keyed on
  `requested_gene` and nothing else, so calling
  `cpg_expression_direction(my_cpgs)` — the mode most people use — left the key
  `NA` and the 413,982-pair table contributed nothing. No error, no warning, and
  the printed summary blamed 450K array coverage for the silence.

  The layer is now always reported, with two new columns deciding what may be
  done with it. `smr_gene` names the gene the causal estimate concerns.
  `smr_gene_match` says whether that is the same gene as the rest of the row.

  This matters more than it sounds, because the two layers usually disagree about
  which gene a CpG regulates: of the 85,220 CpGs present in both the SMR table
  and a catalogue, only 58.6% share even one gene, and since the catalogue commits
  to a single target the realised match rate is lower still. `cg00000029` is
  typical — SMR reaches AKTIP and RBL2, the blood catalogue reaches SCP2.

  So the rule is: report always, promote only on match. A non-matching row shows
  its direction, tier and p-value in the `smr_*` columns but cannot supply
  `best_direction`, because that column answers a question about one specific
  pair. `smr_agreement` is `NA` rather than `FALSE` for those rows: two statements
  about different genes cannot contradict each other, and calling it disagreement
  would manufacture a conflict out of compatible facts.

  Where a gene *is* supplied, catalogue-derived fallback keys are not used at all.
  If you ask for CREBBP, a direction for ADCY9 is not an answer to your question.

  No validation figure changes: S1 95.9% / S2 84.9% / S3 70.4% were measured
  against the table directly. What changes is that they now reach the caller.

* **Two-token identifiers lost their gene.** `.cpgd_parse_input()` required more
  than two underscore-separated fields and discarded the first two by position,
  so `"cg00039463_TRAP1"` — the plainest form anyone would type — parsed to a CpG
  with no gene, and the distance layer reported itself unavailable. Gene
  candidates are now taken from whatever follows the token carrying the `cg`
  identifier, with probe-type codes (`TC21`, `BC11`) dropped by name rather than
  by position.

* **The "no SMR evidence" message asserted a cause it had not checked.** It
  attributed every empty result to GoDMC's 450K coverage, including for CpGs
  sitting in the table with four significant rows. It now distinguishes a CpG
  absent from the table from one present but unmatched on gene, and says which.
  New column `smr_in_table` exposes the same distinction.

* Regression tests in `test-smr-layer.R` covering all three. Each would have
  failed on 2.2.3, and none of the existing tests could have caught them: every
  bug returned a plausible answer rather than an error.

* `DESCRIPTION` coverage figures corrected. They still described a single
  gradient-boosted table of 1,603,777 pairs over 604,456 CpGs, which predated the
  nasal and SMR layers. The shipped data are 3,105,697 pairs over 761,445 CpGs
  across four evidence sources. Analysis code and fitted models are archived at
  <https://doi.org/10.17605/OSF.IO/U3VFK>.

# cpgdirection 2.2.3

* `UNMAPPED` was unreachable. An identifier absent from the lookup table came
  back as `NO_EXPRESSION_ANCHOR` instead, because the right join fills
  non-matching rows with `NA` and the code then read a missing `target_gene` as
  "no gene nearby". The two are different facts and a user needs to tell them
  apart: an absent identifier usually means a typo or a probe from a different
  array, while no expression anchor is a real biological answer about a real
  CpG. The distinction is now captured before `status` is overwritten, and
  `UNMAPPED` rows carry a note saying to check the probe ID and array version.

  Found by a test that had been asserting the documented behaviour since 1.0.0
  and had never once been run.

* No change to any prediction. Affects only how absent identifiers are labelled.

# cpgdirection 2.2.2

* Fixes `could not find function "fifelse"`, introduced by me in 2.2.1. Deleting
  the hand-written `NAMESPACE` also deleted its
  `importFrom(data.table, ...)` list, and the replacement declared only `:=`
  and `.SD`. Three functions are called bare inside `[.data.table` expressions,
  where a `::` prefix would be evaluated in the wrong frame: `fifelse`, `fcase`
  and `setorderv`. All three are now declared.

  The scan I used to find them the first time had a hole: it excluded any name
  used both bare and qualified, which is exactly the case here, since
  `data.table::fifelse` appears elsewhere in the package. Corrected and re-run
  against the full data.table export list.

* `cpgd_lookup_info()` gains the `@param tissue` documentation it never had,
  clearing the last `R CMD check` warning.

* No change to any prediction.

# cpgdirection 2.2.1

The first `R CMD check` the package has ever had. Two errors, two warnings and
three notes, all of them real, none of them affecting a prediction.

* **Documentation is now generated by roxygen rather than written by hand.** The
  hand-written `man/*.Rd` and `NAMESPACE` had drifted from the code -- check
  reported codoc mismatches where `cpgd_build_gene_tss()` and `cpgd_selftest()`
  were documented as taking no arguments -- and roxygen refused to touch files it
  had not generated, so every `devtools::document()` was a no-op. Both are
  deleted and regenerated. A new `R/cpgdirection-package.R` carries the imports.

* **Four failing tests, one cause.** They predated the multi-layer default and
  asserted on `gene_source`, `status` and `call`, which exist only when a tissue
  is named; the default returns the `best_*` column set. They now ask for
  `tissue = "blood"` where they mean the per-tissue shape. Two tests added: one
  checking the default column set, one checking that `target_tissue` fills
  blanks without overwriting existing calls.

* **A test of mine was nonsense** -- it called `readLines()` on a directory. Its
  intent, that every label the ladder emits is a known level, is now tested
  directly, along with the ladder ordering.

* **The example errored** because it selected `target_gene` from a result that
  has no such column in the default mode. It now shows both shapes.

* Fixed: non-ASCII characters in `R/report.R` code (permitted in comments, not
  in code); the `cpgd_report()` example writing into the check directory;
  `PUSH_INSTRUCTIONS.md` flagged as a non-standard top-level file; and the
  "declared imports not used" note for the two annotation packages, which are
  reached through `asNamespace()` and so were invisible to check.

* Remaining and expected: installed size 58 Mb, which is the lookup tables.

# cpgdirection 2.2.0

The Mendelian randomisation layer worked but was invisible, which made a working
layer indistinguishable from a broken one. Three changes, no change to any
prediction.

* **New `smr_agreement` column.** SMR frequently loses the evidence ladder to a
  catalogue consensus, and until now it then vanished from the output entirely —
  including where it *disagreed*. In a 65-CpG HPA panel the catalogue reported
  `-1` for cg16302441/POMC in all three tissues while SMR reported `+1` at
  p = 3.8e-8, and nothing in the output said so. A correlational prediction and
  a causal estimate pointing opposite ways is a finding, not a tie to be settled
  silently by precedence. Disagreements are now flagged in the column, in `note`
  and in the printed summary.

* **The printout reports the SMR layer even when it wins nothing**, giving the
  number of CpGs with SMR evidence, how many supplied `best_direction`, and how
  many were outranked. Where a set has no SMR evidence at all it says so and
  says why: GoDMC is a 450K meta-analysis, so EPIC-only probes were never
  assayed and cannot have an mQTL.

* **When `best_direction` contains blanks, the printout now gives the exact call
  that fills them** rather than leaving the reader to find `target_tissue` in the
  help page.

# cpgdirection 2.1.1

* Fixes `factor level [11] is duplicated`, which broke `cpg_expression_direction()`
  in 2.1.0. The print method held the evidence levels as an inline vector that
  had been hand-edited each time a layer was added; the 1.9.0 and 2.1.0 patches
  each appended `distance_targeted` and `catalogue_targeted`, and `factor()`
  rejects a duplicated level.

* The levels now live in a single constant, `CPGD_EVIDENCE`, in ladder order,
  with a `stopifnot(!anyDuplicated(.))` beside it and a regression test. The
  print method also unions in any level it sees but does not recognise, so a
  future layer cannot be silently dropped from the tally. The underlying cause
  was one list maintained in two places; the fix is that there is now one.

* No change to any prediction. 2.1.0 computed results correctly and failed only
  when printing them.

# cpgdirection 2.0.0

**A causal evidence layer, from summary-data Mendelian randomisation.**

The package no longer depends only on what someone happened to measure. A new
layer asks whether a genetic instrument for methylation at a CpG moves
expression of a gene -- so its coverage is not bounded by which pairs a study
chose to test.

* **413,982 CpG-gene pairs, 103,285 CpGs, 14,056 genes.** Exposure: GoDMC blood
  mQTL meta-analysis (~32,000 individuals). Outcome: eQTLGen blood cis-eQTL
  (31,684 individuals). Exposed as `cpgd_smr_directions()` and folded into
  `best_direction`.

* **Validated, not asserted.** Against the directly measured eQTMs already
  shipped, on the pairs where both exist:

  | tier | pairs | validated on | concordance |
  |---|---:|---:|---:|
  | S1 (>1 instrument, all agreeing) | 29,456 | 2,141 | **95.9% ± 0.8** |
  | S2 (single instrument) | 366,838 | 6,008 | 84.9% ± 0.9 |
  | S3 (instruments disagree) | 17,688 | 456 | 70.4% ± 4.2 |

  Majority-class baseline on the same overlap: 59.8%. S1 beats tier A of the
  catalogue models and approaches directly measured evidence. S3 scoring *below*
  S2 is the expected ordering, and a small check on the method: instrument
  disagreement is a warning, not extra evidence.

* **The ladder is reordered by measured accuracy, not by recency**: measured,
  `smr_high`, `catalogue_consensus`, `smr_moderate`, `catalogue_single`,
  `smr_weak`, then the distance layers.

* **It reaches where the distance curves cannot.** Median CpG-gene distance in
  the SMR table is 135 kb, 90th percentile 657 kb -- precisely the range where
  the distance curves flatten and the catalogue models have no support.

* The direction distribution is **49.2% positive**, against 26.2% from the
  catalogue models. SMR is not thresholded on an imbalanced classification
  problem, so it does not shrink towards the majority class.

* **What SMR does not establish.** A significant result is consistent with
  methylation affecting expression and equally with **linkage**: two distinct
  causal variants in LD, one acting on each trait. HEIDI is not run. Multiple
  concordant instruments (S1) argue against linkage without excluding it. Both
  datasets are whole blood. Palindromic SNPs and unresolvable allele mismatches
  were dropped rather than assumed -- 21,698 comparisons refused on those
  grounds, 9 on allele mismatch.

# cpgdirection 1.9.0

* New `target_tissue` argument. Where the tissues disagree, conflicts are
  resolved using a named tissue rather than abstaining, producing
  `distance_targeted` and `catalogue_targeted` rows. This exists because blood
  holds an effective veto over distance calls -- unanimity is required and blood
  never votes `+1` -- which is correct when the question is about blood and
  wrong when it is not. A saliva panel of brain-expressed genes is not asking
  about blood.

  It is tissue **selection, not voting**. Name the tissue on biological grounds
  before looking at the votes. A majority rule tuned after seeing which way each
  tissue fell is selecting on the outcome, and would be indefensible at review;
  a tissue named in advance is a stated assumption that can be argued with.

  The price is carried in the output: targeted rows report the declared tissue's
  accuracy, not the consensus figure, and stay flagged in
  `direction_uncertain` so the sensitivity analysis still covers them. For
  `solid_tissue` that accuracy is 0.55-0.70, tumour tissue, provisional, with
  brain and kidney below their own majority baselines under leave-one-cohort-out
  -- the clause to quote in a limitations section if solid tissue is being used
  as a brain proxy.

# cpgdirection 1.8.1

* States plainly what 1.8.0 made structural: the blood distance curve peaks at
  **0.449 and never reaches 0.5 at any distance**, and a distance call requires
  all three tissue curves to agree. `distance_only` therefore **cannot return
  `+1` for any CpG, ever**. Every direction that layer contributes is `-1` by
  construction. That is a property of the fitted curves, not of any user's
  panel, and it is the whole explanation for a distance-dominated panel looking
  uniformly negative. Now said in the help page and in
  `cpgd_direction_balance()`.

* Corrected the `cpgd_direction_balance()` footer, which still described the
  pre-1.8.0 pooled-curve rule.

# cpgdirection 1.8.0

**The distance layer no longer averages away tissue disagreement.**

Prompted by the observation that a real HPA panel returned 43 directions of
which every single one was `-1`. Investigating that produced a genuine defect,
though not the one expected.

* The three distance curves are not parallel. Blood asymptotes at 0.449 and
  **never crosses 0.5 at any distance**; nasal epithelium reaches 0.594 and
  solid tissue 0.861. Beyond roughly 20 kb they disagree, often sharply.
  Averaging them yielded a pooled probability near 0.50, which the previous
  version read as "no information" and reported as `distance_uninformative`.
  It was not ignorance: blood was saying `-1` and solid tissue `+1`, both
  clearly.

* A distance call now enters `best_direction` only when **all three curves agree
  in sign**. Where they disagree the row is `distance_tissue_conflict`, with new
  `dist_dir_blood`, `dist_dir_nasal` and `dist_dir_solid` columns so the tissue
  relevant to the question can be read directly. This is the rule the catalogue
  layer already used, now applied consistently.

* On the 65-CpG HPA panel: 45 unanimous `-1` calls and 20 tissue conflicts,
  against the previous 43 calls and 22 uninformative. **`best_direction` gains
  no new `+1`** — the change buys transparency about where the tissues part
  company, not extra positives. In 19 of the 20 conflicts solid tissue says
  `+1` while blood says `-1`.

* `min_distance_info` is demoted from the main gate to a floor on the pooled
  probability, default lowered from 0.10 to 0.02. Unanimity across three
  independently derived curves is evidence in its own right and is no longer
  suppressed for having a modest magnitude.

* **Known limitation, documented not silently fixed.** The isotonic fit is
  monotone by construction and in three bins it pulls a raw estimate across 0.5:
  nasal at 41,827 bp (raw 0.5404, fitted 0.4884), nasal at 69,506 bp (raw
  0.5053, fitted 0.4901) and solid tissue at 25,171 bp (raw 0.5216, fitted
  0.4979). The fitted value also falls outside the raw Wilson interval in nine
  bins overall. The fit is the pre-registered method and has not been changed;
  `distance_curves.csv` ships `p_positive_raw` alongside `p_positive` so the
  discrepancy can be inspected.

# cpgdirection 1.7.1

* `cpgd_selftest()` now passes. `fread` returns an empty character field as `""`
  rather than `NA`, while the package returns `NA_character_`, so every
  abstaining row was scored a mismatch on `target_gene` and `evidence_tier` --
  exactly the 12 of 40 that produced the 70% figures. Numeric fields were
  unaffected, which is why `probability_plus1` and `confidence` already read
  100%. Blanks are now normalised before comparison.

* Removed the remaining data.table shallow-copy warning. It came from
  `cpg_direction_universal()`, which built its per-tissue probability columns
  with `[[<-`; the 1.6.1 fix had only covered the equivalent code in the
  multi-tissue assembly.

# cpgdirection 1.7.0

* Two alternative codings for rows whose sign is nearly free to flip.
  `best_direction` still abstains on them, which is honest but leaves nothing to
  analyse. `best_direction_filled` now takes the point estimate for those rows
  and `best_direction_flipped` takes its opposite, with `direction_uncertain`
  flagging which rows differ.

  **They are a sensitivity analysis, not a menu.** The question is whether a
  downstream result survives both codings, not which coding produces one.
  Choosing the column that fits the data better is selecting on the outcome: it
  doubles the effective number of tests and biases towards significance. The
  printout, the HTML report and the help page all say so.

  A row counts as uncertain when the distance prior sits within
  `min_distance_info` of 0.5, or when a call was made below that confidence.

# cpgdirection 1.6.1

* `cpgd_selftest()` no longer reports `NA%` and fails. It compared directions
  with `==` and let `NA` propagate, so a single row missing on one side turned
  the whole score into `NA`. Missing-on-one-side now counts as a mismatch.

* The test now gates on `probability_plus1`, the model's actual output, together
  with gene assignment and evidence tier, rather than on direction. Direction is
  the sign of the probability, so testing it added nothing and was the one field
  whose *reporting* convention had changed. Verified against the shipped tables:
  probability, confidence, tier and TSS distance agree exactly on all 28
  reference CpGs present in the blood table.

* `cpgd_selftest(regenerate = "file.csv")` writes a fresh reference from the
  installed package, for use after a deliberate change of behaviour.

* Removed a spurious data.table shallow-copy warning. The per-tissue assembly
  used `[[<-`, which is base-R assignment and copies the table, so every later
  `:=` warned about a copy the caller had nothing to do with.

# cpgdirection 1.6.0

* `cpg_expression_direction(..., html_out = "results.html")` writes a formatted
  HTML table on completion, into `getwd()` unless given an absolute path, so
  `setwd()` controls where results land. `cpgd_report()` does the same for a
  result you already have. Uses **gt** when installed and falls back to a plain
  styled table when it is not, so a long run never fails at the last step for
  want of an optional package. Direction cells are shaded by sign and
  `best_evidence` by strength, because that is the column readers skip.

* New `cpgd_direction_balance()`. Reports the +1 / -1 balance in the shipped
  tables and, given a result, tests its catalogue-derived calls against the
  relevant base rate with an exact binomial test. Distance-only calls are
  excluded on purpose: the pooled distance curve sits below 0.5 out to ~19 kb
  and above it only past ~300 kb, so a panel of promoter-proximal CpGs cannot
  receive a +1 from that layer at all, and including those calls would make an
  entirely expected result look anomalous.

* Documented a real property of the models rather than a bug: they predict +1
  less often than it occurs. In blood, 25.8% of predictions are +1 against
  39.5% among directly measured pairs, and the same shrinkage towards the
  majority class appears in nasal epithelium (21.3% vs 38.5%) and solid tissue
  (23.8% vs 46.1%). A +1 call is therefore conservative; a -1 call carries less
  weight than its frequency suggests.

# cpgdirection 1.5.1

Three corrections found by running a real 68-CpG HPA panel through 1.5.0.

* The distance layer no longer requires the three tissue curves to agree on a
  side of 0.5. The blood curve asymptotes at 0.449 and never crosses it, while
  nasal and solid rise past it, so beyond roughly 30 kb they disagree by
  construction. Requiring consensus discarded a perfectly well-determined
  distance and reported `no_evidence` for pairs whose TSS distance was known to
  the base pair. The pooled probability is used instead, with a new
  `agree_universal` column recording whether the curves concurred.

* New `best_evidence` level `distance_uninformative`, and a
  `min_distance_info` argument (default 0.10 on the `|p - 0.5| * 2` scale).
  Where the distance prior sits at 0.47-0.51 it is not a weak answer, it is no
  answer, and labelling it `no_evidence` was as misleading as giving it a
  direction would have been.

* `best_confidence` for catalogue rows now averages over the tissues that
  actually **called**, not over every tissue holding a number. Including
  abstainers pulled single confident calls down towards the abstention
  threshold: `catalogue_single` rows were reporting a median confidence of 0.20
  when the calling tissue was well above its tier threshold.

* `best_tier` is now `U` for distance rows and `M` for measured rows rather than
  `NA`, and the printout warns when most directions rest on distance alone,
  since inverse associations are the majority class in every catalogue and those
  calls sit close to the base rate.

# cpgdirection 1.5.0

`cpg_expression_direction()` now does everything by default.

* `tissue = "all"` is the new default. One call runs the measured catalogue, all
  three tissue models, and the distance fallback, and returns them side by side
  with a single **`best_direction`** column holding the answer from the
  strongest evidence available for that CpG. `best_evidence` says which layer it
  came from, and `best_expected_accuracy` what it is worth — those two must be
  read together, since `measured` is worth ~1.00 and `distance_only` ~0.62.

* The evidence ladder, in order: a measured eQTM for the requested gene; a
  cross-tissue consensus; a single calling tissue; distance alone.

* Where tissues return opposite directions, `best_evidence` is
  `tissue_conflict` and `best_direction` is `NA`. The distance curve is
  deliberately not used to break the tie, because it is weaker than either
  tissue it would be adjudicating.

* Passing a specific tissue returns the previous single-tissue result unchanged.

* `cpg_direction_all_tissues()` is now a thin wrapper on the same engine with
  `universal = FALSE`. Its output gains `tier_*` and `best_*` columns; nothing
  was removed.

# cpgdirection 1.4.1

Two bug fixes, both affecting real panels rather than the examples.

* `annotation_mismatch` was TRUE for every CpG with **no** measured record.
  `strsplit(NA, ";")` returns a length-one element holding `NA` rather than an
  empty vector, so the emptiness guard never fired and absence of evidence was
  reported as evidence of conflict. On a 65-CpG HPA panel this flagged 63 CpGs;
  it now flags 5, all real. Absence of a measured record now returns `NA`.

* The mismatch comparison now accepts compound panel labels. A requested gene of
  `LINC02210_CRHR1` or `ENSG00000253356_STAR` is compared on the whole label,
  both separator forms, and its last token, which is where the informative
  symbol sits in Illumina-style annotations.

* `cpg_direction_universal()` no longer requires `genes`. If omitted, the symbol
  is taken from the identifier, as every other function in the package already
  did. Inputs whose symbol matches no gene in the annotation are dropped with a
  count, not an error.

* Duplicate probe identifiers (`cg09527270_TC11` / `_TC12` / `_TC13`) are still
  collapsed, but the printout now says how many were collapsed instead of
  silently reporting a smaller total.

* Removed `inst/extdata/eqtm_direction_lookup_hg19.csv.gz`, an unreferenced
  18 MB holdover from the single-tissue version, superseded by
  `lookup_blood_hg19.csv.gz`. Installed size drops by 18 MB, to about 52 MB.

# cpgdirection 1.4.0

* The four Bioconductor annotation packages moved from `Suggests` to `Imports`,
  with `biocViews:`, so that a single `BiocManager::install()` resolves
  everything. This is what makes the tarball route require the dependencies to
  be present first.

* New `cpgd_gene_tss(refresh = FALSE)`. Universal mode now builds the TSS
  annotation from TxDb on first use and caches it under
  `tools::R_user_dir("cpgdirection", "cache")` — once per machine rather than
  once per session. Previously `cpgd_build_gene_tss()` had to be run by hand and
  its CSV passed in.

No model, lookup table or calibration changed in 1.4.x. Predictions are
identical to 1.3.0.

# cpgdirection 1.3.0

* Universal mode: `cpg_direction_universal()`, `cpgd_build_gene_tss()`,
  `cpgd_cpg_positions()`. Predicts direction from distance alone, for any
  CpG-gene pair, at evidence tier U. Weaker by construction — AUC 0.616 against
  0.793 for the full blood model — and labelled as such in every row.

* Bundled 930,178 EPIC v2 probe positions in hg19.

# cpgdirection 1.2.0

* Multi-tissue: `cpg_direction_all_tissues()` returns one column per tissue plus
  `consensus_direction`, for use where the sampled tissue is a proxy rather than
  the tissue of interest.

* `cpgd_measured_eqtms()` exposes the 25,056 directly measured CpG-gene pairs,
  and results now carry `annotation_mismatch`: a flag for CpGs whose measured
  association is with a gene other than the array annotation. 23.1% of measured
  pairs involve a gene that is not the nearest measured gene.

* Separate lookup tables, models and calibration for blood, nasal epithelium and
  solid tissue.

# cpgdirection 1.1.1

* `cpgd_lookup_info()` used `tissue` in its body without accepting it as an
  argument.

# cpgdirection 1.1.0

* The identifier parser now extracts `cg[0-9]{6,}` from anywhere in a string.
  Prefixed names such as `Zcg00335286_TC21_MC2R`, produced by z-scoring a
  column, were previously dropped without warning.

* New `REQUESTED_GENE_UNAVAILABLE` status and `allow_nearest_gene` argument. When
  a gene is requested and the table holds no pairing for it, the function no
  longer silently returns a direction for the nearest catalogued gene, which is a
  different gene. `cg02079741` annotated to POMC was returning a confident
  direction for DNMT3A.

* `cpgd_parse_check()` shows what the parser extracted from each identifier.

# cpgdirection 1.0.0

* First release. `cpg_expression_direction()`, `cpgd_selftest()`,
  `cpgd_accuracy_table()`, `cpgd_lookup_info()`.
