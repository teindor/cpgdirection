# How the files in `inst/extdata` were made

Every file the package ships in `inst/extdata` is listed here with its
sources, the rule that produced it, and the script where one exists. The ten
large layers (catalogue lookups, SMR directions, probe positions, brain
directions, saliva-bridge scores, the EPIC v2 probe-gene annotation, distal
links and the Onco-eQTM consensus) are not in this directory: they are served
by the `cpgdirectionData` ExperimentHub package, whose own `inst/scripts`
documents them. Genome build is hg19 throughout unless stated.

All source data are open resources; nothing here derives from controlled-access
individual-level data.

## `measured_eqtms.csv.gz`

Directly measured cis-eQTM pairs with a replicated sign, 75,062 rows
(blood 26,555; nasal epithelium 11,045; solid tissue 37,462).

* **Sources.** The eQTM Atlas (Sriram et al. 2026, bioRxiv
  doi:10.64898/2026.06.07.730721; browser
  <https://shiny.crc.pitt.edu/eqtm_browser/>, code
  <https://github.com/ads303/eQTM-Atlas>, MIT). Per-cohort CSV exports from
  the Downloads tab, downloaded 5-6 October 2026: EVA-PR (four nasal subsets
  and white blood cells), FUSION (muscle), Framingham Heart Study (whole blood;
  the export is capped at 10,000 rows), Colonomics (colon), EvoImmunoPop,
  HELIX, and GTEx (eight tissues: adipose, artery, colon, lung, muscle, nerve,
  ovary, whole blood). Coordinates are hg19 (EnsDb.Hsapiens.v75). The blood
  pairs of the first package build (`source = "packaged_v1"`) came from the
  Atlas exports of the Framingham and HELIX cohorts downloaded 5 August 2026;
  HELIX and EvoImmunoPop exports no longer carry a sign, so the v1 blood pairs
  are kept where the October exports cannot replace them.
* **Rule** (`inst/scripts/build_measured_eqtms.R`). A pair is cis when CpG and
  TSS are on the same chromosome within 250 kb (decided from the coordinates,
  not from the Atlas `eQTM_label`, which is unreliable for FUSION). The sign
  is the sign of `beta`, else of the test statistic; rows without either are
  dropped. Significance follows each cohort's own rule: FDR < 0.05 where the
  export carries an FDR column (FUSION and EVA-PR "All Subjects" at 0.01), else
  the p threshold the Atlas lists for the cohort (Colonomics 1.7e-8,
  EvoImmunoPop 4.99e-5, HELIX 0.05, Framingham 1e-7). Atlas tissues map to
  three classes: blood, monocyte, WBC and leukocyte to `blood`; nasal to
  `nasal_epithelium`; everything else to `solid_tissue`. Within a class the
  sign is the cohort majority (ties dropped); `n_cohorts`, `cohort_agreement`
  and `cohorts` record the provenance. Tumour-derived pairs from the earlier
  solid-tissue catalogue are not in this file (tumour evidence lives in the
  Onco-eQTM layer of `cpgdirectionData`).
* **Columns.** `cpg_id`, `target_gene`, `tissue`, `direction` (+1/-1),
  `tss_dist` (signed CpG-TSS distance, bp), `n_cohorts`, `cohort_agreement`,
  `cohorts`, `source` (`eqtm_atlas` or `packaged_v1`).
* **Reproduce.** From the package root, with the Atlas exports in a directory:
  `Rscript inst/scripts/build_measured_eqtms.R --atlas-dir <dir> --out tools/out`,
  then copy `tools/out/measured_eqtms_v2.csv.gz` to
  `inst/extdata/measured_eqtms.csv.gz`.

## `distance_curves.csv`

The distance prior: the empirical probability that a measured cis association
is positive, as a function of CpG-to-TSS distance, one curve per tissue class
(40 bins each, 120 rows).

* **Sources.** The measured cis-eQTM pairs available to the package at
  version 1.8.0 (NEWS): the blood and nasal-epithelium pairs described above
  (eQTM Atlas exports of Framingham, HELIX and EVA-PR, August 2026) and, for
  solid tissue, the tumour-derived TCGA pairs that trained the solid-tissue
  catalogue model. The solid curve is therefore a tumour curve and the help
  pages say so.
* **Rule.** Pairs were binned by absolute CpG-TSS distance into 40
  equal-count bins per tissue; `dist_mid` is the median distance of the bin,
  `n` its size, `p_positive_raw` the fraction of pairs with a positive sign,
  `ci_lo`/`ci_hi` the Wilson 95% interval, and `p_positive` the isotonic
  (monotone non-decreasing) fit of `p_positive_raw` against distance, which is
  the pre-registered form of the prior. Both the raw and the fitted value are
  shipped so that the three bins where the fit crosses 0.5 against the raw
  estimate can be inspected (NEWS 1.8.0).
* **Script.** Fitted during the 1.x development series; the fitting script is
  not part of the package. The procedure above is complete: `stats::isoreg()`
  on the binned fractions, Wilson intervals from `stats::prop.test()`-style
  arithmetic.

## `epicv2_probe_qc.csv.gz`

Probe reliability flags for 926,858 EPIC v2 CpGs, used by
`cpg_gene_pairs(probe_qc = "exclude")`.

* **Sources.** (1) The Zhou laboratory InfiniumAnnotation manifest for EPIC
  v2, `EPICv2.hg38.manifest.tsv.gz` (<https://zwdzwd.github.io/InfiniumAnnotation>,
  downloaded 23 August 2026), for the general mask (`MASK_general`:
  degenerate or non-unique mapping, SNP-contaminated extension base) and the
  mask reasons; (2) the Zhou EPIC v2 hg19 manifest from the same site, for an
  independent check of the packaged hg19 probe coordinates; (3) the Garvan
  experimental evidence of Peters et al. 2024, BMC Genomics 25:251
  (doi:10.1186/s12864-024-10027-5): WGBS-confirmed cross-hybridising probes,
  BLAT-predicted off-targets and alignment mismatches, from the paper's
  supplementary tables.
* **Rule.** Flags are collapsed from probe to CpG level over replicate
  probes. Exclusion-grade flags use the every-replicate rule:
  `masked_general` (Zhou) and `cross_hybridizing` (Garvan, WGBS-evidenced) are
  TRUE only when every replicate probe of the CpG is flagged; a CpG with at
  least one clean replicate is `masked_partial`. `qc_exclude` is their union.
  Informational flags use the any-replicate rule: `ch_blat_predicted`,
  `mapping_flagged`, `n_offtargets_max`. `pos_hg19_verified` is TRUE where the
  packaged hg19 position equals the Zhou hg19 manifest position (99.5% exact;
  3,632 CpGs flagged unverified; NA where the CpG is absent from that
  manifest). `qc_sources` names the tracks consulted.
* **Script.** Built at version 2.3.0; the build script is not part of the
  package. The rule above is the complete specification.

## `hic_cotargets.csv.gz`

Genes physically contacted by the enhancer a CpG sits in, for 23,194 CpGs
(`cpg_id`, `hic_genes` separated by `;`, `n_hic_genes`). Reported by
`cpg_brain_bridge()` as candidate co-regulation; a contact carries no
direction.

* **Source.** PsychENCODE integrative analysis resource INT-16, adult
  dorsolateral prefrontal cortex Hi-C enhancer-promoter linkages (Wang et al.
  2018, Science 362:eaat8464; <http://resource.psychencode.org>): 90,015
  enhancer-to-gene links over 6,533 genes, hg19 as released.
* **Rule.** A CpG whose hg19 position falls inside a linked enhancer receives
  every gene that enhancer is linked to; no padding was applied to the
  shipped file (a 2 kb padding was used only as a sensitivity check during
  validation). CpGs outside every linked enhancer are absent.
* **Script.** Built at version 2.3.0; not part of the package.

## `bridge_deliverable.csv`

The graded 728-CpG periphery-to-brain deliverable behind
`cpgd_bridge_deliverable()` and `cpg_brain_bridge()`.

* **Sources.** The brain causal-direction layer (tier T1 of the SMR analysis
  Brain-mMeta mQTL x BrainMeta v2 cis-eQTL, served by `cpgdirectionData` as
  `brain_directions`); IMAGE-CpG (Braun et al. 2019, Transl Psychiatry 9:47;
  27 living donors, saliva, blood and buccal against brain, EPIC and 450K
  estimates); AMAZE-CpG (Nishitani et al. 2023, Transl Psychiatry 13:72;
  19 donors, live brain against blood, saliva and buccal, raw and
  cell-composition-adjusted); GSE59685 (Hannon et al. 2015, Epigenetics
  10:1024; blood against four brain regions, n = 71-75, recomputed raw and
  cell-adjusted); Sommerer et al. 2022, Clin Epigenetics 14:139 (buccal
  against prefrontal cortex, n = 120); and GSE15745 (Gibbs et al. 2010, PLoS
  Genet 6:e1000952) as the brain eQTM reference used to validate the
  direction layer. A pediatric coupling table under a CC BY-NC licence was
  deliberately left out.
* **Rule.** A CpG enters with a tier-T1 brain direction and "robust"
  periphery-brain concordance in IMAGE-CpG (Spearman rho >= 0.5 on both the
  EPIC and the 450K estimate, in at least one of saliva, blood or buccal).
  Grade A requires replication in AMAZE in the same tissue (raw or
  cell-adjusted), grade B is IMAGE-CpG only; A+ and A++ add the Hannon and
  Sommerer confirmations (`confirmed_tissues`, `hannon_cortical_*`). `strict_qc`
  is FALSE for the 118 CpGs on probes discontinued in EPIC v2; `open_sea`
  flags isolated probes with lower prior reliability. The headline set is the
  205 A++ CpGs, 184 of them saliva-usable.
* **Script.** Built at version 2.3.0 from the analyses documented in the
  package NEWS and the OSF deposit (doi:10.17605/OSF.IO/U3VFK); not part of
  the package.

## `selftest_reference.csv`

Forty CpGs with the output `cpg_expression_direction()` produced for them
when the reference was last regenerated, used by `cpgd_selftest()` to verify
an installation (`input`, `cpg_id`, `target_gene`, `tss_dist`, `status`,
`direction`, `probability_plus1`, `confidence`, `evidence_tier`).

* **Source and rule.** The package's own blood catalogue lookup; the file is
  written by `cpgd_selftest(regenerate = "<path>")` and copied into
  `inst/extdata`. Last regenerated at version 1.7.1 (NEWS). Regenerate only
  after a deliberate change of behaviour, never to make a failing self-test
  pass.
