# Gene-class covariate: stimulus-responsive vs lineage vs other genes, in the package's own evidence

Run 2026-11-03 20:17. Classes: GO cytokine activity, chemokine activity, response to type I interferon, an immediate-early list and checkpoint ligands (406 genes); 48 lineage/identity genes; everything else 'other'.

## 1. Measured blood vs blood SMR sign agreement

| class | pairs | genes | agree |
|---|---|---|---|
| lineage_identity |    41 |   14 | 0.878 [0.745, 0.947] |
| other | 10963 | 2533 | 0.938 [0.933, 0.943] |
| stimulus_responsive |   111 |   41 | 0.928 [0.864, 0.963] |

## 2. Measured blood vs measured other-tissue sign agreement

| class | pairs | genes | agree |
|---|---|---|---|
| lineage_identity |   41 |   14 | 0.951 [0.839, 0.987] |
| other | 4819 | 1278 | 0.925 [0.918, 0.932] |
| stimulus_responsive |   82 |   24 | 1.000 [0.955, 1.000] |

## 3. Cohort agreement of multi-cohort measured pairs

| tissue | class | pairs | mean_cohort_agreement |
|---|---|---|---|
| blood | lineage_identity |   39 | 1.000 |
| blood | other | 2860 | 1.000 |
| blood | stimulus_responsive |   36 | 1.000 |
| nasal_epithelium | lineage_identity |  158 | 1.000 |
| nasal_epithelium | other | 3838 | 1.000 |
| nasal_epithelium | stimulus_responsive |  111 | 1.000 |
| solid_tissue | lineage_identity |    2 | 1.000 |
| solid_tissue | other | 3894 | 0.997 |
| solid_tissue | stimulus_responsive |   50 | 1.000 |

## 4. Position and sign of measured pairs

| tissue | class | pairs | share_beyond_10kb | share_negative |
|---|---|---|---|---|
| blood | lineage_identity |   183 | 0.678 | 0.689 |
| blood | other | 26091 | 0.689 | 0.562 |
| blood | stimulus_responsive |   281 | 0.534 | 0.744 |
| nasal_epithelium | lineage_identity |   309 | 0.767 | 0.421 |
| nasal_epithelium | other | 10353 | 0.708 | 0.614 |
| nasal_epithelium | stimulus_responsive |   383 | 0.692 | 0.556 |
| solid_tissue | lineage_identity |    35 | 0.800 | 0.743 |
| solid_tissue | other | 37055 | 0.676 | 0.613 |
| solid_tissue | stimulus_responsive |   372 | 0.581 | 0.742 |

## 5. SMR rung mix

| class | pairs | genes | share_replicated_or_high |
|---|---|---|---|
| lineage_identity |   1225 |    41 | 0.076 |
| other | 540224 | 16786 | 0.162 |
| stimulus_responsive |   7867 |   243 | 0.152 |

## Reading

A covariate is justified if stimulus-responsive genes are clearly lower on 1 AND 2 (non-overlapping Wilson intervals
with 'other'), in which case it enters as a reported column (gene_class) and a confidence note, not a rung change;
the stress-response pathways DMSA tests are exactly where it would bite. Lower coverage alone (6) is not evidence of
wrong signs and changes nothing.

