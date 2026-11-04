# November 2026 data build

Run 2026-11-04 08:59.

## measured_eqtms

| tissue | before | after | multi_cohort |
|---|---|---|---|
| blood | 26555 | 45311 | 6997 |
| nasal_epithelium | 11045 | 11045 | 4107 |
| solid_tissue | 37462 | 37462 | 3946 |

- duplicated rows removed: 10; packaged_v1 blood rows replaced by signed HELIX rows: 8,432
- pairs dropped as cohort ties: 118; existing pairs whose sign changed: 13
- liver included: FALSE

## epicv2_probe_qc (flags informational; qc_exclude unchanged)

| cpgs | t2t_ambiguous | pangenome_xpop_ambiguous | of_which_not_qc_exclude |
|---|---|---|---|
| 926858 | 73943 | 110649 | 77802 |

## brain_measured_eqtms (new)

8,815 ROSMAP aged-DLPFC pairs over 8,150 CpGs and 3,416 genes; negative 0.627.

