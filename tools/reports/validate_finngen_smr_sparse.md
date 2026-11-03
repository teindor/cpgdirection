# FinnGen cell-type SMR validation (2026-11-01)

Files: 9; context-pairs: 463048; contexts: blood_B, blood_CD4_T, blood_CD8_T, blood_DC, blood_Mono, blood_NK, blood_other, blood_other_T, blood_PBMC

## 1. Concordance with measured blood eQTMs

Overlap 17692 context-pairs (5047 unique pairs); majority-class baseline 60.0%.

| context | tier | n | concordance (95% CI) |
|---|---|---:|---|
| blood_B | S1 | 78 | 73.1% (62.3-81.7) |
| blood_B | S2 | 1329 | 88.8% (87.0-90.4) |
| blood_CD4_T | S1 | 255 | 93.7% (90.1-96.1) |
| blood_CD4_T | S2 | 2878 | 88.8% (87.6-89.9) |
| blood_CD4_T | S3 | 1 | 0.0% (0.0-79.3) |
| blood_CD8_T | S1 | 155 | 94.8% (90.1-97.4) |
| blood_CD8_T | S2 | 2430 | 89.7% (88.4-90.8) |
| blood_DC | S1 | 38 | 65.8% (49.9-78.8) |
| blood_DC | S2 | 1059 | 86.8% (84.6-88.7) |
| blood_DC | S3 | 2 | 50.0% (9.5-90.5) |
| blood_Mono | S1 | 162 | 92.6% (87.5-95.7) |
| blood_Mono | S2 | 2294 | 88.4% (87.1-89.7) |
| blood_Mono | S3 | 2 | 100.0% (34.2-100.0) |
| blood_NK | S1 | 88 | 93.2% (85.9-96.8) |
| blood_NK | S2 | 1943 | 88.2% (86.7-89.6) |
| blood_PBMC | S1 | 335 | 94.6% (91.7-96.6) |
| blood_PBMC | S2 | 3667 | 90.3% (89.3-91.3) |
| blood_other | S1 | 2 | 100.0% (34.2-100.0) |
| blood_other | S2 | 114 | 93.0% (86.8-96.4) |
| blood_other | S3 | 1 | 100.0% (20.7-100.0) |
| blood_other_T | S1 | 34 | 82.4% (66.5-91.7) |
| blood_other_T | S2 | 825 | 88.2% (85.9-90.3) |

All contexts pooled:

| tier | n | concordance (95% CI) |
|---|---:|---|
| S1 | 1147 | 91.3% (89.5-92.8) |
| S2 | 16539 | 89.0% (88.5-89.5) |
| S3 | 6 | 66.7% (30.0-90.3) |

## 2. Agreement with the bulk GoDMC x eQTLGen layer

Shared pairs: 266682 of 463048 context-pairs (57.6%).

| context | n shared | sign agreement (95% CI) | S1 both | S1 agreement |
|---|---:|---|---:|---|
| blood_B | 20072 | 82.1% (81.6-82.6) | 217 | 83.4% (77.9-87.8) |
| blood_CD4_T | 50616 | 84.2% (83.8-84.5) | 645 | 93.6% (91.5-95.3) |
| blood_CD8_T | 36960 | 85.5% (85.1-85.9) | 420 | 95.0% (92.5-96.7) |
| blood_DC | 13313 | 82.1% (81.5-82.8) | 105 | 79.0% (70.3-85.7) |
| blood_Mono | 40127 | 86.0% (85.6-86.3) | 483 | 95.9% (93.7-97.3) |
| blood_NK | 26815 | 83.9% (83.4-84.3) | 265 | 92.5% (88.6-95.1) |
| blood_PBMC | 67313 | 87.7% (87.5-88.0) | 958 | 96.1% (94.7-97.2) |
| blood_other | 1315 | 89.8% (88.1-91.3) | 7 | 100.0% (64.6-100.0) |
| blood_other_T | 10151 | 81.9% (81.1-82.6) | 74 | 86.5% (76.9-92.5) |

## 3. Cross-context consistency (cell types, PBMC excluded)

Unique pairs: 154680; called in >= 2 cell types: 78952; of those, unanimous: 94.6%, any disagreement: 5.4%.

| n cell types calling | pairs | unanimous |
|---:|---:|---|
| 1 | 75728 | 100.0% (100.0-100.0) |
| 2 | 31968 | 94.2% (94.0-94.5) |
| 3 | 18513 | 93.6% (93.3-94.0) |
| 4 | 11892 | 94.4% (94.0-94.8) |
| 5 | 7579 | 94.5% (93.9-95.0) |
| 6 | 4551 | 96.5% (95.9-97.0) |
| 7 | 3555 | 99.2% (98.8-99.4) |
| 8 | 894 | 100.0% (99.6-100.0) |

PBMC pseudobulk vs each cell type on shared pairs:

| cell type | n shared | agreement (95% CI) |
|---|---:|---|
| blood_B | 26005 | 99.1% (98.9-99.2) |
| blood_CD4_T | 74516 | 99.7% (99.7-99.8) |
| blood_CD8_T | 54911 | 99.8% (99.8-99.9) |
| blood_DC | 14482 | 99.1% (99.0-99.3) |
| blood_Mono | 46535 | 98.9% (98.8-99.0) |
| blood_NK | 36793 | 99.2% (99.1-99.3) |
| blood_other | 1411 | 98.4% (97.7-99.0) |
| blood_other_T | 14411 | 99.8% (99.7-99.8) |

## 4. HEIDI

| context | pass | fail | not tested |
|---|---:|---:|---:|
| blood_B | 0 | 0 | 34482 |
| blood_CD4_T | 0 | 0 | 92183 |
| blood_CD8_T | 0 | 0 | 61660 |
| blood_DC | 0 | 0 | 20073 |
| blood_Mono | 0 | 0 | 69355 |
| blood_NK | 0 | 0 | 44578 |
| blood_PBMC | 0 | 0 | 123039 |
| blood_other | 0 | 0 | 1718 |
| blood_other_T | 0 | 0 | 15960 |

