# SMR layer with outcome-study replication -- build report (v2)

Built 2026-11-03 15:25 from 413,982 eQTLGen pairs, 547,603 FinnGen SMR-tested pairs (p < 5e-08) and 463,048 FinnGen instrument-vote pairs; union 246,608 FinnGen pairs.
Output: `tools/reports/smr_replication_build/smr_directions.csv.gz`, 549,316 rows.

## sign(b_SMR) vs instrument vote, pairs both cover

| ctx | n | agree |
|---|---|---|
| PBMC | 90359 | 0.991 [0.990, 0.991] |
| CD4_T | 67712 | 0.990 [0.990, 0.991] |
| CD8_T | 44791 | 0.991 [0.990, 0.992] |
| Mono | 50443 | 0.991 [0.990, 0.992] |
| NK | 32290 | 0.992 [0.991, 0.993] |
| B | 24914 | 0.988 [0.987, 0.989] |
| DC | 14105 | 0.990 [0.989, 0.992] |
| other_T | 11246 | 0.991 [0.989, 0.993] |
| other |  1135 | 0.997 [0.992, 0.999] |

## Rows by outcome, replication, FinnGen evidence and rung

| outcome | replication | finngen_evidence | smr_rung | N |
|---|---|---|---|---|
| FinnGen | finngen_only | both | moderate |  39442 |
| FinnGen | finngen_only | both | reported |  15734 |
| FinnGen | finngen_only | instrument_vote | reported |  23118 |
| FinnGen | finngen_only | smr_top_snp | moderate |  38786 |
| FinnGen | finngen_only | smr_top_snp | reported |  18254 |
| eQTLGen | concordant | both | replicated |  43854 |
| eQTLGen | concordant | both | weak |   1184 |
| eQTLGen | concordant | instrument_vote | replicated |  13609 |
| eQTLGen | concordant | instrument_vote | weak |   1341 |
| eQTLGen | concordant | smr_top_snp | replicated |  11073 |
| eQTLGen | concordant | smr_top_snp | weak |    573 |
| eQTLGen | concordant_celltype | both | high |    522 |
| eQTLGen | concordant_celltype | both | moderate |  10069 |
| eQTLGen | concordant_celltype | both | weak |    267 |
| eQTLGen | concordant_celltype | instrument_vote | high |    735 |
| eQTLGen | concordant_celltype | instrument_vote | moderate |   4341 |
| eQTLGen | concordant_celltype | instrument_vote | weak |    563 |
| eQTLGen | concordant_celltype | smr_top_snp | high |    586 |
| eQTLGen | concordant_celltype | smr_top_snp | moderate |   3191 |
| eQTLGen | concordant_celltype | smr_top_snp | weak |    207 |
| eQTLGen | discordant | both | weak |   5873 |
| eQTLGen | discordant | instrument_vote | weak |   2746 |
| eQTLGen | discordant | smr_top_snp | weak |   3232 |
| eQTLGen | discordant_celltype | both | high |    160 |
| eQTLGen | discordant_celltype | both | moderate |   3028 |
| eQTLGen | discordant_celltype | both | weak |    125 |
| eQTLGen | discordant_celltype | instrument_vote | high |    252 |
| eQTLGen | discordant_celltype | instrument_vote | moderate |   1355 |
| eQTLGen | discordant_celltype | instrument_vote | weak |    358 |
| eQTLGen | discordant_celltype | smr_top_snp | high |    261 |
| eQTLGen | discordant_celltype | smr_top_snp | moderate |   1571 |
| eQTLGen | discordant_celltype | smr_top_snp | weak |    198 |
| eQTLGen | eqtlgen_only | NA | high |  17572 |
| eQTLGen | eqtlgen_only | NA | moderate | 273602 |
| eQTLGen | eqtlgen_only | NA | weak |  11534 |

## Rows by rung

| smr_rung | rows | cpgs |
|---|---|---|
| moderate | 375385 | 108988 |
| replicated |  68536 |  41759 |
| reported |  57106 |  35811 |
| weak |  28201 |  16923 |
| high |  20088 |  10802 |

## Validation against measured blood eQTMs (Wilson 95% CI)

### By rung

| smr_rung | n | agree |
|---|---|---|
| high | 1094 | 0.957 [0.943, 0.968] |
| moderate | 4617 | 0.894 [0.885, 0.902] |
| replicated | 5182 | 0.976 [0.972, 0.980] |
| reported |  519 | 0.744 [0.704, 0.779] |
| weak |  999 | 0.704 [0.675, 0.731] |

### By rung x replication x tier

| smr_rung | replication | smr_tier | n | agree | baseline |
|---|---|---|---|---|---|
| high | concordant_celltype | S1 |  208 | 0.976 [0.945, 0.990] | 0.630 |
| high | discordant_celltype | S1 |   59 | 0.932 [0.838, 0.973] | 0.661 |
| high | eqtlgen_only | S1 |  827 | 0.954 [0.938, 0.966] | 0.586 |
| moderate | concordant_celltype | S2 |  486 | 0.905 [0.876, 0.928] | 0.626 |
| moderate | discordant_celltype | S2 |  142 | 0.803 [0.730, 0.860] | 0.683 |
| moderate | eqtlgen_only | S2 | 2191 | 0.889 [0.875, 0.901] | 0.616 |
| moderate | finngen_only | S2 | 1798 | 0.904 [0.890, 0.917] | 0.582 |
| replicated | concordant | S1 | 1795 | 0.986 [0.979, 0.990] | 0.605 |
| replicated | concordant | S2 | 3387 | 0.971 [0.965, 0.976] | 0.597 |
| reported | finngen_only | S1 |    2 | 0.500 [0.095, 0.905] | 0.500 |
| reported | finngen_only | S2 |  517 | 0.745 [0.705, 0.780] | 0.621 |
| weak | concordant | S3 |  214 | 0.855 [0.802, 0.896] | 0.547 |
| weak | concordant_celltype | S3 |   48 | 0.750 [0.612, 0.851] | 0.500 |
| weak | discordant | S1 |  111 | 0.919 [0.853, 0.957] | 0.604 |
| weak | discordant | S2 |  302 | 0.593 [0.536, 0.647] | 0.566 |
| weak | discordant | S3 |   68 | 0.382 [0.276, 0.501] | 0.647 |
| weak | discordant_celltype | S3 |   29 | 0.655 [0.473, 0.801] | 0.552 |
| weak | eqtlgen_only | S3 |  227 | 0.696 [0.633, 0.752] | 0.621 |

### FinnGen-only rows by evidence x HEIDI status

| outcome | finngen_evidence | heidi_status | smr_rung | n | agree |
|---|---|---|---|---|---|
| FinnGen | both | fail | moderate | 721 | 0.892 [0.867, 0.912] |
| FinnGen | both | fail | reported |  96 | 0.781 [0.689, 0.852] |
| FinnGen | both | not_tested | moderate |  15 | 1.000 [0.796, 1.000] |
| FinnGen | both | not_tested | reported |   5 | 0.800 [0.376, 0.964] |
| FinnGen | both | pass | moderate |  63 | 0.968 [0.891, 0.991] |
| FinnGen | both | pass | reported |  14 | 0.786 [0.524, 0.924] |
| FinnGen | instrument_vote | fail | reported | 129 | 0.605 [0.518, 0.685] |
| FinnGen | instrument_vote | pass | reported |  12 | 0.917 [0.646, 0.985] |
| FinnGen | smr_top_snp | fail | moderate | 908 | 0.903 [0.882, 0.921] |
| FinnGen | smr_top_snp | fail | reported | 225 | 0.773 [0.714, 0.823] |
| FinnGen | smr_top_snp | not_tested | moderate |   4 | 1.000 [0.510, 1.000] |
| FinnGen | smr_top_snp | pass | reported |  38 | 0.868 [0.727, 0.942] |
| FinnGen | smr_top_snp | pass | moderate |  87 | 0.954 [0.888, 0.982] |

Rung rules are in the header of tools/build_smr_replication.R; the accuracy strings in
R/utils.R (CPGD_SMR_ACCURACY) must match the 'By rung' table above.

