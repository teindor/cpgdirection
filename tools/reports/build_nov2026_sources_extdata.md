# Kidney, placenta, adipose and GTP whole-blood eQTMs folded into the measured layer

Run 2026-11-12 15:16. Add file: `tools/out/measured_eqtms_add_nov2026_sources.csv.gz` (36,464 rows, 31,190 CpGs, 8,357 genes).

| tissue | before | after | multi_cohort | two_groups_agree_before | two_groups_agree_after |
|---|---|---|---|---|---|
| blood | 45311 | 46790 | 7818 | 6333 | 7159 |
| nasal_epithelium | 11045 | 11045 | 4107 | 0 | 0 |
| solid_tissue | 49832 | 76958 | 9875 | 3730 | 6891 |

`two_groups_agree`: pairs seen in two or more donor groups with every cohort agreeing.

| source | tissue | pairs | already_measured_in_tissue | same_sign | opposite | new_to_tissue |
|---|---|---|---|---|---|---|
| liu2022_kidney | solid_tissue | 29,988 | 4,685 | 4,254 | 431 | 25,303 |
| delahaye2018_placenta | solid_tissue | 2,314 | 231 | 211 | 20 | 2,083 |
| voisin2016_adipose | solid_tissue | 427 | 168 | 159 | 9 | 259 |
| kennedy2018_gtp | blood | 3,735 | 2,156 | 2,024 | 132 | 1,579 |

- pairs dropped as cohort ties: blood 100; solid_tissue 412
- existing pairs whose sign changed: solid_tissue 2
- sources after: delahaye2018_placenta 2,075; eqtm_atlas 63,196; helix_dryad_2022 21,820; kennedy2018_gtp 1,579; kreitmaier2026_msk 12,386; liu2022_kidney 25,215; packaged_v1 2,885; voisin2016_adipose 239; wang2021_S2 5,398

