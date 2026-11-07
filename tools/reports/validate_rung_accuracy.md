# Every ladder source, scored in each target tissue against an independent yardstick

Run 2026-11-04 23:01. Agreement with the yardstick's measured sign, Wilson 95% CI. Sources are scored standalone, on every yardstick pair where they call (catalogue sources: formed from predicted calls only, as the ladders form them). `n_cpg`: distinct CpGs.

Yardsticks: adult blood = MESA monocytes + T cells; child blood = HELIX Dryad pairs outside the blood catalogue's DIRECT set; nasal = EVA-PR nasal pairs outside the nasal catalogue's DIRECT set (the nasal catalogue is trained on EVA-PR: its rows are same-study); solid = FUSION and GTEx solid, each held out in turn, pooled, each pair once. Measured cohorts are counted by donor group. Distance is not scored against HELIX or EVA-PR (its blood and nasal curves were fitted on them).

## blood (adult): MESA

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_model_blood_B |     7 |     7 | 1.000 | [0.646, 1.000] |
| 2 | catalogue_single_blood_B |     3 |     3 | 1.000 | [0.439, 1.000] |
| 3 | measured_blood_replicated |   435 |   432 | 0.970 | [0.950, 0.982] |
| 4 | catalogue_consensus_B |    18 |    18 | 0.944 | [0.742, 0.990] |
| 5 | measured_solid_replicated |    78 |    78 | 0.936 | [0.859, 0.972] |
| 6 | smr_replicated_S1 |   244 |   240 | 0.926 | [0.886, 0.953] |
| 7 | smr_weak_discordant |    23 |    23 | 0.913 | [0.732, 0.976] |
| 8 | smr_high |   239 |   239 | 0.908 | [0.865, 0.938] |
| 9 | measured_blood_single |   903 |   880 | 0.888 | [0.866, 0.907] |
| 10 | onco_consensus |   766 |   763 | 0.850 | [0.823, 0.873] |
| 11 | measured_nasal_single |   326 |   313 | 0.847 | [0.803, 0.882] |
| 12 | smr_replicated_S2 |   733 |   721 | 0.846 | [0.818, 0.870] |
| 13 | catalogue_model_nasal_B |    44 |    44 | 0.841 | [0.706, 0.921] |
| 14 | catalogue_single_nasal_B |    30 |    30 | 0.800 | [0.627, 0.905] |
| 15 | catalogue_consensus_A |   534 |   515 | 0.800 | [0.764, 0.831] |
| 16 | smr_moderate | 1,179 | 1,147 | 0.797 | [0.773, 0.819] |
| 17 | catalogue_model_blood_A |   995 |   947 | 0.755 | [0.727, 0.780] |
| 18 | catalogue_single_blood_A |   370 |   363 | 0.751 | [0.705, 0.793] |
| 19 | measured_solid_single |   379 |   370 | 0.749 | [0.703, 0.790] |
| 20 | smr_moderate_s3_replicated |    43 |    41 | 0.721 | [0.573, 0.833] |
| 21 | smr_moderate_discordant_celltype |    70 |    70 | 0.714 | [0.599, 0.807] |
| 22 | smr_moderate_finngen |   424 |   410 | 0.710 | [0.665, 0.751] |
| 23 | catalogue_model_solid_A |   297 |   295 | 0.694 | [0.639, 0.743] |
| 24 | distance_only | 5,196 | 4,884 | 0.668 | [0.655, 0.681] |
| 25 | catalogue_model_nasal_A | 1,430 | 1,376 | 0.659 | [0.634, 0.684] |
| 26 | distance_curve_nasal | 7,025 | 6,584 | 0.640 | [0.628, 0.651] |
| 27 | distance_curve_blood | 7,025 | 6,584 | 0.639 | [0.628, 0.650] |
| 28 | catalogue_single_nasal_A |   756 |   748 | 0.626 | [0.591, 0.659] |
| 29 | distance_curve_solid | 7,025 | 6,584 | 0.610 | [0.598, 0.621] |
| 30 | smr_weak |    73 |    72 | 0.575 | [0.461, 0.682] |
| 31 | catalogue_single_solid_A |    57 |    57 | 0.386 | [0.271, 0.516] |

## blood (child): HELIX (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | measured_solid_replicated |   163 |   155 | 0.994 | [0.966, 0.999] |
| 2 | measured_blood_replicated |   620 |   583 | 0.990 | [0.979, 0.996] |
| 3 | smr_high |   770 |   711 | 0.966 | [0.951, 0.977] |
| 4 | measured_blood_single | 1,550 | 1,416 | 0.964 | [0.953, 0.972] |
| 5 | smr_replicated_S1 |   887 |   828 | 0.963 | [0.948, 0.973] |
| 6 | smr_replicated_S2 | 2,188 | 2,055 | 0.960 | [0.951, 0.967] |
| 7 | smr_moderate | 2,304 | 2,107 | 0.902 | [0.889, 0.913] |
| 8 | measured_nasal_single |   230 |   201 | 0.896 | [0.849, 0.929] |
| 9 | catalogue_single_solid_A |     9 |     9 | 0.889 | [0.565, 0.980] |
| 10 | smr_moderate_s3_replicated |   145 |   145 | 0.855 | [0.789, 0.903] |
| 11 | smr_moderate_finngen | 1,205 | 1,127 | 0.848 | [0.827, 0.867] |
| 12 | measured_solid_single |   820 |   740 | 0.835 | [0.808, 0.859] |
| 13 | onco_consensus |   974 |   942 | 0.830 | [0.805, 0.852] |
| 14 | catalogue_model_solid_A |    11 |    11 | 0.818 | [0.523, 0.949] |
| 15 | smr_moderate_discordant_celltype |   119 |   119 | 0.790 | [0.708, 0.853] |
| 16 | smr_weak_discordant |    72 |    72 | 0.722 | [0.610, 0.812] |
| 17 | catalogue_single_nasal_A |   198 |   176 | 0.672 | [0.604, 0.733] |
| 18 | catalogue_model_nasal_A |   200 |   178 | 0.670 | [0.602, 0.731] |
| 19 | smr_weak |   192 |   177 | 0.635 | [0.565, 0.700] |
| 20 | catalogue_consensus_A |     2 |     2 | 0.500 | [0.095, 0.905] |

## nasal_epithelium (any): EVA-PR nasal (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_single_nasal_A |   6 |   6 | 1.000 | [0.610, 1.000] |
| 2 | measured_blood_replicated |  87 |  80 | 0.977 | [0.920, 0.994] |
| 3 | measured_solid_replicated |  29 |  25 | 0.966 | [0.828, 0.994] |
| 4 | onco_consensus | 279 | 215 | 0.957 | [0.926, 0.975] |
| 5 | measured_solid_single | 182 | 132 | 0.912 | [0.862, 0.945] |
| 6 | measured_blood_single | 228 | 198 | 0.882 | [0.833, 0.917] |
| 7 | smr_replicated_S1 |   8 |   8 | 0.875 | [0.529, 0.978] |
| 8 | smr_high |  74 |  66 | 0.824 | [0.722, 0.894] |
| 9 | smr_weak_discordant |   4 |   4 | 0.750 | [0.301, 0.954] |
| 10 | catalogue_model_nasal_A |  11 |  11 | 0.727 | [0.434, 0.903] |
| 11 | smr_weak |  26 |  25 | 0.692 | [0.500, 0.835] |
| 12 | smr_moderate | 295 | 254 | 0.614 | [0.557, 0.667] |
| 13 | smr_replicated_S2 |  29 |  28 | 0.586 | [0.407, 0.745] |
| 14 | catalogue_single_solid_A |  37 |  37 | 0.568 | [0.409, 0.713] |
| 15 | catalogue_model_solid_A |  39 |  39 | 0.564 | [0.410, 0.707] |
| 16 | catalogue_single_blood_A |  56 |  54 | 0.554 | [0.424, 0.676] |
| 17 | catalogue_model_blood_A |  63 |  61 | 0.524 | [0.403, 0.642] |
| 18 | smr_moderate_finngen |  16 |  16 | 0.500 | [0.280, 0.720] |
| 19 | catalogue_consensus_A |   4 |   4 | 0.250 | [0.046, 0.699] |
| 20 | smr_moderate_s3_replicated |   4 |   4 | 0.250 | [0.046, 0.699] |
| 21 | smr_moderate_discordant_celltype |   2 |   2 | 0.000 | [0.000, 0.658] |

## solid_tissue (any): FUSION + GTEx solid

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | measured_solid_replicated |     16 |     16 | 1.000 | [0.806, 1.000] |
| 2 | measured_blood_replicated |  1,377 |  1,203 | 0.952 | [0.939, 0.962] |
| 3 | measured_solid_single |  1,644 |  1,395 | 0.922 | [0.908, 0.934] |
| 4 | measured_nasal_single |    970 |    817 | 0.901 | [0.881, 0.918] |
| 5 | onco_consensus |  1,977 |  1,849 | 0.874 | [0.859, 0.888] |
| 6 | smr_replicated_S1 |    576 |    528 | 0.873 | [0.844, 0.898] |
| 7 | measured_blood_single |  2,626 |  2,219 | 0.865 | [0.852, 0.878] |
| 8 | catalogue_model_blood_B |     20 |     20 | 0.850 | [0.640, 0.948] |
| 9 | catalogue_model_nasal_B |    182 |    182 | 0.830 | [0.768, 0.877] |
| 10 | catalogue_single_nasal_B |    105 |    105 | 0.829 | [0.745, 0.889] |
| 11 | catalogue_consensus_B |     79 |     79 | 0.823 | [0.724, 0.891] |
| 12 | smr_high |    406 |    359 | 0.815 | [0.775, 0.850] |
| 13 | smr_weak_discordant |     55 |     54 | 0.800 | [0.676, 0.884] |
| 14 | catalogue_single_solid_A |    325 |    325 | 0.778 | [0.730, 0.820] |
| 15 | catalogue_single_blood_B |      9 |      9 | 0.778 | [0.453, 0.937] |
| 16 | smr_replicated_S2 |  1,434 |  1,322 | 0.770 | [0.747, 0.791] |
| 17 | smr_moderate_s3_replicated |     84 |     82 | 0.762 | [0.661, 0.840] |
| 18 | catalogue_consensus_A |  1,301 |  1,245 | 0.751 | [0.727, 0.774] |
| 19 | catalogue_single_nasal_A |  2,795 |  2,627 | 0.707 | [0.690, 0.724] |
| 20 | catalogue_model_solid_A |    941 |    934 | 0.707 | [0.677, 0.735] |
| 21 | smr_moderate_finngen |  1,130 |  1,045 | 0.705 | [0.678, 0.731] |
| 22 | catalogue_model_nasal_A |  4,399 |  4,009 | 0.700 | [0.687, 0.714] |
| 23 | smr_moderate |  1,800 |  1,610 | 0.694 | [0.672, 0.715] |
| 24 | catalogue_model_blood_A |  3,117 |  2,836 | 0.686 | [0.670, 0.702] |
| 25 | catalogue_single_blood_A |  1,644 |  1,517 | 0.675 | [0.652, 0.697] |
| 26 | distance_only | 20,312 | 15,989 | 0.643 | [0.637, 0.650] |
| 27 | smr_moderate_discordant_celltype |     97 |     95 | 0.629 | [0.530, 0.718] |
| 28 | distance_curve_blood | 37,489 | 25,837 | 0.613 | [0.608, 0.618] |
| 29 | smr_weak |    113 |    108 | 0.602 | [0.510, 0.687] |
| 30 | distance_curve_nasal | 37,489 | 25,837 | 0.595 | [0.590, 0.600] |
| 31 | distance_curve_solid | 37,489 | 25,837 | 0.542 | [0.537, 0.547] |

## By fold

| target | age | yardstick | source | n | agree | 95% CI |
|---|---|---|---|---|---|---|
| blood | adult | MESA | catalogue_model_blood_B |      7 | 1.000 | [0.646, 1.000] |
| blood | adult | MESA | catalogue_single_blood_B |      3 | 1.000 | [0.439, 1.000] |
| blood | adult | MESA | measured_blood_replicated |    435 | 0.970 | [0.950, 0.982] |
| blood | adult | MESA | catalogue_consensus_B |     18 | 0.944 | [0.742, 0.990] |
| blood | adult | MESA | measured_solid_replicated |     78 | 0.936 | [0.859, 0.972] |
| blood | adult | MESA | smr_replicated_S1 |    244 | 0.926 | [0.886, 0.953] |
| blood | adult | MESA | smr_weak_discordant |     23 | 0.913 | [0.732, 0.976] |
| blood | adult | MESA | smr_high |    239 | 0.908 | [0.865, 0.938] |
| blood | adult | MESA | measured_blood_single |    903 | 0.888 | [0.866, 0.907] |
| blood | adult | MESA | onco_consensus |    766 | 0.850 | [0.823, 0.873] |
| blood | adult | MESA | measured_nasal_single |    326 | 0.847 | [0.803, 0.882] |
| blood | adult | MESA | smr_replicated_S2 |    733 | 0.846 | [0.818, 0.870] |
| blood | adult | MESA | catalogue_model_nasal_B |     44 | 0.841 | [0.706, 0.921] |
| blood | adult | MESA | catalogue_single_nasal_B |     30 | 0.800 | [0.627, 0.905] |
| blood | adult | MESA | catalogue_consensus_A |    534 | 0.800 | [0.764, 0.831] |
| blood | adult | MESA | smr_moderate |  1,179 | 0.797 | [0.773, 0.819] |
| blood | adult | MESA | catalogue_model_blood_A |    995 | 0.755 | [0.727, 0.780] |
| blood | adult | MESA | catalogue_single_blood_A |    370 | 0.751 | [0.705, 0.793] |
| blood | adult | MESA | measured_solid_single |    379 | 0.749 | [0.703, 0.790] |
| blood | adult | MESA | smr_moderate_s3_replicated |     43 | 0.721 | [0.573, 0.833] |
| blood | adult | MESA | smr_moderate_discordant_celltype |     70 | 0.714 | [0.599, 0.807] |
| blood | adult | MESA | smr_moderate_finngen |    424 | 0.710 | [0.665, 0.751] |
| blood | adult | MESA | catalogue_model_solid_A |    297 | 0.694 | [0.639, 0.743] |
| blood | adult | MESA | distance_only |  5,196 | 0.668 | [0.655, 0.681] |
| blood | adult | MESA | catalogue_model_nasal_A |  1,430 | 0.659 | [0.634, 0.684] |
| blood | adult | MESA | distance_curve_nasal |  7,025 | 0.640 | [0.628, 0.651] |
| blood | adult | MESA | distance_curve_blood |  7,025 | 0.639 | [0.628, 0.650] |
| blood | adult | MESA | catalogue_single_nasal_A |    756 | 0.626 | [0.591, 0.659] |
| blood | adult | MESA | distance_curve_solid |  7,025 | 0.610 | [0.598, 0.621] |
| blood | adult | MESA | smr_weak |     73 | 0.575 | [0.461, 0.682] |
| blood | adult | MESA | catalogue_single_solid_A |     57 | 0.386 | [0.271, 0.516] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_replicated |    163 | 0.994 | [0.966, 0.999] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_replicated |    620 | 0.990 | [0.979, 0.996] |
| blood | child | HELIX (non-catalogue pairs) | smr_high |    770 | 0.966 | [0.951, 0.977] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_single |  1,550 | 0.964 | [0.953, 0.972] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S1 |    887 | 0.963 | [0.948, 0.973] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S2 |  2,188 | 0.960 | [0.951, 0.967] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate |  2,304 | 0.902 | [0.889, 0.913] |
| blood | child | HELIX (non-catalogue pairs) | measured_nasal_single |    230 | 0.896 | [0.849, 0.929] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_solid_A |      9 | 0.889 | [0.565, 0.980] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_s3_replicated |    145 | 0.855 | [0.789, 0.903] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_finngen |  1,205 | 0.848 | [0.827, 0.867] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_single |    820 | 0.835 | [0.808, 0.859] |
| blood | child | HELIX (non-catalogue pairs) | onco_consensus |    974 | 0.830 | [0.805, 0.852] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_solid_A |     11 | 0.818 | [0.523, 0.949] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_discordant_celltype |    119 | 0.790 | [0.708, 0.853] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak_discordant |     72 | 0.722 | [0.610, 0.812] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_A |    198 | 0.672 | [0.604, 0.733] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_A |    200 | 0.670 | [0.602, 0.731] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak |    192 | 0.635 | [0.565, 0.700] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_consensus_A |      2 | 0.500 | [0.095, 0.905] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_nasal_A |      6 | 1.000 | [0.610, 1.000] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_replicated |     87 | 0.977 | [0.920, 0.994] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_replicated |     29 | 0.966 | [0.828, 0.994] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | onco_consensus |    279 | 0.957 | [0.926, 0.975] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_single |    182 | 0.912 | [0.862, 0.945] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_single |    228 | 0.882 | [0.833, 0.917] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S1 |      8 | 0.875 | [0.529, 0.978] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_high |     74 | 0.824 | [0.722, 0.894] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak_discordant |      4 | 0.750 | [0.301, 0.954] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_nasal_A |     11 | 0.727 | [0.434, 0.903] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak |     26 | 0.692 | [0.500, 0.835] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate |    295 | 0.614 | [0.557, 0.667] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S2 |     29 | 0.586 | [0.407, 0.745] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_solid_A |     37 | 0.568 | [0.409, 0.713] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_solid_A |     39 | 0.564 | [0.410, 0.707] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_blood_A |     56 | 0.554 | [0.424, 0.676] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_blood_A |     63 | 0.524 | [0.403, 0.642] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_finngen |     16 | 0.500 | [0.280, 0.720] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_consensus_A |      4 | 0.250 | [0.046, 0.699] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_s3_replicated |      4 | 0.250 | [0.046, 0.699] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_discordant_celltype |      2 | 0.000 | [0.000, 0.658] |
| solid_tissue | any | FUSION | measured_solid_replicated |     16 | 1.000 | [0.806, 1.000] |
| solid_tissue | any | FUSION | measured_blood_replicated |    936 | 0.939 | [0.922, 0.953] |
| solid_tissue | any | FUSION | measured_solid_single |  1,634 | 0.921 | [0.907, 0.933] |
| solid_tissue | any | FUSION | measured_nasal_single |    707 | 0.885 | [0.860, 0.907] |
| solid_tissue | any | FUSION | onco_consensus |  1,629 | 0.870 | [0.853, 0.885] |
| solid_tissue | any | FUSION | catalogue_model_blood_B |     15 | 0.867 | [0.621, 0.963] |
| solid_tissue | any | FUSION | smr_replicated_S1 |    424 | 0.851 | [0.814, 0.882] |
| solid_tissue | any | FUSION | measured_blood_single |  1,715 | 0.836 | [0.817, 0.852] |
| solid_tissue | any | FUSION | catalogue_single_blood_B |      6 | 0.833 | [0.436, 0.970] |
| solid_tissue | any | FUSION | catalogue_consensus_B |     64 | 0.828 | [0.718, 0.901] |
| solid_tissue | any | FUSION | catalogue_model_nasal_B |    156 | 0.827 | [0.760, 0.878] |
| solid_tissue | any | FUSION | smr_high |    291 | 0.818 | [0.769, 0.858] |
| solid_tissue | any | FUSION | catalogue_single_nasal_B |     93 | 0.817 | [0.727, 0.883] |
| solid_tissue | any | FUSION | smr_weak_discordant |     45 | 0.800 | [0.662, 0.891] |
| solid_tissue | any | FUSION | catalogue_single_solid_A |    290 | 0.776 | [0.724, 0.820] |
| solid_tissue | any | FUSION | catalogue_consensus_A |  1,153 | 0.754 | [0.728, 0.778] |
| solid_tissue | any | FUSION | smr_replicated_S2 |  1,158 | 0.734 | [0.708, 0.759] |
| solid_tissue | any | FUSION | smr_moderate_s3_replicated |     69 | 0.725 | [0.610, 0.816] |
| solid_tissue | any | FUSION | catalogue_single_nasal_A |  2,353 | 0.705 | [0.686, 0.723] |
| solid_tissue | any | FUSION | catalogue_model_solid_A |    828 | 0.703 | [0.671, 0.733] |
| solid_tissue | any | FUSION | catalogue_model_nasal_A |  3,774 | 0.703 | [0.688, 0.717] |
| solid_tissue | any | FUSION | smr_moderate_finngen |    941 | 0.683 | [0.653, 0.712] |
| solid_tissue | any | FUSION | smr_moderate |  1,436 | 0.681 | [0.656, 0.705] |
| solid_tissue | any | FUSION | catalogue_model_blood_A |  2,720 | 0.679 | [0.661, 0.696] |
| solid_tissue | any | FUSION | distance_only | 16,111 | 0.667 | [0.660, 0.675] |
| solid_tissue | any | FUSION | catalogue_single_blood_A |  1,404 | 0.665 | [0.640, 0.689] |
| solid_tissue | any | FUSION | distance_curve_blood | 29,694 | 0.635 | [0.630, 0.641] |
| solid_tissue | any | FUSION | smr_moderate_discordant_celltype |     87 | 0.632 | [0.527, 0.726] |
| solid_tissue | any | FUSION | distance_curve_nasal | 29,694 | 0.613 | [0.607, 0.619] |
| solid_tissue | any | FUSION | smr_weak |     92 | 0.587 | [0.485, 0.682] |
| solid_tissue | any | FUSION | distance_curve_solid | 29,694 | 0.546 | [0.541, 0.552] |
| solid_tissue | any | GTEx solid | measured_solid_replicated |     16 | 1.000 | [0.806, 1.000] |
| solid_tissue | any | GTEx solid | measured_blood_replicated |    825 | 0.988 | [0.978, 0.993] |
| solid_tissue | any | GTEx solid | measured_nasal_single |    481 | 0.958 | [0.937, 0.973] |
| solid_tissue | any | GTEx solid | smr_replicated_S1 |    265 | 0.955 | [0.923, 0.974] |
| solid_tissue | any | GTEx solid | smr_replicated_S2 |    425 | 0.941 | [0.915, 0.960] |
| solid_tissue | any | GTEx solid | smr_moderate_s3_replicated |     16 | 0.938 | [0.717, 0.989] |
| solid_tissue | any | GTEx solid | measured_blood_single |  1,268 | 0.935 | [0.920, 0.947] |
| solid_tissue | any | GTEx solid | measured_solid_single |  1,634 | 0.922 | [0.908, 0.934] |
| solid_tissue | any | GTEx solid | smr_weak_discordant |     24 | 0.917 | [0.742, 0.977] |
| solid_tissue | any | GTEx solid | onco_consensus |    551 | 0.900 | [0.872, 0.923] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_B |     19 | 0.895 | [0.686, 0.971] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_B |     34 | 0.853 | [0.699, 0.936] |
| solid_tissue | any | GTEx solid | smr_high |    159 | 0.849 | [0.785, 0.896] |
| solid_tissue | any | GTEx solid | catalogue_single_solid_A |     56 | 0.821 | [0.702, 0.900] |
| solid_tissue | any | GTEx solid | smr_moderate_finngen |    302 | 0.815 | [0.767, 0.854] |
| solid_tissue | any | GTEx solid | catalogue_consensus_B |     16 | 0.812 | [0.570, 0.934] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_B |      5 | 0.800 | [0.376, 0.964] |
| solid_tissue | any | GTEx solid | smr_moderate |    460 | 0.761 | [0.720, 0.798] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_A |    633 | 0.749 | [0.714, 0.781] |
| solid_tissue | any | GTEx solid | catalogue_model_solid_A |    156 | 0.737 | [0.663, 0.800] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_A |    297 | 0.734 | [0.681, 0.781] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_A |    491 | 0.723 | [0.682, 0.761] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_A |    866 | 0.715 | [0.684, 0.744] |
| solid_tissue | any | GTEx solid | catalogue_consensus_A |    194 | 0.711 | [0.644, 0.771] |
| solid_tissue | any | GTEx solid | smr_weak |     33 | 0.667 | [0.496, 0.802] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_B |      3 | 0.667 | [0.208, 0.939] |
| solid_tissue | any | GTEx solid | distance_only |  5,315 | 0.575 | [0.562, 0.588] |
| solid_tissue | any | GTEx solid | distance_curve_nasal |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_blood |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_solid |  9,435 | 0.540 | [0.530, 0.550] |
| solid_tissue | any | GTEx solid | smr_moderate_discordant_celltype |     15 | 0.533 | [0.301, 0.752] |
