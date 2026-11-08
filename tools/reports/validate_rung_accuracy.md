# Every ladder source, scored in each target tissue against an independent yardstick

Run 2026-11-08 11:06. Agreement with the yardstick's measured sign, Wilson 95% CI. Sources are scored standalone, on every yardstick pair where they call (catalogue sources: formed from predicted calls only, as the ladders form them). `n_cpg`: distinct CpGs.

Yardsticks: adult blood = MESA monocytes + T cells; child blood = HELIX Dryad pairs outside the blood catalogue's DIRECT set; nasal = EVA-PR nasal pairs outside the nasal catalogue's DIRECT set (the nasal catalogue is trained on EVA-PR: its rows are same-study); solid = FUSION and GTEx solid, each held out in turn, pooled, each pair once. Measured cohorts are counted by donor group. Distance is not scored against HELIX or EVA-PR (its blood and nasal curves were fitted on them).

## blood (adult): MESA

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_model_blood_B |     8 |     8 | 1.000 | [0.676, 1.000] |
| 2 | catalogue_single_blood_B |     4 |     4 | 1.000 | [0.510, 1.000] |
| 3 | measured_blood_replicated |   435 |   432 | 0.970 | [0.950, 0.982] |
| 4 | catalogue_consensus_B |    20 |    20 | 0.950 | [0.764, 0.991] |
| 5 | measured_solid_replicated |    78 |    78 | 0.936 | [0.859, 0.972] |
| 6 | smr_replicated_S1 |   244 |   240 | 0.926 | [0.886, 0.953] |
| 7 | smr_weak_discordant |    23 |    23 | 0.913 | [0.732, 0.976] |
| 8 | smr_high |   239 |   239 | 0.908 | [0.865, 0.938] |
| 9 | measured_blood_single |   903 |   880 | 0.888 | [0.866, 0.907] |
| 10 | onco_consensus |   766 |   763 | 0.850 | [0.823, 0.873] |
| 11 | catalogue_model_nasal_B |    46 |    46 | 0.848 | [0.718, 0.924] |
| 12 | measured_nasal_single |   326 |   313 | 0.847 | [0.803, 0.882] |
| 13 | smr_replicated_S2 |   733 |   721 | 0.846 | [0.818, 0.870] |
| 14 | catalogue_single_nasal_B |    30 |    30 | 0.800 | [0.627, 0.905] |
| 15 | smr_moderate | 1,179 | 1,147 | 0.797 | [0.773, 0.819] |
| 16 | catalogue_consensus_A |   620 |   600 | 0.779 | [0.745, 0.810] |
| 17 | measured_solid_single |   379 |   370 | 0.749 | [0.703, 0.790] |
| 18 | smr_moderate_s3_replicated |    43 |    41 | 0.721 | [0.573, 0.833] |
| 19 | smr_moderate_discordant_celltype |    70 |    70 | 0.714 | [0.599, 0.807] |
| 20 | catalogue_model_blood_A | 1,590 | 1,515 | 0.713 | [0.690, 0.735] |
| 21 | smr_moderate_finngen |   424 |   410 | 0.710 | [0.665, 0.751] |
| 22 | catalogue_single_blood_A |   865 |   849 | 0.692 | [0.661, 0.722] |
| 23 | catalogue_model_solid_A |   549 |   536 | 0.670 | [0.630, 0.708] |
| 24 | distance_only | 5,196 | 4,884 | 0.668 | [0.655, 0.681] |
| 25 | catalogue_model_nasal_A | 1,500 | 1,447 | 0.660 | [0.636, 0.684] |
| 26 | distance_curve_nasal | 7,025 | 6,584 | 0.640 | [0.628, 0.651] |
| 27 | distance_curve_blood | 7,025 | 6,584 | 0.639 | [0.628, 0.650] |
| 28 | catalogue_single_nasal_A |   797 |   789 | 0.629 | [0.595, 0.661] |
| 29 | distance_curve_solid | 7,025 | 6,584 | 0.610 | [0.598, 0.621] |
| 30 | catalogue_single_solid_A |   227 |   221 | 0.577 | [0.512, 0.640] |
| 31 | smr_weak |    73 |    72 | 0.575 | [0.461, 0.682] |

## blood (child): HELIX (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | smr_weak_discordant |    11 |  11 | 1.000 | [0.741, 1.000] |
| 2 | measured_blood_replicated |   271 | 240 | 0.993 | [0.973, 0.998] |
| 3 | measured_solid_replicated |    75 |  72 | 0.987 | [0.928, 0.998] |
| 4 | smr_high |   414 | 362 | 0.969 | [0.947, 0.982] |
| 5 | measured_blood_single |   634 | 533 | 0.964 | [0.946, 0.976] |
| 6 | measured_nasal_single |   241 | 210 | 0.905 | [0.861, 0.936] |
| 7 | smr_replicated_S1 |   164 | 148 | 0.902 | [0.847, 0.939] |
| 8 | smr_moderate | 1,049 | 910 | 0.899 | [0.879, 0.916] |
| 9 | smr_replicated_S2 |   285 | 252 | 0.891 | [0.850, 0.922] |
| 10 | onco_consensus |   277 | 239 | 0.838 | [0.790, 0.876] |
| 11 | measured_solid_single |   356 | 306 | 0.820 | [0.777, 0.857] |
| 12 | smr_moderate_finngen |   191 | 168 | 0.780 | [0.716, 0.833] |
| 13 | smr_moderate_s3_replicated |    39 |  38 | 0.718 | [0.562, 0.835] |
| 14 | catalogue_consensus_A |     3 |   3 | 0.667 | [0.208, 0.939] |
| 15 | catalogue_model_nasal_A |   283 | 271 | 0.597 | [0.539, 0.653] |
| 16 | catalogue_single_nasal_A |   279 | 268 | 0.595 | [0.536, 0.651] |
| 17 | smr_weak |   111 |  99 | 0.577 | [0.484, 0.664] |
| 18 | smr_moderate_discordant_celltype |    10 |  10 | 0.500 | [0.237, 0.763] |
| 19 | catalogue_model_solid_A |    63 |  63 | 0.413 | [0.300, 0.536] |
| 20 | catalogue_single_solid_A |    59 |  59 | 0.407 | [0.291, 0.534] |

## nasal_epithelium (any): EVA-PR nasal (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_single_nasal_A |   7 |   7 | 1.000 | [0.646, 1.000] |
| 2 | measured_blood_replicated |  84 |  77 | 0.976 | [0.917, 0.993] |
| 3 | measured_solid_replicated |  29 |  25 | 0.966 | [0.828, 0.994] |
| 4 | onco_consensus | 282 | 215 | 0.961 | [0.932, 0.978] |
| 5 | measured_solid_single | 180 | 128 | 0.911 | [0.860, 0.945] |
| 6 | measured_blood_single | 224 | 194 | 0.879 | [0.830, 0.916] |
| 7 | smr_high |  73 |  65 | 0.836 | [0.734, 0.903] |
| 8 | smr_weak_discordant |   4 |   4 | 0.750 | [0.301, 0.954] |
| 9 | catalogue_model_nasal_A |  10 |  10 | 0.700 | [0.397, 0.892] |
| 10 | smr_weak |  26 |  25 | 0.692 | [0.500, 0.835] |
| 11 | smr_replicated_S1 |   3 |   3 | 0.667 | [0.208, 0.939] |
| 12 | smr_moderate | 282 | 242 | 0.613 | [0.555, 0.668] |
| 13 | catalogue_single_solid_A |  25 |  25 | 0.600 | [0.407, 0.766] |
| 14 | smr_replicated_S2 |  27 |  27 | 0.593 | [0.407, 0.755] |
| 15 | catalogue_model_solid_A |  26 |  26 | 0.577 | [0.389, 0.745] |
| 16 | catalogue_single_blood_A |  29 |  29 | 0.552 | [0.375, 0.716] |
| 17 | catalogue_model_blood_A |  33 |  33 | 0.515 | [0.352, 0.675] |
| 18 | smr_moderate_finngen |  13 |  13 | 0.385 | [0.177, 0.645] |
| 19 | smr_moderate_s3_replicated |   4 |   4 | 0.250 | [0.046, 0.699] |
| 20 | catalogue_consensus_A |   3 |   3 | 0.000 | [0.000, 0.561] |
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
| 8 | catalogue_model_nasal_B |    198 |    198 | 0.828 | [0.770, 0.874] |
| 9 | catalogue_single_nasal_B |    114 |    114 | 0.825 | [0.744, 0.883] |
| 10 | catalogue_consensus_B |     89 |     89 | 0.820 | [0.728, 0.886] |
| 11 | smr_high |    406 |    359 | 0.815 | [0.775, 0.850] |
| 12 | catalogue_model_blood_B |     43 |     43 | 0.814 | [0.674, 0.903] |
| 13 | smr_weak_discordant |     55 |     54 | 0.800 | [0.676, 0.884] |
| 14 | catalogue_single_blood_B |     28 |     28 | 0.786 | [0.605, 0.898] |
| 15 | smr_replicated_S2 |  1,434 |  1,322 | 0.770 | [0.747, 0.791] |
| 16 | smr_moderate_s3_replicated |     84 |     82 | 0.762 | [0.661, 0.840] |
| 17 | catalogue_consensus_A |  1,432 |  1,366 | 0.758 | [0.735, 0.779] |
| 18 | catalogue_model_blood_A |  5,407 |  4,839 | 0.707 | [0.694, 0.719] |
| 19 | catalogue_single_blood_A |  3,748 |  3,434 | 0.706 | [0.691, 0.721] |
| 20 | catalogue_single_nasal_A |  2,955 |  2,770 | 0.706 | [0.689, 0.722] |
| 21 | smr_moderate_finngen |  1,130 |  1,045 | 0.705 | [0.678, 0.731] |
| 22 | catalogue_model_nasal_A |  4,612 |  4,201 | 0.699 | [0.686, 0.712] |
| 23 | smr_moderate |  1,800 |  1,610 | 0.694 | [0.672, 0.715] |
| 24 | catalogue_model_solid_A |  1,583 |  1,550 | 0.675 | [0.652, 0.698] |
| 25 | catalogue_single_solid_A |    792 |    779 | 0.668 | [0.634, 0.700] |
| 26 | distance_only | 20,312 | 15,989 | 0.643 | [0.637, 0.650] |
| 27 | smr_moderate_discordant_celltype |     97 |     95 | 0.629 | [0.530, 0.718] |
| 28 | distance_curve_blood | 37,489 | 25,837 | 0.613 | [0.608, 0.618] |
| 29 | smr_weak |    113 |    108 | 0.602 | [0.510, 0.687] |
| 30 | distance_curve_nasal | 37,489 | 25,837 | 0.595 | [0.590, 0.600] |
| 31 | distance_curve_solid | 37,489 | 25,837 | 0.542 | [0.537, 0.547] |

## By fold

| target | age | yardstick | source | n | agree | 95% CI |
|---|---|---|---|---|---|---|
| blood | adult | MESA | catalogue_model_blood_B |      8 | 1.000 | [0.676, 1.000] |
| blood | adult | MESA | catalogue_single_blood_B |      4 | 1.000 | [0.510, 1.000] |
| blood | adult | MESA | measured_blood_replicated |    435 | 0.970 | [0.950, 0.982] |
| blood | adult | MESA | catalogue_consensus_B |     20 | 0.950 | [0.764, 0.991] |
| blood | adult | MESA | measured_solid_replicated |     78 | 0.936 | [0.859, 0.972] |
| blood | adult | MESA | smr_replicated_S1 |    244 | 0.926 | [0.886, 0.953] |
| blood | adult | MESA | smr_weak_discordant |     23 | 0.913 | [0.732, 0.976] |
| blood | adult | MESA | smr_high |    239 | 0.908 | [0.865, 0.938] |
| blood | adult | MESA | measured_blood_single |    903 | 0.888 | [0.866, 0.907] |
| blood | adult | MESA | onco_consensus |    766 | 0.850 | [0.823, 0.873] |
| blood | adult | MESA | catalogue_model_nasal_B |     46 | 0.848 | [0.718, 0.924] |
| blood | adult | MESA | measured_nasal_single |    326 | 0.847 | [0.803, 0.882] |
| blood | adult | MESA | smr_replicated_S2 |    733 | 0.846 | [0.818, 0.870] |
| blood | adult | MESA | catalogue_single_nasal_B |     30 | 0.800 | [0.627, 0.905] |
| blood | adult | MESA | smr_moderate |  1,179 | 0.797 | [0.773, 0.819] |
| blood | adult | MESA | catalogue_consensus_A |    620 | 0.779 | [0.745, 0.810] |
| blood | adult | MESA | measured_solid_single |    379 | 0.749 | [0.703, 0.790] |
| blood | adult | MESA | smr_moderate_s3_replicated |     43 | 0.721 | [0.573, 0.833] |
| blood | adult | MESA | smr_moderate_discordant_celltype |     70 | 0.714 | [0.599, 0.807] |
| blood | adult | MESA | catalogue_model_blood_A |  1,590 | 0.713 | [0.690, 0.735] |
| blood | adult | MESA | smr_moderate_finngen |    424 | 0.710 | [0.665, 0.751] |
| blood | adult | MESA | catalogue_single_blood_A |    865 | 0.692 | [0.661, 0.722] |
| blood | adult | MESA | catalogue_model_solid_A |    549 | 0.670 | [0.630, 0.708] |
| blood | adult | MESA | distance_only |  5,196 | 0.668 | [0.655, 0.681] |
| blood | adult | MESA | catalogue_model_nasal_A |  1,500 | 0.660 | [0.636, 0.684] |
| blood | adult | MESA | distance_curve_nasal |  7,025 | 0.640 | [0.628, 0.651] |
| blood | adult | MESA | distance_curve_blood |  7,025 | 0.639 | [0.628, 0.650] |
| blood | adult | MESA | catalogue_single_nasal_A |    797 | 0.629 | [0.595, 0.661] |
| blood | adult | MESA | distance_curve_solid |  7,025 | 0.610 | [0.598, 0.621] |
| blood | adult | MESA | catalogue_single_solid_A |    227 | 0.577 | [0.512, 0.640] |
| blood | adult | MESA | smr_weak |     73 | 0.575 | [0.461, 0.682] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak_discordant |     11 | 1.000 | [0.741, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_replicated |    271 | 0.993 | [0.973, 0.998] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_replicated |     75 | 0.987 | [0.928, 0.998] |
| blood | child | HELIX (non-catalogue pairs) | smr_high |    414 | 0.969 | [0.947, 0.982] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_single |    634 | 0.964 | [0.946, 0.976] |
| blood | child | HELIX (non-catalogue pairs) | measured_nasal_single |    241 | 0.905 | [0.861, 0.936] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S1 |    164 | 0.902 | [0.847, 0.939] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate |  1,049 | 0.899 | [0.879, 0.916] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S2 |    285 | 0.891 | [0.850, 0.922] |
| blood | child | HELIX (non-catalogue pairs) | onco_consensus |    277 | 0.838 | [0.790, 0.876] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_single |    356 | 0.820 | [0.777, 0.857] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_finngen |    191 | 0.780 | [0.716, 0.833] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_s3_replicated |     39 | 0.718 | [0.562, 0.835] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_consensus_A |      3 | 0.667 | [0.208, 0.939] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_A |    283 | 0.597 | [0.539, 0.653] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_A |    279 | 0.595 | [0.536, 0.651] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak |    111 | 0.577 | [0.484, 0.664] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_discordant_celltype |     10 | 0.500 | [0.237, 0.763] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_solid_A |     63 | 0.413 | [0.300, 0.536] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_solid_A |     59 | 0.407 | [0.291, 0.534] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_nasal_A |      7 | 1.000 | [0.646, 1.000] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_replicated |     84 | 0.976 | [0.917, 0.993] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_replicated |     29 | 0.966 | [0.828, 0.994] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | onco_consensus |    282 | 0.961 | [0.932, 0.978] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_single |    180 | 0.911 | [0.860, 0.945] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_single |    224 | 0.879 | [0.830, 0.916] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_high |     73 | 0.836 | [0.734, 0.903] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak_discordant |      4 | 0.750 | [0.301, 0.954] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_nasal_A |     10 | 0.700 | [0.397, 0.892] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak |     26 | 0.692 | [0.500, 0.835] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S1 |      3 | 0.667 | [0.208, 0.939] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate |    282 | 0.613 | [0.555, 0.668] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_solid_A |     25 | 0.600 | [0.407, 0.766] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S2 |     27 | 0.593 | [0.407, 0.755] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_solid_A |     26 | 0.577 | [0.389, 0.745] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_blood_A |     29 | 0.552 | [0.375, 0.716] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_blood_A |     33 | 0.515 | [0.352, 0.675] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_finngen |     13 | 0.385 | [0.177, 0.645] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_s3_replicated |      4 | 0.250 | [0.046, 0.699] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_consensus_A |      3 | 0.000 | [0.000, 0.561] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_discordant_celltype |      2 | 0.000 | [0.000, 0.658] |
| solid_tissue | any | FUSION | measured_solid_replicated |     16 | 1.000 | [0.806, 1.000] |
| solid_tissue | any | FUSION | measured_blood_replicated |    936 | 0.939 | [0.922, 0.953] |
| solid_tissue | any | FUSION | measured_solid_single |  1,634 | 0.921 | [0.907, 0.933] |
| solid_tissue | any | FUSION | measured_nasal_single |    707 | 0.885 | [0.860, 0.907] |
| solid_tissue | any | FUSION | onco_consensus |  1,629 | 0.870 | [0.853, 0.885] |
| solid_tissue | any | FUSION | smr_replicated_S1 |    424 | 0.851 | [0.814, 0.882] |
| solid_tissue | any | FUSION | measured_blood_single |  1,715 | 0.836 | [0.817, 0.852] |
| solid_tissue | any | FUSION | catalogue_model_blood_B |     36 | 0.833 | [0.681, 0.921] |
| solid_tissue | any | FUSION | catalogue_single_blood_B |     23 | 0.826 | [0.629, 0.930] |
| solid_tissue | any | FUSION | catalogue_model_nasal_B |    171 | 0.825 | [0.761, 0.874] |
| solid_tissue | any | FUSION | catalogue_consensus_B |     73 | 0.822 | [0.719, 0.893] |
| solid_tissue | any | FUSION | smr_high |    291 | 0.818 | [0.769, 0.858] |
| solid_tissue | any | FUSION | catalogue_single_nasal_B |    102 | 0.814 | [0.727, 0.877] |
| solid_tissue | any | FUSION | smr_weak_discordant |     45 | 0.800 | [0.662, 0.891] |
| solid_tissue | any | FUSION | catalogue_consensus_A |  1,268 | 0.761 | [0.737, 0.784] |
| solid_tissue | any | FUSION | smr_replicated_S2 |  1,158 | 0.734 | [0.708, 0.759] |
| solid_tissue | any | FUSION | smr_moderate_s3_replicated |     69 | 0.725 | [0.610, 0.816] |
| solid_tissue | any | FUSION | catalogue_single_nasal_A |  2,500 | 0.706 | [0.688, 0.724] |
| solid_tissue | any | FUSION | catalogue_model_nasal_A |  3,971 | 0.703 | [0.688, 0.717] |
| solid_tissue | any | FUSION | catalogue_model_blood_A |  4,616 | 0.698 | [0.685, 0.711] |
| solid_tissue | any | FUSION | catalogue_single_blood_A |  3,132 | 0.695 | [0.679, 0.711] |
| solid_tissue | any | FUSION | smr_moderate_finngen |    941 | 0.683 | [0.653, 0.712] |
| solid_tissue | any | FUSION | smr_moderate |  1,436 | 0.681 | [0.656, 0.705] |
| solid_tissue | any | FUSION | catalogue_model_solid_A |  1,385 | 0.676 | [0.651, 0.700] |
| solid_tissue | any | FUSION | catalogue_single_solid_A |    688 | 0.673 | [0.637, 0.707] |
| solid_tissue | any | FUSION | distance_only | 16,111 | 0.667 | [0.660, 0.675] |
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
| solid_tissue | any | GTEx solid | catalogue_model_nasal_B |     35 | 0.857 | [0.706, 0.937] |
| solid_tissue | any | GTEx solid | smr_high |    159 | 0.849 | [0.785, 0.896] |
| solid_tissue | any | GTEx solid | catalogue_consensus_B |     17 | 0.824 | [0.590, 0.938] |
| solid_tissue | any | GTEx solid | smr_moderate_finngen |    302 | 0.815 | [0.767, 0.854] |
| solid_tissue | any | GTEx solid | smr_moderate |    460 | 0.761 | [0.720, 0.798] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_A |    752 | 0.758 | [0.726, 0.787] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_A |    970 | 0.747 | [0.719, 0.774] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_A |    649 | 0.738 | [0.703, 0.770] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_B |      7 | 0.714 | [0.359, 0.918] |
| solid_tissue | any | GTEx solid | catalogue_consensus_A |    213 | 0.714 | [0.650, 0.770] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_A |    886 | 0.708 | [0.677, 0.737] |
| solid_tissue | any | GTEx solid | catalogue_model_solid_A |    260 | 0.677 | [0.618, 0.731] |
| solid_tissue | any | GTEx solid | catalogue_single_solid_A |    139 | 0.676 | [0.595, 0.748] |
| solid_tissue | any | GTEx solid | smr_weak |     33 | 0.667 | [0.496, 0.802] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_B |      5 | 0.600 | [0.231, 0.882] |
| solid_tissue | any | GTEx solid | distance_only |  5,315 | 0.575 | [0.562, 0.588] |
| solid_tissue | any | GTEx solid | distance_curve_nasal |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_blood |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_solid |  9,435 | 0.540 | [0.530, 0.550] |
| solid_tissue | any | GTEx solid | smr_moderate_discordant_celltype |     15 | 0.533 | [0.301, 0.752] |
