# Every ladder source, scored in each target tissue against an independent yardstick

Run 2026-11-09 07:02. Agreement with the yardstick's measured sign, Wilson 95% CI. Sources are scored standalone, on every yardstick pair where they call (catalogue sources: formed from predicted calls only, as the ladders form them). `n_cpg`: distinct CpGs.

Yardsticks: adult blood = MESA monocytes + T cells; child blood = HELIX Dryad pairs outside the blood catalogue's DIRECT set; nasal = EVA-PR nasal pairs outside the nasal catalogue's DIRECT set (the nasal catalogue is trained on EVA-PR: its rows are same-study); solid = FUSION and GTEx solid, each held out in turn, pooled, each pair once. Measured cohorts are counted by donor group. Distance is not scored against HELIX or EVA-PR (its blood and nasal curves were fitted on them).

## blood (adult): MESA

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_model_blood_B |    10 |    10 | 1.000 | [0.722, 1.000] |
| 2 | catalogue_single_blood_B |     4 |     4 | 1.000 | [0.510, 1.000] |
| 3 | measured_blood_replicated |   435 |   432 | 0.970 | [0.950, 0.982] |
| 4 | measured_solid_replicated |    78 |    78 | 0.936 | [0.859, 0.972] |
| 5 | smr_replicated_S1 |   244 |   240 | 0.926 | [0.886, 0.953] |
| 6 | smr_weak_discordant |    23 |    23 | 0.913 | [0.732, 0.976] |
| 7 | smr_high |   239 |   239 | 0.908 | [0.865, 0.938] |
| 8 | measured_blood_single |   903 |   880 | 0.888 | [0.866, 0.907] |
| 9 | onco_consensus |   766 |   763 | 0.850 | [0.823, 0.873] |
| 10 | measured_nasal_single |   326 |   313 | 0.847 | [0.803, 0.882] |
| 11 | smr_replicated_S2 |   733 |   721 | 0.846 | [0.818, 0.870] |
| 12 | catalogue_single_nasal_B |    48 |    48 | 0.833 | [0.704, 0.913] |
| 13 | catalogue_consensus_B |    30 |    30 | 0.833 | [0.664, 0.927] |
| 14 | smr_moderate | 1,179 | 1,147 | 0.797 | [0.773, 0.819] |
| 15 | catalogue_model_nasal_B |    78 |    78 | 0.795 | [0.692, 0.870] |
| 16 | catalogue_consensus_A |   641 |   617 | 0.793 | [0.759, 0.822] |
| 17 | measured_solid_single |   379 |   370 | 0.749 | [0.703, 0.790] |
| 18 | smr_moderate_s3_replicated |    43 |    41 | 0.721 | [0.573, 0.833] |
| 19 | catalogue_model_blood_A | 1,587 | 1,518 | 0.715 | [0.692, 0.737] |
| 20 | smr_moderate_discordant_celltype |    70 |    70 | 0.714 | [0.599, 0.807] |
| 21 | smr_moderate_finngen |   424 |   410 | 0.710 | [0.665, 0.751] |
| 22 | catalogue_single_blood_A |   880 |   864 | 0.675 | [0.643, 0.705] |
| 23 | catalogue_model_nasal_A | 1,390 | 1,337 | 0.668 | [0.643, 0.693] |
| 24 | distance_only | 5,196 | 4,884 | 0.668 | [0.655, 0.681] |
| 25 | distance_curve_nasal | 7,025 | 6,584 | 0.640 | [0.628, 0.651] |
| 26 | distance_curve_blood | 7,025 | 6,584 | 0.639 | [0.628, 0.650] |
| 27 | catalogue_model_solid_A |   769 |   747 | 0.627 | [0.592, 0.660] |
| 28 | catalogue_single_nasal_A |   728 |   720 | 0.615 | [0.580, 0.650] |
| 29 | distance_curve_solid | 7,025 | 6,584 | 0.610 | [0.598, 0.621] |
| 30 | smr_weak |    73 |    72 | 0.575 | [0.461, 0.682] |
| 31 | catalogue_single_solid_A |   329 |   322 | 0.547 | [0.493, 0.600] |

## blood (child): HELIX (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | smr_weak_discordant |  11 |  11 | 1.000 | [0.741, 1.000] |
| 2 | catalogue_consensus_A |   4 |   4 | 1.000 | [0.510, 1.000] |
| 3 | catalogue_single_nasal_B |   1 |   1 | 1.000 | [0.207, 1.000] |
| 4 | catalogue_model_nasal_B |   1 |   1 | 1.000 | [0.207, 1.000] |
| 5 | measured_blood_replicated | 266 | 236 | 0.992 | [0.973, 0.998] |
| 6 | measured_solid_replicated |  74 |  71 | 0.986 | [0.927, 0.998] |
| 7 | smr_high | 392 | 343 | 0.967 | [0.944, 0.981] |
| 8 | measured_blood_single | 625 | 535 | 0.965 | [0.947, 0.977] |
| 9 | measured_nasal_single | 240 | 212 | 0.904 | [0.860, 0.935] |
| 10 | smr_moderate | 995 | 861 | 0.895 | [0.875, 0.913] |
| 11 | smr_replicated_S2 | 275 | 254 | 0.891 | [0.849, 0.923] |
| 12 | smr_replicated_S1 | 122 | 117 | 0.885 | [0.817, 0.930] |
| 13 | onco_consensus | 286 | 255 | 0.839 | [0.792, 0.877] |
| 14 | measured_solid_single | 327 | 283 | 0.826 | [0.781, 0.863] |
| 15 | smr_moderate_finngen | 172 | 158 | 0.767 | [0.699, 0.824] |
| 16 | smr_moderate_discordant_celltype |   6 |   6 | 0.667 | [0.300, 0.903] |
| 17 | smr_moderate_s3_replicated |  31 |  30 | 0.613 | [0.438, 0.763] |
| 18 | smr_weak | 108 |  97 | 0.583 | [0.489, 0.672] |
| 19 | catalogue_model_nasal_A | 180 | 174 | 0.550 | [0.477, 0.621] |
| 20 | catalogue_single_nasal_A | 175 | 169 | 0.543 | [0.469, 0.615] |
| 21 | catalogue_model_solid_A |  81 |  79 | 0.444 | [0.341, 0.553] |
| 22 | catalogue_single_solid_A |  76 |  74 | 0.408 | [0.304, 0.520] |

## nasal_epithelium (any): EVA-PR nasal (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_single_nasal_A |   7 |   7 | 1.000 | [0.646, 1.000] |
| 2 | measured_blood_replicated |  85 |  78 | 0.976 | [0.918, 0.994] |
| 3 | measured_solid_replicated |  29 |  25 | 0.966 | [0.828, 0.994] |
| 4 | onco_consensus | 278 | 214 | 0.960 | [0.931, 0.978] |
| 5 | measured_solid_single | 179 | 129 | 0.911 | [0.860, 0.944] |
| 6 | measured_blood_single | 236 | 206 | 0.886 | [0.839, 0.920] |
| 7 | smr_high |  73 |  65 | 0.836 | [0.734, 0.903] |
| 8 | catalogue_model_nasal_A |   9 |   9 | 0.778 | [0.453, 0.937] |
| 9 | smr_weak_discordant |   4 |   4 | 0.750 | [0.301, 0.954] |
| 10 | smr_weak |  26 |  25 | 0.692 | [0.500, 0.835] |
| 11 | smr_replicated_S1 |   3 |   3 | 0.667 | [0.208, 0.939] |
| 12 | smr_moderate | 279 | 239 | 0.609 | [0.551, 0.665] |
| 13 | smr_replicated_S2 |  28 |  28 | 0.607 | [0.424, 0.764] |
| 14 | catalogue_single_solid_A |  55 |  53 | 0.545 | [0.415, 0.670] |
| 15 | catalogue_model_solid_A |  57 |  55 | 0.544 | [0.416, 0.666] |
| 16 | catalogue_single_blood_A |  47 |  46 | 0.468 | [0.333, 0.608] |
| 17 | catalogue_model_blood_A |  51 |  50 | 0.431 | [0.305, 0.567] |
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
| 5 | catalogue_single_blood_B |     27 |     27 | 0.889 | [0.719, 0.961] |
| 6 | onco_consensus |  1,977 |  1,849 | 0.874 | [0.859, 0.888] |
| 7 | smr_replicated_S1 |    576 |    528 | 0.873 | [0.844, 0.898] |
| 8 | measured_blood_single |  2,626 |  2,219 | 0.865 | [0.852, 0.878] |
| 9 | catalogue_model_blood_B |     42 |     42 | 0.857 | [0.722, 0.933] |
| 10 | catalogue_consensus_B |    118 |    118 | 0.831 | [0.753, 0.888] |
| 11 | catalogue_model_nasal_B |    284 |    284 | 0.827 | [0.779, 0.867] |
| 12 | smr_high |    406 |    359 | 0.815 | [0.775, 0.850] |
| 13 | catalogue_single_nasal_B |    168 |    168 | 0.810 | [0.743, 0.862] |
| 14 | smr_weak_discordant |     55 |     54 | 0.800 | [0.676, 0.884] |
| 15 | smr_replicated_S2 |  1,434 |  1,322 | 0.770 | [0.747, 0.791] |
| 16 | smr_moderate_s3_replicated |     84 |     82 | 0.762 | [0.661, 0.840] |
| 17 | catalogue_consensus_A |  1,572 |  1,504 | 0.747 | [0.725, 0.768] |
| 18 | catalogue_single_nasal_A |  2,621 |  2,472 | 0.719 | [0.701, 0.736] |
| 19 | smr_moderate_finngen |  1,130 |  1,045 | 0.705 | [0.678, 0.731] |
| 20 | catalogue_model_nasal_A |  4,264 |  3,898 | 0.705 | [0.691, 0.718] |
| 21 | catalogue_model_blood_A |  5,539 |  4,949 | 0.703 | [0.691, 0.715] |
| 22 | catalogue_single_blood_A |  3,725 |  3,401 | 0.697 | [0.682, 0.711] |
| 23 | smr_moderate |  1,800 |  1,610 | 0.694 | [0.672, 0.715] |
| 24 | distance_only | 20,312 | 15,989 | 0.643 | [0.637, 0.650] |
| 25 | catalogue_model_solid_A |  2,298 |  2,222 | 0.639 | [0.619, 0.658] |
| 26 | smr_moderate_discordant_celltype |     97 |     95 | 0.629 | [0.530, 0.718] |
| 27 | catalogue_single_solid_A |  1,213 |  1,193 | 0.627 | [0.599, 0.653] |
| 28 | distance_curve_blood | 37,489 | 25,837 | 0.613 | [0.608, 0.618] |
| 29 | smr_weak |    113 |    108 | 0.602 | [0.510, 0.687] |
| 30 | distance_curve_nasal | 37,489 | 25,837 | 0.595 | [0.590, 0.600] |
| 31 | distance_curve_solid | 37,489 | 25,837 | 0.542 | [0.537, 0.547] |

## By fold

| target | age | yardstick | source | n | agree | 95% CI |
|---|---|---|---|---|---|---|
| blood | adult | MESA | catalogue_model_blood_B |     10 | 1.000 | [0.722, 1.000] |
| blood | adult | MESA | catalogue_single_blood_B |      4 | 1.000 | [0.510, 1.000] |
| blood | adult | MESA | measured_blood_replicated |    435 | 0.970 | [0.950, 0.982] |
| blood | adult | MESA | measured_solid_replicated |     78 | 0.936 | [0.859, 0.972] |
| blood | adult | MESA | smr_replicated_S1 |    244 | 0.926 | [0.886, 0.953] |
| blood | adult | MESA | smr_weak_discordant |     23 | 0.913 | [0.732, 0.976] |
| blood | adult | MESA | smr_high |    239 | 0.908 | [0.865, 0.938] |
| blood | adult | MESA | measured_blood_single |    903 | 0.888 | [0.866, 0.907] |
| blood | adult | MESA | onco_consensus |    766 | 0.850 | [0.823, 0.873] |
| blood | adult | MESA | measured_nasal_single |    326 | 0.847 | [0.803, 0.882] |
| blood | adult | MESA | smr_replicated_S2 |    733 | 0.846 | [0.818, 0.870] |
| blood | adult | MESA | catalogue_single_nasal_B |     48 | 0.833 | [0.704, 0.913] |
| blood | adult | MESA | catalogue_consensus_B |     30 | 0.833 | [0.664, 0.927] |
| blood | adult | MESA | smr_moderate |  1,179 | 0.797 | [0.773, 0.819] |
| blood | adult | MESA | catalogue_model_nasal_B |     78 | 0.795 | [0.692, 0.870] |
| blood | adult | MESA | catalogue_consensus_A |    641 | 0.793 | [0.759, 0.822] |
| blood | adult | MESA | measured_solid_single |    379 | 0.749 | [0.703, 0.790] |
| blood | adult | MESA | smr_moderate_s3_replicated |     43 | 0.721 | [0.573, 0.833] |
| blood | adult | MESA | catalogue_model_blood_A |  1,587 | 0.715 | [0.692, 0.737] |
| blood | adult | MESA | smr_moderate_discordant_celltype |     70 | 0.714 | [0.599, 0.807] |
| blood | adult | MESA | smr_moderate_finngen |    424 | 0.710 | [0.665, 0.751] |
| blood | adult | MESA | catalogue_single_blood_A |    880 | 0.675 | [0.643, 0.705] |
| blood | adult | MESA | catalogue_model_nasal_A |  1,390 | 0.668 | [0.643, 0.693] |
| blood | adult | MESA | distance_only |  5,196 | 0.668 | [0.655, 0.681] |
| blood | adult | MESA | distance_curve_nasal |  7,025 | 0.640 | [0.628, 0.651] |
| blood | adult | MESA | distance_curve_blood |  7,025 | 0.639 | [0.628, 0.650] |
| blood | adult | MESA | catalogue_model_solid_A |    769 | 0.627 | [0.592, 0.660] |
| blood | adult | MESA | catalogue_single_nasal_A |    728 | 0.615 | [0.580, 0.650] |
| blood | adult | MESA | distance_curve_solid |  7,025 | 0.610 | [0.598, 0.621] |
| blood | adult | MESA | smr_weak |     73 | 0.575 | [0.461, 0.682] |
| blood | adult | MESA | catalogue_single_solid_A |    329 | 0.547 | [0.493, 0.600] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak_discordant |     11 | 1.000 | [0.741, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_consensus_A |      4 | 1.000 | [0.510, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_B |      1 | 1.000 | [0.207, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_B |      1 | 1.000 | [0.207, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_replicated |    266 | 0.992 | [0.973, 0.998] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_replicated |     74 | 0.986 | [0.927, 0.998] |
| blood | child | HELIX (non-catalogue pairs) | smr_high |    392 | 0.967 | [0.944, 0.981] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_single |    625 | 0.965 | [0.947, 0.977] |
| blood | child | HELIX (non-catalogue pairs) | measured_nasal_single |    240 | 0.904 | [0.860, 0.935] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate |    995 | 0.895 | [0.875, 0.913] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S2 |    275 | 0.891 | [0.849, 0.923] |
| blood | child | HELIX (non-catalogue pairs) | smr_replicated_S1 |    122 | 0.885 | [0.817, 0.930] |
| blood | child | HELIX (non-catalogue pairs) | onco_consensus |    286 | 0.839 | [0.792, 0.877] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_single |    327 | 0.826 | [0.781, 0.863] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_finngen |    172 | 0.767 | [0.699, 0.824] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_discordant_celltype |      6 | 0.667 | [0.300, 0.903] |
| blood | child | HELIX (non-catalogue pairs) | smr_moderate_s3_replicated |     31 | 0.613 | [0.438, 0.763] |
| blood | child | HELIX (non-catalogue pairs) | smr_weak |    108 | 0.583 | [0.489, 0.672] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_A |    180 | 0.550 | [0.477, 0.621] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_A |    175 | 0.543 | [0.469, 0.615] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_solid_A |     81 | 0.444 | [0.341, 0.553] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_solid_A |     76 | 0.408 | [0.304, 0.520] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_nasal_A |      7 | 1.000 | [0.646, 1.000] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_replicated |     85 | 0.976 | [0.918, 0.994] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_replicated |     29 | 0.966 | [0.828, 0.994] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | onco_consensus |    278 | 0.960 | [0.931, 0.978] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_single |    179 | 0.911 | [0.860, 0.944] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_single |    236 | 0.886 | [0.839, 0.920] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_high |     73 | 0.836 | [0.734, 0.903] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_nasal_A |      9 | 0.778 | [0.453, 0.937] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak_discordant |      4 | 0.750 | [0.301, 0.954] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_weak |     26 | 0.692 | [0.500, 0.835] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S1 |      3 | 0.667 | [0.208, 0.939] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate |    279 | 0.609 | [0.551, 0.665] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_replicated_S2 |     28 | 0.607 | [0.424, 0.764] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_solid_A |     55 | 0.545 | [0.415, 0.670] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_solid_A |     57 | 0.544 | [0.416, 0.666] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_blood_A |     47 | 0.468 | [0.333, 0.608] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_blood_A |     51 | 0.431 | [0.305, 0.567] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_finngen |     13 | 0.385 | [0.177, 0.645] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_s3_replicated |      4 | 0.250 | [0.046, 0.699] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_consensus_A |      3 | 0.000 | [0.000, 0.561] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_moderate_discordant_celltype |      2 | 0.000 | [0.000, 0.658] |
| solid_tissue | any | FUSION | measured_solid_replicated |     16 | 1.000 | [0.806, 1.000] |
| solid_tissue | any | FUSION | measured_blood_replicated |    936 | 0.939 | [0.922, 0.953] |
| solid_tissue | any | FUSION | measured_solid_single |  1,634 | 0.921 | [0.907, 0.933] |
| solid_tissue | any | FUSION | catalogue_single_blood_B |     23 | 0.913 | [0.732, 0.976] |
| solid_tissue | any | FUSION | catalogue_model_blood_B |     35 | 0.886 | [0.740, 0.955] |
| solid_tissue | any | FUSION | measured_nasal_single |    707 | 0.885 | [0.860, 0.907] |
| solid_tissue | any | FUSION | onco_consensus |  1,629 | 0.870 | [0.853, 0.885] |
| solid_tissue | any | FUSION | smr_replicated_S1 |    424 | 0.851 | [0.814, 0.882] |
| solid_tissue | any | FUSION | measured_blood_single |  1,715 | 0.836 | [0.817, 0.852] |
| solid_tissue | any | FUSION | catalogue_model_nasal_B |    251 | 0.829 | [0.777, 0.870] |
| solid_tissue | any | FUSION | catalogue_single_nasal_B |    152 | 0.822 | [0.754, 0.875] |
| solid_tissue | any | FUSION | catalogue_consensus_B |    100 | 0.820 | [0.733, 0.883] |
| solid_tissue | any | FUSION | smr_high |    291 | 0.818 | [0.769, 0.858] |
| solid_tissue | any | FUSION | smr_weak_discordant |     45 | 0.800 | [0.662, 0.891] |
| solid_tissue | any | FUSION | catalogue_consensus_A |  1,411 | 0.751 | [0.728, 0.773] |
| solid_tissue | any | FUSION | smr_replicated_S2 |  1,158 | 0.734 | [0.708, 0.759] |
| solid_tissue | any | FUSION | smr_moderate_s3_replicated |     69 | 0.725 | [0.610, 0.816] |
| solid_tissue | any | FUSION | catalogue_single_nasal_A |  2,222 | 0.717 | [0.698, 0.735] |
| solid_tissue | any | FUSION | catalogue_model_nasal_A |  3,691 | 0.708 | [0.693, 0.722] |
| solid_tissue | any | FUSION | catalogue_model_blood_A |  4,733 | 0.703 | [0.690, 0.716] |
| solid_tissue | any | FUSION | catalogue_single_blood_A |  3,101 | 0.698 | [0.682, 0.714] |
| solid_tissue | any | FUSION | smr_moderate_finngen |    941 | 0.683 | [0.653, 0.712] |
| solid_tissue | any | FUSION | smr_moderate |  1,436 | 0.681 | [0.656, 0.705] |
| solid_tissue | any | FUSION | distance_only | 16,111 | 0.667 | [0.660, 0.675] |
| solid_tissue | any | FUSION | catalogue_model_solid_A |  2,029 | 0.644 | [0.623, 0.665] |
| solid_tissue | any | FUSION | distance_curve_blood | 29,694 | 0.635 | [0.630, 0.641] |
| solid_tissue | any | FUSION | catalogue_single_solid_A |  1,056 | 0.634 | [0.604, 0.662] |
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
| solid_tissue | any | GTEx solid | catalogue_consensus_B |     24 | 0.917 | [0.742, 0.977] |
| solid_tissue | any | GTEx solid | smr_weak_discordant |     24 | 0.917 | [0.742, 0.977] |
| solid_tissue | any | GTEx solid | onco_consensus |    551 | 0.900 | [0.872, 0.923] |
| solid_tissue | any | GTEx solid | smr_high |    159 | 0.849 | [0.785, 0.896] |
| solid_tissue | any | GTEx solid | smr_moderate_finngen |    302 | 0.815 | [0.767, 0.854] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_B |     52 | 0.808 | [0.681, 0.892] |
| solid_tissue | any | GTEx solid | smr_moderate |    460 | 0.761 | [0.720, 0.798] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_A |    569 | 0.756 | [0.719, 0.789] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_B |      4 | 0.750 | [0.301, 0.954] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_B |      7 | 0.714 | [0.359, 0.918] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_A |    793 | 0.711 | [0.679, 0.742] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_A |    973 | 0.703 | [0.674, 0.731] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_A |    747 | 0.696 | [0.662, 0.728] |
| solid_tissue | any | GTEx solid | catalogue_consensus_A |    213 | 0.690 | [0.625, 0.748] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_B |     29 | 0.690 | [0.508, 0.827] |
| solid_tissue | any | GTEx solid | smr_weak |     33 | 0.667 | [0.496, 0.802] |
| solid_tissue | any | GTEx solid | catalogue_model_solid_A |    338 | 0.618 | [0.565, 0.669] |
| solid_tissue | any | GTEx solid | catalogue_single_solid_A |    194 | 0.603 | [0.533, 0.669] |
| solid_tissue | any | GTEx solid | distance_only |  5,315 | 0.575 | [0.562, 0.588] |
| solid_tissue | any | GTEx solid | distance_curve_nasal |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_blood |  9,435 | 0.545 | [0.535, 0.555] |
| solid_tissue | any | GTEx solid | distance_curve_solid |  9,435 | 0.540 | [0.530, 0.550] |
| solid_tissue | any | GTEx solid | smr_moderate_discordant_celltype |     15 | 0.533 | [0.301, 0.752] |
