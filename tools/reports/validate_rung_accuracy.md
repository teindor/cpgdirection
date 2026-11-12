# Every ladder source, scored in each target tissue against an independent yardstick

Run 2026-11-12 16:25. Agreement with the yardstick's measured sign, Wilson 95% CI. Sources are scored standalone, on every yardstick pair where they call (catalogue sources: formed from predicted calls only, as the ladders form them). `n_cpg`: distinct CpGs.

Yardsticks: adult blood = MESA monocytes + T cells; child blood = HELIX Dryad pairs outside the blood catalogue's DIRECT set; nasal = EVA-PR nasal pairs outside the nasal catalogue's DIRECT set (the nasal catalogue is trained on EVA-PR: its rows are same-study); solid = FUSION, GTEx solid, the osteoarthritis joint tissues and kidney, each held out in turn, pooled, each pair once (with --blueprint, adult blood pools MESA and BLUEPRINT the same way). Measured cohorts are counted by donor group. SMR pairs are scored by accuracy band (rung not `reported`). Distance is not scored against HELIX or EVA-PR (its blood and nasal curves were fitted on them).

## blood (adult): MESA + BLUEPRINT

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | measured_blood_replicated |  4,542 |  3,881 | 0.987 | [0.983, 0.990] |
| 2 | smr_acc_95 | 11,549 |  9,206 | 0.961 | [0.958, 0.965] |
| 3 | measured_blood_single |  9,587 |  7,671 | 0.954 | [0.950, 0.958] |
| 4 | measured_solid_replicated |  1,641 |  1,464 | 0.937 | [0.924, 0.948] |
| 5 | smr_acc_90 |  4,424 |  4,099 | 0.904 | [0.895, 0.912] |
| 6 | measured_nasal_single |  1,800 |  1,558 | 0.887 | [0.872, 0.901] |
| 7 | smr_acc_80 |  3,719 |  3,412 | 0.846 | [0.834, 0.857] |
| 8 | onco_consensus |  3,336 |  3,130 | 0.829 | [0.816, 0.842] |
| 9 | measured_solid_single |  4,158 |  3,640 | 0.819 | [0.807, 0.830] |
| 10 | distance_only_b1 |  4,954 |  4,800 | 0.773 | [0.761, 0.784] |
| 11 | distance_curve_blood_b1 |  4,954 |  4,800 | 0.773 | [0.761, 0.784] |
| 12 | distance_curve_nasal_b1 |  4,954 |  4,800 | 0.773 | [0.761, 0.784] |
| 13 | distance_curve_solid_b1 |  4,954 |  4,800 | 0.773 | [0.761, 0.784] |
| 14 | catalogue_single_nasal_B |    253 |    253 | 0.771 | [0.715, 0.818] |
| 15 | catalogue_model_nasal_B |    401 |    400 | 0.756 | [0.711, 0.795] |
| 16 | catalogue_consensus_A |  2,145 |  2,002 | 0.733 | [0.714, 0.751] |
| 17 | smr_acc_65 |  2,368 |  2,236 | 0.730 | [0.711, 0.747] |
| 18 | catalogue_single_blood_B |     40 |     39 | 0.725 | [0.572, 0.839] |
| 19 | catalogue_consensus_B |    149 |    149 | 0.725 | [0.648, 0.790] |
| 20 | distance_only_b2 |  3,796 |  3,625 | 0.722 | [0.708, 0.736] |
| 21 | distance_curve_blood_b2 |  3,796 |  3,625 | 0.722 | [0.708, 0.736] |
| 22 | distance_curve_nasal_b2 |  3,796 |  3,625 | 0.722 | [0.708, 0.736] |
| 23 | distance_curve_solid_b2 |  3,796 |  3,625 | 0.722 | [0.708, 0.736] |
| 24 | catalogue_model_blood_A |  7,566 |  6,498 | 0.720 | [0.710, 0.730] |
| 25 | catalogue_model_blood_B |     62 |     61 | 0.710 | [0.587, 0.808] |
| 26 | catalogue_single_blood_A |  5,005 |  4,462 | 0.704 | [0.691, 0.716] |
| 27 | distance_only_b3 |  3,865 |  3,669 | 0.628 | [0.613, 0.644] |
| 28 | distance_curve_blood_b3 |  3,865 |  3,669 | 0.628 | [0.613, 0.644] |
| 29 | distance_curve_nasal_b3 |  3,865 |  3,669 | 0.628 | [0.613, 0.644] |
| 30 | distance_curve_solid_b3 |  3,865 |  3,669 | 0.628 | [0.613, 0.644] |
| 31 | smr_acc_55 |    607 |    596 | 0.618 | [0.579, 0.656] |
| 32 | catalogue_model_nasal_A |  6,138 |  5,202 | 0.611 | [0.599, 0.624] |
| 33 | catalogue_single_nasal_A |  3,698 |  3,328 | 0.611 | [0.596, 0.627] |
| 34 | catalogue_single_solid_A |  1,689 |  1,634 | 0.596 | [0.572, 0.619] |
| 35 | catalogue_model_solid_A |  3,393 |  3,151 | 0.586 | [0.570, 0.603] |
| 36 | distance_only_b4 |  3,088 |  2,896 | 0.575 | [0.558, 0.592] |
| 37 | distance_curve_blood_b4 |  3,088 |  2,896 | 0.575 | [0.558, 0.592] |
| 38 | distance_curve_nasal_b4 |  3,088 |  2,896 | 0.575 | [0.558, 0.592] |
| 39 | distance_curve_solid_b4 |  3,088 |  2,896 | 0.575 | [0.558, 0.592] |
| 40 | distance_only_b5 |  9,188 |  7,730 | 0.568 | [0.558, 0.578] |
| 41 | distance_curve_blood_b6 | 18,945 | 13,606 | 0.565 | [0.558, 0.572] |
| 42 | distance_curve_blood_b5 | 13,692 | 10,894 | 0.561 | [0.553, 0.569] |
| 43 | distance_curve_nasal_b5 | 13,692 | 10,894 | 0.561 | [0.553, 0.569] |
| 44 | distance_curve_blood_b7 |    280 |    279 | 0.532 | [0.474, 0.590] |
| 45 | distance_curve_solid_b5 | 13,692 | 10,894 | 0.530 | [0.521, 0.538] |
| 46 | distance_curve_nasal_b6 | 18,945 | 13,606 | 0.522 | [0.514, 0.529] |
| 47 | distance_curve_nasal_b7 |    280 |    279 | 0.468 | [0.410, 0.526] |
| 48 | distance_curve_solid_b7 |    280 |    279 | 0.468 | [0.410, 0.526] |
| 49 | distance_curve_solid_b6 | 18,945 | 13,606 | 0.435 | [0.428, 0.442] |

## blood (child): HELIX (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_consensus_A |   4 |   4 | 1.000 | [0.510, 1.000] |
| 2 | catalogue_single_nasal_B |   1 |   1 | 1.000 | [0.207, 1.000] |
| 3 | catalogue_model_nasal_B |   1 |   1 | 1.000 | [0.207, 1.000] |
| 4 | measured_blood_replicated | 323 | 288 | 0.988 | [0.969, 0.995] |
| 5 | measured_solid_replicated | 124 | 120 | 0.960 | [0.909, 0.983] |
| 6 | measured_blood_single | 631 | 547 | 0.954 | [0.935, 0.968] |
| 7 | smr_acc_90 | 567 | 533 | 0.928 | [0.903, 0.946] |
| 8 | smr_acc_95 | 917 | 825 | 0.926 | [0.907, 0.941] |
| 9 | measured_nasal_single | 240 | 212 | 0.904 | [0.860, 0.935] |
| 10 | smr_acc_80 | 438 | 405 | 0.893 | [0.860, 0.918] |
| 11 | onco_consensus | 286 | 255 | 0.839 | [0.792, 0.877] |
| 12 | measured_solid_single | 417 | 369 | 0.832 | [0.793, 0.865] |
| 13 | smr_acc_55 |  26 |  26 | 0.808 | [0.621, 0.915] |
| 14 | smr_acc_65 | 288 | 263 | 0.764 | [0.712, 0.809] |
| 15 | catalogue_model_nasal_A | 180 | 174 | 0.550 | [0.477, 0.621] |
| 16 | catalogue_single_nasal_A | 175 | 169 | 0.543 | [0.469, 0.615] |
| 17 | catalogue_model_solid_A |  81 |  79 | 0.444 | [0.341, 0.553] |
| 18 | catalogue_single_solid_A |  76 |  74 | 0.408 | [0.304, 0.520] |

## nasal_epithelium (any): EVA-PR nasal (non-catalogue pairs)

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | catalogue_single_nasal_A |   7 |   7 | 1.000 | [0.646, 1.000] |
| 2 | measured_solid_replicated |  59 |  52 | 0.983 | [0.910, 0.997] |
| 3 | measured_blood_replicated |  99 |  90 | 0.970 | [0.915, 0.990] |
| 4 | onco_consensus | 278 | 214 | 0.960 | [0.931, 0.978] |
| 5 | measured_solid_single | 205 | 159 | 0.898 | [0.848, 0.932] |
| 6 | measured_blood_single | 245 | 215 | 0.857 | [0.808, 0.895] |
| 7 | smr_acc_95 | 131 | 115 | 0.779 | [0.700, 0.841] |
| 8 | catalogue_model_nasal_A |   9 |   9 | 0.778 | [0.453, 0.937] |
| 9 | smr_acc_80 | 130 | 124 | 0.615 | [0.530, 0.695] |
| 10 | smr_acc_90 | 110 | 106 | 0.573 | [0.479, 0.661] |
| 11 | catalogue_single_solid_A |  55 |  53 | 0.545 | [0.415, 0.670] |
| 12 | catalogue_model_solid_A |  57 |  55 | 0.544 | [0.416, 0.666] |
| 13 | smr_acc_65 | 143 | 124 | 0.490 | [0.409, 0.571] |
| 14 | catalogue_single_blood_A |  47 |  46 | 0.468 | [0.333, 0.608] |
| 15 | catalogue_model_blood_A |  51 |  50 | 0.431 | [0.305, 0.567] |
| 16 | smr_acc_55 |   9 |   9 | 0.333 | [0.121, 0.646] |
| 17 | catalogue_consensus_A |   3 |   3 | 0.000 | [0.000, 0.561] |

## solid_tissue (any): FUSION + GTEx solid + osteoarthritis joint tissues + kidney

| rank | source | n | n_cpg | agree | 95% CI |
|---|---|---|---|---|---|
| 1 | measured_solid_replicated |  1,755 |  1,579 | 0.941 | [0.929, 0.951] |
| 2 | measured_blood_replicated |  2,208 |  1,964 | 0.928 | [0.916, 0.938] |
| 3 | onco_consensus |  4,762 |  4,586 | 0.895 | [0.886, 0.903] |
| 4 | measured_nasal_single |  1,648 |  1,436 | 0.888 | [0.872, 0.903] |
| 5 | measured_solid_single |  5,974 |  5,389 | 0.880 | [0.872, 0.888] |
| 6 | catalogue_consensus_B |    241 |    241 | 0.851 | [0.800, 0.890] |
| 7 | smr_acc_95 |  5,707 |  4,920 | 0.850 | [0.841, 0.859] |
| 8 | measured_blood_single |  4,287 |  3,660 | 0.838 | [0.827, 0.849] |
| 9 | catalogue_model_blood_B |    104 |    104 | 0.827 | [0.743, 0.888] |
| 10 | catalogue_model_nasal_B |    558 |    557 | 0.815 | [0.781, 0.845] |
| 11 | distance_only_b1 |  6,105 |  5,915 | 0.811 | [0.801, 0.821] |
| 12 | distance_curve_blood_b1 |  6,105 |  5,915 | 0.811 | [0.801, 0.821] |
| 13 | distance_curve_nasal_b1 |  6,105 |  5,915 | 0.811 | [0.801, 0.821] |
| 14 | distance_curve_solid_b1 |  6,105 |  5,915 | 0.811 | [0.801, 0.821] |
| 15 | catalogue_single_blood_B |     66 |     66 | 0.803 | [0.692, 0.881] |
| 16 | catalogue_single_nasal_B |    328 |    327 | 0.790 | [0.742, 0.830] |
| 17 | catalogue_consensus_A |  3,919 |  3,795 | 0.766 | [0.753, 0.779] |
| 18 | distance_only_b2 |  5,741 |  5,508 | 0.762 | [0.751, 0.773] |
| 19 | distance_curve_blood_b2 |  5,741 |  5,508 | 0.762 | [0.751, 0.773] |
| 20 | distance_curve_nasal_b2 |  5,741 |  5,508 | 0.762 | [0.751, 0.773] |
| 21 | distance_curve_solid_b2 |  5,741 |  5,508 | 0.762 | [0.751, 0.773] |
| 22 | smr_acc_90 |  2,743 |  2,564 | 0.741 | [0.724, 0.757] |
| 23 | catalogue_model_nasal_A | 10,497 |  9,575 | 0.699 | [0.690, 0.707] |
| 24 | catalogue_single_nasal_A |  6,432 |  5,989 | 0.693 | [0.681, 0.704] |
| 25 | catalogue_model_blood_A | 12,509 | 11,347 | 0.684 | [0.676, 0.692] |
| 26 | smr_acc_80 |  2,694 |  2,539 | 0.678 | [0.660, 0.695] |
| 27 | catalogue_single_blood_A |  8,138 |  7,545 | 0.674 | [0.663, 0.684] |
| 28 | distance_only_b3 |  7,218 |  6,804 | 0.673 | [0.662, 0.684] |
| 29 | distance_curve_blood_b3 |  7,218 |  6,804 | 0.673 | [0.662, 0.684] |
| 30 | distance_curve_nasal_b3 |  7,218 |  6,804 | 0.673 | [0.662, 0.684] |
| 31 | distance_curve_solid_b3 |  7,218 |  6,804 | 0.673 | [0.662, 0.684] |
| 32 | catalogue_model_solid_A |  5,517 |  5,323 | 0.667 | [0.654, 0.679] |
| 33 | catalogue_single_solid_A |  2,853 |  2,804 | 0.630 | [0.612, 0.647] |
| 34 | distance_only_b5 | 15,401 | 13,329 | 0.602 | [0.594, 0.609] |
| 35 | distance_curve_blood_b5 | 22,632 | 18,951 | 0.598 | [0.591, 0.604] |
| 36 | distance_curve_nasal_b5 | 22,632 | 18,951 | 0.598 | [0.591, 0.604] |
| 37 | distance_only_b4 |  6,116 |  5,617 | 0.595 | [0.582, 0.607] |
| 38 | distance_curve_blood_b4 |  6,116 |  5,617 | 0.595 | [0.582, 0.607] |
| 39 | distance_curve_nasal_b4 |  6,116 |  5,617 | 0.595 | [0.582, 0.607] |
| 40 | distance_curve_solid_b4 |  6,116 |  5,617 | 0.595 | [0.582, 0.607] |
| 41 | smr_acc_65 |  1,915 |  1,817 | 0.572 | [0.550, 0.594] |
| 42 | distance_curve_blood_b6 | 28,234 | 23,529 | 0.561 | [0.555, 0.566] |
| 43 | distance_curve_solid_b5 | 22,632 | 18,951 | 0.541 | [0.535, 0.547] |
| 44 | distance_curve_nasal_b6 | 28,234 | 23,529 | 0.526 | [0.520, 0.532] |
| 45 | smr_acc_55 |    649 |    643 | 0.448 | [0.411, 0.487] |
| 46 | distance_curve_solid_b6 | 28,234 | 23,529 | 0.439 | [0.434, 0.445] |

## By fold

| target | age | yardstick | source | n | agree | 95% CI |
|---|---|---|---|---|---|---|
| blood | adult | BLUEPRINT | measured_blood_replicated |  4,418 | 0.990 | [0.986, 0.992] |
| blood | adult | BLUEPRINT | smr_acc_95 | 10,740 | 0.969 | [0.966, 0.972] |
| blood | adult | BLUEPRINT | measured_blood_single |  8,679 | 0.960 | [0.955, 0.964] |
| blood | adult | BLUEPRINT | measured_solid_replicated |  1,560 | 0.946 | [0.934, 0.956] |
| blood | adult | BLUEPRINT | smr_acc_90 |  3,829 | 0.913 | [0.903, 0.921] |
| blood | adult | BLUEPRINT | measured_nasal_single |  1,631 | 0.900 | [0.885, 0.914] |
| blood | adult | BLUEPRINT | smr_acc_80 |  3,207 | 0.856 | [0.843, 0.868] |
| blood | adult | BLUEPRINT | measured_solid_single |  3,766 | 0.832 | [0.819, 0.843] |
| blood | adult | BLUEPRINT | onco_consensus |  2,815 | 0.829 | [0.815, 0.843] |
| blood | adult | BLUEPRINT | distance_only_b1 |  4,510 | 0.775 | [0.762, 0.787] |
| blood | adult | BLUEPRINT | distance_curve_blood_b1 |  4,510 | 0.775 | [0.762, 0.787] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b1 |  4,510 | 0.775 | [0.762, 0.787] |
| blood | adult | BLUEPRINT | distance_curve_solid_b1 |  4,510 | 0.775 | [0.762, 0.787] |
| blood | adult | BLUEPRINT | catalogue_single_nasal_B |    224 | 0.763 | [0.704, 0.814] |
| blood | adult | BLUEPRINT | catalogue_model_nasal_B |    348 | 0.747 | [0.699, 0.790] |
| blood | adult | BLUEPRINT | smr_acc_65 |  2,035 | 0.745 | [0.726, 0.764] |
| blood | adult | BLUEPRINT | catalogue_model_blood_A |  6,242 | 0.726 | [0.715, 0.737] |
| blood | adult | BLUEPRINT | distance_only_b2 |  3,224 | 0.725 | [0.710, 0.740] |
| blood | adult | BLUEPRINT | distance_curve_blood_b2 |  3,224 | 0.725 | [0.710, 0.740] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b2 |  3,224 | 0.725 | [0.710, 0.740] |
| blood | adult | BLUEPRINT | distance_curve_solid_b2 |  3,224 | 0.725 | [0.710, 0.740] |
| blood | adult | BLUEPRINT | catalogue_consensus_A |  1,601 | 0.716 | [0.694, 0.738] |
| blood | adult | BLUEPRINT | catalogue_single_blood_A |  4,286 | 0.714 | [0.701, 0.728] |
| blood | adult | BLUEPRINT | catalogue_single_blood_B |     38 | 0.711 | [0.552, 0.830] |
| blood | adult | BLUEPRINT | catalogue_consensus_B |    124 | 0.685 | [0.599, 0.761] |
| blood | adult | BLUEPRINT | catalogue_model_blood_B |     54 | 0.667 | [0.534, 0.778] |
| blood | adult | BLUEPRINT | smr_acc_55 |    482 | 0.637 | [0.593, 0.679] |
| blood | adult | BLUEPRINT | distance_only_b3 |  3,018 | 0.624 | [0.606, 0.641] |
| blood | adult | BLUEPRINT | distance_curve_blood_b3 |  3,018 | 0.624 | [0.606, 0.641] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b3 |  3,018 | 0.624 | [0.606, 0.641] |
| blood | adult | BLUEPRINT | distance_curve_solid_b3 |  3,018 | 0.624 | [0.606, 0.641] |
| blood | adult | BLUEPRINT | catalogue_single_nasal_A |  3,195 | 0.616 | [0.599, 0.633] |
| blood | adult | BLUEPRINT | catalogue_single_solid_A |  1,442 | 0.610 | [0.584, 0.634] |
| blood | adult | BLUEPRINT | catalogue_model_nasal_A |  5,083 | 0.603 | [0.589, 0.616] |
| blood | adult | BLUEPRINT | catalogue_model_solid_A |  2,784 | 0.578 | [0.560, 0.597] |
| blood | adult | BLUEPRINT | distance_curve_blood_b6 | 18,039 | 0.566 | [0.559, 0.573] |
| blood | adult | BLUEPRINT | distance_only_b4 |  2,383 | 0.562 | [0.542, 0.582] |
| blood | adult | BLUEPRINT | distance_curve_blood_b4 |  2,383 | 0.562 | [0.542, 0.582] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b4 |  2,383 | 0.562 | [0.542, 0.582] |
| blood | adult | BLUEPRINT | distance_curve_solid_b4 |  2,383 | 0.562 | [0.542, 0.582] |
| blood | adult | BLUEPRINT | distance_only_b5 |  7,576 | 0.556 | [0.545, 0.567] |
| blood | adult | BLUEPRINT | distance_curve_blood_b5 | 11,645 | 0.553 | [0.543, 0.562] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b5 | 11,645 | 0.553 | [0.543, 0.562] |
| blood | adult | BLUEPRINT | distance_curve_solid_b5 | 11,645 | 0.521 | [0.512, 0.530] |
| blood | adult | BLUEPRINT | distance_curve_nasal_b6 | 18,039 | 0.519 | [0.512, 0.526] |
| blood | adult | BLUEPRINT | distance_curve_solid_b6 | 18,039 | 0.434 | [0.427, 0.441] |
| blood | adult | MESA | catalogue_model_blood_B |     10 | 1.000 | [0.722, 1.000] |
| blood | adult | MESA | catalogue_single_blood_B |      4 | 1.000 | [0.510, 1.000] |
| blood | adult | MESA | measured_blood_replicated |    544 | 0.972 | [0.955, 0.983] |
| blood | adult | MESA | smr_acc_95 |  1,523 | 0.913 | [0.898, 0.926] |
| blood | adult | MESA | measured_solid_replicated |    239 | 0.904 | [0.860, 0.935] |
| blood | adult | MESA | measured_blood_single |    909 | 0.899 | [0.877, 0.917] |
| blood | adult | MESA | smr_acc_90 |    748 | 0.870 | [0.844, 0.893] |
| blood | adult | MESA | onco_consensus |    766 | 0.850 | [0.823, 0.873] |
| blood | adult | MESA | measured_nasal_single |    326 | 0.847 | [0.803, 0.882] |
| blood | adult | MESA | catalogue_single_nasal_B |     48 | 0.833 | [0.704, 0.913] |
| blood | adult | MESA | catalogue_consensus_B |     30 | 0.833 | [0.664, 0.927] |
| blood | adult | MESA | catalogue_model_nasal_B |     78 | 0.795 | [0.692, 0.870] |
| blood | adult | MESA | distance_only_b1 |    653 | 0.793 | [0.761, 0.823] |
| blood | adult | MESA | distance_curve_blood_b1 |    653 | 0.793 | [0.761, 0.823] |
| blood | adult | MESA | distance_curve_nasal_b1 |    653 | 0.793 | [0.761, 0.823] |
| blood | adult | MESA | distance_curve_solid_b1 |    653 | 0.793 | [0.761, 0.823] |
| blood | adult | MESA | catalogue_consensus_A |    641 | 0.793 | [0.759, 0.822] |
| blood | adult | MESA | smr_acc_80 |    601 | 0.779 | [0.744, 0.810] |
| blood | adult | MESA | measured_solid_single |    657 | 0.755 | [0.721, 0.786] |
| blood | adult | MESA | distance_only_b2 |    758 | 0.739 | [0.706, 0.769] |
| blood | adult | MESA | distance_curve_blood_b2 |    758 | 0.739 | [0.706, 0.769] |
| blood | adult | MESA | distance_curve_nasal_b2 |    758 | 0.739 | [0.706, 0.769] |
| blood | adult | MESA | distance_curve_solid_b2 |    758 | 0.739 | [0.706, 0.769] |
| blood | adult | MESA | catalogue_model_blood_A |  1,587 | 0.715 | [0.692, 0.737] |
| blood | adult | MESA | catalogue_single_blood_A |    880 | 0.675 | [0.643, 0.705] |
| blood | adult | MESA | catalogue_model_nasal_A |  1,390 | 0.668 | [0.643, 0.693] |
| blood | adult | MESA | distance_only_b3 |  1,010 | 0.651 | [0.622, 0.680] |
| blood | adult | MESA | distance_curve_blood_b3 |  1,010 | 0.651 | [0.622, 0.680] |
| blood | adult | MESA | distance_curve_nasal_b3 |  1,010 | 0.651 | [0.622, 0.680] |
| blood | adult | MESA | distance_curve_solid_b3 |  1,010 | 0.651 | [0.622, 0.680] |
| blood | adult | MESA | smr_acc_65 |    392 | 0.640 | [0.592, 0.686] |
| blood | adult | MESA | distance_only_b4 |    847 | 0.630 | [0.597, 0.662] |
| blood | adult | MESA | distance_curve_blood_b4 |    847 | 0.630 | [0.597, 0.662] |
| blood | adult | MESA | distance_curve_nasal_b4 |    847 | 0.630 | [0.597, 0.662] |
| blood | adult | MESA | distance_curve_solid_b4 |    847 | 0.630 | [0.597, 0.662] |
| blood | adult | MESA | catalogue_model_solid_A |    769 | 0.627 | [0.592, 0.660] |
| blood | adult | MESA | distance_only_b5 |  1,928 | 0.623 | [0.602, 0.645] |
| blood | adult | MESA | catalogue_single_nasal_A |    728 | 0.615 | [0.580, 0.650] |
| blood | adult | MESA | distance_curve_blood_b5 |  2,478 | 0.614 | [0.595, 0.633] |
| blood | adult | MESA | distance_curve_nasal_b5 |  2,478 | 0.614 | [0.595, 0.633] |
| blood | adult | MESA | distance_curve_solid_b5 |  2,478 | 0.578 | [0.558, 0.597] |
| blood | adult | MESA | distance_curve_nasal_b6 |    999 | 0.571 | [0.540, 0.601] |
| blood | adult | MESA | distance_curve_blood_b6 |    999 | 0.551 | [0.520, 0.581] |
| blood | adult | MESA | catalogue_single_solid_A |    329 | 0.547 | [0.493, 0.600] |
| blood | adult | MESA | distance_curve_blood_b7 |    280 | 0.532 | [0.474, 0.590] |
| blood | adult | MESA | smr_acc_55 |    140 | 0.529 | [0.446, 0.609] |
| blood | adult | MESA | distance_curve_nasal_b7 |    280 | 0.468 | [0.410, 0.526] |
| blood | adult | MESA | distance_curve_solid_b7 |    280 | 0.468 | [0.410, 0.526] |
| blood | adult | MESA | distance_curve_solid_b6 |    999 | 0.449 | [0.419, 0.480] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_consensus_A |      4 | 1.000 | [0.510, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_B |      1 | 1.000 | [0.207, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_B |      1 | 1.000 | [0.207, 1.000] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_replicated |    323 | 0.988 | [0.969, 0.995] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_replicated |    124 | 0.960 | [0.909, 0.983] |
| blood | child | HELIX (non-catalogue pairs) | measured_blood_single |    631 | 0.954 | [0.935, 0.968] |
| blood | child | HELIX (non-catalogue pairs) | smr_acc_90 |    567 | 0.928 | [0.903, 0.946] |
| blood | child | HELIX (non-catalogue pairs) | smr_acc_95 |    917 | 0.926 | [0.907, 0.941] |
| blood | child | HELIX (non-catalogue pairs) | measured_nasal_single |    240 | 0.904 | [0.860, 0.935] |
| blood | child | HELIX (non-catalogue pairs) | smr_acc_80 |    438 | 0.893 | [0.860, 0.918] |
| blood | child | HELIX (non-catalogue pairs) | onco_consensus |    286 | 0.839 | [0.792, 0.877] |
| blood | child | HELIX (non-catalogue pairs) | measured_solid_single |    417 | 0.832 | [0.793, 0.865] |
| blood | child | HELIX (non-catalogue pairs) | smr_acc_55 |     26 | 0.808 | [0.621, 0.915] |
| blood | child | HELIX (non-catalogue pairs) | smr_acc_65 |    288 | 0.764 | [0.712, 0.809] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_nasal_A |    180 | 0.550 | [0.477, 0.621] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_nasal_A |    175 | 0.543 | [0.469, 0.615] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_model_solid_A |     81 | 0.444 | [0.341, 0.553] |
| blood | child | HELIX (non-catalogue pairs) | catalogue_single_solid_A |     76 | 0.408 | [0.304, 0.520] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_nasal_A |      7 | 1.000 | [0.646, 1.000] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_replicated |     59 | 0.983 | [0.910, 0.997] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_replicated |     99 | 0.970 | [0.915, 0.990] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | onco_consensus |    278 | 0.960 | [0.931, 0.978] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_solid_single |    205 | 0.898 | [0.848, 0.932] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | measured_blood_single |    245 | 0.857 | [0.808, 0.895] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_acc_95 |    131 | 0.779 | [0.700, 0.841] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_nasal_A |      9 | 0.778 | [0.453, 0.937] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_acc_80 |    130 | 0.615 | [0.530, 0.695] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_acc_90 |    110 | 0.573 | [0.479, 0.661] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_solid_A |     55 | 0.545 | [0.415, 0.670] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_solid_A |     57 | 0.544 | [0.416, 0.666] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_acc_65 |    143 | 0.490 | [0.409, 0.571] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_single_blood_A |     47 | 0.468 | [0.333, 0.608] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_model_blood_A |     51 | 0.431 | [0.305, 0.567] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | smr_acc_55 |      9 | 0.333 | [0.121, 0.646] |
| nasal_epithelium | any | EVA-PR nasal (non-catalogue pairs) | catalogue_consensus_A |      3 | 0.000 | [0.000, 0.561] |
| solid_tissue | any | FUSION | measured_solid_replicated |  1,338 | 0.969 | [0.959, 0.977] |
| solid_tissue | any | FUSION | measured_blood_replicated |    986 | 0.936 | [0.919, 0.950] |
| solid_tissue | any | FUSION | catalogue_single_blood_B |     23 | 0.913 | [0.732, 0.976] |
| solid_tissue | any | FUSION | catalogue_model_blood_B |     35 | 0.886 | [0.740, 0.955] |
| solid_tissue | any | FUSION | measured_nasal_single |    707 | 0.885 | [0.860, 0.907] |
| solid_tissue | any | FUSION | measured_solid_single |  3,537 | 0.873 | [0.862, 0.884] |
| solid_tissue | any | FUSION | onco_consensus |  1,629 | 0.870 | [0.853, 0.885] |
| solid_tissue | any | FUSION | smr_acc_95 |  2,249 | 0.855 | [0.839, 0.869] |
| solid_tissue | any | FUSION | measured_blood_single |  1,697 | 0.835 | [0.817, 0.852] |
| solid_tissue | any | FUSION | catalogue_model_nasal_B |    251 | 0.829 | [0.777, 0.870] |
| solid_tissue | any | FUSION | catalogue_single_nasal_B |    152 | 0.822 | [0.754, 0.875] |
| solid_tissue | any | FUSION | catalogue_consensus_B |    100 | 0.820 | [0.733, 0.883] |
| solid_tissue | any | FUSION | distance_only_b1 |  2,043 | 0.792 | [0.774, 0.809] |
| solid_tissue | any | FUSION | distance_curve_blood_b1 |  2,043 | 0.792 | [0.774, 0.809] |
| solid_tissue | any | FUSION | distance_curve_nasal_b1 |  2,043 | 0.792 | [0.774, 0.809] |
| solid_tissue | any | FUSION | distance_curve_solid_b1 |  2,043 | 0.792 | [0.774, 0.809] |
| solid_tissue | any | FUSION | distance_only_b2 |  2,191 | 0.754 | [0.736, 0.772] |
| solid_tissue | any | FUSION | distance_curve_blood_b2 |  2,191 | 0.754 | [0.736, 0.772] |
| solid_tissue | any | FUSION | distance_curve_nasal_b2 |  2,191 | 0.754 | [0.736, 0.772] |
| solid_tissue | any | FUSION | distance_curve_solid_b2 |  2,191 | 0.754 | [0.736, 0.772] |
| solid_tissue | any | FUSION | catalogue_consensus_A |  1,411 | 0.751 | [0.728, 0.773] |
| solid_tissue | any | FUSION | smr_acc_90 |  1,079 | 0.748 | [0.721, 0.773] |
| solid_tissue | any | FUSION | catalogue_single_nasal_A |  2,222 | 0.717 | [0.698, 0.735] |
| solid_tissue | any | FUSION | catalogue_model_nasal_A |  3,691 | 0.708 | [0.693, 0.722] |
| solid_tissue | any | FUSION | catalogue_model_blood_A |  4,733 | 0.703 | [0.690, 0.716] |
| solid_tissue | any | FUSION | distance_only_b3 |  2,868 | 0.700 | [0.683, 0.717] |
| solid_tissue | any | FUSION | distance_curve_blood_b3 |  2,868 | 0.700 | [0.683, 0.717] |
| solid_tissue | any | FUSION | distance_curve_nasal_b3 |  2,868 | 0.700 | [0.683, 0.717] |
| solid_tissue | any | FUSION | distance_curve_solid_b3 |  2,868 | 0.700 | [0.683, 0.717] |
| solid_tissue | any | FUSION | catalogue_single_blood_A |  3,101 | 0.698 | [0.682, 0.714] |
| solid_tissue | any | FUSION | smr_acc_80 |  1,059 | 0.687 | [0.659, 0.715] |
| solid_tissue | any | FUSION | catalogue_model_solid_A |  2,029 | 0.644 | [0.623, 0.665] |
| solid_tissue | any | FUSION | catalogue_single_solid_A |  1,056 | 0.634 | [0.604, 0.662] |
| solid_tissue | any | FUSION | distance_only_b4 |  2,531 | 0.608 | [0.589, 0.627] |
| solid_tissue | any | FUSION | distance_curve_blood_b4 |  2,531 | 0.608 | [0.589, 0.627] |
| solid_tissue | any | FUSION | distance_curve_nasal_b4 |  2,531 | 0.608 | [0.589, 0.627] |
| solid_tissue | any | FUSION | distance_curve_solid_b4 |  2,531 | 0.608 | [0.589, 0.627] |
| solid_tissue | any | FUSION | distance_only_b5 |  6,478 | 0.607 | [0.595, 0.619] |
| solid_tissue | any | FUSION | distance_curve_blood_b6 | 10,549 | 0.602 | [0.592, 0.611] |
| solid_tissue | any | FUSION | distance_curve_blood_b5 |  9,512 | 0.599 | [0.589, 0.609] |
| solid_tissue | any | FUSION | distance_curve_nasal_b5 |  9,512 | 0.599 | [0.589, 0.609] |
| solid_tissue | any | FUSION | smr_acc_65 |    717 | 0.568 | [0.531, 0.603] |
| solid_tissue | any | FUSION | distance_curve_solid_b5 |  9,512 | 0.547 | [0.537, 0.557] |
| solid_tissue | any | FUSION | distance_curve_nasal_b6 | 10,549 | 0.539 | [0.530, 0.549] |
| solid_tissue | any | FUSION | smr_acc_55 |    239 | 0.452 | [0.390, 0.515] |
| solid_tissue | any | FUSION | distance_curve_solid_b6 | 10,549 | 0.398 | [0.389, 0.408] |
| solid_tissue | any | GTEx solid | measured_blood_replicated |    872 | 0.985 | [0.975, 0.991] |
| solid_tissue | any | GTEx solid | measured_solid_replicated |  1,214 | 0.984 | [0.976, 0.990] |
| solid_tissue | any | GTEx solid | measured_nasal_single |    481 | 0.958 | [0.937, 0.973] |
| solid_tissue | any | GTEx solid | smr_acc_95 |  1,144 | 0.957 | [0.944, 0.967] |
| solid_tissue | any | GTEx solid | measured_blood_single |  1,229 | 0.934 | [0.919, 0.947] |
| solid_tissue | any | GTEx solid | catalogue_consensus_B |     24 | 0.917 | [0.742, 0.977] |
| solid_tissue | any | GTEx solid | measured_solid_single |  1,718 | 0.907 | [0.892, 0.920] |
| solid_tissue | any | GTEx solid | onco_consensus |    551 | 0.900 | [0.872, 0.923] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_B |     52 | 0.808 | [0.681, 0.892] |
| solid_tissue | any | GTEx solid | distance_only_b1 |  1,097 | 0.787 | [0.761, 0.810] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b1 |  1,097 | 0.787 | [0.761, 0.810] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b1 |  1,097 | 0.787 | [0.761, 0.810] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b1 |  1,097 | 0.787 | [0.761, 0.810] |
| solid_tissue | any | GTEx solid | smr_acc_90 |    326 | 0.782 | [0.734, 0.824] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_A |    569 | 0.756 | [0.719, 0.789] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_B |      4 | 0.750 | [0.301, 0.954] |
| solid_tissue | any | GTEx solid | smr_acc_80 |    301 | 0.741 | [0.689, 0.787] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_B |      7 | 0.714 | [0.359, 0.918] |
| solid_tissue | any | GTEx solid | catalogue_model_nasal_A |    793 | 0.711 | [0.679, 0.742] |
| solid_tissue | any | GTEx solid | catalogue_model_blood_A |    973 | 0.703 | [0.674, 0.731] |
| solid_tissue | any | GTEx solid | catalogue_single_blood_A |    747 | 0.696 | [0.662, 0.728] |
| solid_tissue | any | GTEx solid | catalogue_consensus_A |    213 | 0.690 | [0.625, 0.748] |
| solid_tissue | any | GTEx solid | catalogue_single_nasal_B |     29 | 0.690 | [0.508, 0.827] |
| solid_tissue | any | GTEx solid | distance_only_b2 |    701 | 0.683 | [0.648, 0.717] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b2 |    701 | 0.683 | [0.648, 0.717] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b2 |    701 | 0.683 | [0.648, 0.717] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b2 |    701 | 0.683 | [0.648, 0.717] |
| solid_tissue | any | GTEx solid | smr_acc_65 |    217 | 0.631 | [0.565, 0.693] |
| solid_tissue | any | GTEx solid | smr_acc_55 |     45 | 0.622 | [0.476, 0.749] |
| solid_tissue | any | GTEx solid | catalogue_model_solid_A |    338 | 0.618 | [0.565, 0.669] |
| solid_tissue | any | GTEx solid | catalogue_single_solid_A |    194 | 0.603 | [0.533, 0.669] |
| solid_tissue | any | GTEx solid | distance_only_b3 |    718 | 0.503 | [0.466, 0.539] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b3 |    718 | 0.503 | [0.466, 0.539] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b3 |    718 | 0.503 | [0.466, 0.539] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b3 |    718 | 0.503 | [0.466, 0.539] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b5 |  3,023 | 0.501 | [0.483, 0.519] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b5 |  3,023 | 0.501 | [0.483, 0.519] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b6 |  3,077 | 0.500 | [0.483, 0.518] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b6 |  3,077 | 0.500 | [0.483, 0.518] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b6 |  3,077 | 0.500 | [0.482, 0.517] |
| solid_tissue | any | GTEx solid | distance_only_b5 |  1,980 | 0.489 | [0.467, 0.511] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b5 |  3,023 | 0.485 | [0.467, 0.503] |
| solid_tissue | any | GTEx solid | distance_only_b4 |    819 | 0.469 | [0.435, 0.503] |
| solid_tissue | any | GTEx solid | distance_curve_blood_b4 |    819 | 0.469 | [0.435, 0.503] |
| solid_tissue | any | GTEx solid | distance_curve_nasal_b4 |    819 | 0.469 | [0.435, 0.503] |
| solid_tissue | any | GTEx solid | distance_curve_solid_b4 |    819 | 0.469 | [0.435, 0.503] |
| solid_tissue | any | kidney | measured_solid_replicated |  1,368 | 0.978 | [0.969, 0.985] |
| solid_tissue | any | kidney | measured_nasal_single |    885 | 0.938 | [0.920, 0.952] |
| solid_tissue | any | kidney | measured_blood_replicated |  1,058 | 0.937 | [0.920, 0.950] |
| solid_tissue | any | kidney | onco_consensus |  2,479 | 0.929 | [0.918, 0.938] |
| solid_tissue | any | kidney | catalogue_consensus_B |    103 | 0.893 | [0.819, 0.939] |
| solid_tissue | any | kidney | measured_solid_single |  3,413 | 0.880 | [0.869, 0.891] |
| solid_tissue | any | kidney | distance_only_b1 |  2,939 | 0.860 | [0.847, 0.873] |
| solid_tissue | any | kidney | distance_curve_blood_b1 |  2,939 | 0.860 | [0.847, 0.873] |
| solid_tissue | any | kidney | distance_curve_nasal_b1 |  2,939 | 0.860 | [0.847, 0.873] |
| solid_tissue | any | kidney | distance_curve_solid_b1 |  2,939 | 0.860 | [0.847, 0.873] |
| solid_tissue | any | kidney | catalogue_model_nasal_B |    243 | 0.856 | [0.806, 0.895] |
| solid_tissue | any | kidney | smr_acc_95 |  2,408 | 0.855 | [0.840, 0.868] |
| solid_tissue | any | kidney | measured_blood_single |  1,593 | 0.853 | [0.835, 0.870] |
| solid_tissue | any | kidney | catalogue_single_nasal_B |    146 | 0.836 | [0.767, 0.887] |
| solid_tissue | any | kidney | distance_only_b2 |  2,483 | 0.822 | [0.806, 0.836] |
| solid_tissue | any | kidney | distance_curve_blood_b2 |  2,483 | 0.822 | [0.806, 0.836] |
| solid_tissue | any | kidney | distance_curve_nasal_b2 |  2,483 | 0.822 | [0.806, 0.836] |
| solid_tissue | any | kidney | distance_curve_solid_b2 |  2,483 | 0.822 | [0.806, 0.836] |
| solid_tissue | any | kidney | catalogue_consensus_A |  1,906 | 0.817 | [0.799, 0.834] |
| solid_tissue | any | kidney | catalogue_model_blood_B |     51 | 0.804 | [0.675, 0.890] |
| solid_tissue | any | kidney | distance_only_b3 |  2,825 | 0.770 | [0.754, 0.785] |
| solid_tissue | any | kidney | distance_curve_blood_b3 |  2,825 | 0.770 | [0.754, 0.785] |
| solid_tissue | any | kidney | distance_curve_nasal_b3 |  2,825 | 0.770 | [0.754, 0.785] |
| solid_tissue | any | kidney | distance_curve_solid_b3 |  2,825 | 0.770 | [0.754, 0.785] |
| solid_tissue | any | kidney | smr_acc_90 |  1,099 | 0.760 | [0.734, 0.784] |
| solid_tissue | any | kidney | catalogue_single_blood_B |     36 | 0.750 | [0.589, 0.862] |
| solid_tissue | any | kidney | catalogue_model_nasal_A |  5,185 | 0.731 | [0.719, 0.743] |
| solid_tissue | any | kidney | catalogue_model_blood_A |  5,272 | 0.719 | [0.707, 0.731] |
| solid_tissue | any | kidney | catalogue_model_solid_A |  2,618 | 0.718 | [0.701, 0.735] |
| solid_tissue | any | kidney | catalogue_single_nasal_A |  3,209 | 0.711 | [0.696, 0.727] |
| solid_tissue | any | kidney | catalogue_single_blood_A |  3,279 | 0.710 | [0.695, 0.726] |
| solid_tissue | any | kidney | distance_only_b4 |  2,179 | 0.695 | [0.676, 0.714] |
| solid_tissue | any | kidney | distance_curve_blood_b4 |  2,179 | 0.695 | [0.676, 0.714] |
| solid_tissue | any | kidney | distance_curve_nasal_b4 |  2,179 | 0.695 | [0.676, 0.714] |
| solid_tissue | any | kidney | distance_curve_solid_b4 |  2,179 | 0.695 | [0.676, 0.714] |
| solid_tissue | any | kidney | distance_only_b5 |  5,377 | 0.686 | [0.674, 0.698] |
| solid_tissue | any | kidney | distance_curve_blood_b5 |  7,829 | 0.671 | [0.660, 0.681] |
| solid_tissue | any | kidney | distance_curve_nasal_b5 |  7,829 | 0.671 | [0.660, 0.681] |
| solid_tissue | any | kidney | catalogue_single_solid_A |  1,365 | 0.667 | [0.642, 0.692] |
| solid_tissue | any | kidney | smr_acc_80 |  1,063 | 0.664 | [0.635, 0.692] |
| solid_tissue | any | kidney | distance_curve_solid_b5 |  7,829 | 0.585 | [0.574, 0.596] |
| solid_tissue | any | kidney | distance_curve_blood_b6 | 11,733 | 0.558 | [0.549, 0.567] |
| solid_tissue | any | kidney | smr_acc_65 |    794 | 0.530 | [0.495, 0.565] |
| solid_tissue | any | kidney | distance_curve_nasal_b6 | 11,733 | 0.528 | [0.519, 0.537] |
| solid_tissue | any | kidney | distance_curve_solid_b6 | 11,733 | 0.442 | [0.433, 0.451] |
| solid_tissue | any | kidney | smr_acc_55 |    288 | 0.434 | [0.378, 0.492] |
| solid_tissue | any | osteoarthritis joint tissues | measured_solid_replicated |  1,352 | 0.990 | [0.983, 0.994] |
| solid_tissue | any | osteoarthritis joint tissues | measured_blood_replicated |    845 | 0.940 | [0.922, 0.954] |
| solid_tissue | any | osteoarthritis joint tissues | measured_nasal_single |    590 | 0.898 | [0.871, 0.920] |
| solid_tissue | any | osteoarthritis joint tissues | onco_consensus |  1,304 | 0.891 | [0.873, 0.907] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_single_blood_B |      9 | 0.889 | [0.565, 0.980] |
| solid_tissue | any | osteoarthritis joint tissues | smr_acc_95 |  2,209 | 0.887 | [0.873, 0.899] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_model_blood_B |     17 | 0.882 | [0.657, 0.967] |
| solid_tissue | any | osteoarthritis joint tissues | measured_solid_single |  2,991 | 0.875 | [0.862, 0.886] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_consensus_B |     66 | 0.864 | [0.761, 0.927] |
| solid_tissue | any | osteoarthritis joint tissues | measured_blood_single |  1,331 | 0.846 | [0.826, 0.864] |
| solid_tissue | any | osteoarthritis joint tissues | distance_only_b1 |  1,500 | 0.795 | [0.773, 0.814] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b1 |  1,500 | 0.795 | [0.773, 0.814] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b1 |  1,500 | 0.795 | [0.773, 0.814] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b1 |  1,500 | 0.795 | [0.773, 0.814] |
| solid_tissue | any | osteoarthritis joint tissues | smr_acc_90 |    731 | 0.766 | [0.734, 0.795] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_model_nasal_B |    131 | 0.763 | [0.684, 0.828] |
| solid_tissue | any | osteoarthritis joint tissues | distance_only_b2 |  1,453 | 0.723 | [0.699, 0.745] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b2 |  1,453 | 0.723 | [0.699, 0.745] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b2 |  1,453 | 0.723 | [0.699, 0.745] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b2 |  1,453 | 0.723 | [0.699, 0.745] |
| solid_tissue | any | osteoarthritis joint tissues | smr_acc_80 |    673 | 0.713 | [0.678, 0.746] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_consensus_A |  1,035 | 0.712 | [0.684, 0.739] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_single_nasal_B |     66 | 0.682 | [0.562, 0.782] |
| solid_tissue | any | osteoarthritis joint tissues | smr_acc_65 |    423 | 0.667 | [0.620, 0.710] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_model_nasal_A |  2,776 | 0.663 | [0.645, 0.680] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_model_solid_A |  1,295 | 0.659 | [0.632, 0.684] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_single_nasal_A |  1,693 | 0.657 | [0.634, 0.680] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_model_blood_A |  3,317 | 0.622 | [0.606, 0.639] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_single_blood_A |  2,087 | 0.612 | [0.591, 0.633] |
| solid_tissue | any | osteoarthritis joint tissues | catalogue_single_solid_A |    633 | 0.608 | [0.570, 0.645] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b5 |  4,838 | 0.512 | [0.498, 0.526] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b5 |  4,838 | 0.512 | [0.498, 0.526] |
| solid_tissue | any | osteoarthritis joint tissues | distance_only_b3 |  1,748 | 0.509 | [0.486, 0.533] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b3 |  1,748 | 0.509 | [0.486, 0.533] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b3 |  1,748 | 0.509 | [0.486, 0.533] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b3 |  1,748 | 0.509 | [0.486, 0.533] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b6 |  4,622 | 0.506 | [0.492, 0.521] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b6 |  4,622 | 0.506 | [0.491, 0.520] |
| solid_tissue | any | osteoarthritis joint tissues | distance_only_b5 |  3,290 | 0.505 | [0.487, 0.522] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b6 |  4,622 | 0.494 | [0.480, 0.509] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b5 |  4,838 | 0.494 | [0.480, 0.508] |
| solid_tissue | any | osteoarthritis joint tissues | distance_only_b4 |  1,435 | 0.457 | [0.432, 0.483] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_blood_b4 |  1,435 | 0.457 | [0.432, 0.483] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_nasal_b4 |  1,435 | 0.457 | [0.432, 0.483] |
| solid_tissue | any | osteoarthritis joint tissues | distance_curve_solid_b4 |  1,435 | 0.457 | [0.432, 0.483] |
| solid_tissue | any | osteoarthritis joint tissues | smr_acc_55 |    165 | 0.424 | [0.351, 0.501] |
