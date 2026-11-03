# Cell-type SMR layer (GoDMC x FinnGen immune multiome L1) vs references

Layer: 463,048 pairs over 9 contexts, files `/Users/teindor/finngen_multiome/smr/smrd_blood_*.directions.csv.gz`.
Agreement = share of overlapping pairs where the layer sign equals the reference sign; Wilson 95% CI.
Baseline = majority-class share of the reference sign within the same overlap.
Ladder rungs for comparison: catalogue_single 0.62-0.87, catalogue_consensus 0.77-0.87,
smr_moderate 0.84-0.86 (GoDMC x eQTLGen S2), smr_high 0.95-0.97 (S1).
HEIDI is REPORTED, never applied (2.2.5 precedent); the by-stratum tables below are the record for this layer.

## Layer composition

| context | pairs | cpgs | genes | S1 | S2 | S3 | heidi_pass | heidi_fail | not_tested | pos_share |
|---|---|---|---|---|---|---|---|---|---|---|
| blood_PBMC | 123039 | 58444 | 8827 | 2373 | 120636 | 30 | 23825 | 93809 | 5405 | 0.49 |
| blood_CD4_T |  92183 | 49744 | 6989 | 1685 |  90484 | 14 | 18207 | 69896 | 4080 | 0.49 |
| blood_CD8_T |  61660 | 37659 | 5157 |  949 |  60696 | 15 | 13667 | 45389 | 2604 | 0.50 |
| blood_Mono |  69355 | 41487 | 5529 | 1179 |  68152 | 24 | 13889 | 52438 | 3028 | 0.49 |
| blood_NK |  44578 | 29699 | 3817 |  667 |  43900 | 11 | 10138 | 32576 | 1864 | 0.49 |
| blood_B |  34482 | 22807 | 2887 |  615 |  33851 | 16 |  7300 | 25844 | 1338 | 0.50 |
| blood_DC |  20073 | 14385 | 1732 |  275 |  19787 | 11 |  5043 | 14233 |  797 | 0.50 |
| blood_other_T |  15960 | 12701 | 1606 |  225 |  15729 |  6 |  4125 | 11169 |  666 | 0.49 |
| blood_other |   1718 |  1587 |  168 |   12 |   1705 |  1 |   683 |   942 |   93 | 0.48 |

## Reference: measured blood eQTMs

Reference pairs: 16,560; overlapping layer rows: 19,264 (5,474 distinct pairs).

### By context and tier

| context | smr_tier | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | S1 |  352 | 0.943 [0.914, 0.963] | 0.602 | 0.43 | 0.40 |
| blood_PBMC | S2 | 4001 | 0.900 [0.891, 0.909] | 0.588 | 0.43 | 0.41 |
| blood_CD4_T | S1 |  268 | 0.922 [0.883, 0.948] | 0.578 | 0.45 | 0.42 |
| blood_CD4_T | S2 | 3182 | 0.882 [0.871, 0.893] | 0.596 | 0.42 | 0.40 |
| blood_CD4_T | S3 |    1 | 0.000 [0.000, 0.793] | 1.000 | 1.00 | 0.00 |
| blood_CD8_T | S1 |  159 | 0.931 [0.880, 0.961] | 0.623 | 0.41 | 0.38 |
| blood_CD8_T | S2 | 2651 | 0.891 [0.879, 0.903] | 0.592 | 0.43 | 0.41 |
| blood_Mono | S1 |  169 | 0.923 [0.873, 0.954] | 0.580 | 0.43 | 0.42 |
| blood_Mono | S2 | 2462 | 0.883 [0.870, 0.896] | 0.593 | 0.43 | 0.41 |
| blood_Mono | S3 |    2 | 1.000 [0.342, 1.000] | 0.500 | 0.50 | 0.50 |
| blood_NK | S1 |   92 | 0.913 [0.838, 0.955] | 0.663 | 0.38 | 0.34 |
| blood_NK | S2 | 2109 | 0.880 [0.865, 0.893] | 0.583 | 0.44 | 0.42 |
| blood_B | S1 |   81 | 0.716 [0.610, 0.803] | 0.543 | 0.54 | 0.46 |
| blood_B | S2 | 1456 | 0.882 [0.864, 0.897] | 0.592 | 0.43 | 0.41 |
| blood_DC | S1 |   41 | 0.634 [0.481, 0.764] | 0.659 | 0.37 | 0.34 |
| blood_DC | S2 | 1153 | 0.865 [0.844, 0.883] | 0.595 | 0.43 | 0.41 |
| blood_DC | S3 |    2 | 0.500 [0.095, 0.905] | 0.500 | 1.00 | 0.50 |
| blood_other_T | S1 |   34 | 0.824 [0.665, 0.917] | 0.765 | 0.41 | 0.24 |
| blood_other_T | S2 |  909 | 0.882 [0.860, 0.902] | 0.617 | 0.41 | 0.38 |
| blood_other | S1 |    2 | 1.000 [0.342, 1.000] | 0.500 | 0.50 | 0.50 |
| blood_other | S2 |  137 | 0.912 [0.853, 0.949] | 0.533 | 0.51 | 0.47 |
| blood_other | S3 |    1 | 1.000 [0.207, 1.000] | 1.000 | 1.00 | 1.00 |

### By context and HEIDI status

| context | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | fail | 3755 | 0.899 [0.888, 0.908] | 0.591 | 0.43 | 0.41 |
| blood_PBMC | not_tested |   50 | 0.940 [0.838, 0.979] | 0.540 | 0.48 | 0.46 |
| blood_PBMC | pass |  548 | 0.936 [0.912, 0.954] | 0.584 | 0.40 | 0.42 |
| blood_CD4_T | fail | 2915 | 0.880 [0.867, 0.891] | 0.598 | 0.43 | 0.40 |
| blood_CD4_T | not_tested |   40 | 0.925 [0.801, 0.974] | 0.550 | 0.53 | 0.55 |
| blood_CD4_T | pass |  496 | 0.915 [0.888, 0.937] | 0.589 | 0.41 | 0.41 |
| blood_CD8_T | fail | 2277 | 0.886 [0.872, 0.898] | 0.600 | 0.43 | 0.40 |
| blood_CD8_T | not_tested |   37 | 0.892 [0.753, 0.957] | 0.514 | 0.51 | 0.51 |
| blood_CD8_T | pass |  496 | 0.929 [0.903, 0.949] | 0.575 | 0.42 | 0.43 |
| blood_Mono | fail | 2212 | 0.878 [0.864, 0.891] | 0.594 | 0.43 | 0.41 |
| blood_Mono | not_tested |   25 | 0.960 [0.805, 0.993] | 0.520 | 0.56 | 0.52 |
| blood_Mono | pass |  396 | 0.927 [0.897, 0.949] | 0.588 | 0.42 | 0.41 |
| blood_NK | fail | 1715 | 0.865 [0.848, 0.881] | 0.605 | 0.43 | 0.39 |
| blood_NK | not_tested |   26 | 1.000 [0.871, 1.000] | 0.615 | 0.62 | 0.62 |
| blood_NK | pass |  460 | 0.933 [0.906, 0.952] | 0.526 | 0.45 | 0.47 |
| blood_B | fail | 1218 | 0.856 [0.836, 0.875] | 0.590 | 0.43 | 0.41 |
| blood_B | not_tested |   17 | 0.882 [0.657, 0.967] | 0.706 | 0.59 | 0.71 |
| blood_B | pass |  302 | 0.940 [0.908, 0.962] | 0.603 | 0.42 | 0.40 |
| blood_DC | fail |  855 | 0.833 [0.806, 0.856] | 0.593 | 0.45 | 0.41 |
| blood_DC | not_tested |   11 | 0.909 [0.623, 0.984] | 0.636 | 0.73 | 0.64 |
| blood_DC | pass |  330 | 0.915 [0.880, 0.941] | 0.615 | 0.38 | 0.38 |
| blood_other_T | fail |  692 | 0.861 [0.834, 0.885] | 0.614 | 0.42 | 0.39 |
| blood_other_T | not_tested |   14 | 0.929 [0.685, 0.987] | 0.643 | 0.57 | 0.64 |
| blood_other_T | pass |  237 | 0.932 [0.893, 0.958] | 0.662 | 0.36 | 0.34 |
| blood_other | fail |   75 | 0.880 [0.787, 0.936] | 0.560 | 0.51 | 0.44 |
| blood_other | not_tested |    3 | 1.000 [0.438, 1.000] | 0.667 | 0.67 | 0.67 |
| blood_other | pass |   62 | 0.952 [0.867, 0.983] | 0.500 | 0.52 | 0.50 |

### All contexts pooled: tier x HEIDI status

| smr_tier | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | fail |  1084 | 0.894 [0.874, 0.911] | 0.606 | 0.43 | 0.39 |
| S1 | pass |   114 | 0.982 [0.938, 0.995] | 0.579 | 0.44 | 0.42 |
| S2 | fail | 14624 | 0.877 [0.872, 0.882] | 0.596 | 0.43 | 0.40 |
| S2 | not_tested |   223 | 0.933 [0.892, 0.959] | 0.552 | 0.55 | 0.55 |
| S2 | pass |  3213 | 0.927 [0.917, 0.935] | 0.585 | 0.41 | 0.42 |
| S3 | fail |     6 | 0.667 [0.300, 0.903] | 0.500 | 0.83 | 0.50 |

### All contexts pooled: by nsnp_HEIDI (tested pairs)

| nsnp_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|
| 3-5 |   355 | 0.915 [0.882, 0.940] | 0.589 | 0.42 | 0.41 |
| 6-10 |   526 | 0.914 [0.887, 0.935] | 0.587 | 0.40 | 0.41 |
| 11-20 | 18160 | 0.886 [0.881, 0.890] | 0.595 | 0.43 | 0.41 |

### All contexts pooled: tier x instruments

| smr_tier | inst_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | 2 |  1159 | 0.903 [0.884, 0.918] | 0.607 | 0.43 | 0.39 |
| S1 | 3-4 |    39 | 0.897 [0.764, 0.959] | 0.513 | 0.38 | 0.49 |
| S2 | 1 | 18060 | 0.887 [0.882, 0.891] | 0.592 | 0.43 | 0.41 |
| S3 | 3-4 |     6 | 0.667 [0.300, 0.903] | 0.500 | 0.83 | 0.50 |

### By context and |CpG - gene| distance

| context | tss_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | <=1.5kb | 1022 | 0.925 [0.907, 0.939] | 0.761 | 0.27 | 0.24 |
| blood_PBMC | 1.5-10kb |  635 | 0.917 [0.892, 0.936] | 0.638 | 0.39 | 0.36 |
| blood_PBMC | 10-100kb | 2057 | 0.897 [0.883, 0.909] | 0.514 | 0.49 | 0.49 |
| blood_PBMC | >100kb |  639 | 0.879 [0.852, 0.902] | 0.509 | 0.50 | 0.49 |
| blood_CD4_T | <=1.5kb |  801 | 0.903 [0.880, 0.921] | 0.758 | 0.26 | 0.24 |
| blood_CD4_T | 1.5-10kb |  512 | 0.871 [0.839, 0.897] | 0.648 | 0.38 | 0.35 |
| blood_CD4_T | 10-100kb | 1617 | 0.883 [0.867, 0.898] | 0.522 | 0.50 | 0.48 |
| blood_CD4_T | >100kb |  521 | 0.879 [0.848, 0.904] | 0.520 | 0.50 | 0.48 |
| blood_CD8_T | <=1.5kb |  664 | 0.920 [0.897, 0.938] | 0.748 | 0.27 | 0.25 |
| blood_CD8_T | 1.5-10kb |  390 | 0.890 [0.855, 0.917] | 0.667 | 0.36 | 0.33 |
| blood_CD8_T | 10-100kb | 1312 | 0.893 [0.875, 0.908] | 0.511 | 0.51 | 0.49 |
| blood_CD8_T | >100kb |  444 | 0.860 [0.825, 0.890] | 0.545 | 0.51 | 0.45 |
| blood_Mono | <=1.5kb |  618 | 0.895 [0.868, 0.917] | 0.715 | 0.33 | 0.28 |
| blood_Mono | 1.5-10kb |  367 | 0.899 [0.864, 0.926] | 0.632 | 0.40 | 0.37 |
| blood_Mono | 10-100kb | 1263 | 0.884 [0.866, 0.901] | 0.546 | 0.46 | 0.45 |
| blood_Mono | >100kb |  385 | 0.865 [0.827, 0.895] | 0.506 | 0.50 | 0.49 |
| blood_NK | <=1.5kb |  492 | 0.886 [0.855, 0.911] | 0.736 | 0.31 | 0.26 |
| blood_NK | 1.5-10kb |  311 | 0.913 [0.877, 0.940] | 0.704 | 0.31 | 0.30 |
| blood_NK | 10-100kb | 1079 | 0.880 [0.859, 0.898] | 0.501 | 0.50 | 0.50 |
| blood_NK | >100kb |  319 | 0.846 [0.803, 0.882] | 0.536 | 0.54 | 0.46 |
| blood_B | <=1.5kb |  327 | 0.884 [0.845, 0.914] | 0.749 | 0.31 | 0.25 |
| blood_B | 1.5-10kb |  209 | 0.876 [0.824, 0.914] | 0.636 | 0.40 | 0.36 |
| blood_B | 10-100kb |  732 | 0.872 [0.845, 0.894] | 0.544 | 0.46 | 0.46 |
| blood_B | >100kb |  269 | 0.862 [0.816, 0.899] | 0.517 | 0.53 | 0.52 |
| blood_DC | <=1.5kb |  277 | 0.881 [0.837, 0.914] | 0.744 | 0.28 | 0.26 |
| blood_DC | 1.5-10kb |  168 | 0.899 [0.844, 0.936] | 0.649 | 0.38 | 0.35 |
| blood_DC | 10-100kb |  566 | 0.825 [0.792, 0.854] | 0.532 | 0.49 | 0.47 |
| blood_DC | >100kb |  185 | 0.876 [0.820, 0.916] | 0.530 | 0.51 | 0.47 |
| blood_other_T | <=1.5kb |  198 | 0.914 [0.867, 0.946] | 0.783 | 0.24 | 0.22 |
| blood_other_T | 1.5-10kb |  125 | 0.896 [0.830, 0.938] | 0.752 | 0.30 | 0.25 |
| blood_other_T | 10-100kb |  460 | 0.872 [0.838, 0.899] | 0.533 | 0.49 | 0.47 |
| blood_other_T | >100kb |  160 | 0.850 [0.787, 0.897] | 0.581 | 0.46 | 0.42 |
| blood_other | <=1.5kb |   23 | 0.957 [0.790, 0.992] | 0.870 | 0.09 | 0.13 |
| blood_other | 1.5-10kb |   14 | 1.000 [0.785, 1.000] | 0.786 | 0.21 | 0.21 |
| blood_other | 10-100kb |   73 | 0.918 [0.832, 0.962] | 0.630 | 0.71 | 0.63 |
| blood_other | >100kb |   30 | 0.833 [0.664, 0.927] | 0.533 | 0.50 | 0.47 |

HEIDI pass vs fail concordance (pooled): pass 0.929 [0.920, 0.937], fail 0.878 [0.873, 0.883], Fisher p = 1.68e-18.

## Reference: measured nasal eQTMs

Reference pairs: 8,319; overlapping layer rows: 3,536 (1,209 distinct pairs).

### By context and tier

| context | smr_tier | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | S1 |  60 | 0.967 [0.886, 0.991] | 0.750 | 0.28 | 0.25 |
| blood_PBMC | S2 | 722 | 0.688 [0.654, 0.721] | 0.683 | 0.38 | 0.32 |
| blood_CD4_T | S1 |  48 | 0.958 [0.860, 0.988] | 0.708 | 0.29 | 0.29 |
| blood_CD4_T | S2 | 620 | 0.732 [0.696, 0.766] | 0.684 | 0.39 | 0.32 |
| blood_CD4_T | S3 |   1 | 0.000 [0.000, 0.793] | 1.000 | 1.00 | 0.00 |
| blood_CD8_T | S1 |  28 | 0.964 [0.823, 0.994] | 0.786 | 0.25 | 0.21 |
| blood_CD8_T | S2 | 437 | 0.746 [0.703, 0.785] | 0.664 | 0.43 | 0.34 |
| blood_Mono | S1 |  27 | 0.963 [0.817, 0.993] | 0.667 | 0.37 | 0.33 |
| blood_Mono | S2 | 511 | 0.663 [0.621, 0.703] | 0.642 | 0.45 | 0.36 |
| blood_Mono | S3 |   1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |
| blood_NK | S1 |  15 | 0.933 [0.702, 0.988] | 0.733 | 0.33 | 0.27 |
| blood_NK | S2 | 391 | 0.714 [0.667, 0.756] | 0.657 | 0.41 | 0.34 |
| blood_B | S1 |  16 | 0.875 [0.640, 0.965] | 0.500 | 0.62 | 0.50 |
| blood_B | S2 | 270 | 0.748 [0.693, 0.796] | 0.659 | 0.44 | 0.34 |
| blood_DC | S1 |   5 | 1.000 [0.566, 1.000] | 0.800 | 0.20 | 0.20 |
| blood_DC | S2 | 209 | 0.689 [0.623, 0.748] | 0.694 | 0.43 | 0.31 |
| blood_DC | S3 |   1 | 0.000 [0.000, 0.793] | 1.000 | 1.00 | 0.00 |
| blood_other_T | S1 |   6 | 0.833 [0.436, 0.970] | 0.833 | 0.33 | 0.17 |
| blood_other_T | S2 | 142 | 0.683 [0.603, 0.754] | 0.711 | 0.41 | 0.29 |
| blood_other | S2 |  26 | 0.808 [0.621, 0.915] | 0.731 | 0.46 | 0.27 |

### By context and HEIDI status

| context | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | fail | 627 | 0.729 [0.693, 0.762] | 0.689 | 0.38 | 0.31 |
| blood_PBMC | not_tested |  20 | 0.600 [0.387, 0.781] | 0.650 | 0.55 | 0.35 |
| blood_PBMC | pass | 135 | 0.637 [0.553, 0.713] | 0.689 | 0.32 | 0.31 |
| blood_CD4_T | fail | 517 | 0.750 [0.711, 0.786] | 0.694 | 0.38 | 0.31 |
| blood_CD4_T | not_tested |  14 | 0.714 [0.454, 0.883] | 0.500 | 0.50 | 0.50 |
| blood_CD4_T | pass | 138 | 0.739 [0.660, 0.805] | 0.674 | 0.40 | 0.33 |
| blood_CD8_T | fail | 361 | 0.756 [0.709, 0.798] | 0.668 | 0.43 | 0.33 |
| blood_CD8_T | not_tested |  10 | 0.700 [0.397, 0.892] | 0.600 | 0.70 | 0.40 |
| blood_CD8_T | pass |  94 | 0.777 [0.682, 0.849] | 0.691 | 0.36 | 0.31 |
| blood_Mono | fail | 426 | 0.685 [0.640, 0.728] | 0.657 | 0.44 | 0.34 |
| blood_Mono | not_tested |  16 | 0.688 [0.444, 0.858] | 0.500 | 0.44 | 0.50 |
| blood_Mono | pass |  97 | 0.649 [0.550, 0.737] | 0.608 | 0.45 | 0.39 |
| blood_NK | fail | 289 | 0.727 [0.673, 0.775] | 0.689 | 0.40 | 0.31 |
| blood_NK | not_tested |   9 | 0.778 [0.453, 0.937] | 0.667 | 0.56 | 0.33 |
| blood_NK | pass | 108 | 0.704 [0.612, 0.782] | 0.583 | 0.40 | 0.42 |
| blood_B | fail | 220 | 0.750 [0.689, 0.803] | 0.645 | 0.45 | 0.35 |
| blood_B | not_tested |   7 | 0.857 [0.487, 0.974] | 0.714 | 0.43 | 0.29 |
| blood_B | pass |  59 | 0.763 [0.640, 0.853] | 0.661 | 0.44 | 0.34 |
| blood_DC | fail | 153 | 0.673 [0.595, 0.742] | 0.654 | 0.46 | 0.35 |
| blood_DC | not_tested |   1 | 0.000 [0.000, 0.793] | 1.000 | 0.00 | 1.00 |
| blood_DC | pass |  61 | 0.754 [0.633, 0.845] | 0.820 | 0.33 | 0.18 |
| blood_other_T | fail | 108 | 0.704 [0.612, 0.782] | 0.731 | 0.40 | 0.27 |
| blood_other_T | not_tested |   4 | 0.750 [0.301, 0.954] | 0.750 | 0.50 | 0.25 |
| blood_other_T | pass |  36 | 0.639 [0.476, 0.775] | 0.667 | 0.42 | 0.33 |
| blood_other | fail |  14 | 0.786 [0.524, 0.924] | 0.571 | 0.64 | 0.43 |
| blood_other | pass |  12 | 0.833 [0.552, 0.953] | 0.917 | 0.25 | 0.08 |

### All contexts pooled: tier x HEIDI status

| smr_tier | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | fail |  174 | 0.943 [0.897, 0.968] | 0.736 | 0.31 | 0.26 |
| S1 | pass |   31 | 1.000 [0.890, 1.000] | 0.613 | 0.39 | 0.39 |
| S2 | fail | 2538 | 0.713 [0.695, 0.730] | 0.673 | 0.42 | 0.33 |
| S2 | not_tested |   81 | 0.691 [0.584, 0.781] | 0.593 | 0.52 | 0.41 |
| S2 | pass |  709 | 0.695 [0.660, 0.728] | 0.674 | 0.38 | 0.33 |
| S3 | fail |    3 | 0.333 [0.061, 0.792] | 1.000 | 0.67 | 0.00 |

### All contexts pooled: by nsnp_HEIDI (tested pairs)

| nsnp_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|
| 3-5 |  115 | 0.600 [0.509, 0.685] | 0.557 | 0.32 | 0.44 |
| 6-10 |  158 | 0.627 [0.549, 0.698] | 0.791 | 0.51 | 0.21 |
| 11-20 | 3182 | 0.733 [0.717, 0.748] | 0.675 | 0.40 | 0.32 |

### All contexts pooled: tier x instruments

| smr_tier | inst_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | 2 |  191 | 0.948 [0.906, 0.971] | 0.712 | 0.33 | 0.29 |
| S1 | 3-4 |   14 | 1.000 [0.785, 1.000] | 0.786 | 0.21 | 0.21 |
| S2 | 1 | 3328 | 0.709 [0.693, 0.724] | 0.672 | 0.41 | 0.33 |
| S3 | 3-4 |    3 | 0.333 [0.061, 0.792] | 1.000 | 0.67 | 0.00 |

### By context and |CpG - gene| distance

| context | tss_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | <=1.5kb | 199 | 0.784 [0.722, 0.835] | 0.844 | 0.24 | 0.16 |
| blood_PBMC | 1.5-10kb | 170 | 0.712 [0.640, 0.775] | 0.776 | 0.31 | 0.22 |
| blood_PBMC | 10-100kb | 327 | 0.706 [0.655, 0.753] | 0.596 | 0.45 | 0.40 |
| blood_PBMC | >100kb |  86 | 0.547 [0.442, 0.647] | 0.500 | 0.51 | 0.50 |
| blood_CD4_T | <=1.5kb | 159 | 0.818 [0.750, 0.870] | 0.855 | 0.24 | 0.14 |
| blood_CD4_T | 1.5-10kb | 127 | 0.787 [0.708, 0.850] | 0.787 | 0.28 | 0.21 |
| blood_CD4_T | 10-100kb | 289 | 0.713 [0.658, 0.762] | 0.609 | 0.50 | 0.39 |
| blood_CD4_T | >100kb |  94 | 0.681 [0.581, 0.766] | 0.500 | 0.41 | 0.50 |
| blood_CD8_T | <=1.5kb | 126 | 0.825 [0.750, 0.882] | 0.833 | 0.25 | 0.17 |
| blood_CD8_T | 1.5-10kb |  81 | 0.778 [0.676, 0.855] | 0.815 | 0.28 | 0.19 |
| blood_CD8_T | 10-100kb | 204 | 0.740 [0.676, 0.796] | 0.569 | 0.54 | 0.43 |
| blood_CD8_T | >100kb |  54 | 0.648 [0.515, 0.762] | 0.537 | 0.56 | 0.54 |
| blood_Mono | <=1.5kb | 132 | 0.727 [0.646, 0.796] | 0.773 | 0.36 | 0.23 |
| blood_Mono | 1.5-10kb |  90 | 0.622 [0.519, 0.715] | 0.767 | 0.43 | 0.23 |
| blood_Mono | 10-100kb | 238 | 0.672 [0.610, 0.729] | 0.555 | 0.45 | 0.45 |
| blood_Mono | >100kb |  79 | 0.684 [0.575, 0.776] | 0.557 | 0.56 | 0.44 |
| blood_NK | <=1.5kb |  92 | 0.772 [0.676, 0.846] | 0.859 | 0.30 | 0.14 |
| blood_NK | 1.5-10kb |  78 | 0.641 [0.530, 0.739] | 0.769 | 0.28 | 0.23 |
| blood_NK | 10-100kb | 181 | 0.746 [0.678, 0.804] | 0.569 | 0.46 | 0.43 |
| blood_NK | >100kb |  55 | 0.673 [0.541, 0.782] | 0.527 | 0.56 | 0.53 |
| blood_B | <=1.5kb |  70 | 0.771 [0.660, 0.854] | 0.886 | 0.29 | 0.11 |
| blood_B | 1.5-10kb |  55 | 0.745 [0.617, 0.842] | 0.764 | 0.42 | 0.24 |
| blood_B | 10-100kb | 106 | 0.802 [0.716, 0.867] | 0.566 | 0.44 | 0.43 |
| blood_B | >100kb |  55 | 0.655 [0.523, 0.766] | 0.600 | 0.69 | 0.60 |
| blood_DC | <=1.5kb |  59 | 0.746 [0.622, 0.839] | 0.847 | 0.41 | 0.15 |
| blood_DC | 1.5-10kb |  43 | 0.674 [0.525, 0.795] | 0.791 | 0.40 | 0.21 |
| blood_DC | 10-100kb |  89 | 0.719 [0.618, 0.802] | 0.607 | 0.40 | 0.39 |
| blood_DC | >100kb |  24 | 0.500 [0.314, 0.686] | 0.500 | 0.58 | 0.50 |
| blood_other_T | <=1.5kb |  40 | 0.750 [0.598, 0.858] | 0.825 | 0.28 | 0.17 |
| blood_other_T | 1.5-10kb |  17 | 0.765 [0.527, 0.904] | 0.824 | 0.29 | 0.18 |
| blood_other_T | 10-100kb |  67 | 0.701 [0.583, 0.798] | 0.657 | 0.49 | 0.34 |
| blood_other_T | >100kb |  24 | 0.500 [0.314, 0.686] | 0.625 | 0.46 | 0.38 |
| blood_other | <=1.5kb |   7 | 1.000 [0.646, 1.000] | 1.000 | 0.00 | 0.00 |
| blood_other | 1.5-10kb |   5 | 0.800 [0.376, 0.964] | 0.800 | 0.40 | 0.20 |
| blood_other | 10-100kb |  12 | 0.750 [0.468, 0.911] | 0.500 | 0.75 | 0.50 |
| blood_other | >100kb |   2 | 0.500 [0.095, 0.905] | 1.000 | 0.50 | 0.00 |

HEIDI pass vs fail concordance (pooled): pass 0.708 [0.674, 0.740], fail 0.727 [0.710, 0.744], Fisher p = 0.308.

## Reference: measured solid eQTMs

Reference pairs: 1,754; overlapping layer rows: 450 (168 distinct pairs).

### By context and tier

| context | smr_tier | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | S1 |   3 | 0.333 [0.061, 0.792] | 0.667 | 0.67 | 0.67 |
| blood_PBMC | S2 | 107 | 0.673 [0.579, 0.754] | 0.505 | 0.45 | 0.50 |
| blood_CD4_T | S1 |   1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |
| blood_CD4_T | S2 |  91 | 0.615 [0.513, 0.709] | 0.516 | 0.45 | 0.48 |
| blood_CD8_T | S2 |  64 | 0.594 [0.471, 0.705] | 0.547 | 0.42 | 0.55 |
| blood_Mono | S1 |   1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |
| blood_Mono | S2 |  75 | 0.613 [0.500, 0.715] | 0.587 | 0.37 | 0.41 |
| blood_NK | S2 |  43 | 0.581 [0.433, 0.716] | 0.605 | 0.40 | 0.40 |
| blood_B | S1 |   1 | 0.000 [0.000, 0.793] | 1.000 | 0.00 | 1.00 |
| blood_B | S2 |  30 | 0.600 [0.423, 0.754] | 0.533 | 0.40 | 0.47 |
| blood_DC | S2 |  17 | 0.529 [0.310, 0.738] | 0.588 | 0.41 | 0.41 |
| blood_other_T | S2 |  16 | 0.438 [0.231, 0.668] | 0.562 | 0.38 | 0.56 |
| blood_other | S2 |   1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |

### By context and HEIDI status

| context | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | fail | 80 | 0.662 [0.554, 0.757] | 0.512 | 0.45 | 0.51 |
| blood_PBMC | not_tested |  6 | 0.667 [0.300, 0.903] | 0.667 | 0.33 | 0.33 |
| blood_PBMC | pass | 24 | 0.667 [0.467, 0.820] | 0.500 | 0.50 | 0.50 |
| blood_CD4_T | fail | 63 | 0.635 [0.511, 0.743] | 0.524 | 0.41 | 0.52 |
| blood_CD4_T | not_tested |  5 | 0.600 [0.231, 0.882] | 0.600 | 0.40 | 0.40 |
| blood_CD4_T | pass | 24 | 0.583 [0.388, 0.755] | 0.625 | 0.54 | 0.38 |
| blood_CD8_T | fail | 36 | 0.500 [0.345, 0.655] | 0.583 | 0.36 | 0.58 |
| blood_CD8_T | not_tested |  2 | 0.500 [0.095, 0.905] | 1.000 | 0.50 | 0.00 |
| blood_CD8_T | pass | 26 | 0.731 [0.539, 0.863] | 0.538 | 0.50 | 0.54 |
| blood_Mono | fail | 57 | 0.579 [0.450, 0.698] | 0.509 | 0.42 | 0.49 |
| blood_Mono | not_tested |  6 | 0.667 [0.300, 0.903] | 0.667 | 0.33 | 0.33 |
| blood_Mono | pass | 13 | 0.769 [0.497, 0.918] | 0.923 | 0.15 | 0.08 |
| blood_NK | fail | 24 | 0.542 [0.351, 0.721] | 0.625 | 0.33 | 0.38 |
| blood_NK | not_tested |  3 | 0.667 [0.208, 0.939] | 0.667 | 0.67 | 0.33 |
| blood_NK | pass | 16 | 0.625 [0.386, 0.815] | 0.562 | 0.44 | 0.44 |
| blood_B | fail | 19 | 0.526 [0.317, 0.727] | 0.632 | 0.37 | 0.63 |
| blood_B | not_tested |  2 | 1.000 [0.342, 1.000] | 1.000 | 0.00 | 0.00 |
| blood_B | pass | 10 | 0.600 [0.313, 0.832] | 0.700 | 0.50 | 0.30 |
| blood_DC | fail |  9 | 0.556 [0.267, 0.811] | 0.667 | 0.33 | 0.33 |
| blood_DC | not_tested |  1 | 0.000 [0.000, 0.793] | 1.000 | 1.00 | 0.00 |
| blood_DC | pass |  7 | 0.571 [0.250, 0.842] | 0.571 | 0.43 | 0.57 |
| blood_other_T | fail |  6 | 0.333 [0.097, 0.700] | 0.833 | 0.17 | 0.83 |
| blood_other_T | not_tested |  2 | 0.500 [0.095, 0.905] | 1.000 | 0.50 | 0.00 |
| blood_other_T | pass |  8 | 0.500 [0.215, 0.785] | 0.500 | 0.50 | 0.50 |
| blood_other | pass |  1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |

### All contexts pooled: tier x HEIDI status

| smr_tier | heidi_status | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | fail |   6 | 0.500 [0.188, 0.812] | 0.500 | 0.33 | 0.50 |
| S2 | fail | 288 | 0.594 [0.536, 0.649] | 0.517 | 0.40 | 0.52 |
| S2 | not_tested |  27 | 0.630 [0.442, 0.785] | 0.741 | 0.41 | 0.26 |
| S2 | pass | 129 | 0.651 [0.566, 0.728] | 0.581 | 0.46 | 0.42 |

### All contexts pooled: by nsnp_HEIDI (tested pairs)

| nsnp_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|
| 3-5 |  33 | 0.394 [0.247, 0.563] | 0.515 | 0.70 | 0.52 |
| 6-10 |  52 | 0.731 [0.597, 0.832] | 0.769 | 0.15 | 0.23 |
| 11-20 | 338 | 0.612 [0.559, 0.663] | 0.524 | 0.43 | 0.52 |

### All contexts pooled: tier x instruments

| smr_tier | inst_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | 2 |   6 | 0.500 [0.188, 0.812] | 0.500 | 0.33 | 0.50 |
| S2 | 1 | 444 | 0.613 [0.567, 0.657] | 0.527 | 0.42 | 0.47 |

### By context and |CpG - gene| distance

| context | tss_bin | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | <=1.5kb | 38 | 0.605 [0.447, 0.744] | 0.711 | 0.32 | 0.29 |
| blood_PBMC | 1.5-10kb | 18 | 0.778 [0.548, 0.910] | 0.556 | 0.56 | 0.44 |
| blood_PBMC | 10-100kb | 49 | 0.633 [0.493, 0.753] | 0.633 | 0.47 | 0.63 |
| blood_PBMC | >100kb |  5 | 1.000 [0.566, 1.000] | 1.000 | 1.00 | 1.00 |
| blood_CD4_T | <=1.5kb | 33 | 0.636 [0.466, 0.778] | 0.727 | 0.33 | 0.27 |
| blood_CD4_T | 1.5-10kb | 12 | 0.750 [0.468, 0.911] | 0.583 | 0.50 | 0.42 |
| blood_CD4_T | 10-100kb | 41 | 0.561 [0.410, 0.701] | 0.634 | 0.44 | 0.63 |
| blood_CD4_T | >100kb |  6 | 0.667 [0.300, 0.903] | 0.667 | 1.00 | 0.67 |
| blood_CD8_T | <=1.5kb | 24 | 0.625 [0.427, 0.788] | 0.625 | 0.33 | 0.38 |
| blood_CD8_T | 1.5-10kb | 11 | 0.727 [0.434, 0.903] | 0.545 | 0.55 | 0.45 |
| blood_CD8_T | 10-100kb | 27 | 0.481 [0.307, 0.660] | 0.704 | 0.41 | 0.70 |
| blood_CD8_T | >100kb |  2 | 1.000 [0.342, 1.000] | 1.000 | 1.00 | 1.00 |
| blood_Mono | <=1.5kb | 24 | 0.583 [0.388, 0.755] | 0.792 | 0.21 | 0.21 |
| blood_Mono | 1.5-10kb | 12 | 0.750 [0.468, 0.911] | 0.750 | 0.50 | 0.25 |
| blood_Mono | 10-100kb | 35 | 0.571 [0.409, 0.720] | 0.543 | 0.40 | 0.54 |
| blood_Mono | >100kb |  5 | 0.800 [0.376, 0.964] | 0.800 | 0.60 | 0.80 |
| blood_NK | <=1.5kb | 14 | 0.571 [0.326, 0.786] | 0.786 | 0.21 | 0.21 |
| blood_NK | 1.5-10kb |  7 | 0.571 [0.250, 0.842] | 0.714 | 0.43 | 0.29 |
| blood_NK | 10-100kb | 20 | 0.550 [0.342, 0.742] | 0.550 | 0.50 | 0.55 |
| blood_NK | >100kb |  2 | 1.000 [0.342, 1.000] | 0.500 | 0.50 | 0.50 |
| blood_B | <=1.5kb |  8 | 0.500 [0.215, 0.785] | 0.750 | 0.25 | 0.25 |
| blood_B | 1.5-10kb |  4 | 0.750 [0.301, 0.954] | 1.000 | 0.25 | 0.00 |
| blood_B | 10-100kb | 16 | 0.625 [0.386, 0.815] | 0.688 | 0.44 | 0.69 |
| blood_B | >100kb |  3 | 0.333 [0.061, 0.792] | 0.667 | 0.67 | 0.67 |
| blood_DC | <=1.5kb |  4 | 0.750 [0.301, 0.954] | 0.750 | 0.00 | 0.25 |
| blood_DC | 1.5-10kb |  2 | 0.500 [0.095, 0.905] | 1.000 | 0.50 | 0.00 |
| blood_DC | 10-100kb |  9 | 0.444 [0.189, 0.733] | 0.556 | 0.44 | 0.56 |
| blood_DC | >100kb |  2 | 0.500 [0.095, 0.905] | 0.500 | 1.00 | 0.50 |
| blood_other_T | <=1.5kb |  4 | 0.500 [0.150, 0.850] | 0.500 | 0.00 | 0.50 |
| blood_other_T | 1.5-10kb |  4 | 0.500 [0.150, 0.850] | 0.750 | 0.75 | 0.25 |
| blood_other_T | 10-100kb |  8 | 0.375 [0.137, 0.694] | 0.750 | 0.38 | 0.75 |
| blood_other | <=1.5kb |  1 | 1.000 [0.207, 1.000] | 1.000 | 0.00 | 0.00 |

HEIDI pass vs fail concordance (pooled): pass 0.651 [0.566, 0.728], fail 0.592 [0.535, 0.646], Fisher p = 0.279.

## Consistency: existing SMR layer (GoDMC x eQTLGen) -- same exposure, independent outcome

Overlapping rows: 271,440. This is a pipeline check, not independent validation: the instruments are shared.

### By context and layer tier

| context | smr_tier | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| blood_PBMC | S1 |  1704 | 0.913 [0.899, 0.926] | 0.535 | 0.48 | 0.46 |
| blood_PBMC | S2 | 66891 | 0.877 [0.875, 0.880] | 0.516 | 0.48 | 0.48 |
| blood_PBMC | S3 |    18 | 0.611 [0.386, 0.797] | 0.722 | 0.44 | 0.28 |
| blood_CD4_T | S1 |  1190 | 0.875 [0.855, 0.892] | 0.526 | 0.47 | 0.47 |
| blood_CD4_T | S2 | 50400 | 0.842 [0.839, 0.845] | 0.515 | 0.49 | 0.49 |
| blood_CD4_T | S3 |    10 | 0.900 [0.596, 0.982] | 0.800 | 0.30 | 0.20 |
| blood_CD8_T | S1 |   717 | 0.895 [0.871, 0.916] | 0.537 | 0.48 | 0.46 |
| blood_CD8_T | S2 | 36887 | 0.856 [0.852, 0.859] | 0.512 | 0.49 | 0.49 |
| blood_CD8_T | S3 |    11 | 0.818 [0.523, 0.949] | 0.818 | 0.36 | 0.18 |
| blood_Mono | S1 |   886 | 0.886 [0.863, 0.905] | 0.524 | 0.49 | 0.48 |
| blood_Mono | S2 | 39813 | 0.860 [0.857, 0.864] | 0.516 | 0.48 | 0.48 |
| blood_Mono | S3 |    20 | 0.700 [0.481, 0.855] | 0.550 | 0.45 | 0.55 |
| blood_NK | S1 |   474 | 0.861 [0.827, 0.889] | 0.555 | 0.47 | 0.45 |
| blood_NK | S2 | 26800 | 0.840 [0.836, 0.845] | 0.516 | 0.49 | 0.48 |
| blood_NK | S3 |     8 | 0.750 [0.409, 0.929] | 0.875 | 0.38 | 0.12 |
| blood_B | S1 |   431 | 0.784 [0.743, 0.820] | 0.531 | 0.55 | 0.47 |
| blood_B | S2 | 19950 | 0.824 [0.819, 0.829] | 0.512 | 0.49 | 0.49 |
| blood_B | S3 |    12 | 0.417 [0.193, 0.680] | 0.583 | 0.67 | 0.58 |
| blood_DC | S1 |   208 | 0.812 [0.754, 0.860] | 0.625 | 0.50 | 0.38 |
| blood_DC | S2 | 13316 | 0.824 [0.817, 0.830] | 0.518 | 0.49 | 0.48 |
| blood_DC | S3 |    11 | 0.727 [0.434, 0.903] | 0.545 | 0.64 | 0.55 |
| blood_other_T | S1 |   156 | 0.782 [0.711, 0.840] | 0.603 | 0.47 | 0.40 |
| blood_other_T | S2 | 10152 | 0.822 [0.814, 0.829] | 0.525 | 0.48 | 0.47 |
| blood_other_T | S3 |     6 | 0.833 [0.436, 0.970] | 1.000 | 0.17 | 0.00 |
| blood_other | S1 |    11 | 1.000 [0.741, 1.000] | 0.545 | 0.45 | 0.45 |
| blood_other | S2 |  1357 | 0.901 [0.883, 0.915] | 0.541 | 0.47 | 0.46 |
| blood_other | S3 |     1 | 1.000 [0.207, 1.000] | 1.000 | 1.00 | 1.00 |

### Layer tier x eQTLGen tier (pooled)

| smr_tier | ref_tier | n | agree | baseline | layer_pos | ref_pos |
|---|---|---|---|---|---|---|
| S1 | S1 |   3200 | 0.935 [0.926, 0.943] | 0.556 | 0.47 | 0.44 |
| S1 | S2 |   2058 | 0.856 [0.840, 0.870] | 0.521 | 0.48 | 0.48 |
| S1 | S3 |    519 | 0.613 [0.570, 0.654] | 0.503 | 0.57 | 0.50 |
| S2 | S1 |  24942 | 0.866 [0.862, 0.870] | 0.534 | 0.48 | 0.47 |
| S2 | S2 | 226326 | 0.860 [0.858, 0.861] | 0.515 | 0.49 | 0.48 |
| S2 | S3 |  14298 | 0.719 [0.712, 0.726] | 0.513 | 0.51 | 0.51 |
| S3 | S1 |     39 | 0.821 [0.673, 0.910] | 0.872 | 0.31 | 0.13 |
| S3 | S2 |     17 | 0.706 [0.469, 0.867] | 0.882 | 0.29 | 0.12 |
| S3 | S3 |     41 | 0.585 [0.434, 0.722] | 0.683 | 0.66 | 0.68 |

### Coverage: pairs already in the eQTLGen layer vs new

| context | pairs | in_eqtlgen_layer | new_pairs |
|---|---|---|---|
| blood_PBMC | 123039 | 68613 | 54426 |
| blood_CD4_T |  92183 | 51600 | 40583 |
| blood_CD8_T |  61660 | 37615 | 24045 |
| blood_Mono |  69355 | 40719 | 28636 |
| blood_NK |  44578 | 27282 | 17296 |
| blood_B |  34482 | 20393 | 14089 |
| blood_DC |  20073 | 13535 |  6538 |
| blood_other_T |  15960 | 10314 |  5646 |
| blood_other |   1718 |  1369 |   349 |

## Cell types vs pooled PBMC on shared pairs

Disagreement here is either cell-type-specific regulation or noise; the measured-reference tables say which.

| context | shared | agree | only_in_celltype |
|---|---|---|---|
| blood_CD4_T | 74529 | 0.997 [0.997, 0.998] | 17654 |
| blood_CD8_T | 54917 | 0.998 [0.998, 0.999] |  6743 |
| blood_Mono | 46535 | 0.989 [0.988, 0.990] | 22820 |
| blood_NK | 36793 | 0.992 [0.991, 0.993] |  7785 |
| blood_B | 26005 | 0.991 [0.989, 0.992] |  8477 |
| blood_DC | 14482 | 0.991 [0.990, 0.993] |  5591 |
| blood_other_T | 14411 | 0.998 [0.997, 0.998] |  1549 |
| blood_other |  1411 | 0.984 [0.977, 0.990] |   307 |

### Disagreements with PBMC that have a measured blood sign

| context | n | celltype_right | pbmc_right |
|---|---|---|---|
| blood_CD4_T |  7 | 0.571 [0.250, 0.842] | 0.429 [0.158, 0.750] |
| blood_CD8_T |  4 | 1.000 [0.510, 1.000] | 0.000 [0.000, 0.490] |
| blood_Mono | 37 | 0.189 [0.095, 0.342] | 0.811 [0.658, 0.905] |
| blood_NK | 18 | 0.056 [0.010, 0.258] | 0.944 [0.742, 0.990] |
| blood_B | 12 | 0.083 [0.015, 0.354] | 0.917 [0.646, 0.985] |
| blood_DC |  6 | 0.167 [0.030, 0.564] | 0.833 [0.436, 0.970] |
| blood_other_T |  3 | 0.667 [0.208, 0.939] | 0.333 [0.061, 0.792] |
| blood_other |  1 | 0.000 [0.000, 0.793] | 1.000 [0.207, 1.000] |

## Summary

| reference | context | tier | n | agreement | baseline | verdict |
|---|---|---|---|---|---|---|
| blood | blood_B | S1 |   81 | 0.716 | 0.543 | too few pairs |
| blood | blood_B | S2 | 1456 | 0.882 | 0.592 | at smr_moderate level |
| blood | blood_CD4_T | S1 |  268 | 0.922 | 0.578 | at smr_moderate level |
| blood | blood_CD4_T | S2 | 3182 | 0.882 | 0.596 | at smr_moderate level |
| blood | blood_CD4_T | S3 |    1 | 0.000 | 1.000 | too few pairs |
| blood | blood_CD8_T | S1 |  159 | 0.931 | 0.623 | too few pairs |
| blood | blood_CD8_T | S2 | 2651 | 0.891 | 0.592 | at smr_moderate level |
| blood | blood_DC | S1 |   41 | 0.634 | 0.659 | too few pairs |
| blood | blood_DC | S2 | 1153 | 0.865 | 0.595 | at smr_moderate level |
| blood | blood_DC | S3 |    2 | 0.500 | 0.500 | too few pairs |
| blood | blood_Mono | S1 |  169 | 0.923 | 0.580 | too few pairs |
| blood | blood_Mono | S2 | 2462 | 0.883 | 0.593 | at smr_moderate level |
| blood | blood_Mono | S3 |    2 | 1.000 | 0.500 | too few pairs |
| blood | blood_NK | S1 |   92 | 0.913 | 0.663 | too few pairs |
| blood | blood_NK | S2 | 2109 | 0.880 | 0.583 | at smr_moderate level |
| blood | blood_PBMC | S1 |  352 | 0.943 | 0.602 | at smr_moderate level |
| blood | blood_PBMC | S2 | 4001 | 0.900 | 0.588 | at smr_moderate level |
| blood | blood_other | S1 |    2 | 1.000 | 0.500 | too few pairs |
| blood | blood_other | S2 |  137 | 0.912 | 0.533 | too few pairs |
| blood | blood_other | S3 |    1 | 1.000 | 1.000 | too few pairs |
| blood | blood_other_T | S1 |   34 | 0.824 | 0.765 | too few pairs |
| blood | blood_other_T | S2 |  909 | 0.882 | 0.617 | at smr_moderate level |
| nasal | blood_B | S1 |   16 | 0.875 | 0.500 | too few pairs |
| nasal | blood_B | S2 |  270 | 0.748 | 0.659 | at catalogue_single level |
| nasal | blood_CD4_T | S1 |   48 | 0.958 | 0.708 | too few pairs |
| nasal | blood_CD4_T | S2 |  620 | 0.732 | 0.684 | below ladder / not above baseline |
| nasal | blood_CD4_T | S3 |    1 | 0.000 | 1.000 | too few pairs |
| nasal | blood_CD8_T | S1 |   28 | 0.964 | 0.786 | too few pairs |
| nasal | blood_CD8_T | S2 |  437 | 0.746 | 0.664 | at catalogue_single level |
| nasal | blood_DC | S1 |    5 | 1.000 | 0.800 | too few pairs |
| nasal | blood_DC | S2 |  209 | 0.689 | 0.694 | below ladder / not above baseline |
| nasal | blood_DC | S3 |    1 | 0.000 | 1.000 | too few pairs |
| nasal | blood_Mono | S1 |   27 | 0.963 | 0.667 | too few pairs |
| nasal | blood_Mono | S2 |  511 | 0.663 | 0.642 | below ladder / not above baseline |
| nasal | blood_Mono | S3 |    1 | 1.000 | 1.000 | too few pairs |
| nasal | blood_NK | S1 |   15 | 0.933 | 0.733 | too few pairs |
| nasal | blood_NK | S2 |  391 | 0.714 | 0.657 | at catalogue_single level |
| nasal | blood_PBMC | S1 |   60 | 0.967 | 0.750 | too few pairs |
| nasal | blood_PBMC | S2 |  722 | 0.688 | 0.683 | below ladder / not above baseline |
| nasal | blood_other | S2 |   26 | 0.808 | 0.731 | too few pairs |
| nasal | blood_other_T | S1 |    6 | 0.833 | 0.833 | too few pairs |
| nasal | blood_other_T | S2 |  142 | 0.683 | 0.711 | too few pairs |
| solid | blood_B | S1 |    1 | 0.000 | 1.000 | too few pairs |
| solid | blood_B | S2 |   30 | 0.600 | 0.533 | too few pairs |
| solid | blood_CD4_T | S1 |    1 | 1.000 | 1.000 | too few pairs |
| solid | blood_CD4_T | S2 |   91 | 0.615 | 0.516 | too few pairs |
| solid | blood_CD8_T | S2 |   64 | 0.594 | 0.547 | too few pairs |
| solid | blood_DC | S2 |   17 | 0.529 | 0.588 | too few pairs |
| solid | blood_Mono | S1 |    1 | 1.000 | 1.000 | too few pairs |
| solid | blood_Mono | S2 |   75 | 0.613 | 0.587 | too few pairs |
| solid | blood_NK | S2 |   43 | 0.581 | 0.605 | too few pairs |
| solid | blood_PBMC | S1 |    3 | 0.333 | 0.667 | too few pairs |
| solid | blood_PBMC | S2 |  107 | 0.673 | 0.505 | too few pairs |
| solid | blood_other | S2 |    1 | 1.000 | 1.000 | too few pairs |
| solid | blood_other_T | S2 |   16 | 0.438 | 0.562 | too few pairs |

A rung is earned only if the agreement beats the majority baseline AND reaches an existing rung's validated range on at least two independent references with n >= 200.
Expected if the layer behaves like its parent: S1 ~0.96, S2 ~0.85, S3 ~0.70 against measured blood.
The nasal and solid references test transfer out of blood; the brain bridge is not addressed by this layer.

