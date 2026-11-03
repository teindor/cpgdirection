# SMR-tested pairs as a second pair source -- exploration

Source: `/Users/teindor/finngen_multiome/smr/bychr/smrd_<ctx>_chr<c>.smr` (SMR v1.03 two-BESD mode; exposure = GoDMC dense cis mQTL, outcome = FinnGen L1 cis-eQTL; HEIDI vs 1000G EUR).
Direction = sign(b_SMR) at the top cis-mQTL SNP. HEIDI pass = p >= 0.05 (table contract). 'instrument' = pair present in the context's clumped-instrument directions table; 'eQTLGen' = pair present in the shipped GoDMC x eQTLGen layer.
Rungs for comparison: smr_moderate 0.84-0.86, smr_high 0.95-0.97, replicated S2 0.958.

## PBMC

Tested pairs 6,001,027; p_SMR < 5e-8: 144,017; NEW at p < 5e-8: 38,786 (HEIDI pass 4,973), over 20,090 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.991 [0.990, 0.991] of 90,353 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |  2013 | 0.996 [0.992, 0.998] |
| S2 | 88324 | 0.991 [0.990, 0.991] |
| S3 |    16 | 0.875 [0.640, 0.965] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  853 | 0.897 [0.875, 0.916] | 0.552 |
| NEW | 5e-8..1e-5 |  188 | 0.835 [0.775, 0.881] | 0.569 |
| NEW | 1e-5..1e-3 |  328 | 0.720 [0.669, 0.765] | 0.634 |
| NEW | >1e-3 | 1874 | 0.631 [0.609, 0.652] | 0.543 |
| both | <5e-8 | 3193 | 0.937 [0.928, 0.945] | 0.586 |
| both | 5e-8..1e-5 |   85 | 0.882 [0.797, 0.935] | 0.588 |
| both | 1e-5..1e-3 |   84 | 0.762 [0.661, 0.840] | 0.643 |
| both | >1e-3 |  233 | 0.678 [0.616, 0.735] | 0.601 |
| eqtlgen_only | <5e-8 | 1249 | 0.915 [0.898, 0.929] | 0.600 |
| eqtlgen_only | 5e-8..1e-5 |  363 | 0.832 [0.790, 0.867] | 0.667 |
| eqtlgen_only | 1e-5..1e-3 |  471 | 0.781 [0.742, 0.816] | 0.635 |
| eqtlgen_only | >1e-3 | 1692 | 0.688 [0.665, 0.710] | 0.608 |
| instrument_only | <5e-8 |  675 | 0.889 [0.863, 0.910] | 0.584 |
| instrument_only | 5e-8..1e-5 |   23 | 0.783 [0.581, 0.903] | 0.739 |
| instrument_only | 1e-5..1e-3 |   17 | 0.412 [0.216, 0.640] | 0.706 |
| instrument_only | >1e-3 |   43 | 0.605 [0.456, 0.736] | 0.605 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail |  782 | 0.893 [0.869, 0.912] |
| NEW | not_tested |    4 | 1.000 [0.510, 1.000] |
| NEW | pass |   67 | 0.940 [0.856, 0.977] |
| both | fail | 2940 | 0.936 [0.927, 0.944] |
| both | not_tested |   34 | 0.912 [0.770, 0.970] |
| both | pass |  219 | 0.950 [0.912, 0.972] |
| eqtlgen_only | fail | 1075 | 0.915 [0.897, 0.931] |
| eqtlgen_only | not_tested |   12 | 0.917 [0.646, 0.985] |
| eqtlgen_only | pass |  162 | 0.914 [0.860, 0.948] |
| instrument_only | fail |  608 | 0.882 [0.853, 0.905] |
| instrument_only | not_tested |   14 | 1.000 [0.785, 1.000] |
| instrument_only | pass |   53 | 0.943 [0.846, 0.981] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail | 782 | 0.893 [0.869, 0.912] |
| <5e-8 | not_tested |   4 | 1.000 [0.510, 1.000] |
| <5e-8 | pass |  67 | 0.940 [0.856, 0.977] |
| 5e-8..1e-5 | fail | 132 | 0.803 [0.727, 0.862] |
| 5e-8..1e-5 | not_tested |   3 | 0.333 [0.061, 0.792] |
| 5e-8..1e-5 | pass |  53 | 0.943 [0.846, 0.981] |
| 1e-5..1e-3 | fail | 221 | 0.674 [0.610, 0.733] |
| 1e-5..1e-3 | not_tested |  13 | 0.846 [0.578, 0.957] |
| 1e-5..1e-3 | pass |  94 | 0.809 [0.717, 0.875] |
| >1e-3 | fail | 875 | 0.618 [0.586, 0.650] |
| >1e-3 | not_tested | 118 | 0.686 [0.598, 0.763] |
| >1e-3 | pass | 881 | 0.636 [0.603, 0.667] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 845 | 0.901 [0.879, 0.919] |
| 1e-30..1e-15 |   7 | 0.571 [0.250, 0.842] |
| 1e-15..1e-10 |   1 | 0.000 [0.000, 0.793] |

## CD4_T

Tested pairs 5,166,662; p_SMR < 5e-8: 107,483; NEW at p < 5e-8: 28,269 (HEIDI pass 3,798), over 16,609 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.990 [0.990, 0.991] of 67,708 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |  1439 | 0.994 [0.989, 0.997] |
| S2 | 66262 | 0.990 [0.989, 0.991] |
| S3 |     7 | 0.571 [0.250, 0.842] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  706 | 0.907 [0.883, 0.926] | 0.576 |
| NEW | 5e-8..1e-5 |  173 | 0.780 [0.713, 0.836] | 0.578 |
| NEW | 1e-5..1e-3 |  289 | 0.716 [0.662, 0.765] | 0.599 |
| NEW | >1e-3 | 1822 | 0.602 [0.579, 0.624] | 0.583 |
| both | <5e-8 | 2502 | 0.926 [0.915, 0.935] | 0.588 |
| both | 5e-8..1e-5 |   71 | 0.789 [0.680, 0.868] | 0.620 |
| both | 1e-5..1e-3 |   73 | 0.575 [0.461, 0.682] | 0.562 |
| both | >1e-3 |  179 | 0.620 [0.547, 0.688] | 0.648 |
| eqtlgen_only | <5e-8 | 1105 | 0.890 [0.870, 0.907] | 0.605 |
| eqtlgen_only | 5e-8..1e-5 |  323 | 0.845 [0.802, 0.881] | 0.647 |
| eqtlgen_only | 1e-5..1e-3 |  513 | 0.733 [0.693, 0.769] | 0.569 |
| eqtlgen_only | >1e-3 | 2245 | 0.646 [0.626, 0.665] | 0.612 |
| instrument_only | <5e-8 |  548 | 0.891 [0.862, 0.914] | 0.613 |
| instrument_only | 5e-8..1e-5 |   28 | 0.750 [0.566, 0.873] | 0.536 |
| instrument_only | 1e-5..1e-3 |   19 | 0.421 [0.231, 0.637] | 0.684 |
| instrument_only | >1e-3 |   31 | 0.516 [0.348, 0.680] | 0.613 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail |  636 | 0.904 [0.879, 0.925] |
| NEW | not_tested |    5 | 1.000 [0.566, 1.000] |
| NEW | pass |   65 | 0.923 [0.832, 0.967] |
| both | fail | 2293 | 0.925 [0.913, 0.935] |
| both | not_tested |   27 | 0.889 [0.719, 0.961] |
| both | pass |  182 | 0.940 [0.895, 0.966] |
| eqtlgen_only | fail |  950 | 0.889 [0.868, 0.908] |
| eqtlgen_only | not_tested |   11 | 0.909 [0.623, 0.984] |
| eqtlgen_only | pass |  144 | 0.889 [0.827, 0.930] |
| instrument_only | fail |  497 | 0.889 [0.859, 0.914] |
| instrument_only | not_tested |   11 | 1.000 [0.741, 1.000] |
| instrument_only | pass |   40 | 0.875 [0.739, 0.945] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail | 636 | 0.904 [0.879, 0.925] |
| <5e-8 | not_tested |   5 | 1.000 [0.566, 1.000] |
| <5e-8 | pass |  65 | 0.923 [0.832, 0.967] |
| 5e-8..1e-5 | fail | 129 | 0.775 [0.696, 0.839] |
| 5e-8..1e-5 | not_tested |   9 | 0.556 [0.267, 0.811] |
| 5e-8..1e-5 | pass |  35 | 0.857 [0.706, 0.937] |
| 1e-5..1e-3 | fail | 212 | 0.684 [0.619, 0.743] |
| 1e-5..1e-3 | not_tested |   5 | 1.000 [0.566, 1.000] |
| 1e-5..1e-3 | pass |  72 | 0.792 [0.684, 0.869] |
| >1e-3 | fail | 792 | 0.558 [0.523, 0.592] |
| >1e-3 | not_tested |  95 | 0.653 [0.553, 0.741] |
| >1e-3 | pass | 935 | 0.633 [0.602, 0.663] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 698 | 0.911 [0.888, 0.930] |
| 1e-30..1e-15 |   7 | 0.571 [0.250, 0.842] |
| 1e-15..1e-10 |   1 | 0.000 [0.000, 0.793] |

## CD8_T

Tested pairs 4,636,981; p_SMR < 5e-8: 72,908; NEW at p < 5e-8: 18,562 (HEIDI pass 2,924), over 12,080 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.991 [0.990, 0.992] of 44,791 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   809 | 0.988 [0.977, 0.993] |
| S2 | 43975 | 0.991 [0.990, 0.992] |
| S3 |     7 | 0.857 [0.487, 0.974] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  607 | 0.895 [0.868, 0.917] | 0.550 |
| NEW | 5e-8..1e-5 |  189 | 0.831 [0.771, 0.877] | 0.603 |
| NEW | 1e-5..1e-3 |  271 | 0.742 [0.686, 0.790] | 0.502 |
| NEW | >1e-3 | 1815 | 0.632 [0.610, 0.654] | 0.615 |
| both | <5e-8 | 2075 | 0.927 [0.915, 0.937] | 0.582 |
| both | 5e-8..1e-5 |   56 | 0.875 [0.764, 0.938] | 0.554 |
| both | 1e-5..1e-3 |   46 | 0.826 [0.693, 0.909] | 0.717 |
| both | >1e-3 |  155 | 0.671 [0.594, 0.740] | 0.594 |
| eqtlgen_only | <5e-8 | 1029 | 0.878 [0.856, 0.896] | 0.613 |
| eqtlgen_only | 5e-8..1e-5 |  450 | 0.816 [0.777, 0.849] | 0.569 |
| eqtlgen_only | 1e-5..1e-3 |  535 | 0.845 [0.812, 0.873] | 0.617 |
| eqtlgen_only | >1e-3 | 2459 | 0.656 [0.637, 0.675] | 0.617 |
| instrument_only | <5e-8 |  419 | 0.881 [0.846, 0.908] | 0.623 |
| instrument_only | 5e-8..1e-5 |   15 | 0.667 [0.417, 0.848] | 0.800 |
| instrument_only | 1e-5..1e-3 |    9 | 0.667 [0.354, 0.879] | 0.667 |
| instrument_only | >1e-3 |   35 | 0.600 [0.436, 0.744] | 0.771 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail |  542 | 0.886 [0.856, 0.910] |
| NEW | not_tested |    2 | 1.000 [0.342, 1.000] |
| NEW | pass |   63 | 0.968 [0.891, 0.991] |
| both | fail | 1859 | 0.925 [0.912, 0.936] |
| both | not_tested |   27 | 0.889 [0.719, 0.961] |
| both | pass |  189 | 0.952 [0.912, 0.975] |
| eqtlgen_only | fail |  854 | 0.874 [0.850, 0.894] |
| eqtlgen_only | not_tested |   13 | 0.923 [0.667, 0.986] |
| eqtlgen_only | pass |  162 | 0.895 [0.838, 0.933] |
| instrument_only | fail |  350 | 0.869 [0.829, 0.900] |
| instrument_only | not_tested |    9 | 0.889 [0.565, 0.980] |
| instrument_only | pass |   60 | 0.950 [0.863, 0.983] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail | 542 | 0.886 [0.856, 0.910] |
| <5e-8 | not_tested |   2 | 1.000 [0.342, 1.000] |
| <5e-8 | pass |  63 | 0.968 [0.891, 0.991] |
| 5e-8..1e-5 | fail | 135 | 0.822 [0.749, 0.878] |
| 5e-8..1e-5 | not_tested |   7 | 0.857 [0.487, 0.974] |
| 5e-8..1e-5 | pass |  47 | 0.851 [0.723, 0.926] |
| 1e-5..1e-3 | fail | 190 | 0.726 [0.659, 0.785] |
| 1e-5..1e-3 | not_tested |   8 | 0.750 [0.409, 0.929] |
| 1e-5..1e-3 | pass |  73 | 0.781 [0.673, 0.860] |
| >1e-3 | fail | 770 | 0.608 [0.573, 0.642] |
| >1e-3 | not_tested |  89 | 0.708 [0.606, 0.792] |
| >1e-3 | pass | 956 | 0.644 [0.613, 0.674] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 599 | 0.902 [0.875, 0.923] |
| 1e-30..1e-15 |   7 | 0.429 [0.158, 0.750] |
| 1e-15..1e-10 |   1 | 0.000 [0.000, 0.793] |

## Mono

Tested pairs 4,600,907; p_SMR < 5e-8: 80,898; NEW at p < 5e-8: 21,007 (HEIDI pass 2,684), over 12,922 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.991 [0.990, 0.992] of 50,443 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   982 | 0.994 [0.987, 0.997] |
| S2 | 49447 | 0.991 [0.990, 0.992] |
| S3 |    14 | 0.929 [0.685, 0.987] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  524 | 0.865 [0.833, 0.891] | 0.565 |
| NEW | 5e-8..1e-5 |  135 | 0.719 [0.637, 0.788] | 0.637 |
| NEW | 1e-5..1e-3 |  306 | 0.716 [0.663, 0.763] | 0.572 |
| NEW | >1e-3 | 1707 | 0.616 [0.593, 0.639] | 0.640 |
| both | <5e-8 | 1930 | 0.921 [0.908, 0.932] | 0.601 |
| both | 5e-8..1e-5 |   61 | 0.803 [0.687, 0.884] | 0.574 |
| both | 1e-5..1e-3 |   62 | 0.871 [0.766, 0.933] | 0.581 |
| both | >1e-3 |  178 | 0.556 [0.483, 0.627] | 0.539 |
| eqtlgen_only | <5e-8 |  931 | 0.887 [0.865, 0.906] | 0.584 |
| eqtlgen_only | 5e-8..1e-5 |  368 | 0.774 [0.729, 0.814] | 0.579 |
| eqtlgen_only | 1e-5..1e-3 |  511 | 0.787 [0.749, 0.820] | 0.626 |
| eqtlgen_only | >1e-3 | 2458 | 0.643 [0.624, 0.662] | 0.630 |
| instrument_only | <5e-8 |  351 | 0.877 [0.839, 0.908] | 0.613 |
| instrument_only | 5e-8..1e-5 |    9 | 0.667 [0.354, 0.879] | 0.556 |
| instrument_only | 1e-5..1e-3 |    9 | 0.778 [0.453, 0.937] | 0.778 |
| instrument_only | >1e-3 |   33 | 0.697 [0.527, 0.826] | 0.545 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail |  480 | 0.865 [0.831, 0.892] |
| NEW | not_tested |    5 | 1.000 [0.566, 1.000] |
| NEW | pass |   39 | 0.846 [0.703, 0.928] |
| both | fail | 1756 | 0.918 [0.904, 0.930] |
| both | not_tested |   18 | 0.944 [0.742, 0.990] |
| both | pass |  156 | 0.955 [0.910, 0.978] |
| eqtlgen_only | fail |  808 | 0.892 [0.869, 0.912] |
| eqtlgen_only | not_tested |   15 | 1.000 [0.796, 1.000] |
| eqtlgen_only | pass |  108 | 0.833 [0.752, 0.892] |
| instrument_only | fail |  305 | 0.866 [0.823, 0.899] |
| instrument_only | not_tested |    7 | 1.000 [0.646, 1.000] |
| instrument_only | pass |   39 | 0.949 [0.831, 0.986] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail | 480 | 0.865 [0.831, 0.892] |
| <5e-8 | not_tested |   5 | 1.000 [0.566, 1.000] |
| <5e-8 | pass |  39 | 0.846 [0.703, 0.928] |
| 5e-8..1e-5 | fail | 101 | 0.663 [0.567, 0.748] |
| 5e-8..1e-5 | not_tested |   1 | 1.000 [0.207, 1.000] |
| 5e-8..1e-5 | pass |  33 | 0.879 [0.727, 0.952] |
| 1e-5..1e-3 | fail | 209 | 0.632 [0.564, 0.694] |
| 1e-5..1e-3 | not_tested |   7 | 0.714 [0.359, 0.918] |
| 1e-5..1e-3 | pass |  90 | 0.911 [0.834, 0.954] |
| >1e-3 | fail | 647 | 0.614 [0.576, 0.650] |
| >1e-3 | not_tested |  78 | 0.641 [0.530, 0.739] |
| >1e-3 | pass | 982 | 0.616 [0.585, 0.646] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 517 | 0.868 [0.837, 0.895] |
| 1e-30..1e-15 |   6 | 0.500 [0.188, 0.812] |
| 1e-15..1e-10 |   1 | 1.000 [0.207, 1.000] |

## NK

Tested pairs 4,009,017; p_SMR < 5e-8: 52,600; NEW at p < 5e-8: 12,916 (HEIDI pass 1,965), over 9,185 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.992 [0.991, 0.993] of 32,290 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   549 | 0.987 [0.974, 0.994] |
| S2 | 31735 | 0.992 [0.991, 0.993] |
| S3 |     6 | 0.833 [0.436, 0.970] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  454 | 0.894 [0.863, 0.919] | 0.579 |
| NEW | 5e-8..1e-5 |  179 | 0.799 [0.734, 0.851] | 0.570 |
| NEW | 1e-5..1e-3 |  305 | 0.764 [0.713, 0.808] | 0.626 |
| NEW | >1e-3 | 1773 | 0.611 [0.588, 0.633] | 0.586 |
| both | <5e-8 | 1617 | 0.917 [0.902, 0.929] | 0.582 |
| both | 5e-8..1e-5 |   64 | 0.859 [0.754, 0.924] | 0.531 |
| both | 1e-5..1e-3 |   41 | 0.707 [0.555, 0.824] | 0.610 |
| both | >1e-3 |  147 | 0.667 [0.587, 0.738] | 0.571 |
| eqtlgen_only | <5e-8 |  881 | 0.882 [0.859, 0.902] | 0.598 |
| eqtlgen_only | 5e-8..1e-5 |  358 | 0.824 [0.781, 0.860] | 0.615 |
| eqtlgen_only | 1e-5..1e-3 |  534 | 0.749 [0.711, 0.784] | 0.597 |
| eqtlgen_only | >1e-3 | 2610 | 0.672 [0.654, 0.690] | 0.620 |
| instrument_only | <5e-8 |  275 | 0.858 [0.812, 0.894] | 0.618 |
| instrument_only | 5e-8..1e-5 |   12 | 0.833 [0.552, 0.953] | 0.833 |
| instrument_only | 1e-5..1e-3 |   19 | 0.895 [0.686, 0.971] | 0.789 |
| instrument_only | >1e-3 |   26 | 0.500 [0.321, 0.679] | 0.577 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail |  381 | 0.879 [0.843, 0.908] |
| NEW | not_tested |    3 | 1.000 [0.438, 1.000] |
| NEW | pass |   70 | 0.971 [0.902, 0.992] |
| both | fail | 1400 | 0.912 [0.896, 0.926] |
| both | not_tested |   18 | 1.000 [0.824, 1.000] |
| both | pass |  199 | 0.940 [0.898, 0.965] |
| eqtlgen_only | fail |  741 | 0.879 [0.853, 0.900] |
| eqtlgen_only | not_tested |   10 | 1.000 [0.722, 1.000] |
| eqtlgen_only | pass |  130 | 0.892 [0.827, 0.935] |
| instrument_only | fail |  221 | 0.842 [0.788, 0.884] |
| instrument_only | not_tested |    7 | 1.000 [0.646, 1.000] |
| instrument_only | pass |   47 | 0.915 [0.801, 0.966] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail |  381 | 0.879 [0.843, 0.908] |
| <5e-8 | not_tested |    3 | 1.000 [0.438, 1.000] |
| <5e-8 | pass |   70 | 0.971 [0.902, 0.992] |
| 5e-8..1e-5 | fail |  126 | 0.786 [0.706, 0.848] |
| 5e-8..1e-5 | pass |   53 | 0.830 [0.708, 0.908] |
| 1e-5..1e-3 | fail |  161 | 0.733 [0.660, 0.795] |
| 1e-5..1e-3 | not_tested |    5 | 0.600 [0.231, 0.882] |
| 1e-5..1e-3 | pass |  139 | 0.806 [0.732, 0.863] |
| >1e-3 | fail |  674 | 0.631 [0.593, 0.666] |
| >1e-3 | not_tested |   88 | 0.602 [0.498, 0.698] |
| >1e-3 | pass | 1011 | 0.598 [0.568, 0.628] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 446 | 0.901 [0.870, 0.926] |
| 1e-30..1e-15 |   7 | 0.571 [0.250, 0.842] |
| 1e-15..1e-10 |   1 | 0.000 [0.000, 0.793] |

## B

Tested pairs 3,685,964; p_SMR < 5e-8: 43,134; NEW at p < 5e-8: 12,406 (HEIDI pass 1,756), over 7,793 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.988 [0.987, 0.989] of 24,914 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   505 | 0.988 [0.974, 0.995] |
| S2 | 24400 | 0.988 [0.987, 0.990] |
| S3 |     9 | 0.778 [0.453, 0.937] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  376 | 0.902 [0.867, 0.928] | 0.585 |
| NEW | 5e-8..1e-5 |  149 | 0.826 [0.757, 0.878] | 0.685 |
| NEW | 1e-5..1e-3 |  202 | 0.807 [0.747, 0.855] | 0.644 |
| NEW | >1e-3 | 1799 | 0.609 [0.586, 0.632] | 0.615 |
| both | <5e-8 | 1105 | 0.918 [0.900, 0.932] | 0.583 |
| both | 5e-8..1e-5 |   48 | 0.667 [0.525, 0.783] | 0.646 |
| both | 1e-5..1e-3 |   33 | 0.788 [0.622, 0.893] | 0.636 |
| both | >1e-3 |   96 | 0.688 [0.589, 0.771] | 0.646 |
| eqtlgen_only | <5e-8 |  702 | 0.913 [0.890, 0.932] | 0.604 |
| eqtlgen_only | 5e-8..1e-5 |  351 | 0.815 [0.771, 0.852] | 0.630 |
| eqtlgen_only | 1e-5..1e-3 |  539 | 0.863 [0.831, 0.889] | 0.551 |
| eqtlgen_only | >1e-3 | 2882 | 0.650 [0.632, 0.667] | 0.617 |
| instrument_only | <5e-8 |  216 | 0.875 [0.824, 0.913] | 0.583 |
| instrument_only | 5e-8..1e-5 |    7 | 0.714 [0.359, 0.918] | 0.714 |
| instrument_only | 1e-5..1e-3 |    5 | 1.000 [0.566, 1.000] | 0.800 |
| instrument_only | >1e-3 |   27 | 0.407 [0.245, 0.593] | 0.593 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail | 300 | 0.893 [0.853, 0.923] |
| NEW | not_tested |   4 | 1.000 [0.510, 1.000] |
| NEW | pass |  72 | 0.931 [0.848, 0.970] |
| both | fail | 956 | 0.911 [0.891, 0.928] |
| both | not_tested |  10 | 0.900 [0.596, 0.982] |
| both | pass | 139 | 0.964 [0.919, 0.985] |
| eqtlgen_only | fail | 527 | 0.907 [0.879, 0.929] |
| eqtlgen_only | not_tested |  13 | 0.923 [0.667, 0.986] |
| eqtlgen_only | pass | 162 | 0.932 [0.883, 0.962] |
| instrument_only | fail | 182 | 0.857 [0.799, 0.901] |
| instrument_only | not_tested |   6 | 1.000 [0.610, 1.000] |
| instrument_only | pass |  28 | 0.964 [0.823, 0.994] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail |  300 | 0.893 [0.853, 0.923] |
| <5e-8 | not_tested |    4 | 1.000 [0.510, 1.000] |
| <5e-8 | pass |   72 | 0.931 [0.848, 0.970] |
| 5e-8..1e-5 | fail |  104 | 0.837 [0.754, 0.895] |
| 5e-8..1e-5 | not_tested |    9 | 0.778 [0.453, 0.937] |
| 5e-8..1e-5 | pass |   36 | 0.806 [0.650, 0.902] |
| 1e-5..1e-3 | fail |  132 | 0.833 [0.761, 0.887] |
| 1e-5..1e-3 | not_tested |    3 | 0.667 [0.208, 0.939] |
| 1e-5..1e-3 | pass |   67 | 0.761 [0.647, 0.847] |
| >1e-3 | fail |  582 | 0.567 [0.526, 0.607] |
| >1e-3 | not_tested |   65 | 0.569 [0.448, 0.682] |
| >1e-3 | pass | 1152 | 0.633 [0.605, 0.660] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 369 | 0.908 [0.874, 0.933] |
| 1e-30..1e-15 |   7 | 0.571 [0.250, 0.842] |

## DC

Tested pairs 3,059,776; p_SMR < 5e-8: 24,963; NEW at p < 5e-8: 6,639 (HEIDI pass 891), over 4,669 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.990 [0.989, 0.992] of 14,105 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   216 | 0.991 [0.967, 0.997] |
| S2 | 13883 | 0.990 [0.989, 0.992] |
| S3 |     6 | 1.000 [0.610, 1.000] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  279 | 0.853 [0.807, 0.890] | 0.595 |
| NEW | 5e-8..1e-5 |  130 | 0.908 [0.846, 0.946] | 0.531 |
| NEW | 1e-5..1e-3 |  167 | 0.886 [0.829, 0.926] | 0.611 |
| NEW | >1e-3 | 1541 | 0.681 [0.658, 0.704] | 0.622 |
| both | <5e-8 |  854 | 0.893 [0.871, 0.912] | 0.601 |
| both | 5e-8..1e-5 |   31 | 0.935 [0.793, 0.982] | 0.548 |
| both | 1e-5..1e-3 |   44 | 0.818 [0.680, 0.905] | 0.659 |
| both | >1e-3 |   88 | 0.648 [0.544, 0.739] | 0.545 |
| eqtlgen_only | <5e-8 |  610 | 0.941 [0.919, 0.957] | 0.641 |
| eqtlgen_only | 5e-8..1e-5 |  370 | 0.873 [0.835, 0.903] | 0.573 |
| eqtlgen_only | 1e-5..1e-3 |  493 | 0.844 [0.809, 0.873] | 0.643 |
| eqtlgen_only | >1e-3 | 2788 | 0.689 [0.672, 0.706] | 0.598 |
| instrument_only | <5e-8 |  146 | 0.808 [0.737, 0.864] | 0.610 |
| instrument_only | 5e-8..1e-5 |   14 | 0.929 [0.685, 0.987] | 0.571 |
| instrument_only | 1e-5..1e-3 |    4 | 0.750 [0.301, 0.954] | 0.500 |
| instrument_only | >1e-3 |   15 | 0.533 [0.301, 0.752] | 0.533 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail | 229 | 0.843 [0.790, 0.884] |
| NEW | not_tested |   1 | 1.000 [0.207, 1.000] |
| NEW | pass |  49 | 0.898 [0.782, 0.956] |
| both | fail | 699 | 0.880 [0.854, 0.902] |
| both | not_tested |   7 | 0.857 [0.487, 0.974] |
| both | pass | 148 | 0.959 [0.914, 0.981] |
| eqtlgen_only | fail | 476 | 0.941 [0.916, 0.959] |
| eqtlgen_only | not_tested |   1 | 1.000 [0.207, 1.000] |
| eqtlgen_only | pass | 133 | 0.940 [0.886, 0.969] |
| instrument_only | fail | 112 | 0.795 [0.711, 0.859] |
| instrument_only | not_tested |   4 | 1.000 [0.510, 1.000] |
| instrument_only | pass |  30 | 0.833 [0.664, 0.927] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail | 229 | 0.843 [0.790, 0.884] |
| <5e-8 | not_tested |   1 | 1.000 [0.207, 1.000] |
| <5e-8 | pass |  49 | 0.898 [0.782, 0.956] |
| 5e-8..1e-5 | fail |  85 | 0.906 [0.825, 0.952] |
| 5e-8..1e-5 | not_tested |   2 | 1.000 [0.342, 1.000] |
| 5e-8..1e-5 | pass |  43 | 0.907 [0.784, 0.963] |
| 1e-5..1e-3 | fail |  88 | 0.898 [0.817, 0.945] |
| 1e-5..1e-3 | not_tested |   5 | 1.000 [0.566, 1.000] |
| 1e-5..1e-3 | pass |  74 | 0.865 [0.769, 0.925] |
| >1e-3 | fail | 553 | 0.678 [0.638, 0.716] |
| >1e-3 | not_tested |  53 | 0.585 [0.451, 0.707] |
| >1e-3 | pass | 935 | 0.689 [0.658, 0.718] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 275 | 0.855 [0.808, 0.891] |
| 1e-30..1e-15 |   4 | 0.750 [0.301, 0.954] |

## other_T

Tested pairs 2,942,667; p_SMR < 5e-8: 19,174; NEW at p < 5e-8: 4,556 (HEIDI pass 873), over 3,818 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.991 [0.989, 0.993] of 11,246 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   176 | 0.989 [0.960, 0.997] |
| S2 | 11067 | 0.991 [0.989, 0.993] |
| S3 |     3 | 0.667 [0.208, 0.939] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |  224 | 0.866 [0.815, 0.905] | 0.567 |
| NEW | 5e-8..1e-5 |  135 | 0.800 [0.725, 0.859] | 0.504 |
| NEW | 1e-5..1e-3 |  176 | 0.778 [0.711, 0.833] | 0.619 |
| NEW | >1e-3 | 1618 | 0.632 [0.608, 0.655] | 0.609 |
| both | <5e-8 |  681 | 0.912 [0.888, 0.931] | 0.605 |
| both | 5e-8..1e-5 |   25 | 0.800 [0.609, 0.911] | 0.760 |
| both | 1e-5..1e-3 |   36 | 0.889 [0.747, 0.956] | 0.750 |
| both | >1e-3 |   60 | 0.583 [0.457, 0.699] | 0.533 |
| eqtlgen_only | <5e-8 |  561 | 0.873 [0.843, 0.898] | 0.647 |
| eqtlgen_only | 5e-8..1e-5 |  302 | 0.854 [0.810, 0.890] | 0.613 |
| eqtlgen_only | 1e-5..1e-3 |  512 | 0.850 [0.816, 0.878] | 0.596 |
| eqtlgen_only | >1e-3 | 2749 | 0.642 [0.624, 0.659] | 0.590 |
| instrument_only | <5e-8 |  114 | 0.877 [0.804, 0.925] | 0.693 |
| instrument_only | 5e-8..1e-5 |    5 | 0.800 [0.376, 0.964] | 0.600 |
| instrument_only | 1e-5..1e-3 |    7 | 0.571 [0.250, 0.842] | 0.714 |
| instrument_only | >1e-3 |   15 | 0.800 [0.548, 0.930] | 0.667 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail | 178 | 0.837 [0.776, 0.884] |
| NEW | not_tested |   5 | 1.000 [0.566, 1.000] |
| NEW | pass |  41 | 0.976 [0.874, 0.996] |
| both | fail | 544 | 0.901 [0.873, 0.923] |
| both | not_tested |  12 | 0.917 [0.646, 0.985] |
| both | pass | 125 | 0.960 [0.910, 0.983] |
| eqtlgen_only | fail | 399 | 0.845 [0.806, 0.877] |
| eqtlgen_only | not_tested |  14 | 1.000 [0.785, 1.000] |
| eqtlgen_only | pass | 148 | 0.939 [0.888, 0.968] |
| instrument_only | fail |  94 | 0.872 [0.790, 0.925] |
| instrument_only | not_tested |   2 | 1.000 [0.342, 1.000] |
| instrument_only | pass |  18 | 0.889 [0.672, 0.969] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail |  178 | 0.837 [0.776, 0.884] |
| <5e-8 | not_tested |    5 | 1.000 [0.566, 1.000] |
| <5e-8 | pass |   41 | 0.976 [0.874, 0.996] |
| 5e-8..1e-5 | fail |   82 | 0.768 [0.666, 0.846] |
| 5e-8..1e-5 | not_tested |    3 | 1.000 [0.438, 1.000] |
| 5e-8..1e-5 | pass |   50 | 0.840 [0.715, 0.917] |
| 1e-5..1e-3 | fail |  106 | 0.698 [0.605, 0.777] |
| 1e-5..1e-3 | not_tested |    4 | 0.750 [0.301, 0.954] |
| 1e-5..1e-3 | pass |   66 | 0.909 [0.816, 0.958] |
| >1e-3 | fail |  496 | 0.647 [0.604, 0.688] |
| >1e-3 | not_tested |   65 | 0.615 [0.494, 0.724] |
| >1e-3 | pass | 1057 | 0.626 [0.597, 0.655] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 220 | 0.877 [0.827, 0.914] |
| 1e-30..1e-15 |   3 | 0.333 [0.061, 0.792] |
| 1e-15..1e-10 |   1 | 0.000 [0.000, 0.793] |

## other

Tested pairs 1,588,753; p_SMR < 5e-8: 2,426; NEW at p < 5e-8: 629 (HEIDI pass 126), over 587 CpGs.
Sanity: sign(b_SMR) agrees with the instrument-based direction on 0.997 [0.992, 0.999] of 1,135 shared p<5e-8 pairs.

### sanity by instrument tier

| inst_tier | n | agree |
|---|---|---|
| S1 |   10 | 1.000 [0.722, 1.000] |
| S2 | 1125 | 0.997 [0.992, 0.999] |

### measured blood: by novelty x p_SMR

| novelty | p_bin | n | agree | baseline |
|---|---|---|---|---|
| NEW | <5e-8 |   49 | 0.898 [0.782, 0.956] | 0.531 |
| NEW | 5e-8..1e-5 |   37 | 0.865 [0.720, 0.941] | 0.568 |
| NEW | 1e-5..1e-3 |   67 | 0.821 [0.713, 0.894] | 0.642 |
| NEW | >1e-3 | 1143 | 0.591 [0.562, 0.619] | 0.603 |
| both | <5e-8 |  102 | 0.931 [0.865, 0.966] | 0.559 |
| both | 5e-8..1e-5 |    5 | 0.800 [0.376, 0.964] | 0.600 |
| both | 1e-5..1e-3 |    9 | 1.000 [0.701, 1.000] | 0.778 |
| both | >1e-3 |   11 | 0.636 [0.354, 0.848] | 0.545 |
| eqtlgen_only | <5e-8 |  104 | 0.933 [0.868, 0.967] | 0.712 |
| eqtlgen_only | 5e-8..1e-5 |  109 | 0.844 [0.764, 0.900] | 0.624 |
| eqtlgen_only | 1e-5..1e-3 |  221 | 0.932 [0.891, 0.958] | 0.561 |
| eqtlgen_only | >1e-3 | 2271 | 0.668 [0.649, 0.687] | 0.609 |
| instrument_only | <5e-8 |    9 | 1.000 [0.701, 1.000] | 0.556 |
| instrument_only | 5e-8..1e-5 |    2 | 1.000 [0.342, 1.000] | 1.000 |
| instrument_only | 1e-5..1e-3 |    1 | 1.000 [0.207, 1.000] | 1.000 |
| instrument_only | >1e-3 |    1 | 1.000 [0.207, 1.000] | 1.000 |

### measured blood: p<5e-8 by novelty x HEIDI

| novelty | heidi_status | n | agree |
|---|---|---|---|
| NEW | fail | 23 | 0.870 [0.679, 0.955] |
| NEW | not_tested |  2 | 1.000 [0.342, 1.000] |
| NEW | pass | 24 | 0.917 [0.742, 0.977] |
| both | fail | 57 | 0.930 [0.833, 0.972] |
| both | not_tested |  2 | 1.000 [0.342, 1.000] |
| both | pass | 43 | 0.930 [0.814, 0.976] |
| eqtlgen_only | fail | 68 | 0.897 [0.802, 0.949] |
| eqtlgen_only | pass | 36 | 1.000 [0.904, 1.000] |
| instrument_only | fail |  7 | 1.000 [0.646, 1.000] |
| instrument_only | not_tested |  1 | 1.000 [0.207, 1.000] |
| instrument_only | pass |  1 | 1.000 [0.207, 1.000] |

### measured blood: NEW pairs by p_SMR x HEIDI

| p_bin | heidi_status | n | agree |
|---|---|---|---|
| <5e-8 | fail |  23 | 0.870 [0.679, 0.955] |
| <5e-8 | not_tested |   2 | 1.000 [0.342, 1.000] |
| <5e-8 | pass |  24 | 0.917 [0.742, 0.977] |
| 5e-8..1e-5 | fail |  22 | 0.909 [0.722, 0.975] |
| 5e-8..1e-5 | not_tested |   1 | 1.000 [0.207, 1.000] |
| 5e-8..1e-5 | pass |  14 | 0.786 [0.524, 0.924] |
| 1e-5..1e-3 | fail |  22 | 0.773 [0.566, 0.899] |
| 1e-5..1e-3 | not_tested |   1 | 1.000 [0.207, 1.000] |
| 1e-5..1e-3 | pass |  44 | 0.841 [0.706, 0.921] |
| >1e-3 | fail | 202 | 0.649 [0.580, 0.711] |
| >1e-3 | not_tested |  33 | 0.545 [0.380, 0.702] |
| >1e-3 | pass | 908 | 0.579 [0.547, 0.611] |

### measured blood: NEW p<5e-8 pairs by exposure (top mQTL) p

| expo_p | n | agree |
|---|---|---|
| <1e-30 | 49 | 0.898 [0.782, 0.956] |

## Coverage summary

| context | tested | p5e8 | p1e5 | new_p5e8 | new_p5e8_heidi_pass | new_p1e5 | cpgs_new_p5e8 |
|---|---|---|---|---|---|---|---|
| PBMC | 6001027 | 144017 | 218929 | 38786 | 4973 | 87643 | 20090 |
| CD4_T | 5166662 | 107483 | 163675 | 28269 | 3798 | 64650 | 16609 |
| CD8_T | 4636981 |  72908 | 116512 | 18562 | 2924 | 44510 | 12080 |
| Mono | 4600907 |  80898 | 127797 | 21007 | 2684 | 50067 | 12922 |
| NK | 4009017 |  52600 |  84294 | 12916 | 1965 | 31127 |  9185 |
| B | 3685964 |  43134 |  69336 | 12406 | 1756 | 28049 |  7793 |
| DC | 3059776 |  24963 |  42497 |  6639 |  891 | 16087 |  4669 |
| other_T | 2942667 |  19174 |  33640 |  4556 |  873 | 12078 |  3818 |
| other | 1588753 |   2426 |   5070 |   629 |  126 |  1885 |   587 |

## Verdict: NEW p<5e-8 pairs vs measured blood

| context | stratum | n | agreement | baseline | verdict |
|---|---|---|---|---|---|
| PBMC | NEW p<5e-8, HEIDI pass |  67 | 0.940 | 0.612 | too few pairs |
| PBMC | NEW p<5e-8, HEIDI fail | 782 | 0.893 | 0.546 | at smr_moderate level |
| PBMC | NEW p<5e-8, HEIDI not_tested |   4 | 1.000 | 0.750 | too few pairs |
| PBMC | NEW p<5e-8, HEIDI all | 853 | 0.897 | 0.552 | at smr_moderate level |
| CD4_T | NEW p<5e-8, HEIDI pass |  65 | 0.923 | 0.631 | too few pairs |
| CD4_T | NEW p<5e-8, HEIDI fail | 636 | 0.904 | 0.571 | at smr_moderate level |
| CD4_T | NEW p<5e-8, HEIDI not_tested |   5 | 1.000 | 0.600 | too few pairs |
| CD4_T | NEW p<5e-8, HEIDI all | 706 | 0.907 | 0.576 | at smr_moderate level |
| CD8_T | NEW p<5e-8, HEIDI pass |  63 | 0.968 | 0.746 | too few pairs |
| CD8_T | NEW p<5e-8, HEIDI fail | 542 | 0.886 | 0.528 | at smr_moderate level |
| CD8_T | NEW p<5e-8, HEIDI not_tested |   2 | 1.000 | 0.500 | too few pairs |
| CD8_T | NEW p<5e-8, HEIDI all | 607 | 0.895 | 0.550 | at smr_moderate level |
| Mono | NEW p<5e-8, HEIDI pass |  39 | 0.846 | 0.615 | too few pairs |
| Mono | NEW p<5e-8, HEIDI fail | 480 | 0.865 | 0.560 | at smr_moderate level |
| Mono | NEW p<5e-8, HEIDI not_tested |   5 | 1.000 | 0.600 | too few pairs |
| Mono | NEW p<5e-8, HEIDI all | 524 | 0.865 | 0.565 | at smr_moderate level |
| NK | NEW p<5e-8, HEIDI pass |  70 | 0.971 | 0.543 | too few pairs |
| NK | NEW p<5e-8, HEIDI fail | 381 | 0.879 | 0.585 | at smr_moderate level |
| NK | NEW p<5e-8, HEIDI not_tested |   3 | 1.000 | 0.667 | too few pairs |
| NK | NEW p<5e-8, HEIDI all | 454 | 0.894 | 0.579 | at smr_moderate level |
| B | NEW p<5e-8, HEIDI pass |  72 | 0.931 | 0.667 | too few pairs |
| B | NEW p<5e-8, HEIDI fail | 300 | 0.893 | 0.563 | at smr_moderate level |
| B | NEW p<5e-8, HEIDI not_tested |   4 | 1.000 | 0.750 | too few pairs |
| B | NEW p<5e-8, HEIDI all | 376 | 0.902 | 0.585 | at smr_moderate level |
| DC | NEW p<5e-8, HEIDI pass |  49 | 0.898 | 0.653 | too few pairs |
| DC | NEW p<5e-8, HEIDI fail | 229 | 0.843 | 0.581 | at smr_moderate level |
| DC | NEW p<5e-8, HEIDI not_tested |   1 | 1.000 | 1.000 | too few pairs |
| DC | NEW p<5e-8, HEIDI all | 279 | 0.853 | 0.595 | at smr_moderate level |
| other_T | NEW p<5e-8, HEIDI pass |  41 | 0.976 | 0.707 | too few pairs |
| other_T | NEW p<5e-8, HEIDI fail | 178 | 0.837 | 0.545 | too few pairs |
| other_T | NEW p<5e-8, HEIDI not_tested |   5 | 1.000 | 0.800 | too few pairs |
| other_T | NEW p<5e-8, HEIDI all | 224 | 0.866 | 0.567 | at smr_moderate level |
| other | NEW p<5e-8, HEIDI pass |  24 | 0.917 | 0.625 | too few pairs |
| other | NEW p<5e-8, HEIDI fail |  23 | 0.870 | 0.609 | too few pairs |
| other | NEW p<5e-8, HEIDI not_tested |   2 | 1.000 | 1.000 | too few pairs |
| other | NEW p<5e-8, HEIDI all |  49 | 0.898 | 0.531 | too few pairs |

A NEW stratum earns the moderate rung only if its lower Wilson bound clears 0.84 on the measured blood reference;
pairs also found by the instrument table or the eQTLGen layer are not new coverage and are scored only as sanity.

