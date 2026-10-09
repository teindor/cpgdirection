# distal_links (ENCODE-rE2G, blood class) vs peripheral references

Layer: 1,038,305 blood-class CpG-gene links over 444,471 CpGs (26 biosamples; brain class not used here).
Precision = linked genes that are a reference gene for that CpG; recall = reference pairs that are linked;
both on the CpGs the two sources share. Wilson 95% CI. The manifest annotation is the comparator.

## Reference: blood (26,545 pairs, 16,796 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 61.1%
- reference pairs that are linked: 28.0% (7,429)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 10261 | 24688 | 0.301 [0.295, 0.307] | 15604 | 0.476 [0.468, 0.484] |
| rE2G_score_0.2-0.5 |  5866 | 11048 | 0.164 [0.157, 0.171] |  8973 | 0.201 [0.193, 0.210] |
| rE2G_score_0.5-0.9 |  3591 |  6070 | 0.246 [0.235, 0.257] |  5627 | 0.265 [0.254, 0.277] |
| rE2G_score_>=0.9 |  6570 |  7570 | 0.546 [0.534, 0.557] |  9920 | 0.416 [0.407, 0.426] |
| rE2G_biosamples_1-2 |  5481 |  9858 | 0.183 [0.175, 0.191] |  8444 | 0.214 [0.205, 0.223] |
| rE2G_biosamples_3-9 |  5104 |  7876 | 0.288 [0.278, 0.298] |  7807 | 0.291 [0.281, 0.301] |
| rE2G_biosamples_>=10 |  5085 |  6954 | 0.483 [0.471, 0.494] |  7464 | 0.450 [0.438, 0.461] |
| manifest_all | 12834 | 18105 | 0.452 [0.444, 0.459] | 19294 | 0.424 [0.417, 0.431] |
| manifest_on_linked_cpgs |  9832 | 14024 | 0.466 [0.458, 0.474] | 14964 | 0.437 [0.429, 0.445] |
| rE2G_on_manifest_cpgs |  9832 | 23334 | 0.306 [0.300, 0.312] | 14964 | 0.478 [0.470, 0.486] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 16196 | 44902 | 172420 |
| FALSE | TRUE |  2174 | 26781 | 113599 |
| TRUE | FALSE |  2920 | 19518 | 107678 |
| TRUE | TRUE |  5255 | 2613 | 50171 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.565)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 3592 | 0.615 [0.599, 0.630] |
| TRUE | 3837 | 0.768 [0.755, 0.781] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 1807 | 0.553 [0.530, 0.576] |
| 0.5-0.9 | 1491 | 0.634 [0.609, 0.658] |
| >=0.9 | 4131 | 0.778 [0.765, 0.790] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 1804 | 0.578 [0.555, 0.600] |
| 3-9 | 2269 | 0.683 [0.664, 0.702] |
| >=10 | 3356 | 0.764 [0.749, 0.778] |

## Reference: nasal (11,045 pairs, 7,682 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 67.5%
- reference pairs that are linked: 29.3% (3,236)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 5186 | 12659 | 0.256 [0.248, 0.263] | 7377 | 0.439 [0.427, 0.450] |
| rE2G_score_0.2-0.5 | 2911 |  5703 | 0.147 [0.139, 0.157] | 4241 | 0.198 [0.187, 0.211] |
| rE2G_score_0.5-0.9 | 1914 |  3253 | 0.220 [0.207, 0.235] | 2756 | 0.260 [0.244, 0.277] |
| rE2G_score_>=0.9 | 3255 |  3703 | 0.453 [0.437, 0.469] | 4710 | 0.356 [0.343, 0.370] |
| rE2G_biosamples_1-2 | 2833 |  5309 | 0.177 [0.167, 0.188] | 4068 | 0.231 [0.219, 0.245] |
| rE2G_biosamples_3-9 | 2655 |  4030 | 0.269 [0.255, 0.283] | 3818 | 0.284 [0.270, 0.298] |
| rE2G_biosamples_>=10 | 2439 |  3320 | 0.365 [0.349, 0.382] | 3577 | 0.339 [0.324, 0.355] |
| manifest_all | 6373 |  8594 | 0.359 [0.349, 0.369] | 9155 | 0.337 [0.327, 0.347] |
| manifest_on_linked_cpgs | 4948 |  6716 | 0.361 [0.350, 0.373] | 7082 | 0.343 [0.332, 0.354] |
| rE2G_on_manifest_cpgs | 4948 | 12010 | 0.257 [0.249, 0.265] | 7082 | 0.435 [0.424, 0.447] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 6906 | 108453 | 221346 |
| FALSE | TRUE | 1055 | 32832 | 135769 |
| TRUE | FALSE |  903 | 4024 | 90534 |
| TRUE | TRUE | 2181 | 1017 | 29182 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.606)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 1771 | 0.652 [0.630, 0.674] |
| TRUE | 1465 | 0.848 [0.828, 0.865] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 |  841 | 0.561 [0.527, 0.594] |
| 0.5-0.9 |  717 | 0.653 [0.617, 0.687] |
| >=0.9 | 1678 | 0.868 [0.851, 0.884] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 |  941 | 0.708 [0.678, 0.736] |
| 3-9 | 1083 | 0.726 [0.698, 0.751] |
| >=10 | 1212 | 0.780 [0.756, 0.802] |

## Reference: solid (37,462 pairs, 25,891 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 52.3%
- reference pairs that are linked: 25.5% (9,560)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 13544 | 30175 | 0.317 [0.312, 0.322] | 19867 | 0.481 [0.474, 0.488] |
| rE2G_score_0.2-0.5 |  7912 | 14058 | 0.179 [0.173, 0.186] | 11619 | 0.217 [0.209, 0.224] |
| rE2G_score_0.5-0.9 |  4888 |  7970 | 0.295 [0.285, 0.305] |  7667 | 0.307 [0.296, 0.317] |
| rE2G_score_>=0.9 |  7111 |  8147 | 0.576 [0.565, 0.586] | 10587 | 0.443 [0.434, 0.453] |
| rE2G_biosamples_1-2 |  7960 | 13776 | 0.230 [0.223, 0.237] | 11736 | 0.270 [0.262, 0.278] |
| rE2G_biosamples_3-9 |  6612 |  9796 | 0.342 [0.332, 0.351] |  9683 | 0.346 [0.336, 0.355] |
| rE2G_biosamples_>=10 |  4657 |  6603 | 0.461 [0.449, 0.473] |  6966 | 0.437 [0.425, 0.449] |
| manifest_all | 21492 | 29409 | 0.494 [0.488, 0.500] | 30449 | 0.477 [0.471, 0.483] |
| manifest_on_linked_cpgs | 12826 | 18296 | 0.490 [0.483, 0.498] | 18741 | 0.479 [0.472, 0.486] |
| rE2G_on_manifest_cpgs | 12826 | 28480 | 0.320 [0.315, 0.326] | 18741 | 0.487 [0.479, 0.494] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 20555 | 48020 | 176280 |
| FALSE | TRUE |  2380 | 19094 | 86138 |
| TRUE | FALSE |  7347 | 19279 | 117152 |
| TRUE | TRUE |  7180 | 2392 | 47186 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.614)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 5671 | 0.620 [0.607, 0.632] |
| TRUE | 3889 | 0.762 [0.749, 0.776] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 2519 | 0.618 [0.599, 0.637] |
| 0.5-0.9 | 2350 | 0.639 [0.619, 0.658] |
| >=0.9 | 4691 | 0.729 [0.717, 0.742] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 3170 | 0.649 [0.632, 0.665] |
| 3-9 | 3346 | 0.662 [0.646, 0.678] |
| >=10 | 3044 | 0.726 [0.710, 0.741] |

## Reference: smr_S1S2 (396,294 pairs, 102,379 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 60.1%
- reference pairs that are linked: 11.3% (44,779)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 61574 | 141212 | 0.317 [0.315, 0.320] | 246292 | 0.182 [0.180, 0.183] |
| rE2G_score_0.2-0.5 | 35694 |  65265 | 0.252 [0.249, 0.256] | 145458 | 0.113 [0.112, 0.115] |
| rE2G_score_0.5-0.9 | 20483 |  34532 | 0.283 [0.279, 0.288] |  86678 | 0.113 [0.111, 0.115] |
| rE2G_score_>=0.9 | 36137 |  41415 | 0.448 [0.443, 0.452] | 146167 | 0.127 [0.125, 0.129] |
| rE2G_biosamples_1-2 | 33213 |  57534 | 0.229 [0.226, 0.233] | 131801 | 0.100 [0.099, 0.102] |
| rE2G_biosamples_3-9 | 30738 |  45891 | 0.321 [0.317, 0.325] | 126239 | 0.117 [0.115, 0.119] |
| rE2G_biosamples_>=10 | 27353 |  37787 | 0.446 [0.441, 0.451] | 113875 | 0.148 [0.146, 0.150] |
| manifest_all | 85219 | 112715 | 0.373 [0.370, 0.376] | 333747 | 0.126 [0.125, 0.127] |
| manifest_on_linked_cpgs | 57920 |  78046 | 0.378 [0.374, 0.381] | 234381 | 0.126 [0.124, 0.127] |
| rE2G_on_manifest_cpgs | 57920 | 131375 | 0.324 [0.321, 0.326] | 234381 | 0.182 [0.180, 0.183] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 334361 | 178930 | 724166 |
| FALSE | TRUE |  19862 | 42296 | 154664 |
| TRUE | FALSE |  17154 | 17844 | 99771 |
| TRUE | TRUE |  24917 | 12779 | 68766 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.509)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 24552 | 0.531 [0.524, 0.537] |
| TRUE | 20227 | 0.580 [0.573, 0.587] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 16457 | 0.510 [0.502, 0.517] |
| 0.5-0.9 |  9784 | 0.547 [0.537, 0.557] |
| >=0.9 | 18538 | 0.595 [0.588, 0.602] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 13201 | 0.515 [0.507, 0.524] |
| 3-9 | 14735 | 0.546 [0.538, 0.554] |
| >=10 | 16843 | 0.589 [0.582, 0.597] |

## Reading the result

- The layer earns its keep as a TARGET source if rE2G precision on the shared CpGs is at or above the manifest's
  and its recall adds pairs the manifest lacks (Q2, Q3).
- It must NOT be promoted to a direction rung on the strength of Q4: a negative-share above base rate in
  promoter-class links is the TSS-distance effect, already on the ladder as distance_only.
- Score bins should show a precision gradient; if they do not, `distal_score_max` is provenance, not a confidence.

