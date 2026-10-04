# distal_links (ENCODE-rE2G, blood class) vs peripheral references

Layer: 889,779 blood-class CpG-gene links over 429,647 CpGs (26 biosamples; brain class not used here).
Precision = linked genes that are a reference gene for that CpG; recall = reference pairs that are linked;
both on the CpGs the two sources share. Wilson 95% CI. The manifest annotation is the comparator.

## Reference: blood (16,858 pairs, 11,860 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 74.4%
- reference pairs that are linked: 34.6% (5,839)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all |  8818 | 19403 | 0.301 [0.295, 0.307] | 12654 | 0.461 [0.453, 0.470] |
| rE2G_score_0.2-0.5 |  3974 |  7964 | 0.158 [0.150, 0.166] |  5765 | 0.218 [0.207, 0.229] |
| rE2G_score_0.5-0.9 |  2825 |  4905 | 0.241 [0.229, 0.253] |  4173 | 0.283 [0.270, 0.297] |
| rE2G_score_>=0.9 |  5700 |  6534 | 0.521 [0.509, 0.533] |  8176 | 0.416 [0.405, 0.427] |
| rE2G_biosamples_1-2 |  4268 |  8189 | 0.195 [0.187, 0.204] |  6315 | 0.253 [0.242, 0.264] |
| rE2G_biosamples_3-9 |  3966 |  6232 | 0.309 [0.297, 0.320] |  5621 | 0.342 [0.330, 0.355] |
| rE2G_biosamples_>=10 |  3698 |  4982 | 0.465 [0.451, 0.479] |  5184 | 0.447 [0.434, 0.461] |
| manifest_all | 11215 | 15912 | 0.419 [0.411, 0.426] | 15976 | 0.417 [0.409, 0.425] |
| manifest_on_linked_cpgs |  8437 | 12157 | 0.430 [0.422, 0.439] | 12136 | 0.431 [0.422, 0.440] |
| rE2G_on_manifest_cpgs |  8437 | 18211 | 0.307 [0.301, 0.314] | 12136 | 0.461 [0.452, 0.470] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 8628 | 48660 | 172779 |
| FALSE | TRUE | 1567 | 30360 | 115432 |
| TRUE | FALSE | 2391 | 14519 | 103982 |
| TRUE | TRUE | 4272 | 933 | 35285 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.577)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 3050 | 0.617 [0.600, 0.634] |
| TRUE | 2789 | 0.790 [0.775, 0.805] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 1256 | 0.541 [0.513, 0.568] |
| 0.5-0.9 | 1181 | 0.621 [0.593, 0.648] |
| >=0.9 | 3402 | 0.786 [0.772, 0.799] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 1598 | 0.585 [0.561, 0.609] |
| 3-9 | 1923 | 0.687 [0.666, 0.708] |
| >=10 | 2318 | 0.789 [0.772, 0.805] |

## Reference: nasal (8,319 pairs, 6,289 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 75.8%
- reference pairs that are linked: 35.2% (2,932)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 4769 | 10584 | 0.277 [0.269, 0.286] | 6351 | 0.462 [0.449, 0.474] |
| rE2G_score_0.2-0.5 | 2162 |  4385 | 0.144 [0.134, 0.155] | 2821 | 0.224 [0.209, 0.240] |
| rE2G_score_0.5-0.9 | 1622 |  2780 | 0.237 [0.222, 0.254] | 2136 | 0.309 [0.290, 0.329] |
| rE2G_score_>=0.9 | 3014 |  3419 | 0.479 [0.463, 0.496] | 4046 | 0.405 [0.390, 0.420] |
| rE2G_biosamples_1-2 | 2362 |  4662 | 0.188 [0.177, 0.200] | 3069 | 0.286 [0.270, 0.302] |
| rE2G_biosamples_3-9 | 2181 |  3413 | 0.313 [0.298, 0.329] | 2849 | 0.375 [0.357, 0.393] |
| rE2G_biosamples_>=10 | 1907 |  2509 | 0.393 [0.374, 0.413] | 2583 | 0.382 [0.364, 0.401] |
| manifest_all | 5934 |  7927 | 0.382 [0.371, 0.392] | 7892 | 0.383 [0.373, 0.394] |
| manifest_on_linked_cpgs | 4553 |  6125 | 0.384 [0.372, 0.396] | 6085 | 0.387 [0.374, 0.399] |
| rE2G_on_manifest_cpgs | 4553 | 10001 | 0.279 [0.270, 0.287] | 6085 | 0.458 [0.445, 0.470] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 4456 | 99353 | 212268 |
| FALSE | TRUE |  838 | 34358 | 138088 |
| TRUE | FALSE |  931 | 4074 | 76923 |
| TRUE | TRUE | 2094 | 943 | 27007 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.615)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 1719 | 0.654 [0.631, 0.676] |
| TRUE | 1213 | 0.891 [0.872, 0.907] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 |  633 | 0.553 [0.514, 0.591] |
| 0.5-0.9 |  660 | 0.653 [0.616, 0.688] |
| >=0.9 | 1639 | 0.869 [0.852, 0.884] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 |  877 | 0.709 [0.678, 0.738] |
| 3-9 | 1068 | 0.734 [0.707, 0.760] |
| >=10 |  987 | 0.810 [0.784, 0.833] |

## Reference: solid (1,754 pairs, 1,697 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 70.2%
- reference pairs that are linked: 59.5% (1,044)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 1192 | 2424 | 0.431 [0.411, 0.451] | 1239 | 0.843 [0.821, 0.862] |
| rE2G_score_0.2-0.5 |  474 |  935 | 0.171 [0.148, 0.197] |  485 | 0.330 [0.290, 0.373] |
| rE2G_score_0.5-0.9 |  269 |  558 | 0.260 [0.225, 0.298] |  274 | 0.529 [0.470, 0.587] |
| rE2G_score_>=0.9 |  793 |  931 | 0.794 [0.767, 0.819] |  830 | 0.890 [0.867, 0.910] |
| rE2G_biosamples_1-2 |  485 |  935 | 0.232 [0.206, 0.260] |  497 | 0.437 [0.394, 0.481] |
| rE2G_biosamples_3-9 |  425 |  677 | 0.399 [0.363, 0.436] |  436 | 0.619 [0.573, 0.664] |
| rE2G_biosamples_>=10 |  604 |  812 | 0.686 [0.653, 0.717] |  636 | 0.876 [0.848, 0.899] |
| manifest_all | 1696 | 2273 | 0.749 [0.731, 0.766] | 1753 | 0.971 [0.962, 0.978] |
| manifest_on_linked_cpgs | 1192 | 1670 | 0.718 [0.696, 0.739] | 1239 | 0.968 [0.956, 0.976] |
| rE2G_on_manifest_cpgs | 1192 | 2424 | 0.431 [0.411, 0.451] | 1239 | 0.843 [0.821, 0.862] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE |   27 | 11339 | 82934 |
| FALSE | TRUE |   25 | 742 | 55731 |
| TRUE | FALSE |  683 | 21189 | 124774 |
| TRUE | TRUE | 1019 | 777 | 39370 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.539)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 389 | 0.391 [0.344, 0.440] |
| TRUE | 655 | 0.724 [0.688, 0.757] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 160 | 0.312 [0.246, 0.388] |
| 0.5-0.9 | 145 | 0.386 [0.311, 0.467] |
| >=0.9 | 739 | 0.704 [0.670, 0.735] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 217 | 0.429 [0.365, 0.495] |
| 3-9 | 270 | 0.522 [0.463, 0.581] |
| >=10 | 557 | 0.704 [0.665, 0.740] |

## Reference: smr_S1S2 (396,294 pairs, 102,379 CpGs)

### Q1 coverage

- reference CpGs with >= 1 blood-class link: 58.7%
- reference pairs that are linked: 9.8% (38,641)

### Q2 target proposal quality

| proposal | n_cpg | proposals | precision | ref_pairs | recall |
|---|---|---|---|---|---|
| rE2G_all | 60142 | 125198 | 0.309 [0.306, 0.311] | 240222 | 0.161 [0.159, 0.162] |
| rE2G_score_0.2-0.5 | 28027 |  53525 | 0.225 [0.222, 0.229] | 110862 | 0.109 [0.107, 0.111] |
| rE2G_score_0.5-0.9 | 17744 |  30900 | 0.270 [0.265, 0.275] |  74066 | 0.113 [0.111, 0.115] |
| rE2G_score_>=0.9 | 35735 |  40773 | 0.447 [0.442, 0.452] | 143830 | 0.127 [0.125, 0.128] |
| rE2G_biosamples_1-2 | 29082 |  53106 | 0.226 [0.222, 0.230] | 114164 | 0.105 [0.103, 0.107] |
| rE2G_biosamples_3-9 | 26468 |  39750 | 0.315 [0.310, 0.319] | 105937 | 0.118 [0.116, 0.120] |
| rE2G_biosamples_>=10 | 23717 |  32342 | 0.437 [0.431, 0.442] |  98096 | 0.144 [0.142, 0.146] |
| manifest_all | 85219 | 112715 | 0.373 [0.370, 0.376] | 333747 | 0.126 [0.125, 0.127] |
| manifest_on_linked_cpgs | 56608 |  76370 | 0.377 [0.373, 0.380] | 228654 | 0.126 [0.125, 0.127] |
| rE2G_on_manifest_cpgs | 56608 | 115661 | 0.316 [0.313, 0.318] | 228654 | 0.160 [0.158, 0.161] |

### Q3 distance profile of reference pairs by manifest / link status

| in_manifest | linked | n | median_dist | p90_dist |
|---|---|---|---|---|
| FALSE | FALSE | 339374 | 174980 | 718660 |
| FALSE | TRUE |  14849 | 45486 | 164244 |
| TRUE | FALSE |  18279 | 17887 | 100466 |
| TRUE | TRUE |  23792 | 12578 | 67009 |

Pairs linked but NOT in the manifest are the layer's contribution; their distances show whether it reaches beyond the annotation window.

### Q4 sign among linked pairs (reference base rate of -1: 0.509)

By element class (promoter-class elements are near-TSS, i.e. distance information the ladder already uses):

| promoter_class | n | neg_share |
|---|---|---|
| FALSE | 23799 | 0.531 [0.524, 0.537] |
| TRUE | 14842 | 0.598 [0.590, 0.606] |

By rE2G score:

| score_bin | n | neg_share |
|---|---|---|
| 0.2-0.5 | 12063 | 0.505 [0.497, 0.514] |
| 0.5-0.9 |  8356 | 0.546 [0.536, 0.557] |
| >=0.9 | 18222 | 0.595 [0.588, 0.602] |

By number of biosamples:

| bs_bin | n | neg_share |
|---|---|---|
| 1-2 | 12003 | 0.516 [0.507, 0.525] |
| 3-9 | 12517 | 0.557 [0.548, 0.565] |
| >=10 | 14121 | 0.590 [0.582, 0.598] |

## Reading the result

- The layer earns its keep as a TARGET source if rE2G precision on the shared CpGs is at or above the manifest's
  and its recall adds pairs the manifest lacks (Q2, Q3).
- It must NOT be promoted to a direction rung on the strength of Q4: a negative-share above base rate in
  promoter-class links is the TSS-distance effect, already on the ladder as distance_only.
- Score bins should show a precision gradient; if they do not, `distal_score_max` is provenance, not a confidence.

