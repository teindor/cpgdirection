# SMR layer with INTERVAL, rule D (gate 1e-40) and a held-out accuracy per pair

Run 2026-11-12 16:46.

- pairs 768,862: 549,316 of the current layer, 219,546 seen in INTERVAL alone
- measured blood sign available for 25,877 pairs (validation); measured outside blood only for 13,215
- chromosome: package position table 725,502; 450K manifest 42,827; INTERVAL run 0; none 533

## Rungs

| rung | pairs | validated n | SMR sign = measured blood sign |
|---|---|---|---|
| replicated | 108,135 | 13,283 | 0.963 [0.960, 0.966] |
| high | 14,284 | 942 | 0.949 [0.933, 0.961] |
| moderate | 528,204 | 8,608 | 0.859 [0.852, 0.867] |
| weak | 52,878 | 1,377 | 0.755 [0.731, 0.777] |
| reported | 65,361 | 1,667 | 0.412 [0.389, 0.436] |
| all but reported | 703,501 | 24,210 | 0.914 [0.910, 0.917] |

Moves (rows: rung in the current layer; columns: rung now):

```
col_0                 replicated   high  moderate   weak  reported
row_0                                                             
high                        7918  11789         0      0       381
moderate                   33259      0    333638    400     11186
new (INTERVAL alone)           0   2495    183076  33975         0
replicated                 66958      0         0   1578         0
reported                       0      0     11490   4422     52061
weak                           0      0         0  12503      1733
```

## Per-pair accuracy (leave one chromosome out)

- validated pairs 25,875 over 22 chromosomes; AUC 0.850; Brier 0.0770; mean calibration gap over deciles 0.009
- pairs measured outside blood only (n 13,215; never used for fitting): AUC 0.654
- pairs with no chromosome (figure from the all-chromosome fit, which contains none of them): 533

Calibration by decile of the held-out figure:

```
                   n  predicted  observed
q                                        
(0.038, 0.67]   2588      0.425     0.458
(0.67, 0.828]   2588      0.759     0.761
(0.828, 0.896]  2587      0.867     0.872
(0.896, 0.932]  2587      0.915     0.903
(0.932, 0.954]  2588      0.944     0.933
(0.954, 0.971]  2587      0.963     0.957
(0.971, 0.983]  2587      0.978     0.966
(0.983, 0.99]   2588      0.987     0.980
(0.99, 0.994]   2587      0.992     0.987
(0.994, 0.999]  2588      0.996     0.996
```

## Accuracy bands

All pairs:

| band | pairs | validated n | SMR sign = measured blood sign |
|---|---|---|---|
| >= 0.95 | 79,216 | 13,600 | 0.975 [0.972, 0.978] |
| 0.90-0.95 | 130,870 | 4,296 | 0.917 [0.909, 0.925] |
| 0.80-0.90 | 198,683 | 3,437 | 0.865 [0.853, 0.876] |
| 0.65-0.80 | 255,410 | 2,160 | 0.738 [0.719, 0.756] |
| 0.55-0.65 | 65,907 | 658 | 0.611 [0.573, 0.647] |
| 0.50-0.55 | 9,728 | 265 | 0.604 [0.544, 0.661] |
| < 0.50 | 29,048 | 1,461 | 0.325 [0.302, 0.350] |

Pairs that can supply a direction (rung not `reported`):

| band | pairs | validated n | SMR sign = measured blood sign |
|---|---|---|---|
| >= 0.95 | 79,162 | 13,587 | 0.975 [0.973, 0.978] |
| 0.90-0.95 | 130,193 | 4,269 | 0.917 [0.908, 0.925] |
| 0.80-0.90 | 194,037 | 3,352 | 0.866 [0.854, 0.877] |
| 0.65-0.80 | 234,639 | 1,970 | 0.743 [0.723, 0.762] |
| 0.55-0.65 | 50,663 | 502 | 0.643 [0.601, 0.684] |
| 0.50-0.55 | 5,513 | 165 | 0.606 [0.530, 0.677] |
| < 0.50 | 9,294 | 365 | 0.452 [0.402, 0.503] |
| supplying (>= 0.55, not reported) | 688,694 | 23,680 | 0.923 [0.919, 0.926] |

Same bands on pairs measured outside blood only (not used for fitting):

| band | pairs | validated n | SMR sign = measured sign |
|---|---|---|---|
| >= 0.95 | 79,162 | 3,042 | 0.782 [0.767, 0.796] |
| 0.90-0.95 | 130,193 | 2,601 | 0.696 [0.678, 0.713] |
| 0.80-0.90 | 194,037 | 2,781 | 0.633 [0.615, 0.651] |
| 0.65-0.80 | 234,639 | 2,237 | 0.539 [0.518, 0.559] |
| 0.55-0.65 | 50,663 | 688 | 0.433 [0.397, 0.470] |
| 0.50-0.55 | 5,513 | 176 | 0.432 [0.361, 0.506] |
| < 0.50 | 9,294 | 356 | 0.424 [0.374, 0.476] |

`reported`-rung pairs by band (excluded in every band):

| band | pairs | validated n | SMR sign = measured blood sign |
|---|---|---|---|
| >= 0.95 | 54 | 13 | 0.923 [0.667, 0.986] |
| 0.90-0.95 | 677 | 27 | 0.963 [0.817, 0.993] |
| 0.80-0.90 | 4,646 | 85 | 0.835 [0.742, 0.899] |
| 0.65-0.80 | 20,771 | 190 | 0.679 [0.610, 0.741] |
| 0.55-0.65 | 15,244 | 156 | 0.506 [0.429, 0.584] |
| 0.50-0.55 | 4,215 | 100 | 0.600 [0.502, 0.691] |
| < 0.50 | 19,754 | 1,096 | 0.283 [0.257, 0.310] |

## Child blood (HELIX) for the pairs seen in INTERVAL alone

looked up in the HELIX full file: {'untested': 152150, 'null': 56157, 'nominal': 9295, 'sig': 1944}; rule check on 148,719 pairs of the current layer that share those CpGs: same class in 1.0000

## Output

`smr_directions_1.1.0.csv.gz`: 768,862 rows; rungs moderate 528,204, replicated 108,135, reported 65,361, weak 52,878, high 14,284; supplying a direction (accuracy >= 0.55, rung not reported) 688,694.
