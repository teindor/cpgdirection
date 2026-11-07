# Same-person overlap: brain SMR layer vs its cortex references

Run 2026-11-04 14:36. Jaffe signs at p < 0.01 (all-ages age-adjusted model). Cluster bootstrap by CpG, 1000 resamples.

## Overlap-adjusted estimates (no comparison between donor-sharing datasets is used)

| quantity | estimate | lo | hi |
|---|---|---|---|
| ROSMAP_sign_accuracy | 0.981 | 0.962 | 0.998 |
| Jaffe_sign_accuracy | 0.990 | 0.972 | 1.000 |
| T1_ratio_to_ROSMAP | 1.000 | 0.992 | 1.011 |
| T2_ratio_to_ROSMAP | 0.954 | 0.903 | 1.008 |
| T1_accuracy | 0.981 | 0.961 | 1.000 |
| T2_accuracy | 0.960 | 0.929 | 0.991 |
| T1_vs_independent_cortex | 0.964 | 0.926 | 0.998 |
| T2_vs_independent_cortex | 0.943 | 0.899 | 0.984 |
| T1_vs_ROSMAP_insample | 0.995 | 0.991 | 0.998 |
| T2_vs_ROSMAP_insample | 0.939 | 0.930 | 0.947 |
| measured_blood_accuracy | 0.853 | 0.832 | 0.877 |
| measured_blood_vs_ROSMAP | 0.840 | 0.822 | 0.857 |
| blood_SMR_hi_accuracy | 0.776 | 0.752 | 0.804 |
| blood_SMR_lo_accuracy | 0.676 | 0.653 | 0.699 |

## Dose-response by SMR strength

| tier | p_SMR_tertile | n | insample_vs_ROSMAP | ratio_via_independent |
|---|---|---|---|---|
| T1 | 1 | 506 | 0.998 | 1.000 |
| T1 | 2 | 505 | 0.998 | 1.000 |
| T1 | 3 | 506 | 0.988 | 1.000 |
| T2 | 1 | 966 | 0.978 | 1.018 |
| T2 | 2 | 966 | 0.967 | 0.931 |
| T2 | 3 | 966 | 0.872 | 0.858 |

## Log-odds weights for the brain-target ladder

| source | log_odds_weight |
|---|---|
| measured_brain | 3.97 |
| brain_T1 | 3.97 |
| brain_T2 | 3.17 |
| measured_blood | 1.76 |
| blood_SMR_hi | 1.24 |
| blood_SMR_lo | 0.73 |

| judged_by | pairs | ladder_correct | weighted_vote_correct | decisions_differ |
|---|---|---|---|---|
| ROSMAP | 5812 | 0.8999 | 0.9000 | 1 |
| Jaffe |  596 | 0.8087 | 0.8087 | 0 |

T2 opposed by both measured blood and blood SMR (replicated/high), ROSMAP-judged: 43 pairs, ROSMAP sides with T2 on 36. The two blood sources err together on 86 pairs against 14.6 expected if independent, so their weights do not add.
