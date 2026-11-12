# Reported tissue contexts

Run 2026-11-12 17:20. Window 250,000 bp. Rows 32,146. Never a vote and never a source of best_direction.

| context | evidence | pairs | cpgs | genes | negative | two_cohorts_agree | agreement_with_measured |
|---|---|---|---|---|---|---|---|
| fibroblasts | measured_eqtm | 1,344 | 969 | 256 | 0.545 | 0 | 0.721 [0.674, 0.764] (n=380) |
| liver | measured_eqtm | 5,032 | 4,307 | 2,671 | 0.611 | 337 | 0.847 [0.826, 0.865] (n=1,306) |
| pancreatic_islets | measured_eqtm | 19,319 | 15,904 | 3,876 | 0.538 | 0 | 0.805 [0.788, 0.822] (n=2,149) |
| retina | measured_eqtm | 4,486 | 3,930 | 4,485 | 0.579 | 0 | 0.842 [0.792, 0.882] (n=253) |
| saliva | smr | 1,965 | 1,208 | 386 | 0.467 | 0 | 0.762 [0.705, 0.811] (n=248) |

`agreement_with_measured`: sign agreement with every pair of the measured table (any tissue, majority sign); the bar for a voting source is 0.85.

- rows read from the parsed tables: bonder2014_liver  3,788; olsson2014_islets 31,315; advani2024_retina 13,475; wagner2014_fibroblast  1,625
- liver: Bonder 2014 and Horvath liver share 337 pairs, same sign in 337; pairs dropped as cohort ties: 0
- saliva flags: MHC 466, 17q21.31 inversion 134

