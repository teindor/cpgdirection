# CpG coordinates and array membership (HM27, 450K, EPIC v1, EPIC v2)

Run 2026-11-12 18:45. 1,084,297 probe IDs.

## Outside sources against the package, on the probes both hold

| source | shared_with_package | identical | one_base_apart | elsewhere |
|---|---|---|---|---|
| Illumina 450K manifest, hg19 | 391,451 | 0.9985 | 554 | 43 |
| Illumina EPIC v1 coordinates, hg19 | 702,162 | 0.9989 | 789 | 12 |
| Zhou HM450 manifest, hg38 | 335,482 | 0.9999 | 0 | 42 |
| Zhou EPIC manifest, hg38 | 583,214 | 1.0000 | 0 | 17 |
| Zhou EPICv2 manifest, hg38 | 731,983 | 0.9998 | 20 | 115 |
| Zhou HM27 manifest, hg38 | 24,277 | 0.9998 | 0 | 6 |

## Coverage by array

| array | probe_ids | with_hg19 | with_hg38 | no_hg38 | hg38_ambiguous |
|---|---|---|---|---|---|
| HM27 | 27,578 | 26,521 | 27,578 | 0 | 37 |
| 450K | 485,577 | 485,396 | 485,569 | 8 | 944 |
| EPICv1 | 865,918 | 855,290 | 865,904 | 14 | 900 |
| EPICv2 | 930,658 | 930,197 | 930,655 | 3 | 229 |
| all | 1,084,297 | 1,072,123 | 1,084,278 | 19 | 1,094 |

`no_hg38`: no source places the probe on hg38 (the manifests give no mapping). `hg38_ambiguous`: two sources place the probe differently.

## By type of probe ID

| type | probe_ids | no_hg38 | on_alt_contig |
|---|---|---|---|
| cg | 1,080,317 | 19 | 269 |
| ch | 3,091 | 0 | 0 |
| nv | 824 | 0 | 0 |
| rs | 65 | 0 | 0 |

## Source of the hg38 position

| source_hg38 | N |
|---|---|
| cpgdirection (annotation layer; not stored) | 731,824 |
| zhou_hm450_hg38 | 150,129 |
| zhou_epic_hg38 | 145,064 |
| zhou_epicv2_hg38 | 56,206 |
| zhou_hm27_hg38 | 1,055 |
| none | 19 |

hg19 positions supplied here (the rest are the package's): 141,945. hg38 positions stored here: 352,454.

## Reference positions reproduced

`treeB_hm450_hg38_slim.tsv.gz`: 449,659 probes; identical chromosome and position for 449,659; different or missing for 0.

