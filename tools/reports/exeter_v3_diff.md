# Exeter re-annotation swap: Exeter_GENCODEv47 -> Exeter_GENCODEv49

Built 2026-10-02 07:54 from `EPICv2_reannotated_manifest_v3.0.csv.gz`.

| metric | before | after |
|---|---:|---:|
| CpG-gene pairs | 996697 | 1006212 |
| CpGs with >= 1 gene | 791550 | 793233 |
| genes | 38950 | 39289 |
| pairs supported by >= 2 tracks | 84.4% | 83.8% |
| pairs on the Exeter gene track | 851976 | 864733 |
| pairs on the Exeter regulatory track | 416200 | 416200 |

* this Exeter release carries no Promoter_2000bp/Enhancer_5000bp columns; the Exeter_regulatory track was kept unchanged from the previous table
* pairs gained: 9647; pairs lost: 132 (lost pairs were supported only by the old Exeter tracks: 30978)
* pairs whose Exeter feature label changed: 513551
* pairs that gained or lost a TSS200/TSS1500 label: 14214

Examples of changed features (first 20):

| cpg | gene | v47 feature | v49 feature |
|---|---|---|---|
| cg00000236 | VDAC3 | 3UTR | 3UTR_exon |
| cg00000289 | ACTN1 | 3UTR | 3UTR_exon |
| cg00000321 | SFRP1 | Promoter_2000bp;TSS1500 | TSS1500 |
| cg00000622 | NIPA2 | Promoter_2000bp;TSS1500;TSS200 | TSS1500;TSS200 |
| cg00000721 | CARMIL1 | 5UTR;intron | 5UTR_intron;intron |
| cg00000765 | GALNT2 | Enhancer_5000bp;intron | intron |
| cg00000769 | DDX55 | Promoter_2000bp;TSS1500;TSS200 | TSS1500;TSS200 |
| cg00000884 | TLR2 | 5UTR | 5UTR_intron |
| cg00000974 | BMP2 | 5UTR | 5UTR_exon;5UTR_intron;TSS200 |
| cg00001099 | PSKH2 | Promoter_2000bp;intron | intron |
| cg00001103 | PTCHD1-AS | TSS1500;intron | intron;TSS1500 |
| cg00001126 | GBE1 | Enhancer_5000bp;intron | intron |
| cg00001136 | SHANK2 | Promoter_2000bp;intron | intron |
| cg00001224 | ABCB1 | 5UTR | 5UTR_intron |
| cg00001245 | MRPS25 | 5UTR | 5UTR_exon |
| cg00001249 | LRRC9 | 5UTR | 5UTR_intron |
| cg00001349 | MAEL | 5UTR;TSS1500;TSS200 | 5UTR_exon;5UTR_intron;TSS200 |
| cg00001364 | PROX1 | CDS_exon | CDS_exon;intron |
| cg00001446 | ELOVL1 | CDS_exon | CDS_exon;intron |
| cg00001510 | LILRB3 | 5UTR;Promoter_2000bp;TSS1500 | 5UTR_intron;TSS1500 |
