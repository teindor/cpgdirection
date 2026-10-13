# cpgdirection

Predict the direction of a CpG's effect on gene expression, for CpGs with no
measured eQTM record.

```
+1  higher methylation → HIGHER expression
-1  higher methylation → LOWER  expression
```

## Install

From Bioconductor, once the package is accepted (it is under review at
[BiocContributions #174](https://github.com/Bioconductor/BiocContributions/issues/174)):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("cpgdirection")
```

The development version from GitHub, together with the data package it depends
on:

```r
BiocManager::install(c("teindor/cpgdirectionData", "teindor/cpgdirection"))
```

The ten large data layers (catalogue lookups, SMR directions, probe positions,
brain directions, saliva-bridge scores, the EPIC v2 probe-gene annotation,
distal links, Onco-eQTM consensus; about 120 MB) are served by the companion
ExperimentHub package **cpgdirectionData** and fetched on first use, then
cached by ExperimentHub. The small tables ship inside this package. An offline
or pinned copy of the layers can be named with
`options(cpgdirection.data_dir = "<dir>")`; see `?cpgd_has_data`.

`cpg_direction_universal()` (the distance-only layer) builds a gene TSS table
from `TxDb.Hsapiens.UCSC.hg19.knownGene` and `org.Hs.eg.db` on first use and
caches it under `tools::R_user_dir("cpgdirection", "cache")`. Those annotation
packages are Suggests, so install them only if you want that layer
(`BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db",
"GenomicFeatures", "GenomicRanges", "AnnotationDbi"))`, about 400 MB).

Verify the installation with `cpgd_selftest()`.

## Use

```r
library(cpgdirection)

cpgd_selftest()        # run once after installing

res <- cpg_expression_direction(my_cpgs)
```

One call runs everything: the measured catalogue, all three tissue models, and
the distance fallback. Read two columns first.

| column | meaning |
|---|---|
| `best_direction` | +1, −1 or `NA` — the answer from the strongest evidence available |
| `best_evidence` | which layer it came from |

`best_evidence` is the one that stops you over-reading the result:

| value | worth | meaning |
|---|---|---|
| `measured` | 0.96–1.00 | an eQTM cohort measured this CpG against this gene; ~1.00 when 2+ cohorts agree |
| `catalogue_consensus` | 0.64–0.87 | tissues that called it agreed |
| `catalogue_single` | 0.55–0.87 | one tissue only, no corroboration |
| `distance_only` | 0.60–0.65 | distance to the TSS, nothing else |
| `tissue_conflict` | — | tissues disagreed; not resolved |
| `no_evidence` | — | no catalogue holds the pair and no distance available |

A `measured` row and a `distance_only` row are not comparable. Sorting on
`best_direction` without reading `best_evidence` will mix them.

Where tissues disagree the answer is `NA`, deliberately. The distance curve is
weaker than either tissue it would be adjudicating, so it is not used to break
the tie.

For one tissue only, pass it:

```r
cpg_expression_direction(my_cpgs, tissue = "blood")
```

Identifiers are stripped to canonical form automatically, so EPIC v2 suffixes
are fine:

```r
cpg_expression_direction(c("cg00000029_TC21", "cg00050692_TC21_DNMT3A"))
```

In the second case the gene is parsed out of the name and used as the target,
which is usually what you want.

**Supply the gene when you know it.** Without one the function falls back to the
nearest transcription start site, which is frequently not the gene you had in
mind:

```r
cpg_expression_direction(cpgs  = c("cg00000029", "cg00050692"),
                         genes = c("RBL2",       "DNMT3A"))
```

The `gene_source` column records which path was taken — `user_supplied`,
`parsed_from_name`, `catalogue` or `nearest_TSS`. Treat `nearest_TSS` rows as
provisional.

A file path also works, and a `data.frame` with `cpg` and optional `gene`
columns:

```r
cpg_expression_direction("my_cpgs.txt")
```

## Every target gene, one row per pair: `cpg_gene_pairs()`

A CpG is not obliged to regulate one gene, and 20.3% of CpGs with measured
evidence associate with more than one. `cpg_gene_pairs()` is the
pair-preserving interface: give it **bare** CpG identifiers and it discovers
every plausible target automatically -- from the packaged EPIC v2 manifest
annotation, the tissue catalogue, measured eQTMs and SMR -- then resolves
direction and evidence separately for every CpG x gene pair.

```r
pairs <- cpg_gene_pairs(c("cg26261055", "cg06495038", "cg26196496"))

pairs[, .(cpg_id, target_gene, mapping_primary, best_direction, best_evidence)]
#>        cpg_id target_gene mapping_primary best_direction  best_evidence
#>    cg26261055       CRHBP EPICv2_manifest             -1  distance_only
#>    ...
```

You never have to supply a gene. `genes = NULL` means *auto-discover all
supported targets*; writing `"cg26261055_CRHBP"` or passing `genes=` is an
optional way to add or restrict a hypothesis, not a requirement:

```r
# restrict discovered targets to a biological gene set (what DMSA does)
cpg_gene_pairs(my_cpgs, genes = my_gene_set, gene_mode = "filter")

# evaluate explicit CpG-gene hypotheses, element by element
cpg_gene_pairs(my_cpgs, genes = my_genes, gene_mode = "pairwise")
```

The mapping layer itself is multi-source and cross-validated: candidate genes
come from the union of five annotation tracks (Illumina UCSC RefGene and
GencodeV41, the Zhou lab GENCODE v41 transcript annotation, and the Exeter
GENCODE 47 re-annotation including its regulatory-element assignments), and
every pair records which tracks support it (`n_annotation_sources`; 84% of
pairs have two or more). Probe reliability is checked separately: CpGs whose
every replicate probe is flagged as cross-hybridizing, degenerately mapped or
SNP-contaminated are excluded by default (`probe_qc = "exclude"`, 3.4% of the
array), because a probe that reads the wrong locus cannot be rescued by any
direction evidence downstream.

Two facts are kept apart on every row: `mapping_sources` says how the gene
entered the candidate set (manifest annotation is *annotation*, not proof of
regulation), and `best_evidence` says what the direction rests on -- resolved
by the same evidence ladder as `cpg_expression_direction()`, strictly within
the pair. Opposite directions for different genes of one CpG are both
retained, and abstention is pair-specific. For a fixed pair the direction
record is invariant to how the gene was discovered.

`cpg_expression_direction()` is unchanged: it remains the convenience
one-result-per-CpG interface.

## What you get back

One row per CpG:

| column | meaning |
|---|---|
| `direction` | +1, −1, or `NA` |
| `direction_label` | the same in words |
| `probability_plus1` | calibrated probability of a positive association |
| `confidence` | 0–1 |
| `evidence_tier` | A, B or C |
| `expected_accuracy` | what a direction at this tier is worth |
| `call` | +1, −1 or `ABSTAIN` after thresholds |
| `status` | `DIRECT_eQTM`, `PREDICTED`, `ABSTAIN`, `NO_EXPRESSION_ANCHOR`, `UNMAPPED` |
| `gene_source` | how the target gene was chosen |
| `note` | caveat in plain language |

`summary(res)` gives coverage and the direction split; `cpgd_accuracy_table()`
prints the validation figures behind `expected_accuracy`.

## What a direction is worth

| status | expected accuracy | basis |
|---|---|---|
| `DIRECT_eQTM` | 0.96–1.00 | read from a measured catalogue, not predicted; ~1.00 when replicated |
| tier A | 0.77–0.85 | locked test 0.766; independent cohort 0.854 |
| tier B | 0.64–0.85 | locked test 0.642; independent cohort 0.854 |
| baseline | 0.58–0.63 | always guessing −1 |

Tier A means catalogued CpGs lie within 50 kb **and** the target gene is itself
catalogued. Tier B means one of those. Tier C means neither, and the function
always abstains.

## Three things to know before you use this

**It gives direction, not existence.** Every accuracy figure is conditional on
the CpG being an eQTM. The catalogues behind the model contain only significant
associations, so nothing here can tell you whether a CpG affects expression at
all. Feed it a CpG with no real effect and you will get a confident direction
for something that is not there. `ABSTAIN` and `NO_EXPRESSION_ANCHOR` are
correct answers.

**Do not override the abstentions.** Tier C is badly calibrated — that is why it
abstains rather than reporting a number.

**The table is hg19.** EPIC v2 manifests ship GRCh38. Matching by probe
identifier, as this package does, is unaffected. Matching by position is not.

## Tissue

Three tissues ship, each with its own model, calibration and lookup table:

```r
cpg_expression_direction(my_cpgs, tissue = "blood")             # default
cpg_expression_direction(my_cpgs, tissue = "nasal_epithelium")
cpg_expression_direction(my_cpgs, tissue = "solid_tissue")
```

| tissue | cohorts | pairs | CpGs | grouped-CV AUC |
|---|---|---|---|---|
| blood | 4 | 1,663,799 | 641,325 | 0.793 |
| nasal_epithelium | 2 | 1,472,435 | 553,145 | 0.865 |
| solid_tissue | 6 (TCGA) | 849,229 | 454,821 | 0.664 |

**The three are not equally trustworthy, and the package says so in
`expected_accuracy`.**

*Blood* is the only arm with genuine external validation: a pre-registered
locked test set (tier A AUC 0.781, 76.6% accuracy against a 59.2% baseline) and
an independent white-blood-cell cohort (85.4%).

*Nasal epithelium* has the highest cross-validated AUC, but both contributing
cohorts come from one study and share participants, so 0.865 is within-study
consistency, not external validation. Treat it as promising rather than proven.

*Solid tissue* is tumour tissue from TCGA. Two of six tissues (brain, kidney)
did not exceed their own majority baselines under leave-one-cohort-out. It is
provisional, and brain results in particular should not be read as applying to
healthy brain.

The same CpG can return different answers in different tissues. That is the
finding, not a bug: catalogues differ in which CpGs they contain, and base
rates differ from 0.334 to 0.508 across tissues.

## Code provenance

The design, data assembly, validation and all methodological decisions are the
author's. The R code was drafted and edited with the assistance of Claude
(Anthropic) across the development series and was run, tested and checked by
the author; no code was copied from other software projects. Stated here and
in `DESCRIPTION` under Bioconductor's AI and third-party code policy.
Assisted-by: Claude (Anthropic).

## Citation

[Manuscript citation to be added on acceptance.]

Analysis code, fitted models, and the precomputed prediction table are available
at <https://doi.org/10.17605/OSF.IO/U3VFK>.

The package source is released under the MIT licence.
