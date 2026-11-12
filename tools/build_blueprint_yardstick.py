#!/usr/bin/env python3
"""BLUEPRINT cis eQTMs (Chen et al. 2016, Cell 167:1398; monocytes, neutrophils, naive CD4 T cells of about 200
donors; 450K array, RNA-seq) as an adult-blood yardstick for tools/validate_rung_accuracy.R (--blueprint).
PI rulings of 12 Nov 2026: the pairs are a yardstick only in November (they do not vote in the blood table; the
three cell types become separate contexts in December), and the adjustment is the pooled union of two
pre-specified analyses.

Per cell type, tools/blueprint_eqtm_cis.py (matrices prepared by tools/blueprint_convert_matrices.py from the
public BLUEPRINT phenotype matrices) gives the residual Pearson correlation of every CpG-gene pair within
250 kb of the TSS under two adjustments:
  A  sex and age
  B  sex, age, and the first 5 principal components of methylation and of expression
Here: two-sided p from the t distribution (df = n - 4 for A, n - 14 for B), one Benjamini-Hochberg pass at 5%
over both sets of p-values of the cell type, pair kept when either analysis survives, sign from the analysis
with the smaller p. Genes are named as the catalogue names them (cpgdirectionData gene table), else by the
Ensembl GRCh37 symbol.

Usage: python build_blueprint_yardstick.py --npz-dir <dir with pooled450k_<cell>.npz> --gene-table <gene_table.csv.gz>
         --ens <ens_grch37.csv.gz> --measured <inst/extdata/measured_eqtms.csv.gz> --out <dir>
Writes <out>/blueprint_eqtm_pooled.csv.gz (cell, cpg_id, target_gene, ensembl_id, direction, tss_dist, r_A, r_B,
p_min, analysis; tss_dist is the strand-oriented offset of the CpG from the TSS, of which only the size is
used) and <out>/build_blueprint_yardstick.md.
"""
import argparse, datetime
import numpy as np, pandas as pd
from scipy import stats
ap = argparse.ArgumentParser()
for x in ["npz-dir", "gene-table", "ens", "measured", "out"]: ap.add_argument("--" + x, required=True)
a = ap.parse_args()
def bh(p):
    o = np.argsort(p); q = np.empty(len(p)); q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]; return np.minimum(q, 1)
def wil(k, n, z=1.96):
    if n == 0: return (float("nan"), float("nan"))
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h
G = pd.read_csv(a.gene_table).drop_duplicates("gene_id").set_index("gene_id"); E = pd.read_csv(a.ens).drop_duplicates("ensgene").set_index("ensgene")
M = pd.read_csv(a.measured); M["target_gene"] = M.target_gene.str.upper().str.strip(); M = M[M.direction.isin([-1, 1])]
Mb = M[M.tissue == "blood"].drop_duplicates(["cpg_id", "target_gene"])[["cpg_id", "target_gene", "direction"]].rename(columns={"direction": "d_ref"})
names = {"mono": "monocytes", "neut": "neutrophils", "tcel": "T cells"}
L = [f"# BLUEPRINT cis eQTMs, pooled union of two adjustments (adult-blood yardstick)", "", f"Run {datetime.datetime.now():%Y-%m-%d %H:%M}.", "",
     "| cell type | donors | pairs tested | kept | by A only | by B only | by both (sign conflicts) | negative | in measured blood | same sign |", "|---|---|---|---|---|---|---|---|---|---|"]
out = []
for cell in ["mono", "neut", "tcel"]:
    z = np.load(f"{a.npz_dir}/pooled450k_{cell}.npz", allow_pickle=True); n = int(z["n"])
    t = pd.DataFrame({"cell": cell, "cpg_id": z["cpg_ids"][z["cpg"]], "ensembl_id": z["gene_ids"][z["gene"]], "tss_dist": z["dist"]})
    def pv(k):
        r = z[f"r_m{k}_e{k}"].astype(np.float64); df = n - 4 - 2 * k; return r, 2 * stats.t.sf(np.abs(r * np.sqrt(df / np.maximum(1 - r * r, 1e-12))), df)
    (rA, pA), (rB, pB) = pv(0), pv(5)
    q = bh(np.concatenate([pA, pB])); qA, qB = q[:len(t)], q[len(t):]
    keep = (qA < 0.05) | (qB < 0.05); both = (qA < 0.05) & (qB < 0.05)
    t["r_A"], t["r_B"] = np.round(rA, 4), np.round(rB, 4); t["p_min"] = np.minimum(pA, pB); t["analysis"] = np.where(pA <= pB, "A", "B")
    t["direction"] = np.where(pA <= pB, np.sign(rA), np.sign(rB)).astype(int)
    conflicts = int((np.sign(rA) != np.sign(rB))[both].sum())
    k = t[keep & (t.direction != 0)].copy()
    k["target_gene"] = k.ensembl_id.map(G.display).fillna(k.ensembl_id.map(E.symbol)); k = k[k.target_gene.notna() & (k.target_gene.astype(str) != "")]
    k["target_gene"] = k.target_gene.astype(str).str.upper()
    k = k.sort_values("p_min").drop_duplicates(["cpg_id", "target_gene"])
    x = k.merge(Mb, on=["cpg_id", "target_gene"]); kk = int((x.direction == x.d_ref).sum()); lo, hi = wil(kk, len(x))
    L.append(f"| {names[cell]} | {n} | {len(t):,} | {len(k):,} | {int((keep & ~(qB < 0.05)).sum()):,} | {int((keep & ~(qA < 0.05)).sum()):,} | {int(both.sum()):,} ({conflicts}) | {float((k.direction == -1).mean()):.3f} | {len(x):,} | {kk / max(len(x), 1):.3f} [{lo:.3f}, {hi:.3f}] |")
    out.append(k[["cell", "cpg_id", "target_gene", "ensembl_id", "direction", "tss_dist", "r_A", "r_B", "p_min", "analysis"]])
O = pd.concat(out); O.to_csv(f"{a.out}/blueprint_eqtm_pooled.csv.gz", index=False)
u = O.groupby(["cpg_id", "target_gene"]).direction.agg(["sum", "size", "nunique"]).reset_index()
L += ["", f"Rows {len(O):,}; distinct pairs {len(u):,} ({int((u['size'] > 1).sum()):,} in more than one cell type, of which {int((u['nunique'] > 1).sum()):,} with opposite signs between cell types); "
      f"CpGs {O.cpg_id.nunique():,}; genes {O.target_gene.nunique():,}.", ""]
open(f"{a.out}/build_blueprint_yardstick.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
