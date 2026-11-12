# NOTE: exploratory build script kept for provenance; the input paths inside are those of the build machine (12 Nov 2026).
"""Wagner 2014 fibroblasts: CpGs are given by hg18 position. Map to probe IDs through the local hg18->hg19 shift of the paired gene
(hg19 gene start/end minus the hg18 transcript start/end in the file); a shifted position that lands exactly on an array CpG is accepted.
Usage: python -I p2.py <wagner.txt> <out.csv.gz>"""
import sys
import numpy as np
import pandas as pd
w = pd.read_csv(sys.argv[1], sep="\t"); w.columns = ["chr", "pos", "tstart", "tend", "strand", "rho", "p", "target_gene", "probeloc"]
w["c"] = w.chr.astype(str); w["g"] = w.target_gene.astype(str).str.upper()
E = pd.read_csv("/home/claude/w6/inputs/ens_grch37.csv.gz"); E["g"] = E.symbol.astype(str).str.upper(); E["c"] = E.chr.astype(str).str.replace("chr", "")
ann = pd.read_csv("/home/claude/lk/raw/ann_hg19.bed", sep="\t", header=None, comment="t", names=["chr", "s", "e", "name", "sc", "st"])
ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower(); ann["c"] = ann.chr.astype(str).str.replace("chr", "")
ann = ann[ann.cpg_id.str.startswith("cg")]
m = w.reset_index().merge(E[["g", "c", "start", "end"]], on=["g", "c"])
print("rows", len(w), "| rows with an hg19 gene of the same symbol and chromosome", m["index"].nunique())
m["off_a"] = m.start - m.tstart; m["off_b"] = m.end - m.tend
res = {}
for conv, col, add in [("s", "s", 0), ("s+1", "s", 1), ("e", "e", 0), ("e-1", "e", -1)]:
    P = ann.assign(key=ann[col] + add)[["c", "key", "cpg_id"]]
    hits = []
    for oc in ["off_a", "off_b"]:
        for d in [0, 1, -1]:
            h = m.assign(key=m.pos + m[oc] + d).merge(P, on=["c", "key"])[["index", "cpg_id", "key"]].assign(d=d, oc=oc); hits.append(h)
    H = pd.concat(hits)
    exact = H[H.d == 0].drop_duplicates(["index", "cpg_id"])
    print(f"convention {conv:4s}: rows mapped with zero residual {exact['index'].nunique()} | residual counts {H.groupby('d')['index'].nunique().to_dict()}")
    res[conv] = H
H = pd.concat([res["e-1"][res["e-1"].d == 0], res["e"][res["e"].d == 0]])
X = H.drop_duplicates(["index", "cpg_id"])
amb = X.groupby("index").cpg_id.nunique(); X = X[X["index"].isin(amb[amb == 1].index)].drop_duplicates("index")
print("accepted: shifted position on either base of an array CpG | rows mapped to exactly one probe", len(X), "of", len(w), "| ambiguous rows dropped", int((amb > 1).sum()))
P = pd.concat([ann.assign(key=ann.e - 1)[["c", "key", "cpg_id"]], ann.assign(key=ann.e)[["c", "key", "cpg_id"]]])
ch = pd.concat([m.assign(key=m.pos + m[oc] + 37).merge(P, on=["c", "key"])[["index"]] for oc in ["off_a", "off_b"]])["index"].nunique()
print("chance matches after a +37 bp shift:", ch)
o = w.loc[X["index"]].assign(cpg_id=X.cpg_id.values, pos_hg19=X.key.values)
o.to_csv(sys.argv[2], index=False); print("written", len(o), "| CpGs", o.cpg_id.nunique(), "| genes", o.g.nunique(), "| negative", round(float((o.rho < 0).mean()), 3))
