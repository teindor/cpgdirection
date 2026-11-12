# NOTE: exploratory build script kept for provenance; the input paths inside are those of the build machine (12 Nov 2026).
"""Parse the newly downloaded measured-eQTM candidates into one long table. Usage: python -I p1.py <leads_dir> <out.csv.gz>"""
import sys
import numpy as np
import pandas as pd
import openpyxl
U, OUT = sys.argv[1], sys.argv[2]
def sheet(path, name, skip):
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb[name]
    it = ws.iter_rows(values_only=True)
    for _ in range(skip): next(it)
    hdr = [str(h) for h in next(it)]
    rows = [r[:len(hdr)] for r in it]
    df = pd.DataFrame(rows, columns=[h if h != 'None' else f'_x{i}' for i, h in enumerate(hdr)])
    return df
out = []
def add(df, source, tissue, n):
    df = df.copy(); df["source"], df["tissue"], df["n_samples"] = source, tissue, n
    df["target_gene"] = df.target_gene.astype(str).str.strip().str.upper()
    df = df[df.cpg_id.astype(str).str.match(r"^c[gh]\d+") & df.direction.isin([-1, 1]) & ~df.target_gene.isin(["", "NONE", "NAN", "NA"])]
    out.append(df[["source", "tissue", "n_samples", "cpg_id", "target_gene", "direction", "p", "dist"]])
    print(f"{source:22s} rows {len(df):6d} | pairs {df.drop_duplicates(['cpg_id','target_gene']).shape[0]:6d} | CpGs {df.cpg_id.nunique():6d} | genes {df.target_gene.nunique():5d} | negative {float((df.direction == -1).mean()):.3f}", flush=True)

# Bonder 2014 liver
b = sheet(f"{U}/12864_2014_6781_MOESM7_ESM.xlsx", "Supplementary_Table5", 0)
b.columns = ["chr", "probe", "cpg_id", "gene_expr", "gene_cpg", "rho", "p"]
b = b.iloc[1:]; b["rho"] = pd.to_numeric(b.rho, errors="coerce"); b["direction"] = np.sign(b.rho)
b["target_gene"] = b.gene_expr.astype(str).str.split(";"); b = b.explode("target_gene"); b["dist"] = np.nan
add(b, "bonder2014_liver", "liver", 158)
# Delahaye 2018 placenta
d = sheet(f"{U}/pgen.1007785.s032.xlsx", "T20_eQTMs", 1)
d["direction"] = np.sign(pd.to_numeric(d.linear_coefficient, errors="coerce")); d = d.rename(columns={"Gene": "target_gene", "CpG": "cpg_id", "P.value": "p"}); d["dist"] = pd.to_numeric(d.dist, errors="coerce")
add(d, "delahaye2018_placenta", "placenta", 74)
# Olsson 2014 islets
o = sheet(f"{U}/pgen.1004735.s025.xlsx", "a) Direct Correlation CpG-mRNA", 0)
o["direction"] = np.sign(pd.to_numeric(o["T-Stat_CpG vs mRNA"], errors="coerce")); o["cpg_id"] = o["CpG site_Illumina Probe ID"]; o["p"] = o["P-Value CpG vs mRNA"]
o["target_gene"] = o["mRNA_Affymetrix Transcript_Gene Symbol"].astype(str).str.split(r"\s*///\s*|;|,"); o = o.explode("target_gene")
pos, st, s0, s1 = [pd.to_numeric(o[c], errors="coerce") if c != "Strand_Affy Gene" else o[c] for c in ["Mapinfo_CpG site", "Strand_Affy Gene", "Start_Affy Gene", "Stop_Affy Gene"]]
tss = np.where(st == "-", s1, s0); o["dist"] = np.where(st == "-", tss - pos, pos - tss)
add(o, "olsson2014_islets", "pancreatic_islets", 89)
# Kennedy 2018 GTP whole blood
k = pd.read_csv(f"{U}/12864_2018_4842_MOESM1_ESM.txt", sep="\t", low_memory=False)
print("  Kennedy status", k.status.value_counts().to_dict(), "| p max", k["p.val"].max(), "| p <= 1e-11:", int((k["p.val"] <= 1e-11).sum()))
k = k[k.status != "TRANS"].copy(); k["strict"] = k["p.val"] <= 1e-11; print("  Kennedy same-chromosome rows", len(k), "| at p <= 1e-11:", int(k.strict.sum()), "| |distance| <= 1 Mb:", int((pd.to_numeric(k.distance, errors="coerce").abs() <= 1e6).sum())); k["direction"] = np.sign(k.beta); k = k.rename(columns={"CpG.probe": "cpg_id", "annot.gene": "target_gene", "p.val": "p"}); k["dist"] = pd.to_numeric(k.distance, errors="coerce")
add(k, "kennedy2018_gtp_blood", "blood", 333)
# Advani 2024 retina
E = pd.read_csv("/home/claude/w6/inputs/ens_grch37.csv.gz").drop_duplicates("ensgene").set_index("ensgene")
G = pd.read_csv("/home/claude/v2r/final3/pools/gene_table.csv.gz").set_index("gene_id")
r = pd.read_csv(f"{U}/eQTM_significant_pairs_retina.txt", sep="\t")
print("  retina rows", len(r), "| p.adj <= 0.05:", int((r["p.adj"] <= 0.05).sum()), "| p.adj max", r["p.adj"].max(), "| pairs per gene max", int(r.groupby("trait").size().max()))
r["target_gene"] = r.trait.map(G.display).fillna(r.trait.map(E.symbol)); print("  retina genes without a symbol", int(r.target_gene.isna().sum()))
r = r[r.target_gene.notna() & (r["p.adj"] <= 0.05)].copy(); r["direction"] = np.sign(r.beta); r = r.rename(columns={"geno": "cpg_id", "p.nom": "p"}); r["dist"] = np.nan
add(r, "advani2024_retina", "retina", 152)
# Wagner 2014 fibroblasts: CpG given by position
w = pd.read_csv("/home/claude/ln/wagner_mapped.csv.gz"); w["direction"] = np.sign(w.rho); tss = np.where(w.strand == "-", w.tend, w.tstart); w["dist"] = np.where(w.strand == "-", tss - w.pos, w.pos - tss)
add(w, "wagner2014_fibroblast", "fibroblasts", 62)
# Liu 2022 kidney
q = sheet(f"{U}/41588_2022_1097_MOESM4_ESM.xlsx", "Table S12", 1)
q["direction"] = np.sign(pd.to_numeric(q.Beta, errors="coerce")); q = q.rename(columns={"CpG": "cpg_id", "Gene Symbol": "target_gene", "P value": "p"}); q["dist"] = pd.to_numeric(q["Distance (bp)"], errors="coerce")
print("  kidney target types", q["Target Type"].value_counts().to_dict())
add(q, "liu2022_kidney", "kidney", 414)
L = pd.concat(out); L["p"] = pd.to_numeric(L.p, errors="coerce"); L["direction"] = L.direction.astype(int)
L.to_csv(OUT, index=False); print("written", len(L))
