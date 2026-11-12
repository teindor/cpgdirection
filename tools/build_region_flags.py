# NOTE: exploratory build script kept for provenance; the input paths inside are those of the build machine (12 Nov 2026).
"""Lead 5: per-CpG region flags (CoRSIVs 2.0; imprinting and allele-specific methylation) and the test the ruling asked for:
is sign agreement lower inside these regions?  Plus the admission check for Voisin 2016 adipose.
Usage: python -I f5.py <corsiv_media-1.xlsx> <rosenski_MOESM4.xlsx> <voisin_s020.xlsx> <candidate_smr.csv.gz> <out_prefix>"""
import sys, re
import numpy as np, pandas as pd
CX, RX, VX, SMR, OUT = sys.argv[1:6]
def wil(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(float(c - h), 3), round(float(c + h), 3)
ann = pd.read_csv("/home/claude/lk/raw/ann_hg19.bed", sep="\t", header=None, comment="t", names=["chr", "s", "e", "name", "sc", "st"])
ann["cpg_id"] = ann.name.astype(str).str.split("_").str[0].str.lower(); ann["pos"] = ann.s + 1; ann = ann[ann.cpg_id.str.match(r"^c[gh]\d")].drop_duplicates("cpg_id")
# CoRSIVs 2.0: probe-level tables S15 (EPIC v1) and S16 (EPIC v2)
c1 = pd.read_excel(CX, sheet_name="S15"); c2 = pd.read_excel(CX, sheet_name="S16")
cor = pd.concat([c1[["EPIC_probe_ID", "Region_type", "CoRSIV_ID"]], c2[["EPIC_probe_ID", "Region_type", "CoRSIV_ID"]]]); cor["cpg_id"] = cor.EPIC_probe_ID.astype(str).str.split("_").str[0]
print("CoRSIV 2.0 probe rows: EPIC v1", len(c1), "EPIC v2", len(c2), "| unique CpGs", cor.cpg_id.nunique(), "| by screen", cor.drop_duplicates(["cpg_id", "Region_type"]).Region_type.value_counts().to_dict())
in_corsiv = set(cor.cpg_id)
# Rosenski 2025: regions (hg19)
def regions(sheet, col, skip=2):
    d = pd.read_excel(RX, sheet_name=sheet, skiprows=skip); out = []
    for v in d[col].dropna().astype(str):
        for r in v.split(","):
            m = re.match(r"(chr[\w]+):(\d+)-(\d+)", r.strip().replace(",", ""))
            if m: out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return pd.DataFrame(sorted(set(out)), columns=["chr", "s", "e"])
sh = pd.ExcelFile(RX).sheet_names
icr = pd.read_excel(RX, sheet_name=sh[0], skiprows=2).dropna(subset=["chrom"]); icr = icr.rename(columns={"chrom": "chr", "start": "s", "end": "e"})[["chr", "s", "e"]].astype({"s": int, "e": int})
asm = regions(sh[1], "ASM bimodal regions (hg19)"); par = regions(sh[2], "Parental ASM region (hg19)")
print("regions: known iDMRs", len(icr), "| parental-ASM", len(par), "| ASM (sequence-dependent and other)", len(asm), "| median width", int((asm.e - asm.s).median()))
def inside(R):
    hit = set()
    for ch, d in R.groupby("chr"):
        a = ann[ann.chr == ch]; p = a.pos.values; o = np.argsort(p); p = p[o]; ids = a.cpg_id.values[o]
        for s, e in zip(d.s.values, d.e.values):
            i, j = np.searchsorted(p, s), np.searchsorted(p, e, side="right"); hit.update(ids[i:j])
    return hit
in_imp = inside(pd.concat([icr, par])); in_asm = inside(asm)
print("array CpGs flagged: in_corsiv", len(in_corsiv), "| in imprinted / parental-ASM region", len(in_imp), "| in ASM region", len(in_asm), "| ASM not imprinted", len(in_asm - in_imp))
F = ann[["cpg_id"]].copy(); F["in_corsiv"] = F.cpg_id.isin(in_corsiv); F["in_imprinted_or_parental_asm"] = F.cpg_id.isin(in_imp); F["in_asm"] = F.cpg_id.isin(in_asm)
extra = pd.DataFrame({"cpg_id": sorted(in_corsiv - set(F.cpg_id))}); extra["in_corsiv"], extra["in_imprinted_or_parental_asm"], extra["in_asm"] = True, False, False
F = pd.concat([F, extra]); F[F.in_corsiv | F.in_imprinted_or_parental_asm | F.in_asm].to_csv(OUT + "_cpg_region_flags.csv.gz", index=False)
flags = {"CoRSIV 2.0": in_corsiv, "imprinted / parental-ASM": in_imp, "ASM (any)": in_asm}
# test (b): sign agreement inside vs outside
M = pd.read_csv("/home/claude/v2r/smr/inst/extdata/measured_eqtms.csv.gz"); M["target_gene"] = M.target_gene.str.upper().str.strip(); M = M[M.direction.isin([-1, 1])]
a = pd.read_csv("/mnt/user-data/uploads/JASPAR Validation/GITHUB_cpgdirection/tools/out/measured_eqtms_atlas_long.csv.gz", usecols=["cohort", "cpg_id", "target_gene", "direction"]); a["target_gene"] = a.target_gene.astype(str).str.upper()
grp = lambda c: "EVA-PR" if c.startswith("EVA-PR") else "GTEx" if c.startswith("GTEx") else c
a["group"] = a.cohort.map(grp); g = a.groupby(["cpg_id", "target_gene", "group"]).direction.sum().reset_index(); g = g[g.direction != 0]; g["direction"] = np.sign(g.direction)
u = g.groupby(["cpg_id", "target_gene"]).direction.agg(n="size", s="sum").reset_index(); u = u[u.n >= 2]; u["unanimous"] = u.s.abs() == u.n
S = pd.read_csv(SMR, usecols=["cpg_id", "target_gene", "direction", "smr_callable", "smr_rung"], low_memory=False); S["target_gene"] = S.target_gene.astype(str).str.upper(); S = S[S.smr_callable.astype(str) == "True"]
Mb = M[M.tissue == "blood"].drop_duplicates(["cpg_id", "target_gene"])[["cpg_id", "target_gene", "direction"]]
x = S.merge(Mb, on=["cpg_id", "target_gene"], suffixes=("", "_m")); x["ok"] = x.direction == x.direction_m
B = pd.read_csv("/home/claude/bp/blueprint_pairs_fdr05.csv.gz"); B = B[B.adjustment == "B"].drop_duplicates(["cell", "cpg_id", "target_gene"])
y = B.merge(Mb, on=["cpg_id", "target_gene"], suffixes=("", "_m")); y["ok"] = y.direction == y.direction_m
z = S.merge(B.groupby(["cpg_id", "target_gene"]).direction.first().reset_index(), on=["cpg_id", "target_gene"], suffixes=("", "_b")); z["ok"] = z.direction == z.direction_b
rows = []
def cmp(lab, d, col):
    for fn, fs in flags.items():
        i, o = d[d.cpg_id.isin(fs)], d[~d.cpg_id.isin(fs)]
        if len(i) >= 20:
            ki, ko = int(i[col].sum()), int(o[col].sum()); print(f"  {lab:54s} {fn:26s} inside {ki/len(i):.3f} {wil(ki, len(i))} n {len(i):5d} | outside {ko/len(o):.3f} n {len(o)}")
            rows.append(dict(check=lab, flag=fn, inside=ki / len(i), n_inside=len(i), lo=wil(ki, len(i))[0], hi=wil(ki, len(i))[1], outside=ko / len(o), n_outside=len(o)))
print("\nsign agreement inside vs outside the flagged regions:")
cmp("measured pairs seen in 2+ donor groups: unanimous", u, "unanimous"); cmp("SMR callable vs measured blood: same sign", x, "ok"); cmp("BLUEPRINT vs measured blood: same sign", y, "ok"); cmp("SMR callable vs BLUEPRINT: same sign", z, "ok")
pd.DataFrame(rows).to_csv(OUT + "_flag_tests.csv", index=False)
print("\nshare of pairs carrying a flag: measured table", {k: round(float(M.cpg_id.isin(v).mean()), 4) for k, v in flags.items()}, "| SMR callable", {k: round(float(S.cpg_id.isin(v).mean()), 4) for k, v in flags.items()})
# Voisin 2016 adipose
v = pd.read_excel(VX); v["direction"] = np.sign(v.Beta).astype(int); v["target_gene"] = v["Gene Symbol"].astype(str).str.upper().str.split(r"\s*///\s*"); v = v.explode("target_gene").rename(columns={"TargetID": "cpg_id"})
v = v[v.target_gene.ne("NAN") & v.target_gene.ne("---")].drop_duplicates(["cpg_id", "target_gene"])
print(f"\nVoisin 2016 adipose: pairs {len(v)} | CpGs {v.cpg_id.nunique()} | genes {v.target_gene.nunique()} | negative {float((v.direction == -1).mean()):.3f}")
anyM = M.groupby(["cpg_id", "target_gene"]).direction.sum().reset_index(); anyM = anyM[anyM.direction != 0]; anyM["direction"] = np.sign(anyM.direction)
N = pd.read_csv("/home/claude/ln/admission_pairs.csv.gz")
for rn, R in [("all shared measured pairs (any tissue)", anyM), ("measured solid", M[M.tissue == "solid_tissue"]), ("measured blood", M[M.tissue == "blood"]), ("measured solid, 2+ cohorts", M[(M.tissue == "solid_tissue") & (M.n_cohorts >= 2)]),
              ("kidney (Liu 2022)", N[N.source == "liu2022_kidney"]), ("islets (Olsson 2014, same lab)", N[N.source == "olsson2014_islets"]), ("SMR callable", S)]:
    q = v.merge(R[["cpg_id", "target_gene", "direction"]].drop_duplicates(["cpg_id", "target_gene"]), on=["cpg_id", "target_gene"], suffixes=("", "_r")); k = int((q.direction == q.direction_r).sum())
    if len(q): print(f"  {rn:40s} shared {len(q):4d} same sign {k/len(q):.3f} {wil(k, len(q))}")
v[["cpg_id", "target_gene", "direction", "Beta", "P-value", "FDR"]].to_csv(OUT + "_voisin2016_adipose_pairs.csv", index=False)
