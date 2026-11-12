#!/usr/bin/env python3
"""One INTERVAL record per CpG-gene pair for the SMR layer (input of tools/build_smr_layer_interval.py).

GoDMC cis-mQTL instruments (exposure) against the INTERVAL whole-blood RNA-seq cis-eQTL summary statistics
(outcome), SMR 1.03 with HEIDI, run per chromosome as for FinnGen (tools/finngen_smr_celltype.R). Two records
exist per pair and are merged here:
  top-SNP SMR   the SMR test at the top cis-mQTL, pairs at p < 5e-8 (<interval_full_smr_sig.tsv.gz>)
  vote          the sign vote of the pair's clumped genome-wide-significant instruments (<interval_full_directions.csv.gz>)
`evidence` says which exist: both (same sign), smr_top_snp, instrument_vote, conflict (opposite signs).
The pair's sign is the top-SNP SMR sign where that test exists, else the vote; against the sign measured in blood
the top-SNP sign of a conflict is right 0.77 of the time and the vote 0.23, so a conflict keeps the top-SNP sign
at tier S3.

Usage: python collect_interval_smr.py <interval_full_directions.csv.gz> <interval_full_smr_sig.tsv.gz>
         <phenotype_summary.tsv.gz> <out.csv.gz> --ens <ens_grch37.csv.gz> --gene-table <gene_table.csv.gz>
         --measured <inst/extdata/measured_eqtms.csv.gz> --chr <csv.gz with cpg_id, chr> [<more chr files> ...]
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
for x in ["dirf", "smrf", "ph", "out"]: ap.add_argument(x)
for x in ["ens", "gene-table", "measured"]: ap.add_argument("--" + x, required=True)
ap.add_argument("--chr", nargs="+", required=True)
a = ap.parse_args()
DIRF, SMRF, PH, OUT = a.dirf, a.smrf, a.ph, a.out
def wil(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(c - h, 3), round(c + h, 3)
I = pd.read_csv(DIRF)
T = pd.read_csv(SMRF, sep="\t")
print("instrument-vote pairs", len(I), "| SMR top-SNP pairs at p < 5e-8", len(T), "| columns", T.columns.tolist())
T = T[T.b_SMR.notna() & (T.b_SMR != 0)].copy()
T["d_top"] = np.sign(T.b_SMR).astype(int)
T = T.rename(columns={"Expo_ID": "cpg_id", "Outco_ID": "ensembl_id", "Outco_Gene": "gene_run", "topSNP": "top_snp", "p_SMR": "p_top", "p_HEIDI": "p_heidi_top", "nsnp_HEIDI": "nsnp_top"})
T["dist_top"] = T.Outco_bp - T.Expo_bp
T = T.sort_values("p_top").drop_duplicates(["cpg_id", "ensembl_id"])
I = I.rename(columns={"direction": "d_vote", "smr_tier": "tier_vote", "p_SMR": "p_vote_top", "target_gene": "gene_vote"}).drop_duplicates(["cpg_id", "ensembl_id"])
M = T[["cpg_id", "ensembl_id", "gene_run", "d_top", "p_top", "p_heidi_top", "nsnp_top", "top_snp", "dist_top", "Expo_bp", "Outco_bp"]].merge(
    I[["cpg_id", "ensembl_id", "gene_vote", "d_vote", "tier_vote", "n_instruments", "instrument_agreement", "cpg_gene_dist", "p_vote_top", "p_HEIDI", "nsnp_HEIDI", "top_instrument"]],
    on=["cpg_id", "ensembl_id"], how="outer")
M["evidence"] = np.where(M.d_top.notna() & M.d_vote.notna(), np.where(M.d_top == M.d_vote, "both", "conflict"), np.where(M.d_top.notna(), "smr_top_snp", "instrument_vote"))
print("pairs", len(M), "| evidence", M.evidence.value_counts().to_dict())
ph = pd.read_csv(PH, sep="\t").drop_duplicates("phenotype_id").set_index("phenotype_id")
E = pd.read_csv(a.ens).drop_duplicates("ensgene").set_index("ensgene")
G = pd.read_csv(a.gene_table).set_index("gene_id")
M["target_gene"] = M.ensembl_id.map(G.display).fillna(M.ensembl_id.map(E.symbol)).fillna(M.ensembl_id.map(ph.gene_name)).fillna(M.gene_run).fillna(M.gene_vote).astype(str).str.upper()
M["gene_source"] = np.where(M.ensembl_id.isin(G.index), "gene_table_0.99.7", np.where(M.ensembl_id.isin(E.index), "ensembl_grch37", "interval_name_only"))
# which sign is right when the top SNP and the clumped instruments disagree?
B = pd.read_csv(a.measured); B = B[(B.tissue == "blood") & B.direction.isin([-1, 1])]
B["target_gene"] = B.target_gene.str.upper().str.strip(); B = B.drop_duplicates(["cpg_id", "target_gene"])[["cpg_id", "target_gene", "direction"]].rename(columns={"direction": "d_meas"})
x = M.merge(B, on=["cpg_id", "target_gene"])
for ev in ["both", "smr_top_snp", "instrument_vote", "conflict"]:
    s = x[x.evidence == ev]; n = len(s)
    if ev == "conflict":
        k1 = int((s.d_top == s.d_meas).sum()); k2 = int((s.d_vote == s.d_meas).sum())
        print(f"  conflict: n {n} | top-SNP sign right {k1/n:.3f} {wil(k1, n)} | instrument vote right {k2/n:.3f} {wil(k2, n)} | by vote tier: " + str({t: (len(g), round(float((g.d_vote == g.d_meas).mean()), 3)) for t, g in s.groupby('tier_vote')}))
    else:
        d = s.d_top if ev != "instrument_vote" else s.d_vote
        k = int((d == s.d_meas).sum()); print(f"  {ev}: n {n} | right {k/n:.3f} {wil(k, n)}")
s = x[x.evidence == "both"]
for t, g in s.groupby("tier_vote"):
    print(f"    both, vote tier {t}: n {len(g)} right {float((g.d_top == g.d_meas).mean()):.3f}")
# record: the top-SNP SMR sign where the top-SNP test passes 5e-8 (as for FinnGen in build_smr_replication.R), else the instrument vote; a conflict takes the top-SNP sign (right 0.77 against 0.23) and tier S3
M["direction"] = np.where(M.d_top.notna(), M.d_top, M.d_vote)
M["smr_tier"] = np.where(M.evidence == "smr_top_snp", "S2", np.where(M.evidence == "conflict", "S3", M.tier_vote))
M["p_SMR"] = np.where(M.p_top.notna(), M.p_top, M.p_vote_top)
M["n_instruments"] = M.n_instruments.fillna(1)
M["instrument_agreement"] = M.instrument_agreement.fillna(1.0)
M["tss_dist"] = np.where(M.cpg_gene_dist.notna(), M.cpg_gene_dist, M.dist_top)
def chrs(f):
    d = pd.read_csv(f, low_memory=False); d = d.rename(columns={"cpg": "cpg_id", "cpg_chr": "chr"})
    return d[["cpg_id", "chr"]]
CH = pd.concat([chrs(f) for f in a.chr]).drop_duplicates("cpg_id").set_index("cpg_id").chr
M["chr"] = M.cpg_id.map(CH)
M["p_HEIDI_out"] = np.where(M.p_heidi_top.notna(), M.p_heidi_top, M.p_HEIDI)
M["Ensembl_ID"] = M.ensembl_id
out = M[["cpg_id", "target_gene", "Ensembl_ID", "gene_source", "direction", "smr_tier", "p_SMR", "n_instruments", "instrument_agreement", "tss_dist", "chr", "evidence", "d_top", "d_vote", "p_HEIDI_out", "top_snp", "top_instrument"]].copy()
out["direction"] = out.direction.astype(int)
out = out.sort_values("p_SMR").drop_duplicates(["cpg_id", "target_gene"])
print("final INTERVAL pairs", len(out), "| CpGs", out.cpg_id.nunique(), "| genes", out.Ensembl_ID.nunique(), "| tiers", out.smr_tier.value_counts().to_dict(), "| p_SMR missing", int(out.p_SMR.isna().sum()),
      "| p_SMR >= 5e-8 (vote only, top SNP weaker)", int((out.p_SMR >= 5e-8).sum()), "| no chr", int(out.chr.isna().sum()))
out.to_csv(OUT, index=False)
