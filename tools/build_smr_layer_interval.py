#!/usr/bin/env python3
"""SMR layer of cpgdirection 1.1.0: INTERVAL blood RNA-seq as a third SMR outcome, the rung rule ruled on
11-12 Nov 2026, and a held-out accuracy for every pair.

Inputs
  --smr        the current smr_directions table (GoDMC mQTL x eQTLGen, FinnGen PBMC and cell types; csv.gz)
  --interval   one INTERVAL record per CpG-gene pair (tools/collect_interval_smr.py)
  --measured   inst/extdata/measured_eqtms.csv.gz (its blood pairs are the validation signs)
  --positions  cpg_id, chr of the package position table (csv.gz; cpgd_cpg_positions())
  --hm450      Illumina 450K manifest v1.1 (GEO GPL13534), chromosome of the CpGs the package table lacks (optional)
  --helix-rows rows of the HELIX cell-adjusted full file (Dryad fxpnvx0t0, eQTM_autosome_adj.cells.txt.gz) for the
               CpGs of the pairs seen in INTERVAL alone: columns CpG, p.value, sigPair, CpG_pos, TC_gene_TSS, TC_gene
               (optional; fills `child_blood_helix` for those pairs by the rule of tools/ingest_helix_dryad.R)
  --out        output directory

Rung rule ("rule D", strength gate p < 1e-40)
  A pair is `replicated` when at least two bulk outcomes (eQTLGen, FinnGen PBMC, INTERVAL) give the same sign and
  none dissents, and either the instruments are concordant (tier S1) or one outcome reaches p < 1e-40; an absent
  outcome is neutral; a dissent demotes (two against one: weak; one against one: reported).
  INTERVAL admissibility (12 Nov): the INTERVAL sign is the top-SNP SMR sign where that test passes 5e-8, else the
  vote of the clumped instruments. Only a passing top-SNP record that does not conflict with the instrument vote
  can promote a pair; any INTERVAL dissent demotes; a pair seen in INTERVAL alone enters by its own tier when the
  record is admissible and as weak otherwise.

Per-pair accuracy
  Gradient-boosted classifier of "SMR sign equals the measured blood sign" on the pair's SMR evidence (strength
  and agreement of each outcome, tiers, instruments, HEIDI, distance, called sign), fitted leaving one chromosome
  out: every pair's figure comes from a model that saw no pair of its chromosome. `smr_accuracy` is that figure.
  The package groups it into bands (>= 0.95, 0.90-0.95, 0.80-0.90, 0.65-0.80, 0.55-0.65); pairs below 0.55 and
  pairs on the `reported` rung supply no direction.

Writes <out>/smr_directions_1.1.0.csv.gz and <out>/build_smr_layer_interval.md.
"""
import argparse
import datetime
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss

ap = argparse.ArgumentParser()
for a in ["smr", "interval", "measured", "positions", "out"]: ap.add_argument("--" + a, required=True)
ap.add_argument("--hm450"); ap.add_argument("--helix-rows", dest="helix"); ap.add_argument("--gate", type=float, default=40.0)
a = ap.parse_args()
RUNGS = ["replicated", "high", "moderate", "weak", "reported"]
BANDS = [(0.95, 1.01, ">= 0.95"), (0.90, 0.95, "0.90-0.95"), (0.80, 0.90, "0.80-0.90"), (0.65, 0.80, "0.65-0.80"), (0.55, 0.65, "0.55-0.65"), (0.50, 0.55, "0.50-0.55"), (0, 0.50, "< 0.50")]
L = []
def say(*x):
    for v in x: print(v, flush=True); L.append(str(v))
def wil(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h
def row(lab, m, has, okv):
    mm = m & has; n = int(mm.sum()); k = int(okv[mm].sum()); lo, hi = wil(k, n)
    return f"| {lab} | {int(m.sum()):,} | {n:,} | " + (f"{k / n:.3f} [{lo:.3f}, {hi:.3f}] |" if n else "- |")

S = pd.read_csv(a.smr, low_memory=False, keep_default_na=False, na_values=["", "NA"]); S["tg"] = S.target_gene.astype(str).str.upper()   # "null" is a value of child_blood_helix
assert not S.duplicated(["cpg_id", "tg"]).any()
I = pd.read_csv(a.interval, dtype={"chr": str}, low_memory=False)
I = I.rename(columns={"target_gene": "tg", "direction": "d_iv", "smr_tier": "tier_iv", "p_SMR": "p_iv", "evidence": "ev_iv", "n_instruments": "ni_iv", "instrument_agreement": "ia_iv",
                      "tss_dist": "dist_iv", "chr": "chr_iv", "p_HEIDI_out": "p_heidi_iv"})
I["tg"] = I.tg.astype(str).str.upper()
M = pd.read_csv(a.measured); M["tg"] = M.target_gene.str.upper().str.strip(); M = M[M.direction.isin([-1, 1])]
Mb = M[M.tissue == "blood"].drop_duplicates(["cpg_id", "tg"])[["cpg_id", "tg", "direction"]].rename(columns={"direction": "d_meas"})
other = M[M.tissue != "blood"].groupby(["cpg_id", "tg"]).direction.sum().reset_index().rename(columns={"direction": "d_other"})
other = other[other.d_other != 0]; other["d_other"] = np.sign(other.d_other)
A = S.merge(I[["cpg_id", "tg", "Ensembl_ID", "d_iv", "tier_iv", "p_iv", "ev_iv", "ni_iv", "ia_iv", "dist_iv", "chr_iv", "p_heidi_iv", "top_snp"]], on=["cpg_id", "tg"], how="outer", indicator=True)
A["in_layer"] = A._merge != "right_only"; new = (~A.in_layer).values
A.loc[new, "direction"] = A.loc[new, "d_iv"]; A.loc[new, "smr_tier"] = A.loc[new, "tier_iv"]; A.loc[new, "outcome"] = "INTERVAL"; A.loc[new, "replication"] = "interval_only"; A.loc[new, "smr_rung"] = "none"
A["iv"] = np.where(A.d_iv.isna(), "absent", np.where(A.d_iv == A.direction, "agrees", "disagrees")); A.loc[new, "iv"] = "interval_only"
A["src"] = np.where(A.outcome == "eQTLGen", "eQTLGen", np.where(A.outcome == "FinnGen_PBMC", "FinnGen_PBMC", np.where(A.outcome == "INTERVAL", "INTERVAL", "FinnGen_celltype")))
A = A.merge(Mb, on=["cpg_id", "tg"], how="left").merge(other, on=["cpg_id", "tg"], how="left").drop(columns="_merge").reset_index(drop=True)

# chromosome of each pair (held-out unit): package position table, then the 450K manifest, then the INTERVAL run
P = pd.read_csv(a.positions, usecols=["cpg_id", "chr"]).drop_duplicates("cpg_id").set_index("cpg_id").chr.astype(str).str.replace("chr", "")
A["chrom"] = A.cpg_id.map(P); n_pkg = int(A.chrom.notna().sum())
if a.hm450:
    H = pd.read_csv(a.hm450, skiprows=7, usecols=["Name", "Genome_Build", "CHR"], dtype=str, low_memory=False)
    H = H[H.Name.str.match(r"^c[gh]\d", na=False) & (H.Genome_Build == "37") & H.CHR.notna()].drop_duplicates("Name").set_index("Name").CHR
    A["chrom"] = A.chrom.fillna(A.cpg_id.map(H))
n_hm = int(A.chrom.notna().sum()) - n_pkg
A["chrom"] = A.chrom.fillna(A.chr_iv.astype(str).str.replace(r"\.0$", "", regex=True).where(A.chr_iv.notna()))
n_iv = int(A.chrom.notna().sum()) - n_pkg - n_hm
A["chrom"] = A.chrom.fillna("nan")
say(f"# SMR layer with INTERVAL, rule D (gate 1e-{a.gate:g}) and a held-out accuracy per pair", "")
say(f"Run {datetime.datetime.now():%Y-%m-%d %H:%M}.", "")
say(f"- pairs {len(A):,}: {int(A.in_layer.sum()):,} of the current layer, {int(new.sum()):,} seen in INTERVAL alone")
say(f"- measured blood sign available for {int(A.d_meas.notna().sum()):,} pairs (validation); measured outside blood only for {int((A.d_other.notna() & A.d_meas.isna()).sum()):,}")
say(f"- chromosome: package position table {n_pkg:,}; 450K manifest {n_hm:,}; INTERVAL run {n_iv:,}; none {int((A.chrom == 'nan').sum()):,}", "")

def lg(x): return -np.log10(pd.to_numeric(x, errors="coerce").clip(lower=1e-300)).fillna(0).values
eq = (A.src == "eQTLGen").values
pbo = ((A.src == "FinnGen_PBMC") & (A.smr_rung == "moderate")).values
pbx = ((A.src == "FinnGen_PBMC") & (A.smr_rung != "moderate")).values
ct = (A.src == "FinnGen_celltype").values
ivo = (A.src == "INTERVAL").values
rep = A.replication.astype(str).values
conc, disc = eq & (rep == "concordant"), eq & (rep == "discordant")
cconc, cdisc = eq & (rep == "concordant_celltype"), eq & (rep == "discordant_celltype")
has_iv = A.d_iv.notna().values
adm = A.ev_iv.isin(["both", "smr_top_snp"]).values
iva_any, ivd = (A.iv == "agrees").values, (A.iv == "disagrees").values
iva = iva_any & adm
tier = A.smr_tier.astype(str).values
s12 = np.isin(tier, ["S1", "S2"])
lp, lpi, lpf = lg(A.p_SMR), lg(A.p_iv), lg(A.finngen_p_SMR)
lp[new] = 0
strength = np.maximum(np.maximum(lp, np.where(adm, lpi, 0)), np.where(conc | disc, lpf, 0))
cur = A.smr_rung.astype(str).values
base = np.where(tier == "S1", "high", np.where(tier == "S2", "moderate", "weak"))
rung = cur.copy()
rung[ivo] = np.where(adm[ivo], base[ivo], "weak")
rung[(ct | pbx) & iva] = "moderate"
rung[(ct | pbx) & ivd] = "reported"
ok_gate = (tier == "S1") | (strength >= a.gate)
rung[eq & ~conc & ~disc & iva & s12 & ok_gate] = "replicated"
rung[pbo & iva & ok_gate] = "replicated"
rung[disc & iva] = "weak"
rung[conc & ivd] = "weak"
rung[eq & ~conc & ivd] = "reported"
rung[pbo & ivd] = "reported"
assert set(rung) <= set(RUNGS), set(rung)
dm = A.d_meas.values.astype(float); do = A.d_other.values.astype(float); sign = A.direction.values.astype(float)
has = ~np.isnan(dm); okv = sign == dm
say("## Rungs", "", "| rung | pairs | validated n | SMR sign = measured blood sign |", "|---|---|---|---|")
for r in RUNGS: say(row(r, rung == r, has, okv))
say(row("all but reported", np.isin(rung, RUNGS[:4]), has, okv), "")
say("Moves (rows: rung in the current layer; columns: rung now):", "", "```", pd.crosstab(np.where(ivo, "new (INTERVAL alone)", cur), rung).reindex(columns=RUNGS).fillna(0).astype(int).to_string(), "```", "")

c100 = lambda v: np.clip(v, 0, 100)
F = pd.DataFrame(index=A.index)
F["z_eq"] = np.where(eq, c100(lp), 0)
F["z_pb"] = np.where(conc, c100(lpf), np.where(disc, -c100(lpf), np.where(A.src.values == "FinnGen_PBMC", c100(lp), 0)))
F["z_ct"] = np.where(cconc, c100(lpf), np.where(cdisc, -c100(lpf), np.where(ct, c100(lp), 0)))
F["z_iv"] = np.where(iva_any | ivo, c100(lpi), np.where(ivd, -c100(lpi), 0))
F["has_eq"], F["has_pb"], F["has_ct"], F["has_iv"] = eq, (conc | disc | (A.src.values == "FinnGen_PBMC")), (cconc | cdisc | ct), has_iv
F["pb_subthreshold"] = pbx
F["S1"], F["S3"] = tier == "S1", tier == "S3"
F["iv_S1"] = (A.tier_iv.astype(str) == "S1").values
F["fg_S1"] = (A.finngen_tier.astype(str) == "S1").values
F["n_inst"] = np.where(new, 0, pd.to_numeric(A.n_instruments, errors="coerce").fillna(0).clip(upper=10).values)
F["inst_agree"] = pd.to_numeric(A.instrument_agreement, errors="coerce").fillna(1).values
F["ni_iv"] = A.ni_iv.fillna(0).clip(upper=10).values
F["heidi_pass"] = (A.heidi_status.astype(str) == "pass").values
F["heidi_fail"] = (A.heidi_status.astype(str) == "fail").values
F["lp_heidi"] = c100(lg(A.p_HEIDI))
F["n_ctx"] = A.finngen_contexts.astype(str).str.count(r"\|").where(A.finngen_contexts.notna(), -1).values + 1
dist = pd.to_numeric(A.cpg_gene_dist, errors="coerce").fillna(A.dist_iv)
F["absdist_log"] = np.log10(dist.abs().fillna(1e6).clip(lower=1) + 1).values
F["called_neg"] = sign < 0
F["iv_top"], F["iv_vote_only"], F["iv_conflict"] = adm, (A.ev_iv == "instrument_vote").values, (A.ev_iv == "conflict").values
F["iv_both"] = (A.ev_iv == "both").values
F["lp_heidi_iv"] = c100(lg(A.p_heidi_iv))
X = F.astype(float).values
chrom = A.chrom.values
hasc = has & (chrom != "nan")
y = okv.astype(int)
kw = dict(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=100, l2_regularization=1.0, random_state=1)
ext = ~np.isnan(do) & ~has
pred = np.full(len(A), np.nan); pext = np.zeros(len(A)); folds = sorted(set(chrom[chrom != "nan"]))
for c in folds:
    tr = hasc & (chrom != c)
    m = HistGradientBoostingClassifier(**kw).fit(X[tr], y[tr])
    te = chrom == c; pred[te] = m.predict_proba(X[te])[:, 1]; pext[ext] += m.predict_proba(X[ext])[:, 1] / len(folds)
rest = np.isnan(pred)
if rest.any():      # no chromosome known: such a pair is in no training set, so the all-chromosome fit has not seen it either
    pred[rest] = HistGradientBoostingClassifier(**kw).fit(X[hasc], y[hasc]).predict_proba(X[rest])[:, 1]
p = pred[hasc]; yy = y[hasc]
cal = pd.DataFrame({"q": pd.qcut(p, 10, duplicates="drop"), "p": p, "y": yy}).groupby("q", observed=True).agg(n=("y", "size"), predicted=("p", "mean"), observed=("y", "mean"))
yext = (sign[ext] == do[ext]).astype(int)
say("## Per-pair accuracy (leave one chromosome out)", "")
say(f"- validated pairs {int(hasc.sum()):,} over {len(folds)} chromosomes; AUC {roc_auc_score(yy, p):.3f}; Brier {brier_score_loss(yy, p):.4f}; mean calibration gap over deciles {float((cal.predicted - cal.observed).abs().mean()):.3f}")
say(f"- pairs measured outside blood only (n {int(ext.sum()):,}; never used for fitting): AUC {roc_auc_score(yext, pext[ext]):.3f}")
say(f"- pairs with no chromosome (figure from the all-chromosome fit, which contains none of them): {int(rest.sum()):,}", "")
say("Calibration by decile of the held-out figure:", "", "```", cal.round(3).to_string(), "```", "")
acc = np.round(pred, 3)
called = np.isin(rung, RUNGS[:4])
say("## Accuracy bands", "", "All pairs:", "", "| band | pairs | validated n | SMR sign = measured blood sign |", "|---|---|---|---|")
for lo, hi, lab in BANDS: say(row(lab, (acc >= lo) & (acc < hi), has, okv))
say("", "Pairs that can supply a direction (rung not `reported`):", "", "| band | pairs | validated n | SMR sign = measured blood sign |", "|---|---|---|---|")
for lo, hi, lab in BANDS: say(row(lab, called & (acc >= lo) & (acc < hi), has, okv))
say(row("supplying (>= 0.55, not reported)", called & (acc >= 0.55), has, okv))
say("", "Same bands on pairs measured outside blood only (not used for fitting):", "", "| band | pairs | validated n | SMR sign = measured sign |", "|---|---|---|---|")
hx = ~np.isnan(do) & ~has; okx = sign == do
for lo, hi, lab in BANDS: say(row(lab, called & (acc >= lo) & (acc < hi), hx, okx))
say("", "`reported`-rung pairs by band (excluded in every band):", "", "| band | pairs | validated n | SMR sign = measured blood sign |", "|---|---|---|---|")
for lo, hi, lab in BANDS: say(row(lab, ~called & (acc >= lo) & (acc < hi), has, okv))

keep = ["cpg_id", "target_gene", "direction", "smr_tier", "p_SMR", "n_instruments", "instrument_agreement", "cpg_gene_dist", "top_instrument", "p_HEIDI", "nsnp_HEIDI", "heidi_status", "outcome", "replication",
        "finngen_direction", "finngen_tier", "finngen_p_SMR", "finngen_heidi_status", "finngen_contexts", "finngen_evidence", "child_blood_helix"]
T = A[keep].copy()
T["target_gene"] = T.target_gene.fillna(A.tg)
for c, s in [("p_SMR", "p_iv"), ("n_instruments", "ni_iv"), ("instrument_agreement", "ia_iv"), ("cpg_gene_dist", "dist_iv"), ("top_instrument", "top_snp"), ("p_HEIDI", "p_heidi_iv")]:
    T.loc[new, c] = A.loc[new, s]
T["direction"] = sign.astype(int)
# HEIDI status of the pairs seen in INTERVAL alone, by the rule of the current layer (pass: p_HEIDI >= 0.05)
ph = pd.to_numeric(A.p_heidi_iv, errors="coerce")
T.loc[new, "heidi_status"] = np.where(ph.isna(), "not_tested", np.where(ph >= 0.05, "pass", "fail"))[new]
hx_note = "not looked up (no --helix-rows): left empty"
if a.helix:
    Hx = pd.read_csv(a.helix, sep="\t", usecols=["CpG", "p.value", "sigPair", "CpG_pos", "TC_gene_TSS", "TC_gene"], low_memory=False)
    Hx = Hx[Hx.TC_gene.notna() & ((Hx.TC_gene_TSS - Hx.CpG_pos).abs() <= 250000)].copy()
    Hx["tg"] = Hx.TC_gene.astype(str).str.strip().str.upper(); Hx["sig"] = Hx.sigPair.astype(str).str.upper() == "TRUE"
    Hg = Hx.groupby(["CpG", "tg"]).agg(p=("p.value", "min"), sig=("sig", "any")).reset_index().rename(columns={"CpG": "cpg_id"})
    Hg["cb"] = np.where(Hg.sig, "sig", np.where(Hg.p < 0.05, "nominal", "null"))
    cb = A[["cpg_id", "tg"]].merge(Hg[["cpg_id", "tg", "cb"]], on=["cpg_id", "tg"], how="left").cb.values
    covered = A.cpg_id.isin(set(Hx.CpG)).values
    old = T.child_blood_helix.astype(str).values
    chk = ~new & covered & np.isin(old, ["sig", "nominal", "null"])
    same = float((old[chk] == pd.Series(cb[chk]).fillna("untested").values).mean()) if chk.any() else float("nan")
    T.loc[new, "child_blood_helix"] = pd.Series(cb, index=A.index)[new].fillna("untested")
    vc = T.loc[new, "child_blood_helix"].value_counts().to_dict()
    hx_note = f"looked up in the HELIX full file: {vc}; rule check on {int(chk.sum()):,} pairs of the current layer that share those CpGs: same class in {same:.4f}"
T["interval_direction"] = A.d_iv.astype("Int64"); T["interval_tier"] = A.tier_iv; T["interval_p_SMR"] = A.p_iv; T["interval_evidence"] = A.ev_iv; T["interval_p_HEIDI"] = A.p_heidi_iv
T["interval_status"] = A.iv
T["smr_rung"] = rung; T["smr_accuracy"] = acc
assert not T.duplicated(["cpg_id", "target_gene"]).any() and T.smr_accuracy.notna().all() and T.direction.isin([-1, 1]).all()
T.to_csv(f"{a.out}/smr_directions_1.1.0.csv.gz", index=False)
say("", "## Child blood (HELIX) for the pairs seen in INTERVAL alone", "", hx_note)
say("", "## Output", "", f"`smr_directions_1.1.0.csv.gz`: {len(T):,} rows; rungs " + ", ".join(f"{k} {v:,}" for k, v in T.smr_rung.value_counts().items())
    + f"; supplying a direction (accuracy >= 0.55, rung not reported) {int((called & (acc >= 0.55)).sum()):,}.")
open(f"{a.out}/build_smr_layer_interval.md", "w").write("\n".join(L) + "\n")
