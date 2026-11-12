"""BLUEPRINT cis eQTMs under different numbers of removed principal components. Saves r (float16) for every tested cis pair and configuration.
Configurations: K components removed from both data types for K in 0,1,2,3,5,10,15,20; expression only (10); methylation only (10).
Usage: python3 eqtm_sweep.py <cell> <index file> <out.npz>"""
import sys, csv, time
import numpy as np, pandas as pd
cell, idx, out = sys.argv[1], sys.argv[2], sys.argv[3]
t0 = time.time(); W = 250000; KMAX = 20
r = csv.reader(open(idx), delimiter="\t"); h = next(r); iD, iA, iS = h.index("DONOR_ID"), h.index("DONOR_AGE"), h.index("DONOR_SEX")
cov = {}
for row in r:
    if len(row) > max(iD, iA, iS):
        try: a = np.mean([float(x) for x in row[iA].split("-")])
        except ValueError: continue
        if row[iS] in ("Male", "Female"): cov[row[iD]] = (a, 1.0 if row[iS] == "Male" else 0.0)
gc = open(f"work/{cell}_gene.cols.txt").read().split(); mc = open(f"work/{cell}_meth.cols.txt").read().split()
don = [d for d in gc if d in set(mc) and d in cov]; n = len(don)
G = np.load(f"work/{cell}_gene.b00.npy")[:, [gc.index(d) for d in don]]; M = np.load(f"work/{cell}_meth.b00.npy")[:, [mc.index(d) for d in don]]
gid = pd.read_csv(f"work/{cell}_gene.b00.ids.txt", header=None)[0].str.split(".").str[0].values; cid = pd.read_csv(f"work/{cell}_meth.b00.ids.txt", header=None)[0].values
rg = pd.read_csv("ref_gene_grch37.tsv.gz", sep="\t").drop_duplicates("ensgene").set_index("ensgene"); rc = pd.read_csv("ref_cpg_hg19_plus450k.tsv.gz", sep="\t").set_index("cpg_id")
gmin = G.min(axis=1, keepdims=True); gk = (G.std(axis=1) > 0) & ((G == gmin).mean(axis=1) < 0.2) & pd.Index(gid).isin(rg.index); ck = (M.std(axis=1) > 0) & pd.Index(cid).isin(rc.index)
G, gid, M, cid = G[gk].astype(np.float64), gid[gk], M[ck], cid[ck]
age = np.array([cov[d][0] for d in don]); sex = np.array([cov[d][1] for d in don]); C = np.column_stack([np.ones(n), sex, age]); Q, _ = np.linalg.qr(C)
def z(X):
    X = X - X.mean(axis=1, keepdims=True); s = np.sqrt((X * X).sum(axis=1, keepdims=True)); s[s == 0] = 1; return X / s
GA = G - (G @ Q) @ Q.T; Gc = GA - GA.mean(axis=1, keepdims=True); wg, vg = np.linalg.eigh(Gc.T @ Gc); og = np.argsort(wg)[::-1][:KMAX]; pg = vg[:, og]; varg = wg[og] / wg.sum()
MA = np.empty_like(M, dtype=np.float32); S = np.zeros((n, n))
for i in range(0, len(M), 50000):
    b = M[i:i + 50000].astype(np.float64); b = b - (b @ Q) @ Q.T; MA[i:i + 50000] = b; bc = b - b.mean(axis=1, keepdims=True); S += bc.T @ bc
wm, vm = np.linalg.eigh(S); om = np.argsort(wm)[::-1][:KMAX]; pm = vm[:, om]; varm = wm[om] / wm.sum()
# raw-data components (before removing sex and age) to see what they track
Graw = G - G.mean(axis=1, keepdims=True); w0, v0 = np.linalg.eigh(Graw.T @ Graw); pg0 = v0[:, np.argsort(w0)[::-1][:5]]
cross = np.abs(np.corrcoef(np.column_stack([pm[:, :10], pg[:, :10]]).T)[:10, 10:])
print(f"{cell}: n {n} | variance per component, methylation {np.round(varm[:10], 3).tolist()} | expression {np.round(varg[:10], 3).tolist()}")
print("  |correlation| between methylation components (rows 1-5) and expression components (columns 1-5):"); print(np.round(cross[:5, :5], 2))
print("  largest |correlation| of each of the first 10 methylation components with any of the first 10 expression components:", np.round(cross.max(axis=1), 2).tolist())
print("  raw expression components 1-5 vs sex, age:", np.round([abs(np.corrcoef(pg0[:, k], sex)[0, 1]) for k in range(5)], 2).tolist(), np.round([abs(np.corrcoef(pg0[:, k], age)[0, 1]) for k in range(5)], 2).tolist(), flush=True)
configs = [(0, 0), (5, 5)]     # (methylation K, expression K)
gchr, gtss, gstr = rg.loc[gid, "chr"].values, rg.loc[gid, "tss"].values.astype(np.int64), rg.loc[gid, "strand"].values
cchr, cpos = rc.loc[cid, "chr"].values, rc.loc[cid, "pos"].values.astype(np.int64)
Zg = {}
for ke in sorted(set(c[1] for c in configs)):
    X = GA if ke == 0 else GA - (GA @ pg[:, :ke]) @ pg[:, :ke].T; Zg[ke] = z(X)
CI, GI, DI = [], [], []; R = {c: [] for c in configs}
for ch in sorted(set(gchr) & set(cchr)):
    ci = np.where(cchr == ch)[0]; ci = ci[np.argsort(cpos[ci])]; p = cpos[ci]; b = MA[ci].astype(np.float64)
    Zc = {km: z(b if km == 0 else b - (b @ pm[:, :km]) @ pm[:, :km].T) for km in sorted(set(c[0] for c in configs))}
    for gi in np.where(gchr == ch)[0]:
        lo, hi = np.searchsorted(p, gtss[gi] - W), np.searchsorted(p, gtss[gi] + W, side="right")
        if hi > lo:
            CI.append(ci[lo:hi]); GI.append(np.full(hi - lo, gi, np.int32)); d = p[lo:hi] - gtss[gi]; DI.append(d if gstr[gi] >= 0 else -d)
            for (km, ke) in configs: R[(km, ke)].append(Zc[km][lo:hi] @ Zg[ke][gi])
res = {f"r_m{km}_e{ke}": np.concatenate(R[(km, ke)]).astype(np.float32) for (km, ke) in configs}
np.savez(out, cpg=np.concatenate(CI).astype(np.int32), gene=np.concatenate(GI), dist=np.concatenate(DI).astype(np.int32), cpg_ids=cid, gene_ids=gid, n=n, cross=cross, varm=varm, varg=varg, **res)
print(f"  cis pairs {len(res['r_m0_e0'])} | configurations {len(configs)} | {time.time()-t0:.0f} s")
