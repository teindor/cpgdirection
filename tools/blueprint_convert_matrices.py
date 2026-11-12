"""Convert a BLUEPRINT matrix (txt.gz) to float32 .npy plus row ids, in row blocks so a call stays short.
Usage: python3 conv.py <in.txt.gz> <out_prefix> <block_index> <block_rows>"""
import sys, time, gzip, io
import numpy as np, pandas as pd
src, out, bi, br = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
t0 = time.time()
hdr = gzip.open(src, "rt").readline().rstrip("\n").split("\t")
df = pd.read_csv(src, sep="\t", header=None, skiprows=1 + bi * br, nrows=br, na_values=["NA", "NaN", ""], engine="c", low_memory=False)
if len(df) == 0: print("empty block"); sys.exit(0)
ids = df.iloc[:, 0].astype(str).values; X = df.iloc[:, 1:].to_numpy(dtype=np.float32)
np.save(f"{out}.b{bi:02d}.npy", X); pd.Series(ids).to_csv(f"{out}.b{bi:02d}.ids.txt", index=False, header=False)
if bi == 0: open(f"{out}.cols.txt", "w").write("\n".join(hdr[1:]) + "\n")
print(f"block {bi}: rows {len(df)} cols {X.shape[1]} NaN {int(np.isnan(X).sum())} | {time.time() - t0:.0f} s")
