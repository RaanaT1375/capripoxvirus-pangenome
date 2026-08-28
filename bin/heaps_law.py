#!/usr/bin/env python3
"""Pangenome accumulation curve and Heaps' law fit (n = k * N^gamma).

gamma > 0 formally defines an open pangenome, but its magnitude matters more
than its sign: a small exponent means the accessory repertoire saturates
quickly and new genomes contribute rare elements rather than novel functions.
"""
import argparse
import numpy as np, pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="binary presence/absence CSV")
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--out-curve", required=True)
    ap.add_argument("--out-fit", required=True)
    a = ap.parse_args()

    m = pd.read_csv(a.matrix, index_col=0)
    arr = (m.values > 0)
    n_og, n_gen = arr.shape
    rng = np.random.default_rng(a.seed)

    pan  = np.zeros((a.permutations, n_gen), dtype=int)
    core = np.zeros((a.permutations, n_gen), dtype=int)
    for p in range(a.permutations):
        order = rng.permutation(n_gen)
        seen  = np.zeros(n_og, dtype=bool)
        allg  = np.ones(n_og, dtype=bool)
        for i, g in enumerate(order):
            seen |= arr[:, g]
            allg &= arr[:, g]
            pan[p, i]  = seen.sum()
            core[p, i] = allg.sum()

    N = np.arange(1, n_gen + 1)
    curve = pd.DataFrame({
        "n_genomes": N,
        "pan_mean":  pan.mean(0),  "pan_sd":  pan.std(0),
        "core_mean": core.mean(0), "core_sd": core.std(0),
    })
    curve.to_csv(a.out_curve, sep="\t", index=False)

    # fit on N >= 3, where the log-log relationship is defined and stable
    mask = N >= 3
    slope, intercept = np.polyfit(np.log(N[mask]), np.log(pan.mean(0)[mask]), 1)
    k = float(np.exp(intercept))
    pred = k * N[mask] ** slope
    obs  = pan.mean(0)[mask]
    r2 = 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()

    pd.DataFrame([{"k": round(k, 4), "gamma": round(float(slope), 4),
                   "r_squared": round(float(r2), 4),
                   "verdict": "open" if slope > 0 else "closed",
                   "permutations": a.permutations}]
                 ).to_csv(a.out_fit, sep="\t", index=False)
    print(f"Heaps' law: n = {k:.3f} * N^{slope:.4f}   (R2 = {r2:.4f})")

if __name__ == "__main__":
    main()
