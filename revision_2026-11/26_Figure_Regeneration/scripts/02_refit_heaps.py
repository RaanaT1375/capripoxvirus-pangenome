#!/usr/bin/env python3
"""Stage 26, step 2 — refit Heaps' law on the 178-orthogroup matrix.

Rebuilds the pan-genome and core accumulation curves after excluding the
orthogroups listed in 01_inputs/nonviral_orthogroups.txt, and fits
n = k * N^gamma. This supplies the k that Figure 1B's caption needs.

The fit is reported two ways:
  log-log OLS   log n = log k + gamma log N, ordinary least squares
  nonlinear     direct least squares on n = k N^gamma (scipy, if available)
The log-log fit is the conventional one and is what the manuscript quotes;
the nonlinear fit is printed as a check that the two agree.

Both are fitted to the MEAN curve across permutations, which is what
reproduces the submitted value: the previous fit, k = 140.712 and
gamma = 0.0679 on 206 orthogroups, gives 140.712 * 290^0.0679 = 206.6 at
N = 290. The script re-reports that identity for the new fit so you can see
the curve passes through the observed pan-genome size.

Outputs
  02_results/accumulation_curves.csv      per-N mean/sd for pan and core
  03_summary/heaps_fit.csv                k, gamma, R^2 for both fits
  03_summary/02_refit_heaps.log
"""
import os
import sys

def _locate_common():
    env = os.environ.get("COMMON_DIR")
    if env and os.path.exists(os.path.join(env, "_common.py")):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.exists(os.path.join(d, "_common.py")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    sys.exit(
        "could not find _common.py in any parent of\n  "
        + os.path.abspath(__file__)
        + "\nIt belongs in the revision_2026-11 directory alongside the "
          "numbered stages.\nIf it lives somewhere else, point at it:\n"
          "  export COMMON_DIR=/path/to/the/directory/holding/_common.py")


sys.path.insert(0, _locate_common())
from _common import stage_dir, load_pa                            # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402

PERMS = int(os.environ.get("PERMS", "100"))
SEED = int(os.environ.get("SEED", "20260817"))

D = stage_dir("26_Figure_Regeneration")
NONVIRAL = os.path.join(D, "01_inputs", "nonviral_orthogroups.txt")
LOG = open(os.path.join(D, "03_summary", "02_refit_heaps.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


if not os.path.exists(NONVIRAL):
    sys.exit(f"{NONVIRAL} not found; run 01_nonviral_screen.py first")

drop = {l.strip() for l in open(NONVIRAL) if l.strip()}
pa = load_pa()
say(f"presence/absence matrix: {pa.shape[0]} orthogroups x {pa.shape[1]} genomes")
say(f"excluding {len(drop)} non-viral orthogroups")
pa = pa.loc[[i for i in pa.index if i not in drop]]
say(f"retained: {pa.shape[0]} orthogroups")

M = pa.values.astype(bool)          # orthogroup x genome
n_og, n_gen = M.shape

# ---- partition, for the Figure 1D counts -----------------------------------
freq = M.sum(axis=1) / n_gen
part = pd.Series(index=pa.index, dtype=object)
part[freq >= 0.99] = "core"
part[(freq >= 0.95) & (freq < 0.99)] = "soft-core"
part[(freq >= 0.15) & (freq < 0.95)] = "shell"
part[freq < 0.15] = "cloud"
counts = part.value_counts().reindex(["core", "soft-core", "shell", "cloud"])
say("\npartition of the retained matrix:")
for k_, v in counts.items():
    say(f"  {k_:<10} {int(v):>4}  ({100*v/n_og:.1f}%)")
say(f"  {'total':<10} {n_og:>4}")
say("\nmanuscript states: core 132, soft-core 10, shell 14, cloud 22, "
    "total 178")
if list(counts.values) != [132, 10, 14, 22]:
    say("  *** MISMATCH: the partition above differs from the manuscript. ***")
    say("  Report this before regenerating the figures.")

# ---- accumulation curves ----------------------------------------------------
rng = np.random.default_rng(SEED)
pan = np.zeros((PERMS, n_gen), dtype=int)
core = np.zeros((PERMS, n_gen), dtype=int)
for p in range(PERMS):
    order = rng.permutation(n_gen)
    seen = np.zeros(n_og, dtype=bool)
    allg = np.ones(n_og, dtype=bool)
    for i, g in enumerate(order):
        col = M[:, g]
        seen |= col
        allg &= col
        pan[p, i] = seen.sum()
        core[p, i] = allg.sum()

N = np.arange(1, n_gen + 1)
curves = pd.DataFrame({
    "N": N,
    "pan_mean": pan.mean(axis=0), "pan_sd": pan.std(axis=0),
    "core_mean": core.mean(axis=0), "core_sd": core.std(axis=0),
})
curves.to_csv(os.path.join(D, "02_results", "accumulation_curves.csv"),
              index=False)
say(f"\naccumulation curves over {PERMS} permutations written")
say(f"  pan-genome at N={n_gen}: {curves.pan_mean.iloc[-1]:.1f}")
say(f"  strict core at N={n_gen}: {curves.core_mean.iloc[-1]:.1f}")

# ---- Heaps' law -------------------------------------------------------------
x, y = N.astype(float), curves.pan_mean.values

gam, logk = np.polyfit(np.log(x), np.log(y), 1)
k_ll = float(np.exp(logk))
pred = k_ll * x ** gam
ss = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

rows = [dict(fit="log-log OLS", k=k_ll, gamma=float(gam), R2=float(ss),
             predicted_at_Nmax=float(k_ll * n_gen ** gam))]
say("\nHeaps' law, log-log OLS:")
say(f"  k      = {k_ll:.3f}")
say(f"  gamma  = {gam:.4f}")
say(f"  R^2    = {ss:.4f}")
say(f"  k * {n_gen}^gamma = {k_ll * n_gen ** gam:.1f} "
    f"(observed {y[-1]:.1f})")

try:
    from scipy.optimize import curve_fit
    popt, _ = curve_fit(lambda n, k, g: k * n ** g, x, y,
                        p0=[y[0], 0.05], maxfev=20000)
    k_nl, g_nl = float(popt[0]), float(popt[1])
    pred2 = k_nl * x ** g_nl
    ss2 = 1 - ((y - pred2) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rows.append(dict(fit="nonlinear LS", k=k_nl, gamma=g_nl, R2=float(ss2),
                     predicted_at_Nmax=float(k_nl * n_gen ** g_nl)))
    say("\nHeaps' law, nonlinear least squares:")
    say(f"  k      = {k_nl:.3f}")
    say(f"  gamma  = {g_nl:.4f}")
    say(f"  R^2    = {ss2:.4f}")
except ImportError:
    say("\n(scipy not available; nonlinear cross-check skipped)")

pd.DataFrame(rows).to_csv(
    os.path.join(D, "03_summary", "heaps_fit.csv"), index=False)

say("\nmanuscript states gamma = 0.0351 on this matrix.")
say(f"this run gives gamma = {gam:.4f} by the conventional log-log fit.")
say("\nFor the Figure 1B caption, insert:")
say(f"  ($n = {k_ll:.3f} \\times N^{{{gam:.4f}}}$)")
say("\nwritten: 03_summary/heaps_fit.csv")
LOG.close()
