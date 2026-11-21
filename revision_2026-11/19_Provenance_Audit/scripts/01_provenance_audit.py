#!/usr/bin/env python3
"""Stage 19 — sequence provenance as a confounder.

The dataset has two halves: 237 single-contig GenBank assemblies and 53
genomes assembled de novo here from SRA reads. De novo assembly and Prokka
annotation of a fragmented draft can create ORFs that a curated single-contig
submission would never carry, and those land exclusively in the low-frequency
compartment. This stage quantifies that, so the manuscript can state which
results are robust to provenance and which are not.

Three tests:
  A  per-orthogroup presence in SRA-derived vs GenBank genomes
  B  Heaps' law refitted on GenBank genomes alone (is the pan-genome still
     "open" by the manuscript's own gamma > 0.05 criterion?)
  C  provenance of the robust pan-GWAS hits from stage 11, and whether the
     geographic traits themselves are confounded with provenance

Outputs
  02_results/orthogroup_provenance.csv
  02_results/heaps_fit_comparison.csv
  02_results/pangwas_hit_provenance.csv
  03_summary/provenance_audit.log
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (ANNOT_ALL, TRAITS, ROOT, load_pa, stage_dir,  # noqa: E402
                     is_sra)
import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402
from scipy.stats import fisher_exact, mannwhitneyu                 # noqa: E402

N_PERM = int(os.environ.get("HEAPS_PERM", "100"))
SEED = 20260817                      # same seed as the original Heaps fit
D = stage_dir("19_Provenance_Audit")
LOG = open(os.path.join(D, "03_summary", "provenance_audit.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


PA = load_pa()
ann = pd.read_csv(ANNOT_ALL).set_index("Orthogroup")["product"]
sra = [c for c in PA.columns if is_sra(c)]
gb = [c for c in PA.columns if not is_sra(c)]
say(f"{len(sra)} SRA-derived genomes, {len(gb)} GenBank genomes\n")

# ---------------------------------------------------------------- A ---------
freq = PA.sum(axis=1)
cls = pd.cut(freq / PA.shape[1], [-.01, .15, .95, .99, 1.01],
             labels=["cloud", "shell", "soft_core", "core"])
n_sra = PA[sra].sum(axis=1)
n_gb = PA[gb].sum(axis=1)
tab = pd.DataFrame({"pangenome_class": cls, "n_total": freq,
                    "n_sra": n_sra, "n_genbank": n_gb,
                    "pct_sra": (100 * n_sra / len(sra)).round(1),
                    "pct_genbank": (100 * n_gb / len(gb)).round(1)})
tab["fisher_p"] = [fisher_exact([[int(a), len(sra) - int(a)],
                                 [int(b), len(gb) - int(b)]])[1]
                   for a, b in zip(n_sra, n_gb)]
tab["sra_exclusive"] = (tab.n_genbank == 0)
tab["product"] = [str(ann.get(o, ""))[:80] for o in tab.index]
tab.to_csv(os.path.join(D, "02_results", "orthogroup_provenance.csv"))

say("=== A. orthogroups absent from every GenBank genome, by class ===")
say(tab[tab.sra_exclusive].pangenome_class.value_counts().to_string())
for c in ["cloud", "shell", "soft_core", "core"]:
    s = tab[tab.pangenome_class == c]
    if not len(s):
        continue
    load_s = s.n_sra.sum() / len(sra)
    load_g = s.n_genbank.sum() / len(gb)
    u = mannwhitneyu(PA.loc[s.index, sra].sum(axis=0),
                     PA.loc[s.index, gb].sum(axis=0)).pvalue
    say(f"  {c:<10} n={len(s):<4} per-genome load  SRA {load_s:5.2f}  "
        f"GenBank {load_g:5.2f}  ratio "
        f"{(load_s/load_g if load_g else float('inf')):5.1f}x  MWU p={u:.2e}")

# ---------------------------------------------------------------- B ---------
def heaps(mat, permutations=N_PERM, seed=SEED):
    arr = mat.values > 0
    n_og, n_gen = arr.shape
    rng = np.random.default_rng(seed)
    pan = np.zeros((permutations, n_gen), dtype=int)
    for p in range(permutations):
        order = rng.permutation(n_gen)
        seen = np.zeros(n_og, dtype=bool)
        for i, g in enumerate(order):
            seen |= arr[:, g]
            pan[p, i] = seen.sum()
    N = np.arange(1, n_gen + 1)
    m = N >= 3
    slope, intercept = np.polyfit(np.log(N[m]), np.log(pan.mean(0)[m]), 1)
    obs, pred = pan.mean(0)[m], np.exp(intercept) * N[m] ** slope
    r2 = 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()
    return float(np.exp(intercept)), float(slope), float(r2)

say("\n=== B. Heaps' law, by genome provenance ===")
rows = []
for tag, cols in [("all290", list(PA.columns)), ("genbank237", gb),
                  ("sra53", sra)]:
    sub = PA[cols]
    sub = sub.loc[sub.sum(axis=1) > 0]
    k, gam, r2 = heaps(sub)
    verdict = "open" if gam > 0.05 else "closed"
    rows.append(dict(dataset=tag, n_genomes=len(cols),
                     n_orthogroups=sub.shape[0], k=round(k, 3),
                     gamma=round(gam, 4), r_squared=round(r2, 4),
                     verdict=verdict,
                     strict_core=int((sub.sum(axis=1) == len(cols)).sum())))
    say(f"  {tag:<12} n={len(cols):<4} pan={sub.shape[0]:<4} "
        f"gamma={gam:.4f} (R2={r2:.4f}) -> {verdict}")
pd.DataFrame(rows).to_csv(
    os.path.join(D, "02_results", "heaps_fit_comparison.csv"), index=False)

# ---------------------------------------------------------------- C ---------
say("\n=== C. provenance of the robust pan-GWAS hits (stage 11) ===")
hit_file = f"{ROOT}/11_Scoary/02_summary/scoary_tier1_robust.csv"
if os.path.exists(hit_file):
    hits = pd.read_csv(hit_file)
    out = []
    for _, h in hits.iterrows():
        og = h["Gene"]
        if og not in PA.index:
            continue
        a, b = int(n_sra[og]), int(n_gb[og])
        p = fisher_exact([[a, len(sra) - a], [b, len(gb) - b]])[1]
        out.append(dict(Orthogroup=og, trait=h.get("Trait", ""),
                        n_sra=a, n_genbank=b,
                        pct_sra=round(100 * a / len(sra), 1),
                        pct_genbank=round(100 * b / len(gb), 1),
                        fisher_p=p, sra_exclusive=(b == 0),
                        product=str(ann.get(og, ""))[:60]))
    hp = pd.DataFrame(out)
    hp.to_csv(os.path.join(D, "02_results", "pangwas_hit_provenance.csv"),
              index=False)
    say(hp[["Orthogroup", "trait", "n_sra", "n_genbank", "pct_sra",
            "pct_genbank", "fisher_p", "sra_exclusive"]].to_string(index=False))
else:
    say(f"  skipped: {hit_file} not found")

say("\n=== C2. are the traits themselves confounded with provenance? ===")
tr = pd.read_csv(TRAITS, index_col=0)
s_idx = [i for i in tr.index if is_sra(i)]
g_idx = [i for i in tr.index if not is_sra(i)]
rows = []
for t in [c for c in ["Africa", "Asia", "Europe", "Vaccine", "Recombinant"]
          if c in tr.columns]:
    a, b = int(tr.loc[s_idx, t].sum()), int(tr.loc[g_idx, t].sum())
    p = fisher_exact([[a, len(s_idx) - a], [b, len(g_idx) - b]])[1]
    rows.append(dict(trait=t, n_sra=a, n_genbank=b,
                     pct_sra=round(100 * a / len(s_idx), 1),
                     pct_genbank=round(100 * b / len(g_idx), 1), fisher_p=p))
    say(f"  {t:<12} SRA {a:>3}/{len(s_idx)} ({100*a/len(s_idx):4.1f}%)  "
        f"GenBank {b:>3}/{len(g_idx)} ({100*b/len(g_idx):4.1f}%)  p={p:.2e}")
pd.DataFrame(rows).to_csv(
    os.path.join(D, "02_results", "trait_provenance.csv"), index=False)

say("\nwritten: 02_results/*.csv")
LOG.close()
