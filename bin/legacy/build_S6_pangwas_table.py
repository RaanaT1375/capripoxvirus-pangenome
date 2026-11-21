#!/usr/bin/env python3
"""
Build the per-trait pan-GWAS summary table (Supplementary File S6).

Reads the existing Scoary result files, the trait matrix and the gene
presence/absence matrix, and writes one row per trait -- including traits that
returned no significant association, so that negative results are reported.

Nothing is recomputed and no existing file is modified.

Run from:  /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
Output:    11_Scoary/03_summary/S6_pangwas_per_trait_summary.csv
           11_Scoary/03_summary/S6_pangwas_significant_pairs.csv
"""
import os, glob, sys
import numpy as np
import pandas as pd

BASE = "/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis"
TRAITS = "11_Scoary/01_inputs/traits.csv"
GPA = "11_Scoary/01_inputs/gene_presence_absence.csv"
OUTDIR = "11_Scoary/03_summary"

# Significance criterion used in the manuscript: BOTH filters must be satisfied.
# This conjunction reproduces exactly 92 gene-trait pairs / 16 traits / 18 orthogroups.
EMP_ALPHA = 0.05          # permutation empirical p
BON_ALPHA = 0.05          # per-trait Bonferroni-corrected p returned by Scoary
MIN_PAIRS = 3             # Tier 1 convergence criterion (with net positive support)

os.chdir(BASE)
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------------------------------------------- locate results
cands = [d for d in glob.glob("11_Scoary/**/", recursive=True)
         if "SUPERSEDED" not in d and glob.glob(os.path.join(d, "*.results.csv"))]
if not cands:
    cands = [d for d in glob.glob("11_Scoary/**/", recursive=True)
             if "SUPERSEDED" not in d and glob.glob(os.path.join(d, "*.csv"))
             and any("_" in os.path.basename(f) for f in glob.glob(os.path.join(d, "*.csv")))]
if not cands:
    sys.exit("ERROR: no Scoary result directory found under 11_Scoary/ "
             "(excluding _SUPERSEDED_*). Point RESULTS_DIR at it manually.")

RESULTS_DIR = sorted(cands, key=lambda d: -len(glob.glob(os.path.join(d, "*.csv"))))[0]
files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*.results.csv"))) \
        or sorted(glob.glob(os.path.join(RESULTS_DIR, "*.csv")))
print(f"Scoary results directory : {RESULTS_DIR}")
print(f"Per-trait result files    : {len(files)}")
if files:
    print(f"Columns in first file     : {list(pd.read_csv(files[0], nrows=1).columns)}\n")

# ----------------------------------------------------------------- helpers
def pick(cols, *cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    for c in cands:
        for lc, orig in low.items():
            if c.lower() in lc:
                return orig
    return None

# ----------------------------------------------------------------- inputs
traits = pd.read_csv(TRAITS, index_col=0)
print(f"Trait matrix: {traits.shape[0]} genomes x {traits.shape[1]} traits")

gpa = pd.read_csv(GPA, index_col=0, low_memory=False)
# reduce to the binary presence/absence block: keep only columns that look like genomes
bin_cols = [c for c in gpa.columns if set(pd.unique(gpa[c].astype(str))) <= {"0", "1"}]
if bin_cols:
    mat = gpa[bin_cols].astype(int)
else:
    meta_like = {"Non-unique Gene name", "Annotation", "No. isolates", "No. sequences"}
    mat = gpa.drop(columns=[c for c in gpa.columns if c in meta_like], errors="ignore")
    mat = mat.notna().astype(int)

row_sum = mat.sum(axis=1)
n_gen = mat.shape[1]
n_variable = int(((row_sum > 0) & (row_sum < n_gen)).sum())
print(f"Presence/absence matrix: {mat.shape[0]} orthogroups x {n_gen} genomes")
print(f"Orthogroups with testable variance: {n_variable}\n")

# ----------------------------------------------------------------- per trait
rows, allsig = [], []

for tname in traits.columns:
    col = traits[tname].dropna()
    n_pos = int((col == 1).sum())
    n_neg = int((col == 0).sum())

    match = [f for f in files
             if os.path.basename(f).split(".")[0].lower() == str(tname).lower()]
    if not match:
        match = [f for f in files if str(tname).lower() in os.path.basename(f).lower()]

    rec = dict(trait=tname, n_genomes_positive=n_pos, n_genomes_negative=n_neg,
               n_orthogroups_tested=n_variable,
               statistical_method="Fisher's exact test with phylogeny-aware pairwise comparisons",
               multiple_testing_correction="per-trait Bonferroni < 0.05 AND "
                                           "1,000-permutation empirical p < 0.05",
               n_significant_empirical=0, n_passing_bonferroni=0,
               n_tier1_convergent=0, best_orthogroup=None, best_odds_ratio=None,
               best_empirical_p=None, best_bonferroni_p=None,
               best_max_supporting_pairs=None, best_max_opposing_pairs=None,
               result="no significant association")

    if match:
        d = pd.read_csv(match[0])
        c_gene = pick(d.columns, "Gene", "Orthogroup")
        c_emp = pick(d.columns, "Empirical_p", "empirical")
        c_bon = pick(d.columns, "Bonferroni_p", "bonferroni")
        c_or = pick(d.columns, "Odds_ratio", "odds")
        c_sup = pick(d.columns, "Max_supporting_pairs", "supporting")
        c_opp = pick(d.columns, "Max_opposing_pairs", "opposing")

        mask = pd.Series(True, index=d.index)
        if c_emp: mask &= pd.to_numeric(d[c_emp], errors="coerce") < EMP_ALPHA
        if c_bon: mask &= pd.to_numeric(d[c_bon], errors="coerce") < BON_ALPHA
        sig = d[mask]
        rec["n_significant_empirical"] = int(len(sig))
        if c_bon:
            rec["n_passing_bonferroni"] = int(
                (pd.to_numeric(d[c_bon], errors="coerce") < 0.05).sum())
        if c_sup and c_opp:
            sup = pd.to_numeric(sig[c_sup], errors="coerce")
            opp = pd.to_numeric(sig[c_opp], errors="coerce")
            rec["n_tier1_convergent"] = int(((sup >= MIN_PAIRS) & ((sup - opp) > 0)).sum())

        if len(sig):
            key = c_emp or c_bon
            best = sig.loc[pd.to_numeric(sig[key], errors="coerce").idxmin()]
            rec.update(best_orthogroup=best.get(c_gene),
                       best_odds_ratio=best.get(c_or),
                       best_empirical_p=best.get(c_emp),
                       best_bonferroni_p=best.get(c_bon),
                       best_max_supporting_pairs=best.get(c_sup),
                       best_max_opposing_pairs=best.get(c_opp),
                       result=f"{len(sig)} significant gene-trait pair(s)")
            keep = [c for c in [c_gene, c_or, c_emp, c_bon, c_sup, c_opp] if c]
            tmp = sig[keep].copy(); tmp.insert(0, "trait", tname)
            allsig.append(tmp)
    else:
        rec["result"] = "no Scoary output file located for this trait"

    rows.append(rec)
    print(f"{tname:28s} pos={n_pos:4d} neg={n_neg:4d}  "
          f"sig={rec['n_significant_empirical']:3d}  "
          f"bonf={rec['n_passing_bonferroni']:3d}  tier1={rec['n_tier1_convergent']:3d}")

summary = pd.DataFrame(rows)
summary.to_csv(f"{OUTDIR}/S6_pangwas_per_trait_summary.csv", index=False)
print(f"\nWrote {OUTDIR}/S6_pangwas_per_trait_summary.csv  ({len(summary)} traits)")

if allsig:
    sigdf = pd.concat(allsig, ignore_index=True)
    sigdf.to_csv(f"{OUTDIR}/S6_pangwas_significant_pairs.csv", index=False)
    print(f"Wrote {OUTDIR}/S6_pangwas_significant_pairs.csv  ({len(sigdf)} pairs)")
    print("\nCHECK: total significant pairs should equal 92, "
          "spanning 16 traits and 18 orthogroups.")
    c_gene = [c for c in sigdf.columns if c.lower() in ("gene", "orthogroup")]
    print(f"  pairs={len(sigdf)}  traits={sigdf.trait.nunique()}"
          + (f"  orthogroups={sigdf[c_gene[0]].nunique()}" if c_gene else ""))
else:
    print("No significant pairs collected -- check the results directory path.")
