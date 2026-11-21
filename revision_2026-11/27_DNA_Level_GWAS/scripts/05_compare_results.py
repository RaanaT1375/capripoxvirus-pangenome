#!/usr/bin/env python3
"""Stage 27, step 5 — the DNA-level association result, against the old one.

Applies exactly the filters the manuscript describes: Bonferroni < 0.05 and
empirical p < 0.05 from 1,000 tree-constrained permutations, then the
convergence criterion of at least three supporting pairs with net positive
support. Reports the new convergent set, and what became of each of the
seven associations the annotation-based analysis produced.

Outputs
  03_summary/scoary_dna_significant.csv
  03_summary/scoary_dna_convergent.csv
  03_summary/old_vs_new.csv
  03_summary/05_compare_results.log
"""
import os
import re
import sys
import glob


def _locate_common():
    env = os.environ.get("COMMON_DIR")
    if env and os.path.exists(os.path.join(env, "_common.py")):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.exists(os.path.join(d, "_common.py")):
            return d
        p = os.path.dirname(d)
        if p == d:
            break
        d = p
    sys.exit("could not find _common.py above " + os.path.abspath(__file__)
             + "\nSet COMMON_DIR to the directory that holds it.")


sys.path.insert(0, _locate_common())
from _common import REV, stage_dir                                # noqa: E402

import pandas as pd                                               # noqa: E402

BONF, EMP, MINPAIRS = 0.05, 0.05, 3

# the seven convergent associations from the annotation-based matrix
OLD = [("OG0000157", "Asia"), ("OG0000155", "Asia"),
       ("OG0000154", "Europe"), ("OG0000159", "Africa"),
       ("OG0000161", "Africa"), ("OG0000156", "Africa"),
       ("OG0000152", "Vaccine")]

D = stage_dir("27_DNA_Level_GWAS")
RES = os.path.join(D, "02_results", "scoary_results_dna")
LOG = open(os.path.join(D, "03_summary", "05_compare_results.log"), "w")
STAMP = re.compile(r"_(\d{2}_\d{2}_\d{4}(?:_\d+)?)\.results\.csv$")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


files = sorted(glob.glob(os.path.join(RES, "**", "*.results.csv"),
                         recursive=True))
if not files:
    sys.exit(f"no Scoary results under {RES}; run 04_run_scoary.sbatch first")
say(f"{len(files)} result files")

frames = []
for f in files:
    base = os.path.basename(f)
    m = STAMP.search(base)
    trait = base[:m.start()] if m else base.replace(".results.csv", "")
    t = pd.read_csv(f)
    t["Trait"] = trait
    frames.append(t)
allr = pd.concat(frames, ignore_index=True)
say(f"{len(allr)} rows across {allr.Trait.nunique()} traits")


def col(df, *names):
    for n in names:
        if n in df.columns:
            return n
    sys.exit(f"none of these columns found: {names}\nhave: {list(df.columns)}")


c_bonf = col(allr, "Bonferroni_p", "Bonferroni_p_value")
c_emp = col(allr, "Empirical_p", "Empirical_p_value")
c_or = col(allr, "Odds_ratio", "Odds_Ratio")
c_sup = col(allr, "Max_supporting_pairs", "Best_supporting_pairs")
c_opp = col(allr, "Max_opposing_pairs", "Worst_opposing_pairs")

sig = allr[(allr[c_bonf] < BONF) & (allr[c_emp] < EMP)].copy()
sig.to_csv(os.path.join(D, "03_summary", "scoary_dna_significant.csv"),
           index=False)
say(f"\npassing both filters: {len(sig)} pairs, {sig.Trait.nunique()} traits, "
    f"{sig.Gene.nunique()} orthogroups")
say("annotation-based matrix gave: 92 pairs, 16 traits, 18 orthogroups")

conv = sig[(sig[c_sup] >= MINPAIRS) &
           ((sig[c_sup] - sig[c_opp]) > 0)].sort_values(["Trait", c_emp])
conv.to_csv(os.path.join(D, "03_summary", "scoary_dna_convergent.csv"),
            index=False)
say(f"\nmeeting the convergence criterion: {len(conv)}")
say("annotation-based matrix gave: 7")
if len(conv):
    say(f"\n  {'Trait':<24}{'Orthogroup':<13}{'OR':>10}{'emp p':>11}  pairs")
    for r in conv.itertuples():
        say(f"  {getattr(r, 'Trait'):<24}{getattr(r, 'Gene'):<13}"
            f"{getattr(r, c_or):>10.4g}{getattr(r, c_emp):>11.4g}  "
            f"{int(getattr(r, c_sup))}:{int(getattr(r, c_opp))}")
else:
    say("\n  none. On the DNA-level matrix no association survives both the")
    say("  significance filters and the convergence criterion.")

# ---- what happened to the old seven ----------------------------------------
say("\n" + "=" * 70)
say("THE SEVEN ANNOTATION-BASED ASSOCIATIONS, RE-TESTED")
say("=" * 70)
rows = []
for og, trait in OLD:
    hit = allr[(allr.Gene == og) & (allr.Trait.str.lower() == trait.lower())]
    if hit.empty:
        say(f"{og:<12} {trait:<10} absent from the DNA-level output "
            f"(no variance left to test)")
        rows.append(dict(orthogroup=og, trait=trait, status="not testable"))
        continue
    r = hit.sort_values(c_emp).iloc[0]
    passes = (r[c_bonf] < BONF) and (r[c_emp] < EMP)
    converg = (r[c_sup] >= MINPAIRS) and (r[c_sup] - r[c_opp] > 0)
    verdict = ("still convergent" if (passes and converg)
               else "significant, not convergent" if passes
               else "no longer significant")
    say(f"{og:<12} {trait:<10} OR={r[c_or]:<10.4g} emp_p={r[c_emp]:<9.4g} "
        f"bonf={r[c_bonf]:<10.4g} {int(r[c_sup])}:{int(r[c_opp])}  {verdict}")
    rows.append(dict(orthogroup=og, trait=trait, odds_ratio=r[c_or],
                     empirical_p=r[c_emp], bonferroni_p=r[c_bonf],
                     supporting=int(r[c_sup]), opposing=int(r[c_opp]),
                     status=verdict))
pd.DataFrame(rows).to_csv(os.path.join(D, "03_summary", "old_vs_new.csv"),
                          index=False)

say("\nHow to read this:")
say("  Any association that survives here is supported by sequence rather")
say("  than by gene calls, and can be reported as a biological result.")
say("  Any that disappears was an artefact of annotation. If none survives,")
say("  that is a clean negative result and the stronger claim is the")
say("  methodological one.")
LOG.close()
