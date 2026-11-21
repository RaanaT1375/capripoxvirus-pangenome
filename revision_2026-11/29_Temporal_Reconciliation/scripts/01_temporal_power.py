#!/usr/bin/env python3
"""Stage 29 — reconcile the absent temporal signal with the published genus rate.

The manuscript reports no temporal signal within the post-2000 LSDV population
(slope -5.40e-06 substitutions/site/year, R^2 < 0.001, p = 0.887, n = 185) and
argues that structure dominates time. That conclusion is almost certainly
right, but as submitted it rests on two things that do not hold up:

  1. An arithmetic error. The text states that at mu = 1.0e-5
     substitutions/site/year the 110,321-bp core alignment is expected to
     accumulate "about 5 substitutions over 25 years and about 15 over the
     71-year span". The correct figures are 1.0e-5 x 110,321 x 25 = 27.6 and
     x 71 = 78.3, so the quoted values are low by a factor of 5.5, in the
     direction that flatters the argument. The ratio of structural divergence
     (100-250 substitutions across Node_147) to the temporal expectation
     falls from the implied 20-50x to 4-9x. Still dominant, much less
     dramatic.

  2. An untested assertion. "No signal" and "no detectable signal" are
     different claims, and only the second is supported unless the power to
     detect a slope of mu is actually computed. A reviewer is entitled to
     ask whether a clock of the published strength would have been visible
     in these data at all.

This stage replaces the assertion with a number. Four parts:

  expectation   expected substitutions over any span, at mu and over a range
                of rates, computed correctly, against the observed structural
                divergence

  refit         the root-to-tip regression recomputed from the stored
                per-isolate data, to confirm the published slope and to
                obtain the residual scatter the power calculation needs

  randomisation the accepted test for absence of temporal signal: dates are
                permuted across tips and the observed slope compared with
                the resulting null distribution

  power         the minimum slope detectable at 80% power given the actual
                date spread and the observed residual scatter, and the power
                to detect a slope of exactly mu. This is the part that
                decides whether the negative result is informative or merely
                underpowered.

Outputs
  02_results/temporal_refit.csv          per-isolate values used
  02_results/rate_expectations.csv       expected substitutions by rate and span
  03_summary/power_summary.csv           the power calculation
  03_summary/01_temporal_power.log
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

# ---- parameters -------------------------------------------------------------
MU = float(os.environ.get("MU", "1.0e-5"))          # published genus rate
ALN = int(os.environ.get("ALN", "110321"))          # core alignment length
PI = float(os.environ.get("PI", "4.87e-3"))         # within-LSDV diversity
EPOCH = int(os.environ.get("EPOCH", "2000"))        # post-epoch cutoff
NPERM = int(os.environ.get("NPERM", "10000"))
SEED = int(os.environ.get("SEED", "20261120"))
STRUCT_LO = float(os.environ.get("STRUCT_LO", "100"))   # Node_147 divergence
STRUCT_HI = float(os.environ.get("STRUCT_HI", "250"))

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, ".."))
for sub in ("01_inputs", "02_results", "03_summary"):
    os.makedirs(os.path.join(D, sub), exist_ok=True)
LOG = open(os.path.join(D, "03_summary", "01_temporal_power.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


say("stage 29 - temporal reconciliation and power analysis")
say(f"mu={MU:.2e}  alignment={ALN} bp  pi={PI:.2e}  epoch>={EPOCH}  "
    f"permutations={NPERM}")

# ---- locate the per-isolate temporal table ---------------------------------
# Supplementary File S4 holds genome id, collection year, root-to-tip distance
# and continent for the dated isolates. Find whatever file carries those.
#
# Do NOT walk the project tree with a recursive glob. On a shared parallel
# filesystem holding read sets, assemblies and IQ-TREE run directories, a
# `**/*.csv` pattern stats every directory and takes many minutes with no
# output, which looks exactly like a hang. Known locations are tried first and
# any fallback search is bounded and prunes the heavy directories.
KNOWN = [
    "09_Temporal_Signal/root_to_tip_LSDV_only.csv",
    "09_Temporal_Signal/root_to_tip_data.csv",
    "revision_2026-11/16_Temporal/root_to_tip_LSDV_only.csv",
]
SKIP_DIRS = {"99_Backups", "01_Raw_Data", "02_QC", "03_Assembly", "04_Annotation",
             "asr_runs", "logs", ".git", "__pycache__", "work", ".nextflow"}


def usable(f):
    """True if f is a CSV carrying both a date and a root-to-tip column."""
    try:
        t = pd.read_csv(f, nrows=5)
    except Exception:                                            # noqa: BLE001
        return False
    low = {c.lower(): c for c in t.columns}
    has_date = any(k in low for k in
                   ("year", "collection_year", "date", "decimal_year"))
    has_rtt = any("root" in k or "rtt" in k for k in low)
    return has_date and has_rtt


def find_temporal():
    env = os.environ.get("TEMPORAL")
    if env:
        if os.path.exists(env):
            return env
        say(f"TEMPORAL is set but does not exist:\n  {env}")
        return None

    root = os.path.abspath(os.path.join(D, "..", ".."))
    say(f"looking for the per-isolate temporal table under\n  {root}")

    for rel in KNOWN:
        f = os.path.join(root, rel)
        if os.path.exists(f) and usable(f):
            say(f"  found at the expected location: {rel}")
            return f

    say("  not at a known location; bounded search (depth 3, heavy "
        "directories pruned)")
    for dirpath, dirnames, filenames in os.walk(root):
        depth = dirpath[len(root):].count(os.sep)
        if depth >= 3:
            dirnames[:] = []
            continue
        dirnames[:] = [x for x in dirnames if x not in SKIP_DIRS]
        for fn in filenames:
            n = fn.lower()
            if not n.endswith(".csv"):
                continue
            if not any(k in n for k in
                       ("temporal", "root", "rtt", "clock", "date")):
                continue
            f = os.path.join(dirpath, fn)
            if usable(f):
                say(f"  found: {os.path.relpath(f, root)}")
                return f
    return None


src = find_temporal()
if not src:
    say("Could not find a per-isolate temporal table.")
    say("It needs a collection-year column and a root-to-tip column.")
    say("\nLook for it with:")
    say("  find /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis \\")
    say("       -name '*.csv' | xargs grep -l -i 'root.*tip\\|rtt' 2>/dev/null")
    say("\nThen re-run with:  TEMPORAL=/path/to/file.csv python3 "
        "scripts/01_temporal_power.py")
    LOG.close()
    sys.exit(1)

say(f"temporal table: {src}")
t = pd.read_csv(src)
low = {c.lower(): c for c in t.columns}
ycol = next((low[k] for k in ("decimal_year", "collection_year", "year", "date")
             if k in low), None)
rcol = next((low[k] for k in low if "root" in k or "rtt" in k), None)
ccol = next((low[k] for k in ("continent", "region", "country") if k in low),
            None)
say(f"  date column      : {ycol}")
say(f"  root-to-tip column: {rcol}")
say(f"  {len(t)} rows")
say("  columns present    : " + ", ".join(map(str, t.columns)))

d = t[[ycol, rcol] + ([ccol] if ccol else [])].copy()
d.columns = ["year", "rtt"] + (["continent"] if ccol else [])
for c in ("year", "rtt"):
    d[c] = pd.to_numeric(d[c], errors="coerce")
d = d.dropna(subset=["year", "rtt"])
say(f"  {len(d)} rows with both values, years "
    f"{d.year.min():.0f}-{d.year.max():.0f}")

# ---- part 1: the expectation, computed correctly ---------------------------
say("\n" + "=" * 72)
say("1. EXPECTED SUBSTITUTIONS -- the manuscript's arithmetic corrected")
say("=" * 72)
span_full = float(d.year.max() - d.year.min())
rows = []
for rate in (1.0e-6, 2.0e-6, 5.0e-6, MU, 2.0e-5, 5.0e-5):
    for span, label in ((25.0, "25 yr"), (span_full, f"{span_full:.0f} yr")):
        rows.append(dict(rate=rate, span_years=round(span, 1),
                         expected_substitutions=round(rate * ALN * span, 2),
                         expected_per_site=rate * span))
exp = pd.DataFrame(rows)
exp.to_csv(os.path.join(D, "02_results", "rate_expectations.csv"), index=False)

say(f"alignment {ALN} bp, dated span {span_full:.0f} years")
say(f"\n  at the published mu = {MU:.1e} substitutions/site/year:")
for span, lab in ((25.0, "25 years"), (span_full, f"{span_full:.0f} years")):
    say(f"    over {lab:<12}: {MU * ALN * span:6.1f} substitutions")
say("\n  the manuscript states 'about 5' over 25 years and 'about 15' over 71.")
say(f"  correct values are {MU * ALN * 25:.1f} and {MU * ALN * 71:.1f}: "
    f"low by {(MU * ALN * 25) / 5:.1f}x.")
implied = 5.0 / 25.0 / ALN
say(f"  the quoted pair is internally consistent with a rate of {implied:.2e},")
say(f"  which is {MU / implied:.1f}x below the rate the Methods cite.")

say(f"\n  structural divergence across Node_147: "
    f"{STRUCT_LO:.0f}-{STRUCT_HI:.0f} substitutions")
say(f"  ratio to the 25-year temporal expectation: "
    f"{STRUCT_LO / (MU * ALN * 25):.1f}-{STRUCT_HI / (MU * ALN * 25):.1f}x")
say(f"  (the manuscript's figures imply {STRUCT_LO / 5:.0f}-{STRUCT_HI / 5:.0f}x)")
say("\n  The conclusion that structure dominates time survives. Its margin")
say("  does not: state it as a few-fold excess, not an order of magnitude.")

# ---- part 2: refit ---------------------------------------------------------
say("\n" + "=" * 72)
say("2. ROOT-TO-TIP REGRESSION, REFITTED")
say("=" * 72)


def fit(sub, label):
    if len(sub) < 10:
        say(f"  {label}: only {len(sub)} isolates, skipped")
        return None
    r = stats.linregress(sub.year, sub.rtt)
    say(f"  {label:<26} n={len(sub):<5} slope={r.slope:+.3e}  "
        f"R2={r.rvalue ** 2:.4f}  p={r.pvalue:.3g}")
    return r


full = fit(d, "all dated isolates")
post = d[d.year >= EPOCH]
pfit = fit(post, f"post-{EPOCH} only")
if pfit is None:
    say("\nnot enough post-epoch isolates for the power analysis")
    LOG.close()
    sys.exit(1)

say(f"\n  published: all dated -4.87e-05, R2 0.111; "
    f"post-{EPOCH} -5.40e-06, R2 <0.001, p 0.887")
say("  compare the refit above against those; a mismatch means this table is")
say("  not the one the manuscript used and the rest of this run is void.")

resid = post.rtt - (pfit.intercept + pfit.slope * post.year)
sd_resid = float(np.std(resid, ddof=2))
sd_year = float(np.std(post.year, ddof=1))
n = len(post)
say(f"\n  residual SD (post-{EPOCH})      : {sd_resid:.3e} substitutions/site")
say(f"  SD of collection year          : {sd_year:.2f} years")
say(f"  that residual scatter is {sd_resid * ALN:.1f} substitutions "
    f"across the alignment")

d.to_csv(os.path.join(D, "02_results", "temporal_refit.csv"), index=False)

# ---- part 3: date randomisation -------------------------------------------
say("\n" + "=" * 72)
say("3. DATE RANDOMISATION")
say("=" * 72)
rng = np.random.default_rng(SEED)
yv, rv = post.year.to_numpy(float), post.rtt.to_numpy(float)
null = np.empty(NPERM)
for i in range(NPERM):
    null[i] = stats.linregress(rng.permutation(yv), rv).slope
p_two = float((np.abs(null) >= abs(pfit.slope)).mean())
say(f"  {NPERM} permutations of collection date across tips")
say(f"  observed slope          : {pfit.slope:+.3e}")
say(f"  null slope mean, SD     : {null.mean():+.3e}, {null.std():.3e}")
say(f"  null 95% interval       : {np.percentile(null, 2.5):+.3e} to "
    f"{np.percentile(null, 97.5):+.3e}")
say(f"  empirical two-sided p   : {p_two:.4f}")
say("  the observed slope sits inside the null distribution"
    if p_two >= 0.05 else
    "  the observed slope falls outside the null distribution")

# ---- part 4: power --------------------------------------------------------
say("\n" + "=" * 72)
say("4. POWER -- would a clock of the published strength have been seen?")
say("=" * 72)
se_slope = sd_resid / (sd_year * np.sqrt(n - 1))
tcrit = stats.t.ppf(0.975, n - 2)
mds = tcrit * se_slope + stats.norm.ppf(0.80) * se_slope   # 80% power, two-sided
say(f"  standard error of the slope      : {se_slope:.3e}")
say(f"  minimum slope detectable at 80%  : {mds:.3e} substitutions/site/year")
say(f"  published genus rate mu          : {MU:.3e}")
say(f"  ratio mu / minimum detectable    : {MU / mds:.2f}")

ncp = MU / se_slope
power_mu = float(stats.norm.sf(tcrit - ncp) + stats.norm.cdf(-tcrit - ncp))
say(f"\n  power to detect a slope of exactly mu : {power_mu:.1%}")

rows = [dict(quantity="residual_sd_subs_per_site", value=sd_resid),
        dict(quantity="sd_collection_year", value=sd_year),
        dict(quantity="n_post_epoch", value=n),
        dict(quantity="se_slope", value=se_slope),
        dict(quantity="min_detectable_slope_80pct", value=mds),
        dict(quantity="published_mu", value=MU),
        dict(quantity="power_at_mu", value=power_mu),
        dict(quantity="randomisation_p", value=p_two)]
pd.DataFrame(rows).to_csv(os.path.join(D, "03_summary", "power_summary.csv"),
                          index=False)

say("\n  --- how to read this ---")
if power_mu < 0.5:
    say(f"  Power is {power_mu:.1%}. A clock running at the published genus")
    say("  rate would NOT reliably have been detected in these data, so the")
    say("  absent signal is consistent with that rate rather than evidence")
    say("  against it. This is the result that makes the negative finding")
    say("  informative: say that the data cannot resolve a clock of this")
    say("  strength, not that no clock exists.")
else:
    say(f"  Power is {power_mu:.1%}. A clock at the published rate should")
    say("  have been visible, and it is not. That is a stronger claim than")
    say("  the manuscript currently makes and it needs an explanation:")
    say("  recombination disrupting root-to-tip additivity, imprecise")
    say("  collection dates, or a within-LSDV rate below the genus estimate.")
    say("  Do not report the absence as unsurprising without addressing it.")

say("\n  CAVEAT, which must be stated wherever this number is reported:")
say("  root-to-tip distances are not independent -- tips sharing ancestry")
say("  share branches -- so the effective sample size is below n and the")
say("  true standard error is larger than the one above. This power estimate")
say("  is therefore an UPPER bound, and the real power is lower. That makes")
say("  a low-power verdict safe and a high-power verdict provisional.")

# ---- timescale ------------------------------------------------------------
say("\n" + "=" * 72)
say("5. WHAT THE STANDING DIVERSITY CORRESPONDS TO IN TIME")
say("=" * 72)
tmrca = PI / (2 * MU)
say(f"  within-LSDV diversity pi         : {PI:.3e}")
say(f"  pi / 2mu                         : {tmrca:.0f} years")
say(f"  pi across the alignment          : {PI * ALN:.0f} pairwise differences")
say(f"\n  Two mean pairwise differences of {PI * ALN:.0f} substitutions at "
    f"{MU * ALN:.2f}")
say(f"  substitutions per year place the common ancestor of sampled LSDV")
say(f"  diversity on the order of {tmrca:.0f} years ago.")
say("  This is a back-of-envelope conversion under a strict clock and")
say("  constant population size, neither of which holds here. Report it as")
say("  an order of magnitude, with those assumptions named, or not at all.")

say("\nwritten: 02_results/rate_expectations.csv, temporal_refit.csv")
say("        03_summary/power_summary.csv")
LOG.close()
