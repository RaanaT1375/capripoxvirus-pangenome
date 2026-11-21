#!/usr/bin/env python3
"""Stage 28, step 7 — separate clustered substitutions from distributed ones.

Step 2 named the fixed differences. Reading the names revealed something the
counts hid: in some genes every substitution falls inside a few codons.

  OG0000120  16 changes in codons 177-205, a 29-codon window
  OG0000067   4 changes in codons 60-65, six codons, with Ds = 0

A 29-codon window is exactly the unalignable C-terminal block the manuscript
already attributes OG0000120's aBSREL signal to, recovered here independently
and without being told about it. That is reassuring about the diagnosis and
disqualifying for the substitution list: 16 changes in one window are one
alignment event, not 16 adaptive substitutions, and pooling them inflates
every per-gene count they enter.

This step applies the same test the manuscript already used on OG0000120 to
every gene with at least MIN_N non-synonymous fixed differences, so the
clustered ones can be set aside before any locus is named in the paper.

Diagnostics per gene:

  span         exact probability that n substitutions drawn uniformly from a
               gene of L codons would span no more than the observed number
               of codons. This is the deciding test: it is exact, needs no
               permutations, and is scale-free.

  window       the largest number of substitutions in any W-codon window,
               with a permutation p value. Kept for comparison only. A fixed
               window misleads on short genes: on an 89-codon gene a
               30-codon window is a third of the sequence, so four
               substitutions inside it look unremarkable (p = 0.10) even
               when their actual span is six codons (p = 3.5e-04).

  uniformity   Kolmogorov-Smirnov against a uniform null. Reported for
               completeness; weak at n = 3-5 and not used to decide.

  gaps         gap frequency in the substitution columns against the gene's
               own mean, for tested genes only, and only when the local
               frequency clears GAP_FLOOR. Without that floor a ratio on a
               near-zero baseline flags noise: 0.23% against 0.05% is a
               4.4x ratio and fewer than one gapped cell in three hundred.

A gene is flagged when its Holm-corrected span probability is below ALPHA, or
when the gap rule fires. Flagged genes are not deleted -- they go to their
own file, to be reported as alignment-limited rather than silently dropped.

Do NOT read the gap column as evidence about mechanism. Where clustering is
strongest the gap frequency at the substitution columns is zero while the
gene averages over 10%, and that is circular, not informative: step 2 admits
a site only when every codon position reconstructs at PP >= 0.90 and the
ingroup is strictly monomorphic, which a column gapped in a third of the
taxa cannot satisfy. The filter therefore SELECTS gap-free columns, so
measuring gaps at them is guaranteed to return a low number. Gap evidence
about a block has to be measured over the block, as a sliding window, not at
the filter-selected columns inside it.

The gap column is kept only to catch the opposite case -- a gene whose
substitution columns are themselves materially gappy, which the floor below
makes a meaningful flag rather than a ratio on noise.

Outputs
  02_results/clustering_tests.csv        every gene, every diagnostic
  03_summary/substitutions_clean.csv     substitutions in unflagged genes
  03_summary/substitutions_flagged.csv   substitutions in flagged genes
  03_summary/07_substitution_clustering.log
"""
import os
import sys
import glob
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import ROOT, opener                                     # noqa: E402,F401
from _common import CODON_DIR                                     # noqa: E402

import numpy as np                                                # noqa: E402
import pandas as pd                                               # noqa: E402
from Bio import SeqIO                                             # noqa: E402
from scipy.stats import kstest                                    # noqa: E402

MIN_N = int(os.environ.get("MIN_N", "3"))
W = int(os.environ.get("WINDOW", "30"))
N_PERM = int(os.environ.get("N_PERM", "20000"))
ALPHA = float(os.environ.get("ALPHA", "0.05"))
GAP_MULT = float(os.environ.get("GAP_MULT", "2.0"))
GAP_FLOOR = float(os.environ.get("GAP_FLOOR", "0.05"))
SEED = int(os.environ.get("SEED", "20261120"))


def p_span(n, L, s):
    """Exact P(span <= s) for n positions drawn without replacement from L.

    The number of n-subsets of {1..L} whose span is exactly w is
    (L - w + 1) * C(w - 2, n - 2), so the tail sums in closed form. This is
    scale-free: it compares the observed span with the gene's own length
    instead of with a window width chosen in advance, which is what the
    fixed-window statistic got wrong on short genes.
    """
    if n < 2 or L < n:
        return float("nan")
    s = min(s, L)
    acc = sum((L - w + 1) * comb(w - 2, n - 2) for w in range(n, s + 1))
    return acc / comb(L, n)

D, say, close = opener("07_substitution_clustering.log")
SUB = os.path.join(D, "02_results", "fixed_substitutions.csv")

if not os.path.exists(SUB):
    say(f"{SUB} not found; run 02_ancestral_substitutions.py first")
    close()
    sys.exit(1)

sub = pd.read_csv(SUB)
say(f"{len(sub)} fixed differences read from step 2")

# ---- find the non-synonymous ones ------------------------------------------
# Step 2 may label this column in one of several ways; accept any of them
# rather than guessing one and silently getting an empty frame.
cand = [c for c in sub.columns
        if c.lower() in ("type", "kind", "class", "syn", "nonsyn",
                         "synonymous", "is_nonsyn", "effect")]
nsyn = None
for c in cand:
    v = sub[c].astype(str).str.lower()
    hit = v.str.startswith("non") | v.str.contains("nonsyn") | (v == "false")
    if hit.any() and not hit.all():
        nsyn = sub[hit].copy()
        say(f"non-synonymous sites taken from column '{c}': {len(nsyn)}")
        break
if nsyn is None:
    # fall back to comparing the residues directly
    aa_cols = [c for c in sub.columns if "aa" in c.lower()
               or "residue" in c.lower()]
    if len(aa_cols) >= 2:
        a, b = aa_cols[0], aa_cols[1]
        nsyn = sub[sub[a].astype(str) != sub[b].astype(str)].copy()
        say(f"non-synonymous sites inferred from {a} != {b}: {len(nsyn)}")
    else:
        say("could not identify the non-synonymous sites. Columns present:")
        say("  " + ", ".join(sub.columns))
        say("\nSet the column explicitly and re-run, or report this output.")
        close()
        sys.exit(1)

if "codon" not in nsyn.columns:
    say("no 'codon' column in the step 2 output; cannot test positions")
    close()
    sys.exit(1)

# ---- gene lengths and gap profiles from the codon alignments ---------------
lengths, gapfreq = {}, {}
for f in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
    og = os.path.basename(f).replace("_codon.fasta", "")
    recs = [str(r.seq) for r in SeqIO.parse(f, "fasta")]
    if not recs:
        continue
    L = min(len(r) for r in recs) // 3
    lengths[og] = L
    arr = np.array([list(r[:L * 3]) for r in recs])
    # gap frequency per codon: any gap character in the three positions
    g = (arr == "-") | (arr == "N") | (arr == "n")
    per_codon = g.reshape(len(recs), L, 3).any(axis=2).mean(axis=0)
    gapfreq[og] = per_codon
say(f"{len(lengths)} codon alignments read for gene lengths and gap profiles")

rng = np.random.default_rng(SEED)
rows = []
for og, grp in nsyn.groupby("orthogroup"):
    pos = sorted(int(p) for p in grp.codon)
    n = len(pos)
    L = lengths.get(og)
    if L is None or L < 2:
        rows.append(dict(orthogroup=og, n=n, gene_len=np.nan,
                         span=np.nan, max_in_window=np.nan,
                         p_window=np.nan, ks_D=np.nan, p_ks=np.nan,
                         gap_at_sites=np.nan, gap_gene_mean=np.nan,
                         gap_ratio=np.nan, flag="", reason="no alignment"))
        continue
    span = pos[-1] - pos[0] + 1

    # largest count in any W-codon window, kept for comparison only
    def max_window(p):
        p = np.sort(np.asarray(p))
        best = 0
        for i in range(len(p)):
            best = max(best, int(np.sum((p >= p[i]) & (p < p[i] + W))))
        return best

    obs_w = max_window(pos)

    if n >= MIN_N:
        perm = np.empty(N_PERM, int)
        for i in range(N_PERM):
            perm[i] = max_window(rng.choice(L, size=n, replace=False) + 1)
        p_w = float((perm >= obs_w).sum() + 1) / (N_PERM + 1)
        p_sp = p_span(n, L, span)
        ks_D, p_ks = kstest((np.array(pos) - 0.5) / L, "uniform")
    else:
        p_w = p_sp = np.nan
        ks_D = p_ks = np.nan

    gp = gapfreq.get(og)
    if gp is not None and len(gp) >= L:
        idx = [min(max(p - 1, 0), len(gp) - 1) for p in pos]
        g_at = float(np.mean(gp[idx]))
        g_mean = float(np.mean(gp))
        ratio = g_at / g_mean if g_mean > 0 else np.nan
    else:
        g_at = g_mean = ratio = np.nan

    # The gap rule fires only for a gene that was tested AND whose local gap
    # frequency is materially above zero. Without the floor, a ratio computed
    # on a near-zero baseline flags noise: 0.23% against 0.05% is a 4.4x
    # ratio and under one gapped cell in three hundred.
    reasons = []
    if (n >= MIN_N and not np.isnan(ratio) and ratio > GAP_MULT
            and g_at >= GAP_FLOOR):
        reasons.append(f"substitution columns {ratio:.1f}x gappier than the "
                       f"gene ({100 * g_at:.1f}% vs {100 * g_mean:.1f}%)")
    rows.append(dict(orthogroup=og, n=n, gene_len=L, span=span,
                     max_in_window=obs_w, p_window=p_w, p_span=p_sp,
                     ks_D=ks_D, p_ks=p_ks,
                     gap_at_sites=g_at, gap_gene_mean=g_mean,
                     gap_ratio=ratio,
                     flag="", reason="; ".join(reasons)))

res = pd.DataFrame(rows)

# ---- Holm correction on the span statistic ---------------------------------
# The span probability is exact and scale-free, so it is the deciding test.
# The fixed-width window is kept in the table but no longer decides anything:
# on an 89-codon gene a 30-codon window is a third of the sequence, so four
# substitutions landing inside it is unremarkable (p = 0.10) even when their
# actual span is six codons (p = 3.5e-04).
m = res.p_span.notna()
k = int(m.sum())
res["p_span_holm"] = np.nan
if k:
    order = res.index[m][np.argsort(res.loc[m, "p_span"].values)]
    running = 0.0
    for i, idx in enumerate(order):
        adj = min(1.0, res.at[idx, "p_span"] * (k - i))
        running = max(running, adj)
        res.at[idx, "p_span_holm"] = running

for idx in res.index:
    rs = [res.at[idx, "reason"]] if res.at[idx, "reason"] else []
    ph = res.at[idx, "p_span_holm"]
    if not np.isnan(ph) and ph < ALPHA:
        rs.insert(0, f"{int(res.at[idx, 'n'])} substitutions span only "
                     f"{int(res.at[idx, 'span'])} codons of "
                     f"{int(res.at[idx, 'gene_len'])} "
                     f"(p={res.at[idx, 'p_span']:.2e}, Holm {ph:.2e})")
    res.at[idx, "reason"] = "; ".join(rs)
    res.at[idx, "flag"] = "FLAGGED" if rs else ""

res = res.sort_values(["flag", "n"], ascending=[False, False])
res.to_csv(os.path.join(D, "02_results", "clustering_tests.csv"), index=False)

tested = res[res.n >= MIN_N]
say(f"\n{len(res)} genes carry a non-synonymous fixed difference; "
    f"{len(tested)} have at least {MIN_N} and were tested")
say(f"deciding test: exact span probability, Holm-corrected over {k} genes")
say(f"gap rule: ratio > {GAP_MULT}x AND local gap frequency >= "
    f"{100 * GAP_FLOOR:.0f}%, tested genes only")
say(f"window statistic ({W} codons, {N_PERM} permutations, seed {SEED}) "
    f"reported for comparison only")

say(f"\n=== genes with at least {MIN_N} non-synonymous fixed differences ===")
say(f"  {'gene':<13}{'n':>3}{'len':>6}{'span':>6}"
    f"{'p_span':>10}{'Holm':>10}{'p_win':>8}{'p_KS':>8}"
    f"{'gap@':>7}{'mean':>7}  flag")
for r in tested.itertuples():
    say(f"  {r.orthogroup:<13}{r.n:>3}{int(r.gene_len):>6}{int(r.span):>6}"
        f"{r.p_span:>10.2e}{r.p_span_holm:>10.2e}"
        f"{r.p_window:>8.3f}{r.p_ks:>8.3f}"
        f"{100 * r.gap_at_sites:>6.1f}%{100 * r.gap_gene_mean:>6.1f}%"
        f"  {r.flag}")

flagged = res[res.flag == "FLAGGED"]
say(f"\n=== {len(flagged)} genes flagged as alignment-limited ===")
for r in flagged.itertuples():
    say(f"  {r.orthogroup}: {r.reason}")
if len(flagged):
    say("\n  These are not adaptive substitution lists. A dense run of changes")
    say("  in a few dozen codons is one localized event however many codons")
    say("  it spans, so each flagged gene contributes one uncertain event,")
    say("  not n independent ones. Do not name their residues as candidate")
    say("  adaptive sites.")
    low = [r for r in flagged.itertuples()
           if not np.isnan(r.gap_at_sites)
           and r.gap_at_sites < 0.01 and r.gap_gene_mean > 0.05]
    if low:
        say("\n  Do NOT draw a mechanism from the gap columns below:")
        for r in low:
            say(f"    {r.orthogroup}: gaps at the substitution columns "
                f"{100 * r.gap_at_sites:.1f}% against "
                f"{100 * r.gap_gene_mean:.1f}% for the gene.")
        say("    This is circular. Step 2 admits a site only when every codon")
        say("    position reconstructs at PP >= 0.90 and the ingroup is")
        say("    strictly monomorphic, which a column gapped in a third of")
        say("    the taxa cannot satisfy, so the filter selects gap-free")
        say("    columns and a low number here is guaranteed. It says nothing")
        say("    about whether the surrounding block is well aligned.")
        say("    Measure that over the block with a sliding window instead.")

keep = set(res.orthogroup[res.flag != "FLAGGED"])
clean = nsyn[nsyn.orthogroup.isin(keep)]
bad = nsyn[~nsyn.orthogroup.isin(keep)]
clean.to_csv(os.path.join(D, "03_summary", "substitutions_clean.csv"),
             index=False)
bad.to_csv(os.path.join(D, "03_summary", "substitutions_flagged.csv"),
           index=False)

say(f"\n=== effect on the counts ===")
say(f"  non-synonymous fixed differences, all genes : {len(nsyn)}")
say(f"  in flagged genes                           : {len(bad)} "
    f"({100 * len(bad) / max(len(nsyn), 1):.0f}%)")
say(f"  surviving                                  : {len(clean)} "
    f"across {clean.orthogroup.nunique()} genes")
if len(bad):
    worst = bad.orthogroup.value_counts()
    say(f"\n  largest single contributor removed: {worst.index[0]} "
        f"with {int(worst.iloc[0])} changes")
    say("  Before this screen that one gene supplied "
        f"{100 * int(worst.iloc[0]) / max(len(nsyn), 1):.0f}% of the whole")
    say("  non-synonymous total, which is the dominance step 2 flagged.")

say("\n=== what this does and does not change ===")
say("  The MK test itself is NOT affected: alpha is estimated from counts")
say("  over all sites, and removing sites post hoc because they cluster")
say("  would bias it. This screen governs which individual loci and")
say("  residues may be NAMED as candidates, which is a separate claim.")
say("  Report both numbers and the screen, so a reviewer can see that the")
say("  named candidates were filtered and the test statistic was not.")

say("\nwritten: 02_results/clustering_tests.csv")
say("        03_summary/substitutions_clean.csv, substitutions_flagged.csv")
close()
