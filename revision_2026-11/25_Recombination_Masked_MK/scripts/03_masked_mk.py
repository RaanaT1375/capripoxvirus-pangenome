#!/usr/bin/env python3
"""Stage 25, step 3 — McDonald-Kreitman with recombinant regions masked.

Takes the Gubbins intervals from step 2, which are in supermatrix
coordinates, maps them onto the per-gene codon alignments through the
partition table, and recomputes alpha with those codons excluded. Everything
else is stage 17 unchanged: the ML ancestral states, Nei-Gojobori pathway
averaging, the 5% minor-variant cutoff, the 90% coverage filter, and the
per-gene stratified CMH test.

A codon is masked when any of its three positions falls inside a predicted
recombinant interval on any branch. That is deliberately conservative: it
removes every site recombination is inferred to have touched anywhere in the
tree, not only in the genomes that carry it.

Four ingroups are reported so that the masking can be read against the
partition result that motivated it:
  all          all non-vaccine LSDV genomes
  masked       the same genomes, recombinant codons removed
  norecomb     only genomes with no inferred recombination
  reconly      only genomes descending from the recombination events

Outputs
  02_results/masked_codons.csv          per-orthogroup masked codon counts
  03_summary/mk_recombination_masked.csv
  03_summary/03_masked_mk.log
"""
import os
import sys
import glob
import collections

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (CODON_DIR, MK_SCRIPTS, TRAITS, groups,      # noqa: E402
                     vaccine_ids, stage_dir, REV)
sys.path.insert(0, MK_SCRIPTS)
from mk_test_v2 import path_counts, clean, cmh                   # noqa: E402

import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from Bio import SeqIO, Phylo                                     # noqa: E402

CUT, COV, PP = 0.05, 0.90, 0.90
D = stage_dir("25_Recombination_Masked_MK")
ASR = f"{REV}/17_ASR_MK/02_results/asr_runs"
GFF = os.path.join(D, "02_results", "lsdv_core_gubbins.recombination_predictions.gff")
PARTS = os.path.join(D, "01_inputs", "partitions_nt.txt")
LOG = open(os.path.join(D, "03_summary", "03_masked_mk.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


if not os.path.exists(ASR):
    sys.exit(f"stage 17 ASR output not found at {ASR}; run stage 17 first")
if not os.path.exists(GFF):
    sys.exit(f"{GFF} not found; run 02_run_gubbins.sbatch first")

# ---- masked codons per orthogroup ------------------------------------------
parts = pd.read_csv(PARTS, sep="\t")
total = int(parts.end.max())
mask = np.zeros(total + 1, dtype=bool)
n_iv = 0
for ln in open(GFF):
    if ln.startswith("#") or not ln.strip():
        continue
    f = ln.split("\t")
    if len(f) < 5:
        continue
    try:
        s, e = int(f[3]), int(f[4])
    except ValueError:
        continue
    mask[max(1, s):min(total, e) + 1] = True
    n_iv += 1
say(f"{n_iv} predicted recombinant intervals covering "
    f"{int(mask.sum())} of {total} sites ({100*mask.sum()/total:.1f}%)")

masked_codons = {}
rows = []
for r in parts.itertuples():
    n_cod = r.length // 3
    bad = set()
    for k in range(n_cod):
        a = r.start + 3 * k
        if mask[a] or mask[a + 1] or mask[a + 2]:
            bad.add(k)
    masked_codons[r.orthogroup] = bad
    rows.append(dict(orthogroup=r.orthogroup, codons=n_cod,
                     masked=len(bad), pct=round(100 * len(bad) / n_cod, 1)))
md = pd.DataFrame(rows)
md.to_csv(os.path.join(D, "02_results", "masked_codons.csv"), index=False)
say(f"codons masked: {md.masked.sum()} of {md.codons.sum()} "
    f"({100*md.masked.sum()/md.codons.sum():.1f}%); "
    f"{int((md.masked == 0).sum())} of {len(md)} orthogroups untouched")

# ---- ingroups ---------------------------------------------------------------
G = groups()
VACC = vaccine_ids()
TR = pd.read_csv(TRAITS, index_col=0)
base = [g for g in sorted(G["LSDV"]) if g not in VACC]
norec = [g for g in base if TR.loc[g, "No_recombination"] == 1]
rec = [g for g in base if TR.loc[g, "No_recombination"] == 0]
say(f"\ningroups: all={len(base)}  norecomb={len(norec)}  reconly={len(rec)}")

OUTG = G["SPPV"] | G["GTPV"]


def anc_nodes(tf):
    t = Phylo.read(tf, "newick")
    for cl in t.get_nonterminals():
        if cl.name is None and cl.confidence is not None:
            cl.name = str(cl.confidence)
    best = None
    for cl in t.get_nonterminals():
        tips = {x.name for x in cl.get_terminals()}
        if len(tips) >= 2 and tips <= OUTG:
            if best is None or len(tips) > len({x.name for x in best.get_terminals()}):
                best = cl
    if best is None:
        return None
    return best.name


def anc_seq(stf, node):
    df = pd.read_csv(stf, sep="\t", comment="#",
                     usecols=["Node", "Site", "State", "p_A", "p_C", "p_G", "p_T"])
    df = df[df.Node == node].sort_values("Site")
    return ("".join(df.State.values),
            df[["p_A", "p_C", "p_G", "p_T"]].values.max(axis=1))


def scan(ing, use_mask):
    out = []
    for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
        og = os.path.basename(aln).replace("_codon.fasta", "")
        tf, stf = f"{ASR}/{og}.treefile", f"{ASR}/{og}.state"
        if not (os.path.exists(tf) and os.path.exists(stf)):
            continue
        node = anc_nodes(tf)
        if node is None:
            continue
        o_seq, o_pp = anc_seq(stf, node)
        recs = {x.id: str(x.seq).upper() for x in SeqIO.parse(aln, "fasta")}
        seqs = [recs[i] for i in ing if i in recs]
        if not seqs:
            continue
        L = len(next(iter(recs.values()))) // 3
        bad = masked_codons.get(og, set()) if use_mask else set()
        Pn = Ps = Dn = Ds = 0.0
        for k in range(L):
            if k in bad:
                continue
            sl = slice(3 * k, 3 * k + 3)
            anc = o_seq[sl].upper()
            if not clean(anc) or o_pp[sl].min() < PP:
                continue
            ic = [x for x in (s[sl] for s in seqs) if clean(x)]
            if len(ic) < COV * len(seqs):
                continue
            icc = collections.Counter(ic)
            major = icc.most_common(1)[0][0]
            if len(icc) > 1:
                for v, cnt in icc.items():
                    if v == major or cnt / len(ic) < CUT:
                        continue
                    s_, n_ = path_counts(major, v); Ps += s_; Pn += n_
            elif major != anc:
                s_, n_ = path_counts(anc, major); Ds += s_; Dn += n_
        out.append(dict(Orthogroup=og, Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds))
    return pd.DataFrame(out)


def alpha(df):
    pn, ps, dn, ds = df[["Pn", "Ps", "Dn", "Ds"]].sum()
    if ps == 0 or ds == 0 or dn == 0:
        return np.nan, np.nan, (pn, ps, dn, ds)
    tabs = [[[round(r.Dn), round(r.Ds)], [round(r.Pn), round(r.Ps)]]
            for r in df.itertuples()]
    p, _, _ = cmh(tabs)
    return 1 - (pn / ps) / (dn / ds), p, (pn, ps, dn, ds)


res = []
for label, ing, m in [("all non-vaccine LSDV", base, False),
                      ("recombinant codons masked", base, True),
                      ("non-recombinant genomes only", norec, False),
                      ("recombinant-descended only", rec, False)]:
    df = scan(ing, m)
    a, p, c = alpha(df)
    say(f"\n=== {label} (n={len(ing)}{', masked' if m else ''}) ===")
    say(f"  Pn={c[0]:.1f} Ps={c[1]:.1f} Dn={c[2]:.1f} Ds={c[3]:.1f}")
    say(f"  alpha={a:.3f}  CMH p={p:.3e}")
    res.append(dict(analysis=label, n_ingroup=len(ing), masked=m,
                    Pn=c[0], Ps=c[1], Dn=c[2], Ds=c[3], alpha=a, CMH_p=p))

pd.DataFrame(res).to_csv(
    os.path.join(D, "03_summary", "mk_recombination_masked.csv"), index=False)
say("\nwritten: 03_summary/mk_recombination_masked.csv")
say("\nReading: if alpha stays positive and significant with recombinant codons")
say("masked while every genome is retained, the adaptive signal is not an")
say("artefact of introgressed polymorphism. If it collapses as it does when")
say("the recombinant genomes are dropped, the central claim must be rewritten.")
LOG.close()
