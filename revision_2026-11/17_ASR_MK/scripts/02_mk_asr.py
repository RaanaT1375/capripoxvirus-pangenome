#!/usr/bin/env python3
"""Stage 17, step 2 — McDonald-Kreitman test with ML ancestral states.

Replaces the >=90% outgroup-consensus polarization of 14_MK_Test with the
marginal ML reconstruction from step 1. Everything downstream of the
ancestral call is byte-for-byte the stage-14 logic: Nei-Gojobori pathway
averaging, the 5% minor-variant cutoff, the 90% ingroup coverage filter,
sites that are polymorphic counted as P only, and a per-gene stratified CMH
test. Only the ancestral state changes, so the new alpha is directly
comparable with the published alpha = 0.604.

Two estimators are reported:
  branch : ancestor = ANC_OUT (MRCA of SPPV+GTPV)  - direct swap for consensus
  stem   : divergence counted ANC_OUT -> ANC_LSDV  - substitutions localised
           to the LSDV stem edge, the most conservative formulation

Plus a posterior-probability sweep, a gene bootstrap CI, and a recomputation
with OG0000120 (OX-2, the indel artefact documented in stage 16) removed.

Outputs
  02_results/mk_asr_pp<NN>_allgenes.csv   per-gene counts at each threshold
  03_summary/mk_asr_robustness.csv        the table for the manuscript
  03_summary/02_mk_asr.log                full log
"""
import os
import sys
import glob
import collections

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (CODON_DIR, MK_SCRIPTS, groups, vaccine_ids,  # noqa: E402
                     stage_dir)
sys.path.insert(0, MK_SCRIPTS)
from mk_test_v2 import path_counts, clean, cmh                   # noqa: E402

import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402
from Bio import SeqIO, Phylo                                     # noqa: E402

CUT, COV = 0.05, 0.90                      # identical to 14_MK_Test
PP_GRID = [0.50, 0.90, 0.95, 0.99]
PRIMARY = 0.90
ARTEFACT = "OG0000120"                     # OX-2, see 16_aBSREL_stem
N_BOOT = 10000

D = stage_dir("17_ASR_MK")
RUNS = os.path.join(D, "02_results", "asr_runs")
LOG = open(os.path.join(D, "03_summary", "02_mk_asr.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


G = groups()
OUTG = G["SPPV"] | G["GTPV"]
VACC = vaccine_ids()
ING = [g for g in sorted(G["LSDV"]) if g not in VACC]
say(f"ingroup: {len(ING)} LSDV genomes "
    f"({len(G['LSDV'])} minus {len(G['LSDV'] & VACC)} vaccine-derived)")
say(f"outgroup: {len(G['SPPV'])} SPPV + {len(G['GTPV'])} GTPV "
    f"| {len(VACC)} vaccine genomes in the dataset overall\n")


def anc_nodes(treefile):
    """(ANC_OUT, ANC_LSDV) node labels. ANC_LSDV is the far end of the stem."""
    t = Phylo.read(treefile, "newick")
    for cl in t.get_nonterminals():
        if cl.name is None and cl.confidence is not None:
            cl.name = str(cl.confidence)
    best = None
    for cl in t.get_nonterminals():
        tips = {x.name for x in cl.get_terminals()}
        if len(tips) >= 2 and tips <= OUTG:
            if best is None or len(tips) > len(
                    {x.name for x in best.get_terminals()}):
                best = cl
    if best is None:
        return None, None
    path = t.get_path(best)
    parent = path[-2] if len(path) >= 2 else t.root
    return best.name, (parent.name if parent and parent.name else None)


def anc_seqs(statefile, nodes):
    df = pd.read_csv(statefile, sep="\t", comment="#",
                     usecols=["Node", "Site", "State",
                              "p_A", "p_C", "p_G", "p_T"])
    df = df[df.Node.isin([n for n in nodes if n])]
    out = {}
    for n, g in df.groupby("Node"):
        g = g.sort_values("Site")
        out[n] = ("".join(g.State.values),
                  g[["p_A", "p_C", "p_G", "p_T"]].values.max(axis=1))
    return out


def codon_at(seq, pp, k, pp_min):
    sl = slice(3 * k, 3 * k + 3)
    c = seq[sl].upper()
    return None if (not clean(c) or pp[sl].min() < pp_min) else c


def scan(pp_min):
    rows = []
    for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
        og = os.path.basename(aln).replace("_codon.fasta", "")
        tf, stf = f"{RUNS}/{og}.treefile", f"{RUNS}/{og}.state"
        if not (os.path.exists(tf) and os.path.exists(stf)):
            continue
        n_out, n_lsdv = anc_nodes(tf)
        if n_out is None:
            say(f"  skip {og}: no outgroup clade"); continue
        seqs = anc_seqs(stf, [n_out, n_lsdv])
        if n_out not in seqs:
            say(f"  skip {og}: {n_out} absent from .state"); continue
        o_seq, o_pp = seqs[n_out]
        l_seq, l_pp = seqs.get(n_lsdv, (None, None))

        rec = {x.id: str(x.seq).upper() for x in SeqIO.parse(aln, "fasta")}
        ing = [rec[i] for i in ING if i in rec]
        L = len(next(iter(rec.values()))) // 3
        Pn = Ps = Dn = Ds = sDn = sDs = 0.0
        n_anc = n_stem = 0

        for k in range(L):
            sl = slice(3 * k, 3 * k + 3)
            anc = codon_at(o_seq, o_pp, k, pp_min)
            if anc is None:
                continue
            n_anc += 1
            ic = [x for x in (s[sl] for s in ing) if clean(x)]
            if len(ic) >= COV * len(ing):
                icc = collections.Counter(ic)
                major = icc.most_common(1)[0][0]
                if len(icc) > 1:
                    for v, cnt in icc.items():
                        if v == major or cnt / len(ic) < CUT:
                            continue
                        s_, n_ = path_counts(major, v); Ps += s_; Pn += n_
                elif major != anc:
                    s_, n_ = path_counts(anc, major); Ds += s_; Dn += n_
            if l_seq is not None:
                la = codon_at(l_seq, l_pp, k, pp_min)
                if la is not None:
                    n_stem += 1
                    if la != anc:
                        s_, n_ = path_counts(anc, la); sDs += s_; sDn += n_
        rows.append(dict(Orthogroup=og, Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds,
                         sDn=sDn, sDs=sDs, n_anc_sites=n_anc,
                         n_stem_sites=n_stem, n_codons=L,
                         node_out=n_out, node_lsdv=n_lsdv))
    return pd.DataFrame(rows)


def alpha_of(df, pn, ps, dn, ds):
    P_n, P_s, D_n, D_s = df[[pn, ps, dn, ds]].sum()
    if P_s == 0 or D_s == 0 or D_n == 0:
        return np.nan, np.nan, (P_n, P_s, D_n, D_s)
    tabs = [[[round(r[dn]), round(r[ds])], [round(r[pn]), round(r[ps])]]
            for _, r in df.iterrows()]
    p, _, _ = cmh(tabs)
    return 1 - (P_n / P_s) / (D_n / D_s), p, (P_n, P_s, D_n, D_s)


summary = []
primary_df = None
for pp in PP_GRID:
    df = scan(pp)
    df.to_csv(os.path.join(D, "02_results",
                           f"mk_asr_pp{int(pp*100)}_allgenes.csv"), index=False)
    frac = 100 * df.n_anc_sites.sum() / df.n_codons.sum()
    a_b, p_b, c_b = alpha_of(df, "Pn", "Ps", "Dn", "Ds")
    a_s, p_s, c_s = alpha_of(df, "Pn", "Ps", "sDn", "sDs")
    say(f"\n=== PP >= {pp:.2f} | confident ancestral codons "
        f"{df.n_anc_sites.sum():.0f}/{df.n_codons.sum():.0f} ({frac:.1f}%) ===")
    say(f"  branch : Pn={c_b[0]:.1f} Ps={c_b[1]:.1f} Dn={c_b[2]:.1f} "
        f"Ds={c_b[3]:.1f} -> alpha={a_b:.3f}  CMH p={p_b:.3e}")
    say(f"  stem   : Dn={c_s[2]:.1f} Ds={c_s[3]:.1f} -> alpha={a_s:.3f}  "
        f"CMH p={p_s:.3e}")
    summary += [dict(estimator="ASR branch (ANC_OUT)", pp_min=pp,
                     pct_sites=round(frac, 1), Pn=c_b[0], Ps=c_b[1],
                     Dn=c_b[2], Ds=c_b[3], alpha=a_b, CMH_p=p_b,
                     n_genes=len(df), n_ingroup=len(ING)),
                dict(estimator="ASR stem (ANC_OUT->ANC_LSDV)", pp_min=pp,
                     pct_sites=round(frac, 1), Pn=c_s[0], Ps=c_s[1],
                     Dn=c_s[2], Ds=c_s[3], alpha=a_s, CMH_p=p_s,
                     n_genes=len(df), n_ingroup=len(ING))]
    if pp == PRIMARY:
        primary_df = df

# ---- artefact removal and gene bootstrap, both at the primary threshold -----
d = primary_df
sub = d[d.Orthogroup != ARTEFACT]
a_x, p_x, c_x = alpha_of(sub, "Pn", "Ps", "Dn", "Ds")
say(f"\n=== PP >= {PRIMARY:.2f}, {ARTEFACT} (OX-2) removed ===")
say(f"  alpha={a_x:.3f}  CMH p={p_x:.3e}")
summary.append(dict(estimator=f"ASR branch, {ARTEFACT} removed",
                    pp_min=PRIMARY, pct_sites=np.nan, Pn=c_x[0], Ps=c_x[1],
                    Dn=c_x[2], Ds=c_x[3], alpha=a_x, CMH_p=p_x,
                    n_genes=len(sub), n_ingroup=len(ING)))

rng = np.random.default_rng(20261116)
boot = []
for _ in range(N_BOOT):
    s = d.iloc[rng.integers(0, len(d), len(d))]
    pn, ps, dn, ds = s[["Pn", "Ps", "Dn", "Ds"]].sum()
    if ps > 0 and ds > 0 and dn > 0:
        boot.append(1 - (pn / ps) / (dn / ds))
lo, hi = np.percentile(boot, [2.5, 97.5])
say(f"\ngene bootstrap ({N_BOOT} replicates) 95% CI: [{lo:.2f}, {hi:.2f}]")

pd.DataFrame(summary).to_csv(
    os.path.join(D, "03_summary", "mk_asr_robustness.csv"), index=False)
with open(os.path.join(D, "03_summary", "mk_asr_bootstrap_ci.txt"), "w") as fh:
    fh.write(f"n_replicates\t{N_BOOT}\nlower_2.5\t{lo:.4f}\nupper_97.5\t{hi:.4f}\n")
say("\nwritten: 03_summary/mk_asr_robustness.csv")
LOG.close()
