#!/usr/bin/env python3
"""Stage 28, step 2 — name the substitutions that fixed on the LSDV branch.

Your positive result is currently a single number: about 0.6 of fixed
amino-acid differences on the branch leading to LSDV were fixed by
directional selection. The data to make that concrete are already computed.
Stage 17 reconstructed an ancestral codon at 16,018 of 16,215 sites, so the
specific substitutions can be listed: gene, codon, ancestral residue,
derived residue.

That turns an abstract estimate into a table of candidate loci, which is
the part a reader can follow up and cite. With alpha near 0.6, roughly
three in five of the non-synonymous entries are expected to be adaptive,
though the test cannot say which.

A site is included when the ancestral codon is confidently reconstructed
(every position PP >= 0.90), the LSDV ingroup is monomorphic at it, and the
LSDV codon differs from the ancestral one. These are exactly the sites that
contribute Dn and Ds to the MK test, so the counts here must reproduce it.

Requires stage 17 output. If it is missing, run that stage first:
  python3 ../17_ASR_MK/scripts/01_run_asr.py

Outputs
  02_results/fixed_substitutions.csv     every fixed difference, annotated
  03_summary/substitutions_by_gene.csv   per-gene Dn/Ds and the residues
  03_summary/02_ancestral_substitutions.log
"""
import os
import sys
import glob
import collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _s28 import ROOT, opener                                     # noqa: E402
from _common import CODON_DIR, REV, groups, vaccine_ids           # noqa: E402

import pandas as pd                                               # noqa: E402
from Bio import SeqIO, Phylo                                      # noqa: E402
from Bio.Data import CodonTable                                   # noqa: E402

PP = float(os.environ.get("PP", "0.90"))
ASR = f"{REV}/17_ASR_MK/02_results/asr_runs"
TBL = CodonTable.unambiguous_dna_by_id[1].forward_table

D, say, close = opener("02_ancestral_substitutions.log")

if not os.path.isdir(ASR):
    say(f"stage 17 ASR output not found at\n  {ASR}\n")
    say("Run stage 17 first (it needs IQ-TREE and Biopython):")
    say("  conda activate iqtree_env")
    say("  cd $PANG/revision_2026-11/17_ASR_MK")
    say("  python3 scripts/01_run_asr.py")
    say("\nThen re-run this script.")
    close()
    sys.exit(1)


def aa(c):
    c = c.upper()
    if len(c) != 3 or any(x not in "ACGT" for x in c):
        return None
    return "*" if c in ("TAA", "TAG", "TGA") else TBL.get(c)


G = groups()
VACC = vaccine_ids()
ingroup = [g for g in sorted(G["LSDV"]) if g not in VACC]
OUTG = G["SPPV"] | G["GTPV"]
say(f"ingroup: {len(ingroup)} non-vaccine LSDV genomes")
say(f"outgroup: {len(OUTG)} SPPV + GTPV genomes")
say(f"ancestral posterior probability threshold: {PP}")


def anc_node(tf):
    t = Phylo.read(tf, "newick")
    for cl in t.get_nonterminals():
        if cl.name is None and cl.confidence is not None:
            cl.name = str(cl.confidence)
    best = None
    for cl in t.get_nonterminals():
        tips = {x.name for x in cl.get_terminals()}
        if len(tips) >= 2 and tips <= OUTG:
            if best is None or len(tips) > len({x.name for x
                                                in best.get_terminals()}):
                best = cl
    return best.name if best else None


rows = []
poly = []
n_gene = 0
for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
    og = os.path.basename(aln).replace("_codon.fasta", "")
    tf, stf = f"{ASR}/{og}.treefile", f"{ASR}/{og}.state"
    if not (os.path.exists(tf) and os.path.exists(stf)):
        continue
    node = anc_node(tf)
    if node is None:
        continue
    st = pd.read_csv(stf, sep="\t", comment="#",
                     usecols=["Node", "Site", "State", "p_A", "p_C",
                              "p_G", "p_T"])
    st = st[st.Node == node].sort_values("Site")
    anc_seq = "".join(st.State.values)
    anc_pp = st[["p_A", "p_C", "p_G", "p_T"]].values.max(axis=1)

    recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln, "fasta")}
    seqs = [recs[i] for i in ingroup if i in recs]
    if not seqs:
        continue
    n_gene += 1
    L = len(next(iter(recs.values()))) // 3
    for k in range(L):
        sl = slice(3 * k, 3 * k + 3)
        a_cod = anc_seq[sl].upper()
        if aa(a_cod) is None or anc_pp[sl].min() < PP:
            continue
        ic = [s[sl] for s in seqs if aa(s[sl]) is not None]
        if len(ic) < 0.90 * len(seqs):
            continue
        cnt = collections.Counter(ic)
        if len(cnt) > 1:
            # polymorphic, not a fixed difference, but the amino-acid changes
            # segregating here are the internal null for how radical a
            # non-synonymous change is when selection has NOT fixed it
            major = cnt.most_common(1)[0][0]
            for v, c2 in cnt.items():
                if v == major or c2 / len(ic) < 0.05:
                    continue
                m_aa, v_aa = aa(major), aa(v)
                if m_aa and v_aa and m_aa != v_aa:
                    poly.append((og, k + 1, m_aa, v_aa))
            continue
        d_cod = cnt.most_common(1)[0][0]
        if d_cod == a_cod:
            continue
        a_aa, d_aa = aa(a_cod), aa(d_cod)
        rows.append(dict(
            orthogroup=og, codon=k + 1,
            ancestral_codon=a_cod, lsdv_codon=d_cod,
            ancestral_aa=a_aa, lsdv_aa=d_aa,
            synonymous=(a_aa == d_aa),
            substitution=(f"{a_aa}{k + 1}{d_aa}" if a_aa != d_aa else ""),
            min_anc_pp=round(float(anc_pp[sl].min()), 4),
            n_ingroup_called=len(ic)))

if not rows:
    say("\nno fixed differences recovered; check that the ASR output and the "
        "codon alignments refer to the same orthogroups")
    close()
    sys.exit(1)

sub = pd.DataFrame(rows)
sub.to_csv(os.path.join(D, "02_results", "fixed_substitutions.csv"),
           index=False)
nsyn = sub[~sub.synonymous]
say(f"\n{n_gene} genes analysed")
say(f"{len(sub)} fixed differences on the LSDV branch")
say(f"  non-synonymous {len(nsyn)}")
say(f"  synonymous     {int(sub.synonymous.sum())}")
mk_total = 74.8 + 83.2
say(f"\nAgainst the MK test: Dn 74.8 + Ds 83.2 = {mk_total:.1f} substitutions,")
say(f"versus {len(sub)} changed codons here, a ratio of "
    f"{mk_total / max(len(sub), 1):.2f} substitutions per codon. Nei-Gojobori")
say("averages over pathways, so a codon differing at two or three positions")
say("contributes more than one; a ratio between 1 and 2 is expected and")
say("means the two counts agree. A ratio near 1 or above 2 would not.")

per = (sub.groupby("orthogroup")
       .agg(Dn=("synonymous", lambda s: int((~s).sum())),
            Ds=("synonymous", "sum"),
            n_sites=("codon", "count")).reset_index())
per["Ds"] = per.Ds.astype(int)
# map by orthogroup rather than assigning a differently-indexed Series, which
# silently produced an empty column
subs_by_gene = (nsyn.groupby("orthogroup").substitution
                .apply(lambda x: ";".join(sorted(x))))
per["substitutions"] = per.orthogroup.map(subs_by_gene).fillna("")
per = per.sort_values("Dn", ascending=False)
per.to_csv(os.path.join(D, "03_summary", "substitutions_by_gene.csv"),
           index=False)

say(f"\n=== genes with the most non-synonymous fixed differences ===")
say(f"  {'orthogroup':<13}{'Dn':>4}{'Ds':>4}   substitutions")
for r in per.head(15).itertuples():
    say(f"  {r.orthogroup:<13}{r.Dn:>4}{r.Ds:>4}   "
        f"{r.substitutions[:70]}")

no_fixed = n_gene - len(per)
no_nonsyn = int((per.Dn == 0).sum()) + no_fixed
say(f"\n=== how the {n_gene} genes divide ===")
say(f"  no fixed difference of any kind     : {no_fixed}")
say(f"  fixed differences, none amino-acid  : {int((per.Dn == 0).sum())}")
say(f"  at least one amino-acid change      : {int((per.Dn > 0).sum())}")
say(f"\nThe manuscript reports 17 of 60 genes with no fixed difference at")
say(f"all. This run gives {no_fixed}, because a site is used here only when")
say("the ingroup is strictly monomorphic and every ancestral codon position")
say(f"reaches PP >= {PP}, which is stricter than the MK counting. Reconcile")
say("the two before quoting either number.")

ART = os.environ.get("ARTEFACT_OG", "OG0000120")
if ART in set(per.orthogroup):
    r = per[per.orthogroup == ART].iloc[0]
    say(f"\n=== {ART} dominates the amino-acid count ===")
    say(f"  it contributes Dn = {int(r.Dn)} of {len(nsyn)} "
        f"({100 * r.Dn / len(nsyn):.0f}%)")
    say(f"  excluding it: Dn = {len(nsyn) - int(r.Dn)}, "
        f"Ds = {int(sub.synonymous.sum()) - int(r.Ds)}")
    say(f"  This is the OX-2/CD200 locus whose aBSREL signal the manuscript")
    say(f"  already attributes to a 29-codon unalignable block. Its")
    say(f"  substitutions are in the output table but should be reported")
    say(f"  separately, not pooled with the rest.")

# radical vs conservative, by Grantham-style charge/polarity grouping
CLASS = {}
for grp, res in [("nonpolar", "AVLIMFWPG"), ("polar", "STCYNQ"),
                 ("positive", "KRH"), ("negative", "DE")]:
    for r in res:
        CLASS[r] = grp
nsyn = nsyn.copy()
nsyn["radical"] = [CLASS.get(a) != CLASS.get(b)
                   for a, b in zip(nsyn.ancestral_aa, nsyn.lsdv_aa)]
say(f"\n=== character of the {len(nsyn)} non-synonymous changes ===")
say(f"  change physicochemical class : {int(nsyn.radical.sum())} "
    f"({100 * nsyn.radical.mean():.1f}%)")
say(f"  conservative within class    : {int((~nsyn.radical).sum())}")

# A bare percentage means nothing without a null. The amino-acid changes
# segregating as polymorphism within LSDV are the internal comparison: they
# are changes selection has NOT fixed.
pdf = pd.DataFrame(poly, columns=["orthogroup", "codon", "aa1", "aa2"])
if len(pdf) >= 10:
    pdf["radical"] = [CLASS.get(a) != CLASS.get(b)
                      for a, b in zip(pdf.aa1, pdf.aa2)]
    pdf.to_csv(os.path.join(D, "02_results",
                            "polymorphic_nonsyn.csv"), index=False)
    from scipy.stats import fisher_exact
    tab = [[int(nsyn.radical.sum()), int((~nsyn.radical).sum())],
           [int(pdf.radical.sum()), int((~pdf.radical).sum())]]
    orv, pv = fisher_exact(tab)
    say(f"\n  internal null: {len(pdf)} non-synonymous changes segregating")
    say(f"  as polymorphism within LSDV (>= 5% minor allele), of which")
    say(f"  {int(pdf.radical.sum())} ({100 * pdf.radical.mean():.1f}%) are radical")
    say(f"\n  fixed {100 * nsyn.radical.mean():.1f}% vs polymorphic "
        f"{100 * pdf.radical.mean():.1f}%  OR = {orv:.2f}, "
        f"Fisher p = {pv:.3f}")
    if pv >= 0.05:
        say("  No difference. The fixed changes are no more radical than the")
        say("  ones still segregating, so this gives no additional evidence")
        say("  either way. Report it as such.")
    elif orv > 1:
        say("  Fixed changes are MORE radical than segregating ones, which is")
        say("  what directional selection predicts.")
    else:
        say("  Fixed changes are LESS radical than segregating ones, which is")
        say("  what purifying selection on radical changes predicts, and is")
        say("  evidence against reading these fixations as adaptive.")
else:
    say(f"\n  only {len(pdf)} polymorphic non-synonymous changes recovered,")
    say("  too few for an internal comparison; the percentage above has no")
    say("  null to be read against and should not be interpreted alone.")

say("\nwritten: 02_results/fixed_substitutions.csv")
if len(poly) >= 10:
    say("        02_results/polymorphic_nonsyn.csv")
say("        03_summary/substitutions_by_gene.csv")
close()
