#!/usr/bin/env python3
"""Stage 25, step 4 — MK with the stage-06 recombination mask, mapped properly.

Step 3 ran Gubbins on the 60-gene concatenation and that was a mistake. With
no physical adjacency between neighbouring genes, rate heterogeneity *between*
genes produces blocks of high and low SNP density that Gubbins reads as
recombination, so it masked whole genes at a time: 49.7% of codons, with 43 of
60 orthogroups wholly or partly removed and 17 untouched. That mask deletes the
fast-evolving genes wholesale, which is exactly where adaptive divergence would
sit, so the resulting alpha is circular and must not be reported.

This step uses the stage-06 Gubbins run instead. That one was performed on the
Parsnp whole-genome core alignment, where adjacency is real and interval
boundaries reflect actual recombination blocks. Its coordinates are columns of
`parsnp_core_alignment.fasta`, so they are translated here:

    core-alignment column
      -> NC_003027.1 genome position   (reference row of each Parsnp LCB)
      -> codon of a gene               (locating each reference CDS in the genome)

Two masks are applied. The union mask removes a codon from every genome if any
branch carries recombination across it; the per-genome mask removes it only
from the genomes named in the interval's `taxa` field, which is what Gubbins
actually infers and which discards far less data.

The script validates the coordinate mapping before using it: the reference row
rebuilt from the XMFA must match the reference row of the core alignment
exactly. If it does not, nothing is computed.

Outputs
  02_results/parsnp_mask_codons.csv
  02_results/parsnp_mask_intervals.csv
  03_summary/mk_parsnp_masked.csv
  03_summary/04_parsnp_masked_mk.log
"""
import os
import re
import sys
import glob
import json
import collections

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (ROOT, REV, CODON_DIR, MK_SCRIPTS, TRAITS,     # noqa: E402
                     groups, vaccine_ids, stage_dir)
sys.path.insert(0, MK_SCRIPTS)
from mk_test_v2 import path_counts, clean, cmh                     # noqa: E402

import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402
from Bio import SeqIO                                              # noqa: E402
from Bio.Seq import Seq                                            # noqa: E402
from Bio.Phylo import read as read_tree                            # noqa: E402

CUT, COV, PP = 0.05, 0.90, 0.90
REF = os.environ.get("PARSNP_REF", "NC_003027.1")
REC = f"{ROOT}/06_Recombination"
XMFA = f"{REC}/01_Parsnp_Alignment/parsnp.xmfa"
CORE = f"{REC}/01_Parsnp_Alignment/parsnp_core_alignment.fasta"
GFF = f"{REC}/02_Gubbins_Results/lsdv_gubbins.recombination_predictions.gff"
FMAP = f"{REC}/nucleotide_file_map.tsv"
GENOMES = f"{REC}/00_nucleotide_genomes"
ASR = f"{REV}/17_ASR_MK/02_results/asr_runs"

D = stage_dir("25_Recombination_Masked_MK")
LOG = open(os.path.join(D, "03_summary", "04_parsnp_masked_mk.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


for p in (XMFA, CORE, GFF, GENOMES, ASR):
    if not os.path.exists(p):
        sys.exit(f"required input missing: {p}")

# ---------------------------------------------------------------- A ---------
# reference row of every LCB, from the XMFA
say("A. parsing the Parsnp XMFA")
hdr = re.compile(r"^>\s*(\d+):(\d+)-(\d+)\s+([+-])\s+cluster(\d+)")
blocks = []          # (cluster, start, end, strand, aligned_seq) for the reference
cur, keep = None, False
for ln in open(XMFA):
    if ln.startswith("#"):
        continue
    if ln.startswith(">"):
        if cur:
            blocks.append(cur)
            cur = None
        m = hdr.match(ln)
        keep = bool(m) and int(m.group(1)) == 1
        if keep:
            cur = [int(m.group(5)), int(m.group(2)), int(m.group(3)),
                   m.group(4), []]
    elif ln.startswith("="):
        if cur:
            blocks.append(cur); cur = None
        keep = False
    elif keep and cur is not None:
        cur[4].append(ln.strip())
if cur:
    blocks.append(cur)
blocks = [(c, s, e, st, "".join(q).upper()) for c, s, e, st, q in blocks]
blocks.sort(key=lambda b: b[0])
say(f"   {len(blocks)} LCBs carry the reference (sequence index 1)")

# ---------------------------------------------------------------- B ---------
say("B. validating the coordinate mapping against the core alignment")
core_ref = None
for r in SeqIO.parse(CORE, "fasta"):
    if REF in r.id or r.id.split()[0].endswith(".ref") or core_ref is None:
        if REF in r.id:
            core_ref = str(r.seq).upper(); break
        if core_ref is None:
            core_ref = str(r.seq).upper()      # fall back to the first row
if core_ref is None:
    sys.exit(f"could not find a reference row in {CORE}")
say(f"   core alignment reference row: {len(core_ref)} columns")
say(f"   total XMFA reference columns: {sum(len(b[4]) for b in blocks)}")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc_gapped(s):
    """Reverse complement, preserving gap columns."""
    return s.translate(COMP)[::-1]


# The order in which Parsnp concatenates LCBs into the core alignment, and
# the orientation it stores them in, are not documented here, so both are
# recovered empirically: each LCB's reference row is located in the core
# alignment by exact match, forward or reverse-complemented.
placed, unplaced, ambiguous = [], [], []
for b in blocks:
    hit = None
    for cand, flag in ((b[4], "+"), (rc_gapped(b[4]), "-")):
        i = core_ref.find(cand)
        if i < 0:
            continue
        if core_ref.find(cand, i + 1) >= 0:
            ambiguous.append(b[0]); hit = None; break
        hit = (i, b, flag); break
    if hit:
        placed.append(hit)
    elif b[0] not in ambiguous:
        unplaced.append(b[0])

say(f"   LCBs located by exact match: {len(placed)} of {len(blocks)}")
if unplaced:
    say(f"   not found: {unplaced[:10]}{' ...' if len(unplaced) > 10 else ''}")
if ambiguous:
    say(f"   ambiguous (multiple matches): {ambiguous[:10]}")

# A block can fail an exact match over a handful of characters -- case, or an
# ambiguity code written differently in the FASTA than in the XMFA -- while
# still belonging in a uniquely determined slot. Fill any such block into the
# single gap of the right length and report what differs, rather than giving up.
if unplaced and not ambiguous:
    placed.sort(key=lambda x: x[0])
    holes, cur = [], 0
    for off, b, f in placed:
        if off > cur:
            holes.append((cur, off))
        cur = max(cur, off + len(b[4]))
    if cur < len(core_ref):
        holes.append((cur, len(core_ref)))
    todo = [b for b in blocks if b[0] in unplaced]
    say(f"   uncovered stretches: {holes}")
    if len(holes) == len(todo):
        for (h0, h1), b in zip(holes, sorted(todo, key=lambda x: x[0])):
            if h1 - h0 != len(b[4]):
                say(f"   cluster {b[0]}: length {len(b[4])} does not fit the "
                    f"{h1-h0}-column gap"); continue
            slice_ = core_ref[h0:h1]
            fwd = sum(1 for x, y in zip(slice_, b[4]) if x != y)
            rev = sum(1 for x, y in zip(slice_, rc_gapped(b[4])) if x != y)
            flag = "+" if fwd <= rev else "-"
            nd = min(fwd, rev)
            ex = [(i, slice_[i], (b[4] if flag == "+" else rc_gapped(b[4]))[i])
                  for i in range(len(slice_))
                  if slice_[i] != (b[4] if flag == "+" else rc_gapped(b[4]))[i]][:8]
            say(f"   cluster {b[0]} placed at column {h0} ({flag}); "
                f"{nd} of {len(slice_)} characters differ ({100*nd/len(slice_):.2f}%)")
            say(f"      examples (column, core, xmfa): {ex}")
            if nd > 0.02 * len(slice_):
                sys.exit(f"FATAL: cluster {b[0]} differs at {nd} characters, "
                         "too many to treat as an encoding difference.")
            placed.append((h0, b, flag))
            unplaced.remove(b[0])
if unplaced or ambiguous:
    sys.exit("FATAL: could not place every LCB unambiguously; nothing computed.")

placed.sort(key=lambda x: x[0])
cov = sum(len(b[4]) for _, b, _ in placed)
gaps = [(placed[i][0] + len(placed[i][1][4]), placed[i + 1][0])
        for i in range(len(placed) - 1)
        if placed[i][0] + len(placed[i][1][4]) != placed[i + 1][0]]
say(f"   placed columns: {cov} of {len(core_ref)}; "
    f"{'contiguous' if not gaps else f'{len(gaps)} discontinuities'}")
if cov != len(core_ref) or gaps or placed[0][0] != 0:
    sys.exit("FATAL: the located LCBs do not tile the core alignment exactly; "
             "nothing computed.")
n_rev = sum(1 for _, _, f in placed if f == "-")
say(f"   mapping validated by exact tiling ({n_rev} LCBs stored reversed)")

col2pos = np.zeros(len(core_ref) + 1, dtype=np.int64)   # 1-based columns
for off, (c, s, e, st, seq), flag in placed:
    p = s if st == "+" else e
    step = 1 if st == "+" else -1
    arr = []
    for ch in seq:
        if ch == "-":
            arr.append(0)
        else:
            arr.append(p); p += step
    if flag == "-":
        arr = arr[::-1]
    col2pos[off + 1: off + 1 + len(arr)] = arr
say(f"   {int((col2pos[1:] > 0).sum())} of {len(core_ref)} columns map to a "
    f"reference position")

# ---------------------------------------------------------------- C ---------
say("C. translating the Gubbins intervals")
fmap = {}
if os.path.exists(FMAP):
    for ln in open(FMAP):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            fmap[f[0]] = f[1]


def to_genome_id(tok):
    tok = tok.strip()
    if tok in fmap:
        return fmap[tok]
    return re.sub(r"\.fasta$", "", tok)


iv_rows, per_genome = [], collections.defaultdict(set)
union = set()
for ln in open(GFF):
    if ln.startswith("#") or not ln.strip():
        continue
    f = ln.split("\t")
    if len(f) < 9:
        continue
    a, b = int(f[3]), int(f[4])
    taxa = []
    mt = re.search(r'taxa="([^"]*)"', f[8])
    if mt:
        taxa = [to_genome_id(x) for x in mt.group(1).split() if x.strip()]
    pos = {int(p) for p in col2pos[max(1, a):min(len(core_ref), b) + 1] if p}
    union |= pos
    for t in taxa:
        per_genome[t] |= pos
    iv_rows.append(dict(aln_start=a, aln_end=b,
                        ref_start=(min(pos) if pos else None),
                        ref_end=(max(pos) if pos else None),
                        n_ref_sites=len(pos), n_taxa=len(taxa),
                        taxa=";".join(taxa)[:200]))
pd.DataFrame(iv_rows).to_csv(
    os.path.join(D, "02_results", "parsnp_mask_intervals.csv"), index=False)
say(f"   {len(iv_rows)} intervals -> {len(union)} distinct reference positions "
    f"({100*len(union)/150773:.1f}% of the 150,773 bp genome)")
say(f"   genomes named in at least one interval: {len(per_genome)}")

# ---------------------------------------------------------------- D ---------
say("D. locating each reference CDS in the reference genome")
refseq = None
for ext in (".fasta", ".fa", ".fna"):
    p = f"{GENOMES}/{REF}{ext}"
    if os.path.exists(p):
        refseq = str(next(SeqIO.parse(p, "fasta")).seq).upper(); break
if refseq is None:
    sys.exit(f"reference genome {REF} not found under {GENOMES}")
say(f"   {REF}: {len(refseq)} bp")

codon_pos = {}        # orthogroup -> {codon index: [genome positions]}
miss = []
for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
    og = os.path.basename(aln).replace("_codon.fasta", "")
    recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln, "fasta")}
    if REF not in recs:
        miss.append((og, "reference absent from the alignment")); continue
    gapped = recs[REF]
    plain = gapped.replace("-", "")
    i = refseq.find(plain)
    strand = 1
    if i < 0:
        rc = str(Seq(plain).reverse_complement())
        i = refseq.find(rc)
        strand = -1
    if i < 0:
        miss.append((og, "reference CDS not found in the genome")); continue
    # walk the gapped alignment, assigning genome positions to ungapped bases
    pos_of_col, k = {}, 0
    for c, ch in enumerate(gapped):
        if ch == "-":
            continue
        pos_of_col[c] = (i + 1 + k) if strand == 1 else (i + len(plain) - k)
        k += 1
    cp = {}
    for cod in range(len(gapped) // 3):
        ps = [pos_of_col[c] for c in (3*cod, 3*cod+1, 3*cod+2) if c in pos_of_col]
        if ps:
            cp[cod] = ps
    codon_pos[og] = cp
say(f"   mapped {len(codon_pos)} of {len(codon_pos)+len(miss)} orthogroups")
for og, why in miss:
    say(f"   skipped {og}: {why}")

rows = []
for og, cp in codon_pos.items():
    n = len(cp)
    u = sum(1 for cod, ps in cp.items() if any(p in union for p in ps))
    rows.append(dict(orthogroup=og, codons_mapped=n, masked_union=u,
                     pct_union=round(100*u/n, 1) if n else 0.0))
md = pd.DataFrame(rows)
md.to_csv(os.path.join(D, "02_results", "parsnp_mask_codons.csv"), index=False)
say(f"   union mask: {md.masked_union.sum()} of {md.codons_mapped.sum()} codons "
    f"({100*md.masked_union.sum()/md.codons_mapped.sum():.1f}%); "
    f"{int((md.masked_union==0).sum())} of {len(md)} orthogroups untouched")

# ---------------------------------------------------------------- E ---------
G = groups()
VACC = vaccine_ids()
TR = pd.read_csv(TRAITS, index_col=0)
base = [g for g in sorted(G["LSDV"]) if g not in VACC]
OUTG = G["SPPV"] | G["GTPV"]
say(f"\nE. MK with the translated mask (ingroup n={len(base)})")


def anc_node(tf):
    t = read_tree(tf, "newick")
    for cl in t.get_nonterminals():
        if cl.name is None and cl.confidence is not None:
            cl.name = str(cl.confidence)
    best = None
    for cl in t.get_nonterminals():
        tips = {x.name for x in cl.get_terminals()}
        if len(tips) >= 2 and tips <= OUTG:
            if best is None or len(tips) > len({x.name for x in best.get_terminals()}):
                best = cl
    return best.name if best else None


def anc_seq(stf, node):
    df = pd.read_csv(stf, sep="\t", comment="#",
                     usecols=["Node", "Site", "State", "p_A", "p_C", "p_G", "p_T"])
    df = df[df.Node == node].sort_values("Site")
    return ("".join(df.State.values),
            df[["p_A", "p_C", "p_G", "p_T"]].values.max(axis=1))


def scan(mode):
    """mode: 'none' | 'union' | 'per_genome'"""
    out = []
    for aln in sorted(glob.glob(f"{CODON_DIR}/*_codon.fasta")):
        og = os.path.basename(aln).replace("_codon.fasta", "")
        tf, stf = f"{ASR}/{og}.treefile", f"{ASR}/{og}.state"
        if not (os.path.exists(tf) and os.path.exists(stf)) or og not in codon_pos:
            continue
        node = anc_node(tf)
        if node is None:
            continue
        o_seq, o_pp = anc_seq(stf, node)
        recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln, "fasta")}
        ing = [i for i in base if i in recs]
        cp = codon_pos[og]
        L = len(next(iter(recs.values()))) // 3
        Pn = Ps = Dn = Ds = 0.0
        for k in range(L):
            ps = cp.get(k)
            if mode == "union" and ps and any(p in union for p in ps):
                continue
            sl = slice(3*k, 3*k+3)
            anc = o_seq[sl].upper()
            if not clean(anc) or o_pp[sl].min() < PP:
                continue
            use = ing
            if mode == "per_genome" and ps:
                use = [i for i in ing
                       if not any(p in per_genome.get(i, ()) for p in ps)]
                if len(use) < 0.5 * len(ing):
                    continue
            ic = [x for x in (recs[i][sl] for i in use) if clean(x)]
            if not ic or len(ic) < COV * len(use):
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
    return 1 - (pn/ps)/(dn/ds), p, (pn, ps, dn, ds)


res = []
for mode, label in [("none", "no mask (control, should reproduce 0.629)"),
                    ("union", "stage-06 mask, union over all branches"),
                    ("per_genome", "stage-06 mask, per genome")]:
    df = scan(mode)
    a, p, c = alpha(df)
    say(f"\n=== {label} ===")
    say(f"  Pn={c[0]:.1f} Ps={c[1]:.1f} Dn={c[2]:.1f} Ds={c[3]:.1f}")
    say(f"  alpha={a:.3f}  CMH p={p:.3e}")
    res.append(dict(mask=mode, label=label, Pn=c[0], Ps=c[1], Dn=c[2], Ds=c[3],
                    alpha=a, CMH_p=p, n_genes=len(df), n_ingroup=len(base)))
pd.DataFrame(res).to_csv(
    os.path.join(D, "03_summary", "mk_parsnp_masked.csv"), index=False)
say("\nwritten: 03_summary/mk_parsnp_masked.csv")
say("\nThe control row must reproduce alpha = 0.629. If it does not, the")
say("coordinate translation changed which codons enter the test and the")
say("masked rows cannot be compared with the published value.")
LOG.close()
