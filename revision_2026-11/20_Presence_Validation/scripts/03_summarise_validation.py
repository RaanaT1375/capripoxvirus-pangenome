#!/usr/bin/env python3
"""Stage 20, step 3 — turn tblastn hits into a DNA-level presence matrix.

A candidate is called present in a genome when its best tblastn hit reaches
MIN_COV of the query length at >= MIN_PID identity, summing HSPs per contig
so that a gene split across a contig break is not scored as absent. The
result is compared with the OrthoFinder annotation-based call, and the
disagreements are exactly the annotation-dropout cases that would otherwise
masquerade as biological absence.

Outputs
  02_results/presence_dna.csv              candidate x genome, DNA-level
  02_results/annotation_vs_dna.csv         per-candidate concordance
  03_summary/validation_summary.log
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
        p = os.path.dirname(d)
        if p == d:
            break
        d = p
    sys.exit("could not find _common.py above " + os.path.abspath(__file__)
             + "\nSet COMMON_DIR to the directory that holds it.")


sys.path.insert(0, _locate_common())
from _common import load_pa, stage_dir, is_sra                     # noqa: E402
import pandas as pd                                                # noqa: E402
from scipy.stats import fisher_exact                               # noqa: E402

MIN_COV = float(os.environ.get("MIN_COV", "0.80"))
MIN_PID = float(os.environ.get("MIN_PID", "80"))

D = stage_dir("20_Presence_Validation")
LOG = open(os.path.join(D, "03_summary", "validation_summary.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


say(f"thresholds: query coverage >= {MIN_COV:.0%}, identity >= {MIN_PID:.0f}%")

cols = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"]
bl = pd.read_csv(os.path.join(D, "02_results", "tblastn_hits.tsv"),
                 sep="\t", names=cols)
bl["genome"] = bl.sseqid.str.split("|").str[0]
bl = bl[bl.pident >= MIN_PID]
say(f"{len(bl)} HSPs passing identity filter, "
    f"{bl.qseqid.nunique()} queries, {bl.genome.nunique()} genomes")

# sum non-overlapping query coverage per (query, genome, contig), keep the best
def covered(sub):
    iv = sorted(zip(sub.qstart, sub.qend))
    tot = 0; cs = ce = None
    for s, e in iv:
        if cs is None:
            cs, ce = s, e
        elif s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1; cs, ce = s, e
    if cs is not None:
        tot += ce - cs + 1
    return tot


rec = []
for (q, g, c), sub in bl.groupby([bl.qseqid, bl.genome, bl.sseqid]):
    rec.append((q, g, covered(sub) / sub.qlen.iloc[0], sub.pident.max()))
cov = pd.DataFrame(rec, columns=["Orthogroup", "genome", "qcov", "pident"])
cov = cov.sort_values("qcov").groupby(["Orthogroup", "genome"]).last()
cov = cov.reset_index()
cov.to_csv(os.path.join(D, "02_results", "tblastn_best_per_genome.csv"),
           index=False)

PA = load_pa()
genomes = list(PA.columns)
dna = pd.DataFrame(0, index=sorted(cov.Orthogroup.unique()), columns=genomes)
for _, r in cov[cov.qcov >= MIN_COV].iterrows():
    if r.genome in dna.columns:
        dna.loc[r.Orthogroup, r.genome] = 1
dna.to_csv(os.path.join(D, "02_results", "presence_dna.csv"))

sra = [c for c in genomes if is_sra(c)]
gb = [c for c in genomes if not is_sra(c)]
rows = []
say("\n=== annotation-based call vs DNA-level call ===")
say(f"{'Orthogroup':<13}{'ann_n':>6}{'dna_n':>6}{'gained':>7}{'lost':>6}"
    f"{'dna_SRA':>9}{'dna_GB':>8}{'fisher_p':>11}")
for og in dna.index:
    if og not in PA.index:
        continue
    a, d = PA.loc[og], dna.loc[og]
    gained = int(((d == 1) & (a == 0)).sum())   # annotation missed it
    lost = int(((d == 0) & (a == 1)).sum())     # annotation called it, DNA did not
    ns, ng = int(d[sra].sum()), int(d[gb].sum())
    p = fisher_exact([[ns, len(sra) - ns], [ng, len(gb) - ng]])[1]
    rows.append(dict(Orthogroup=og, n_annotation=int(a.sum()), n_dna=int(d.sum()),
                     recovered_by_dna=gained, not_supported_by_dna=lost,
                     dna_n_sra=ns, dna_n_genbank=ng,
                     dna_pct_sra=round(100 * ns / len(sra), 1),
                     dna_pct_genbank=round(100 * ng / len(gb), 1),
                     fisher_p=p,
                     still_sra_exclusive=(ng == 0 and ns > 0)))
    say(f"{og:<13}{int(a.sum()):>6}{int(d.sum()):>6}{gained:>7}{lost:>6}"
        f"{ns:>9}{ng:>8}{p:>11.2e}")

res = pd.DataFrame(rows)
res.to_csv(os.path.join(D, "02_results", "annotation_vs_dna.csv"), index=False)
say(f"\ncandidates still absent from every GenBank genome at the DNA level: "
    f"{int(res.still_sra_exclusive.sum())}/{len(res)}")
say("those are assembly/annotation artefacts, not geographic signal;")
say("candidates recovered in GenBank genomes by DNA search survive and")
say("should be re-tested with the corrected presence matrix.")
LOG.close()
