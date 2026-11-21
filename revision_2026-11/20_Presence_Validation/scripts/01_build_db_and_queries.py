#!/usr/bin/env python3
"""Stage 20, step 1 — build a tblastn database and the query set.

Reviewer-proofing for the pan-GWAS result. OrthoFinder calls a gene absent
when Prokka did not predict it, which conflates true absence with annotation
dropout in a fragmented de novo assembly. Stage 19 shows the candidate loci
are distributed very unevenly between the two halves of the dataset, so each
candidate is re-tested at the DNA level, directly against the assemblies,
independently of annotation.

Step 1 writes:
  01_inputs/all_genomes.fna   290 assemblies, contigs renamed <genome>|<contig>
  01_inputs/queries.faa       one representative protein per candidate
  01_inputs/query_index.csv   which orthogroup each query came from
and makes the BLAST database.

Candidates default to the 30 provenance-exclusive orthogroups written by
stage 26 step 1, plus every convergent pan-GWAS association so that the
associations that are NOT provenance-exclusive act as controls. Override
with CANDIDATES="OG0000159,OG0000161,...".
"""
import os
import sys
import glob
import subprocess
import shutil

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
from _common import (ROOT, REV, OG_SEQ_DIR, GENOME_DIR, ANNOT_ALL,  # noqa: E402
                     stage_dir)
import pandas as pd                                                 # noqa: E402
from Bio import SeqIO                                               # noqa: E402

D = stage_dir("20_Presence_Validation")
INP = os.path.join(D, "01_inputs")

if not shutil.which("makeblastdb") or not shutil.which("tblastn"):
    sys.exit("BLAST+ not found. Activate an environment that has it, e.g.\n"
             "  conda activate prokka_env\n"
             "or install it:  conda create -n blast_env -c bioconda blast")

# ---- candidate orthogroups --------------------------------------------------
cands = []
if os.environ.get("CANDIDATES"):
    cands = [x.strip() for x in os.environ["CANDIDATES"].split(",") if x.strip()]
else:
    # the provenance-exclusive orthogroups written by stage 26 step 1
    prov = f"{REV}/26_Figure_Regeneration/01_inputs/nonviral_orthogroups.txt"
    if os.path.exists(prov):
        got = [l.strip() for l in open(prov) if l.strip()]
        cands += got
        print(f"{len(got)} provenance-exclusive orthogroups from stage 26")
    # plus every convergent pan-GWAS association, so the controls are tested
    conv = f"{REV}/26_Figure_Regeneration/03_summary/scoary_convergent.csv"
    if os.path.exists(conv):
        got = sorted(set(pd.read_csv(conv)["Gene"]))
        cands += got
        print(f"{len(got)} convergent pan-GWAS orthogroups from stage 26")
    legacy = f"{ROOT}/11_Scoary/02_summary/scoary_tier1_robust.csv"
    if not os.path.exists(conv) and os.path.exists(legacy):
        cands += list(pd.read_csv(legacy)["Gene"])
cands = sorted(set(cands))
if not cands:
    sys.exit(
        "no candidate orthogroups.\n"
        "Run stage 26 step 1 first:\n"
        "  python3 ../26_Figure_Regeneration/scripts/01_nonviral_screen.py\n"
        "or name them yourself:\n"
        '  CANDIDATES="OG0000159,OG0000161,OG0000156" python3 '
        "scripts/01_build_db_and_queries.py")
print(f"{len(cands)} candidate orthogroups: {', '.join(cands[:12])}"
      f"{' ...' if len(cands) > 12 else ''}")

# ---- queries: longest protein per orthogroup -------------------------------
ann = pd.read_csv(ANNOT_ALL).set_index("Orthogroup")["product"]
idx, n_q = [], 0
with open(f"{INP}/queries.faa", "w") as out:
    for og in cands:
        hits = (glob.glob(f"{OG_SEQ_DIR}/{og}.fa")
                + glob.glob(f"{OG_SEQ_DIR}/{og}.faa"))
        if not hits:
            print(f"  WARNING: no sequence file for {og}"); continue
        recs = list(SeqIO.parse(hits[0], "fasta"))
        if not recs:
            print(f"  WARNING: {hits[0]} is empty"); continue
        rep = max(recs, key=lambda r: len(r.seq))
        out.write(f">{og}\n{str(rep.seq).rstrip('*')}\n")
        idx.append(dict(Orthogroup=og, representative=rep.id,
                        length_aa=len(str(rep.seq).rstrip('*')),
                        n_members=len(recs),
                        product=str(ann.get(og, ""))[:80]))
        n_q += 1
pd.DataFrame(idx).to_csv(f"{INP}/query_index.csv", index=False)
print(f"wrote {n_q} queries -> 01_inputs/queries.faa")

# ---- subject: all assemblies in one file, contigs tagged by genome ---------
gen_files = sorted(glob.glob(f"{GENOME_DIR}/*"))
gen_files = [f for f in gen_files
             if os.path.splitext(f)[1].lower() in
             (".fa", ".fna", ".fasta", ".fas", ".seq")]
if not gen_files:
    sys.exit(f"no genome FASTA files under {GENOME_DIR}")
n_g = n_c = 0
with open(f"{INP}/all_genomes.fna", "w") as out:
    for f in gen_files:
        gid = os.path.splitext(os.path.basename(f))[0]
        for rec in SeqIO.parse(f, "fasta"):
            out.write(f">{gid}|{rec.id}\n{str(rec.seq)}\n")
            n_c += 1
        n_g += 1
print(f"wrote {n_g} genomes / {n_c} contigs -> 01_inputs/all_genomes.fna")

subprocess.run(["makeblastdb", "-in", f"{INP}/all_genomes.fna",
                "-dbtype", "nucl", "-out", f"{INP}/genomes_db"], check=True)
print("BLAST database built: 01_inputs/genomes_db")
