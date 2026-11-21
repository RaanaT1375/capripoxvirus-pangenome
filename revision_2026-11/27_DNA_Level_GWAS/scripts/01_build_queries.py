#!/usr/bin/env python3
"""Stage 27, step 1 — queries for every orthogroup, not just the candidates.

Stage 20 tested 34 orthogroups at the DNA level and found that the
annotation-based presence/absence matrix does not describe gene content: the
three strongest pan-GWAS associations are present in 269-289 of 290 genomes
once the sequence itself is searched. This stage rebuilds the whole matrix
the same way and re-runs the association analysis on it.

Step 1 writes one representative protein per orthogroup and makes (or
reuses) the nucleotide BLAST database of all 290 assemblies.

Outputs
  01_inputs/queries.faa        one representative per orthogroup
  01_inputs/query_index.csv    orthogroup, representative, length, product
  01_inputs/genomes_db.*       the BLAST database, or a note that stage 20's
                               is being reused
"""
import os
import sys
import glob
import shutil
import subprocess


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
from _common import (REV, OG_SEQ_DIR, GENOME_DIR, ANNOT_ALL,      # noqa: E402
                     stage_dir, load_pa)

import pandas as pd                                               # noqa: E402
from Bio import SeqIO                                             # noqa: E402

D = stage_dir("27_DNA_Level_GWAS")
INP = os.path.join(D, "01_inputs")

if not shutil.which("makeblastdb") or not shutil.which("tblastn"):
    sys.exit("BLAST+ not found. Activate an environment that has it:\n"
             "  conda activate prokka_env")

pa = load_pa()
ogs = list(pa.index)
print(f"{len(ogs)} orthogroups in the annotation-based matrix")

# ---- queries ----------------------------------------------------------------
ann = {}
if os.path.exists(ANNOT_ALL):
    try:
        a = pd.read_csv(ANNOT_ALL)
        if {"Orthogroup", "product"} <= set(a.columns):
            ann = a.set_index("Orthogroup")["product"].to_dict()
    except Exception as e:                                        # noqa: BLE001
        print(f"  (could not read annotations: {e})")

idx, missing = [], []
with open(f"{INP}/queries.faa", "w") as out:
    for og in ogs:
        hits = (glob.glob(f"{OG_SEQ_DIR}/{og}.fa")
                + glob.glob(f"{OG_SEQ_DIR}/{og}.faa"))
        if not hits:
            missing.append(og)
            continue
        recs = list(SeqIO.parse(hits[0], "fasta"))
        if not recs:
            missing.append(og)
            continue
        rep = max(recs, key=lambda r: len(r.seq))
        seq = str(rep.seq).rstrip("*")
        out.write(f">{og}\n{seq}\n")
        idx.append(dict(Orthogroup=og, representative=rep.id,
                        length_aa=len(seq), n_members=len(recs),
                        product=str(ann.get(og, ""))[:80]))
pd.DataFrame(idx).to_csv(f"{INP}/query_index.csv", index=False)
print(f"wrote {len(idx)} queries -> 01_inputs/queries.faa")
if missing:
    print(f"WARNING: no sequence file for {len(missing)} orthogroups: "
          f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    print("  these cannot be tested at the DNA level and will be carried "
          "over from the annotation matrix in step 3")

# ---- subject database -------------------------------------------------------
reuse = f"{REV}/20_Presence_Validation/01_inputs/genomes_db"
if glob.glob(reuse + ".n*"):
    print(f"reusing the BLAST database from stage 20:\n  {reuse}")
    with open(f"{INP}/genomes_db.reused", "w") as fh:
        fh.write(reuse + "\n")
    sys.exit(0)

gen_files = [f for f in sorted(glob.glob(f"{GENOME_DIR}/*"))
             if os.path.splitext(f)[1].lower() in
             (".fa", ".fna", ".fasta", ".fas", ".seq")]
if not gen_files:
    sys.exit(f"no genome FASTA files under {GENOME_DIR}")
n_g = n_c = 0
with open(f"{INP}/all_genomes.fna", "w") as out:
    for f in gen_files:
        gid = os.path.splitext(os.path.basename(f))[0]
        for rec in SeqIO.parse(f, "fasta"):
            out.write(f">{gid}|{rec.id}\n{rec.seq}\n")
            n_c += 1
        n_g += 1
print(f"wrote {n_g} genomes / {n_c} contigs -> 01_inputs/all_genomes.fna")
subprocess.run(["makeblastdb", "-in", f"{INP}/all_genomes.fna",
                "-dbtype", "nucl", "-out", f"{INP}/genomes_db"], check=True)
print("BLAST database built: 01_inputs/genomes_db")
