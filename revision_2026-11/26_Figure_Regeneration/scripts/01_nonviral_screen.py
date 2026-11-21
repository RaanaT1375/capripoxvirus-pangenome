#!/usr/bin/env python3
"""Stage 26, step 1 — fix the list of non-viral orthogroups.

Everything downstream depends on knowing exactly which orthogroups are
excluded, so this writes that list once, to
01_inputs/nonviral_orthogroups.txt, and the later scripts read it.

Routes, tried in order. The first three need nothing but pandas.

  1. $NONVIRAL_LIST points at a file with one orthogroup ID per line.
  2. 01_inputs/nonviral_orthogroups.txt already exists.
  3. SRA-exclusive orthogroups, computed here from the presence/absence
     matrix. An orthogroup present only in genomes assembled here from SRA
     reads, and in no GenBank genome, is the signature of contamination
     carried through assembly: a real viral gene would appear in at least one
     curated genome. This is the default and needs no other stage to have
     been run. If stage 19 happens to be present its sra_exclusive column is
     read as a cross-check.
  4. DIAMOND against the capripoxvirus reference proteomes, if diamond is on
     the path. Only reached if route 3 cannot split the genomes by name.

Whichever route is used, the annotation of every candidate is printed, so
you can see whether they look like host mitochondrial and retroviral genes
as the manuscript states, or like something else.

This script does NOT stop if the count is not 28. It reports the number it
found in a banner and writes the list, because blocking the pipeline on a
number helps nobody. Read the banner: if it is not 28, the manuscript's
"206 - 28 = 178" is wrong and must be changed to match, not the other way
round.

Outputs
  01_inputs/nonviral_orthogroups.txt
  02_results/nonviral_candidates.csv     every candidate, with annotation
  03_summary/01_nonviral_screen.log
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
from _common import (ROOT, REV, PA_CSV, OG_SEQ_DIR,               # noqa: E402
                     stage_dir, load_pa, is_sra)

import pandas as pd                                               # noqa: E402

D = stage_dir("26_Figure_Regeneration")
INP = os.path.join(D, "01_inputs")
OUT = os.path.join(INP, "nonviral_orthogroups.txt")
LOG = open(os.path.join(D, "03_summary", "01_nonviral_screen.log"), "w")

EXPECTED = 28
PROV = f"{REV}/19_Provenance_Audit/02_results/orthogroup_provenance.csv"

# Annotation words that indicate host mitochondrial, retroviral or bacterial
# origin rather than a poxvirus gene. Used for reporting, and as a
# cross-check on whichever route supplies the list.
NONVIRAL_WORDS = [
    "nadh", "cytochrome", "ubiquinone", "atp synthase", "ribosomal",
    "16s", "12s", "trna", "mitochond", "coi", "cox1", "cox2", "cox3",
    "gag", "pol polyprotein", "pro-pol", "reverse transcriptase",
    "integrase", "retrovir", "endogenous retro", "transposase",
    "phage", "flagell", "elongation factor tu",
]


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


def annotations():
    """{orthogroup: annotation} from the Scoary presence/absence CSV."""
    try:
        g = pd.read_csv(PA_CSV, low_memory=False,
                        usecols=lambda c: c in ("Gene", "Annotation",
                                                "Non-unique Gene name"))
    except Exception:                                             # noqa: BLE001
        return {}
    if "Gene" not in g.columns:
        return {}
    col = "Annotation" if "Annotation" in g.columns else None
    if col is None:
        return {}
    return dict(zip(g.Gene, g[col].astype(object).fillna("")))


ANN = annotations()


def looks_nonviral(og):
    a = str(ANN.get(og, "")).lower()
    return any(w in a for w in NONVIRAL_WORDS)


def finish(ids, source):
    ids = sorted(set(ids))
    pa = load_pa()
    freq = pa.sum(axis=1) / pa.shape[1]

    rows = []
    for og in ids:
        n = int(pa.loc[og].sum()) if og in pa.index else -1
        rows.append(dict(orthogroup=og, in_matrix=og in pa.index,
                         n_genomes=n,
                         pct_genomes=round(100 * n / pa.shape[1], 2)
                         if n >= 0 else None,
                         annotation=ANN.get(og, ""),
                         annotation_looks_nonviral=looks_nonviral(og)))
    cand = pd.DataFrame(rows)
    cand.to_csv(os.path.join(D, "02_results", "nonviral_candidates.csv"),
                index=False)

    say(f"\nsource: {source}")
    say(f"candidates: {len(ids)}")
    say("")
    say(f"  {'orthogroup':<13}{'genomes':>8}  annotation")
    for r in cand.itertuples():
        flag = " *" if r.annotation_looks_nonviral else "  "
        say(f"  {r.orthogroup:<13}{r.n_genomes:>8}{flag}"
            f"{str(r.annotation)[:60]}")
    n_flag = int(cand.annotation_looks_nonviral.sum())
    say(f"\n  * = annotation contains a host, retroviral or bacterial term "
        f"({n_flag} of {len(cand)})")
    if len(cand) and n_flag == 0:
        say("  NOTE: none of the annotations looks non-viral. If they are all")
        say("  'hypothetical protein' that is expected, since contaminant")
        say("  fragments often go unannotated; but it means the annotation")
        say("  gives no independent support for calling them non-viral.")

    missing = [i for i in ids if i not in pa.index]
    if missing:
        say(f"\nWARNING: {len(missing)} not in the presence/absence matrix: "
            f"{', '.join(missing)}")

    kept = [i for i in pa.index if i not in set(ids)]
    say("")
    say("=" * 64)
    say(f"  pan-genome before exclusion : {len(pa)}")
    say(f"  non-viral orthogroups       : {len(ids)}")
    say(f"  pan-genome after exclusion  : {len(kept)}")
    say("=" * 64)
    if len(ids) == EXPECTED and len(kept) == len(pa) - EXPECTED:
        say(f"  matches the manuscript ({len(pa)} - {EXPECTED} = {len(kept)})")
    else:
        say(f"  THE MANUSCRIPT SAYS 206 - {EXPECTED} = 178.")
        say(f"  This run says {len(pa)} - {len(ids)} = {len(kept)}.")
        say("  The manuscript must be changed to match this, not the reverse.")
        say("  Send this output before using the figures.")
    say("=" * 64)

    with open(OUT, "w") as fh:
        fh.write("\n".join(ids) + "\n")
    say(f"\nwritten: {OUT}")
    say(f"written: {os.path.join(D, '02_results', 'nonviral_candidates.csv')}")
    LOG.close()
    sys.exit(0)


# ---- route 1: explicit file -------------------------------------------------
env = os.environ.get("NONVIRAL_LIST")
if env:
    if not os.path.exists(env):
        sys.exit(f"$NONVIRAL_LIST set to {env} but that file does not exist")
    ids = [l.strip() for l in open(env) if l.strip() and not l.startswith("#")]
    finish(ids, f"$NONVIRAL_LIST = {env}")

# ---- route 2: already written -----------------------------------------------
if os.path.exists(OUT) and os.environ.get("FORCE_REDERIVE") != "1":
    ids = [l.strip() for l in open(OUT) if l.strip()]
    say(f"{OUT} already exists; reusing it.")
    say("To re-derive:  export FORCE_REDERIVE=1")
    finish(ids, "existing 01_inputs/nonviral_orthogroups.txt")

# ---- route 3: SRA-exclusive orthogroups, computed here ----------------------
# An orthogroup present only in genomes assembled here from SRA reads, and in
# no GenBank genome, is the signature of contamination carried through
# assembly: a real viral gene would appear in at least one curated genome.
# This is exactly stage 19's sra_exclusive column, and it needs nothing but
# the presence/absence matrix, so it is computed here rather than depending on
# stage 19 having been run.
say("computing SRA-exclusive orthogroups from the presence/absence matrix")

PA = load_pa()
sra_cols = [c for c in PA.columns if is_sra(c)]
gb_cols = [c for c in PA.columns if not is_sra(c)]
say(f"{len(PA.columns)} genomes: {len(sra_cols)} SRA-derived "
    f"(SRR/ERR/DRR), {len(gb_cols)} GenBank")

if not sra_cols or not gb_cols:
    say("\nCannot split the genomes by provenance: one side is empty.")
    say("Genome names do not start with SRR/ERR/DRR, so SRA-derived "
        "assemblies cannot be told apart from GenBank ones by name.")
    say("Falling through to DIAMOND.\n")
else:
    n_sra = PA[sra_cols].sum(axis=1)
    n_gb = PA[gb_cols].sum(axis=1)
    freq = PA.sum(axis=1) / PA.shape[1]
    cls = pd.Series("cloud", index=PA.index)
    cls[freq >= 0.15] = "shell"
    cls[freq >= 0.95] = "soft-core"
    cls[freq >= 0.99] = "core"

    sra_only = PA.index[(n_gb == 0) & (n_sra > 0)].tolist()
    say(f"\n{len(sra_only)} orthogroups are SRA-exclusive "
        f"(present in at least one SRA assembly, absent from every "
        f"GenBank genome)")
    say("\nby pan-genome compartment:")
    for k, v in cls[sra_only].value_counts().items():
        say(f"  {str(k):<12} {int(v)}")
    say("\nThe manuscript states that all 28 fall in the cloud compartment.")

    # cross-check against stage 19 if it happens to be there
    if os.path.exists(PROV):
        try:
            prov = pd.read_csv(PROV, index_col=0)
            if "sra_exclusive" in prov.columns:
                theirs = set(prov.index[prov.sra_exclusive.astype(bool)])
                same = theirs == set(sra_only)
                say(f"\nstage 19 output is present and "
                    f"{'agrees' if same else 'DISAGREES'} "
                    f"({len(theirs)} orthogroups)")
        except Exception as e:                                    # noqa: BLE001
            say(f"\n(could not read stage 19 output for cross-check: {e})")

    rule = os.environ.get("NONVIRAL_RULE", "sra_exclusive")
    if rule == "sra_exclusive_cloud":
        sel = [o for o in sra_only if cls[o] == "cloud"]
        finish(sel, "SRA-exclusive AND cloud "
                    "(NONVIRAL_RULE=sra_exclusive_cloud)")
    say("\nTo restrict to the cloud compartment instead:")
    say("  NONVIRAL_RULE=sra_exclusive_cloud FORCE_REDERIVE=1 \\")
    say("    python3 scripts/01_nonviral_screen.py")
    finish(sra_only, "SRA-exclusive orthogroups, computed from the "
                     "presence/absence matrix")

# ---- route 4: DIAMOND -------------------------------------------------------
if not shutil.which("diamond"):
    sys.exit(
        "diamond is not on the path either.\n\n"
        "Route 3 could not run either, so genome names do not identify\n"
        "which assemblies came from SRA reads. Options:\n"
        "  1. activate an environment that has diamond:\n"
        "       conda activate prokka_env\n"
        "     then re-run this script\n"
        "  2. supply the IDs directly if you have them written down:\n"
        f"       printf 'OG0000xxx\\nOG0000yyy\\n' > {OUT}\n"
        "       python3 scripts/01_nonviral_screen.py")

say("re-deriving with DIAMOND")
refs = []
for d in [os.environ.get("REF_PROTEOMES"),
          f"{ROOT}/00_Metadata/reference_proteomes",
          f"{ROOT}/02_Annotation/reference_proteomes",
          f"{ROOT}/00_Metadata"]:
    if d and os.path.isdir(d):
        refs += glob.glob(os.path.join(d, "*.faa")) + \
                glob.glob(os.path.join(d, "*.fasta"))
refs = [r for r in refs if any(k in os.path.basename(r).upper()
                               for k in ("LSDV", "SPPV", "GTPV", "AF325528",
                                         "NC_004002", "NC_004003"))]
if not refs:
    sys.exit("could not find the LSDV/SPPV/GTPV reference proteomes.\n"
             "Point at the directory holding them:\n"
             "  export REF_PROTEOMES=/path/to/dir\n"
             "or use route 1 or 3 above.")
say("reference proteomes: " + ", ".join(os.path.basename(r) for r in refs))

ref_faa = os.path.join(INP, "capripox_refs.faa")
with open(ref_faa, "w") as out:
    for r in refs:
        out.write(open(r).read())

from Bio import SeqIO                                             # noqa: E402
qry = os.path.join(INP, "og_representatives.faa")
n_og = 0
with open(qry, "w") as out:
    for f in sorted(glob.glob(os.path.join(OG_SEQ_DIR, "OG*.fa"))):
        recs = list(SeqIO.parse(f, "fasta"))
        if not recs:
            continue
        best = max(recs, key=lambda r: len(r.seq))
        out.write(f">{os.path.basename(f).replace('.fa', '')}\n{best.seq}\n")
        n_og += 1
say(f"{n_og} orthogroup representatives")

db = os.path.join(INP, "capripox_refs")
subprocess.run(["diamond", "makedb", "--in", ref_faa, "-d", db],
               check=True, capture_output=True)
hits = os.path.join(D, "02_results", "og_vs_capripox_refs.tsv")
subprocess.run(["diamond", "blastp", "-q", qry, "-d", db, "-o", hits,
                "--evalue", "1e-5", "--max-target-seqs", "1", "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue"],
               check=True)
hit_ogs = {l.split("\t")[0] for l in open(hits) if l.strip()}
all_ogs = {os.path.basename(f).replace(".fa", "")
           for f in glob.glob(os.path.join(OG_SEQ_DIR, "OG*.fa"))}
finish(sorted(all_ogs - hit_ogs),
       "DIAMOND: no hit against the capripoxvirus reference proteomes")
