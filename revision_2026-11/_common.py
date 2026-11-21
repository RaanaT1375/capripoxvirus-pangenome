#!/usr/bin/env python3
"""Shared paths and helpers for the revision analyses (stages 17-21).

All stage scripts import this. The analysis root is taken from $PANG,
defaulting to the cluster path, so the same code runs unchanged on a
laptop by exporting PANG=/path/to/Pangenome_Analysis.
"""
import os
import shutil
import subprocess

ROOT = os.environ.get(
    "PANG", "/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis")
REV = os.path.join(ROOT, "revision_2026-11")

# ---- canonical inputs produced by stages 01-16 ------------------------------
SPECIES = f"{ROOT}/00_Metadata/species_assignment.csv"
TRAITS = f"{ROOT}/11_Scoary/01_inputs/traits.csv"
PA_CSV = f"{ROOT}/11_Scoary/01_inputs/gene_presence_absence.csv"
CODON_DIR = f"{ROOT}/10_Selection_Pressure/01_codon_alignments"
TREE_DIR = f"{ROOT}/10_Selection_Pressure/02_pruned_trees"
CORE_TREE = f"{ROOT}/05_Phylogeny/04_tree/core_genome_ML_tree_v3.treefile"
CORE_CONTREE = f"{ROOT}/05_Phylogeny/04_tree/core_genome_ML_tree_v3.contree"
ANNOT_60 = f"{ROOT}/10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv"
ANNOT_ALL = f"{ROOT}/00_Metadata/species_specific_genes.csv"
OG_SEQ_DIR = (f"{ROOT}/03_OrthoFinder_Results/OrthoFinder_Results_v3/"
              f"Results_Aug04/Orthogroup_Sequences")
GENE_COUNT = (f"{ROOT}/03_OrthoFinder_Results/OrthoFinder_Results_v3/"
              f"Results_Aug04/Orthogroups/Orthogroups.GeneCount.tsv")
GENOME_DIR = f"{ROOT}/06_Recombination/00_nucleotide_genomes"
MK_SCRIPTS = f"{ROOT}/14_MK_Test/scripts"

META_COLS = [
    'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates',
    'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
    'Order within Fragment', 'Accessory Fragment',
    'Accessory Order with Fragment', 'QC', 'Min group size nuc',
    'Max group size nuc', 'Avg group size nuc']


def iqtree_bin():
    """Locate an IQ-TREE 2.x/3.x binary; $IQTREE overrides the search."""
    cand = [os.environ.get("IQTREE")] if os.environ.get("IQTREE") else []
    cand += ["iqtree2", "iqtree3", "iqtree"]
    for c in cand:
        p = shutil.which(c) if c and not os.path.sep in c else c
        if p and os.path.exists(p):
            return p
    raise SystemExit(
        "IQ-TREE not found. Load it first, e.g.\n"
        "  conda activate iqtree_env\n"
        "or point at it explicitly:  export IQTREE=/path/to/iqtree2")


def run_iqtree(args, log=None):
    r = subprocess.run([iqtree_bin()] + args, capture_output=True, text=True)
    if log:
        with open(log, "a") as fh:
            fh.write(r.stdout + "\n" + r.stderr + "\n")
    return r


def load_pa():
    """Orthogroup x genome 0/1 presence matrix from the Scoary input."""
    import pandas as pd
    g = pd.read_csv(PA_CSV, low_memory=False)
    gen = [c for c in g.columns if c not in META_COLS]
    return g.set_index("Gene")[gen].astype(int)


def groups():
    """{'LSDV': set(...), 'SPPV': set(...), 'GTPV': set(...)}"""
    import pandas as pd
    sp = pd.read_csv(SPECIES)
    return {s: set(sp.loc[sp.species == s, "Name"]) for s in
            ["LSDV", "SPPV", "GTPV"]}


def vaccine_ids():
    import pandas as pd
    return set(pd.read_csv(TRAITS, index_col=0).query("Vaccine == 1").index)


def stage_dir(stage):
    d = os.path.join(REV, stage)
    for sub in ("01_inputs", "02_results", "03_summary", "logs"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    return d


def is_sra(name):
    return str(name).startswith(("SRR", "ERR", "DRR"))
