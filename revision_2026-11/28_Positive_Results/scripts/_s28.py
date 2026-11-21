#!/usr/bin/env python3
"""Shared helpers for stage 28."""
import os
import sys
import glob


def locate_common():
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


sys.path.insert(0, locate_common())
from _common import ROOT, stage_dir                               # noqa: E402,F401

STAGE = "28_Positive_Results"


def opener(name):
    """Return (say, close) writing both to stdout and to a log."""
    d = stage_dir(STAGE)
    fh = open(os.path.join(d, "03_summary", name), "w")

    def say(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        fh.write(s + "\n")

    return d, say, fh.close


# ---- functional categories --------------------------------------------------
# Keyword rules applied to the product annotation, in order. Override the whole
# classification by placing a CSV with columns orthogroup,category in
# 01_inputs/functional_categories.csv.
CATEGORY_RULES = [
    ("Transcription",      ["rna polymerase", "transcription", "poly(a)",
                            "capping", "mrna", "elongation", "termination",
                            "nucleoside triphosphatase"]),
    ("DNA replication",    ["dna polymerase", "helicase", "primase", "ligase",
                            "topoisomerase", "single-strand", "ssdna",
                            "uracil", "thymidine", "ribonucleotide reductase",
                            "nucleotide"]),
    ("Virion assembly",    ["imv", "virion", "core protein", "scaffold",
                            "morphogenesis", "maturation", "protease",
                            "myristyl", "packaging", "telomere"]),
    ("Membrane / EEV",     ["eev", "membrane", "envelope", "glycoprotein",
                            "fusion", "entry", "hemagglutinin"]),
    ("Host interaction",   ["host range", "ankyrin", "kelch", "interleukin",
                            "interferon", "apoptosis", "ox-2", "cd200",
                            "cd47", "chemokine", "tnf", "serpin",
                            "immunomod",
                            # abbreviations and families the first pass missed
                            "il-1", "il-10", "il-18", "ifn", "pkr",
                            "eif2a", "n1r", "p28-like", "lap/phd",
                            "growth factor", "egf", "binding protein",
                            "complement", "mhc", "nf-kb", "caspase",
                            "bcl-2", "steroid"]),
    ("Redox / metabolism", ["dismutase", "glutaredoxin", "thioredoxin",
                            "kinase", "phosphatase", "metabolic",
                            "s-s bond", "disulphide", "disulfide",
                            "superoxide"]),
]


def classify(product):
    p = str(product).lower()
    for cat, words in CATEGORY_RULES:
        if any(w in p for w in words):
            return cat
    return "Unassigned"


def find_selection_table():
    """Locate a per-gene selection summary with an omega column."""
    import pandas as pd
    pats = [f"{ROOT}/10_Selection_Pressure/**/*.csv",
            f"{ROOT}/13_Selection_LSDV/**/*.csv",
            f"{ROOT}/99_Supplementary/**/*.csv"]
    seen = []
    for pat in pats:
        for f in sorted(glob.glob(pat, recursive=True)):
            try:
                t = pd.read_csv(f, nrows=5)
            except Exception:                                     # noqa: BLE001
                continue
            cols = {c.lower(): c for c in t.columns}
            og = next((cols[c] for c in cols
                       if c in ("orthogroup", "gene", "og", "locus")), None)
            om = next((cols[c] for c in cols
                       if "mg94" in c or c in ("omega", "omega_mg94",
                                               "global_omega", "dn_ds",
                                               "dnds")), None)
            if og and om:
                return f, og, om
            seen.append(f)
    return None, None, seen
