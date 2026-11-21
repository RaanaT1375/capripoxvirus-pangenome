#!/usr/bin/env python3
"""Stage 21, step 1 — audit of the core orthogroups excluded from stages 05-17.

Reviewer 2, comment 2: restricting phylogeny and selection analyses to the 60
single-copy core orthogroups may bias the result toward conserved genomic
components, because host-interaction, virulence and immunomodulatory genes
are the ones most likely to be paralogous or presence/absence variable and
therefore excluded.

This builds the supplementary table that lets a reader check that directly.
Every core orthogroup is listed with the reason it was or was not carried
into the phylogeny, and the decomposition is exact:

  core (present in >= 99% of genomes)                    132
    - present in all 290 but not single copy everywhere   51  -> multi-copy
    - not present in all 290                              21  -> partial PAV
    = single-copy and universal                           60  -> analysed

Outputs
  02_results/core_orthogroup_audit.csv       all 132, with reason + annotation
  03_summary/S_excluded_core_orthogroups.csv the 72 excluded (supplementary)
  03_summary/core_audit.log
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (ANNOT_ALL, ANNOT_60, GENE_COUNT, load_pa,    # noqa: E402
                     stage_dir)
import pandas as pd                                               # noqa: E402

# functional classes the reviewer explicitly names
KEYWORDS = ("kelch", "ankyrin", "IL-1", "interferon", "IFN", "chemokine",
            "OX-2", "host range", "immun", "semaphorin", "CD47", "EEV",
            "TNF", "serpin", "apoptosis", "virulence")

D = stage_dir("21_Core_Audit")
LOG = open(os.path.join(D, "03_summary", "core_audit.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


PA = load_pa()
N = PA.shape[1]
freq = PA.sum(axis=1)
ann = pd.read_csv(ANNOT_ALL).set_index("Orthogroup")
analysed = set(pd.read_csv(ANNOT_60)["Orthogroup"])

core = sorted(freq[freq >= 0.99 * N].index)
strict = set(freq[freq == N].index)
say(f"genomes: {N}")
say(f"core (>= 99%)            : {len(core)}")
say(f"strict core (all {N})    : {len(strict)}")
say(f"single-copy and universal: {len(analysed)}")
say(f"excluded from phylogeny  : {len(set(core) - analysed)}\n")

# ---- copy numbers, if the OrthoFinder gene-count table is available --------
gc = None
if os.path.exists(GENE_COUNT):
    gc = pd.read_csv(GENE_COUNT, sep="\t", index_col=0)
    gc = gc.drop(columns=[c for c in gc.columns if c.lower() == "total"],
                 errors="ignore")
    say(f"copy numbers from {os.path.basename(GENE_COUNT)}")
else:
    say(f"NOTE: {GENE_COUNT} not found - copy-number columns will be blank")

rows = []
for og in core:
    in_all = og in strict
    kept = og in analysed
    if kept:
        reason = "analysed (single-copy, universal)"
    elif not in_all:
        reason = "excluded: partial presence/absence"
    else:
        reason = "excluded: multi-copy in >= 1 genome"
    r = dict(Orthogroup=og, n_present=int(freq[og]),
             n_absent=int(N - freq[og]),
             present_in_all=in_all, analysed=kept, reason=reason,
             product=str(ann["product"].get(og, "")))
    if gc is not None and og in gc.index:
        v = gc.loc[og]
        r.update(max_copies=int(v.max()),
                 n_genomes_multicopy=int((v > 1).sum()),
                 n_genomes_zero=int((v == 0).sum()))
    r["host_interaction_keyword"] = any(
        k.lower() in r["product"].lower() for k in KEYWORDS)
    rows.append(r)

t = pd.DataFrame(rows).set_index("Orthogroup")
t.to_csv(os.path.join(D, "02_results", "core_orthogroup_audit.csv"))

excl = t[~t.analysed].copy()
excl.to_csv(os.path.join(D, "03_summary", "S_excluded_core_orthogroups.csv"))

say("=== reason breakdown for the excluded orthogroups ===")
say(excl.reason.value_counts().to_string())
say("\n=== partial-PAV orthogroups: how many genomes are they missing from? ===")
pav = excl[~excl.present_in_all]
say(pav.n_absent.value_counts().sort_index().to_string())
say(f"  (max {int(pav.n_absent.max())} of {N} genomes)")

flag = excl[excl.host_interaction_keyword]
say(f"\n=== host-interaction / immunomodulatory genes among the excluded: "
    f"{len(flag)} ===")
show = ["n_present", "reason"] + (
    ["max_copies", "n_genomes_multicopy"] if gc is not None else [])
say(flag[show + ["product"]].assign(
    product=flag["product"].str.slice(0, 55)).to_string())
say("\nThis is the reviewer's concern made explicit: these genes ARE excluded,")
say("and the Discussion must say so rather than imply the 60 are representative.")
say("\nwritten: 03_summary/S_excluded_core_orthogroups.csv")
LOG.close()
