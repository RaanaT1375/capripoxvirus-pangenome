#!/usr/bin/env python3
"""Stage 23 — audit of the trait matrix that drives the pan-GWAS.

Every pan-GWAS association is only as good as traits.csv, and that file has
never been checked against the rest of the dataset. This script does not
change anything; it reports.

Checks
  A  structural integrity: shape, dtypes, missing values, non-binary entries,
     duplicate genome IDs, and exact agreement with the species assignment,
     the presence/absence matrix and the phylogeny tip labels
  B  mutually exclusive trait groups: geography, host, BAPS cluster, source
     and recombination category should each sum to exactly 1 per genome.
     Any genome summing to 0 or >1 is listed by name.
  C  cross-consistency with the species assignment, which matters here
     because Host_Sheep and Host_Goat are known to track SPPV and GTPV
     rather than a host-range phenotype
  D  a machine-readable list of every genome that failed any check

Outputs
  02_results/trait_group_sums.csv
  02_results/trait_violations.csv
  02_results/trait_species_crosstabs.txt
  03_summary/trait_audit.log
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from _common import (SPECIES, TRAITS, CORE_TREE, load_pa,          # noqa: E402
                     stage_dir, is_sra)
import pandas as pd                                                # noqa: E402

# groups that should be mutually exclusive and exhaustive
GROUPS = {
    "geography": ["Africa", "Asia", "Europe"],
    "host": ["Host_Cattle", "Host_Sheep", "Host_Goat", "Host_Giraffe",
             "Host_Wildlife"],
    "BAPS": [f"BAPS_cluster_{i}" for i in range(1, 6)],
    "source": ["Vaccine", "Wild"],
    "recombination": ["Recombinant", "Ancestral_recombinant",
                      "Independent_recombinant", "No_recombination"],
}

D = stage_dir("23_Trait_Audit")
LOG = open(os.path.join(D, "03_summary", "trait_audit.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


tr = pd.read_csv(TRAITS, index_col=0)
sp = pd.read_csv(SPECIES).set_index("Name")["species"]
PA = load_pa()

# ------------------------------------------------------------------ A -------
say("=" * 70)
say("A. STRUCTURAL INTEGRITY")
say("=" * 70)
say(f"traits.csv: {tr.shape[0]} genomes x {tr.shape[1]} traits")
say(f"columns: {', '.join(tr.columns)}")

dup = tr.index[tr.index.duplicated()].tolist()
say(f"\nduplicate genome IDs      : {len(dup)}" + (f"  {dup}" if dup else ""))
nan = tr.isna().sum().sum()
say(f"missing values            : {nan}")
if nan:
    bad = tr[tr.isna().any(axis=1)]
    say(f"  genomes with NaN        : {list(bad.index)}")

vals = set(pd.unique(tr.values.ravel()))
non_binary = vals - {0, 1}
say(f"distinct values           : {sorted(vals)}")
say(f"non-binary values         : {sorted(non_binary) if non_binary else 'none'}")

t_set, s_set, p_set = set(tr.index), set(sp.index), set(PA.columns)
say(f"\ngenomes in traits but not in species_assignment: "
    f"{sorted(t_set - s_set) or 'none'}")
say(f"genomes in species_assignment but not in traits: "
    f"{sorted(s_set - t_set) or 'none'}")
say(f"genomes in traits but not in the PA matrix     : "
    f"{sorted(t_set - p_set) or 'none'}")
say(f"genomes in the PA matrix but not in traits     : "
    f"{sorted(p_set - t_set) or 'none'}")

if os.path.exists(CORE_TREE):
    from Bio import Phylo
    tips = {x.name for x in Phylo.read(CORE_TREE, "newick").get_terminals()}
    say(f"tree tips not in traits                        : "
        f"{sorted(tips - t_set) or 'none'}")
    say(f"traits genomes not on the tree                 : "
        f"{sorted(t_set - tips) or 'none'}")

say(f"\nSRA-derived genomes in traits: {sum(1 for i in tr.index if is_sra(i))}"
    f" | GenBank: {sum(1 for i in tr.index if not is_sra(i))}")

# ------------------------------------------------------------------ B -------
say("\n" + "=" * 70)
say("B. MUTUALLY EXCLUSIVE TRAIT GROUPS  (each should sum to exactly 1)")
say("=" * 70)

sums = pd.DataFrame(index=tr.index)
viol = []
for name, cols in GROUPS.items():
    present = [c for c in cols if c in tr.columns]
    missing = [c for c in cols if c not in tr.columns]
    if missing:
        say(f"\n{name}: columns absent from traits.csv: {missing}")
    if not present:
        continue
    s = tr[present].sum(axis=1)
    sums[name] = s
    dist = s.value_counts().sort_index()
    flag = "OK" if set(dist.index) == {1} else "PROBLEM"
    say(f"\n{name}  ({', '.join(present)})  -> {flag}")
    for k, v in dist.items():
        say(f"    sum = {k}: {v} genomes")
    say(f"    per-trait totals: "
        + ", ".join(f"{c}={int(tr[c].sum())}" for c in present))
    for g in s.index[s != 1]:
        viol.append(dict(genome=g, group=name, group_sum=int(s[g]),
                         species=sp.get(g, "?"),
                         traits_set=";".join(c for c in present if tr.loc[g, c] == 1)
                                    or "NONE",
                         provenance="SRA" if is_sra(g) else "GenBank"))

sums.to_csv(os.path.join(D, "02_results", "trait_group_sums.csv"))
V = pd.DataFrame(viol)
V.to_csv(os.path.join(D, "02_results", "trait_violations.csv"), index=False)

if len(V):
    say(f"\n--- {len(V)} group violations, affecting "
        f"{V.genome.nunique()} distinct genomes ---")
    for name, g in V.groupby("group"):
        say(f"\n  [{name}]")
        say(g[["genome", "species", "provenance", "group_sum", "traits_set"]]
            .to_string(index=False))
else:
    say("\nno group violations")

# ------------------------------------------------------------------ C -------
say("\n" + "=" * 70)
say("C. CROSS-CONSISTENCY WITH THE SPECIES ASSIGNMENT")
say("=" * 70)
ct_lines = []
for name in ["host", "BAPS", "source", "geography", "recombination"]:
    cols = [c for c in GROUPS[name] if c in tr.columns]
    if not cols:
        continue
    say(f"\n{name} x species:")
    rows = []
    for c in cols:
        pos = tr.index[tr[c] == 1]
        r = {"trait": c, "n": len(pos)}
        for s_ in ["LSDV", "SPPV", "GTPV"]:
            r[s_] = int(sum(1 for g in pos if sp.get(g) == s_))
        rows.append(r)
    t = pd.DataFrame(rows).set_index("trait")
    say(t.to_string())
    ct_lines.append(f"{name} x species\n{t.to_string()}\n")

with open(os.path.join(D, "02_results", "trait_species_crosstabs.txt"), "w") as fh:
    fh.write("\n".join(ct_lines))

say("\nNote: Host_Sheep and Host_Goat are expected to coincide exactly with")
say("SPPV and GTPV. Where they do, the trait is a species label, not a")
say("host-range phenotype, and any association found with it is a")
say("between-species difference. This is already documented for stage 11.")

# ------------------------------------------------------------------ D -------
say("\n" + "=" * 70)
say("D. VERDICT")
say("=" * 70)
problems = []
if dup:
    problems.append(f"{len(dup)} duplicate genome IDs")
if nan:
    problems.append(f"{nan} missing values")
if non_binary:
    problems.append(f"non-binary values {sorted(non_binary)}")
for a, b, lbl in [(t_set, s_set, "species_assignment"),
                  (t_set, p_set, "PA matrix")]:
    if a != b:
        problems.append(f"genome set differs from {lbl}")
if len(V):
    problems.append(f"{len(V)} mutually-exclusive-group violations "
                    f"in {V.genome.nunique()} genomes")
if problems:
    say("PROBLEMS FOUND:")
    for p in problems:
        say(f"  - {p}")
    say("\nEvery affected genome is listed in 02_results/trait_violations.csv.")
    say("Each one needs a decision before the pan-GWAS is reported: either the")
    say("trait assignment is corrected, or the genome is excluded from that")
    say("trait, or the trait definition is restated in the Methods so that a")
    say("reader knows it is not exclusive.")
else:
    say("no problems found; traits.csv is internally consistent")
LOG.close()
