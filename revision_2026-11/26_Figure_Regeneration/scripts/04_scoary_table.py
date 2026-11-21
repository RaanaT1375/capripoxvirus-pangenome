#!/usr/bin/env python3
"""Stage 26, step 4 — parse the Scoary output and extract the Table 5 rows.

Two things this has to get right, both of which the first version got wrong.

1. ONE FILE PER TRAIT. Scoary writes <Trait>_<DD_MM_YYYY_HHMMSS>.results.csv,
   so a directory that has been run more than once holds several files per
   trait and naive globbing counts every association twice. By default this
   keeps only the newest file per trait and reports what it dropped.
   Override with SCOARY_KEEP=all to see the undeduplicated numbers.

2. A TABLE ROW IS A (GENE, TRAIT) PAIR, NOT A GENE. The same orthogroup is
   associated with several traits, so looking one up by gene alone returns
   whichever trait happens to have the smallest p. Every row wanted below is
   named by both.

Outputs
  02_results/scoary_all_parsed.csv     every row from the files actually used
  02_results/scoary_files_used.csv     which file was used for each trait
  03_summary/scoary_significant.csv    rows passing both filters
  03_summary/scoary_convergent.csv     rows also meeting the convergence rule
  03_summary/table5_rows.csv           the Table 5 rows
  03_summary/04_scoary_table.log
"""
import os
import re
import sys
import glob


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
from _common import ROOT, REV, stage_dir                          # noqa: E402

import pandas as pd                                               # noqa: E402

D = stage_dir("26_Figure_Regeneration")
LOG = open(os.path.join(D, "03_summary", "04_scoary_table.log"), "w")

BONF, EMP, MINPAIRS = 0.05, 0.05, 3
KEEP = os.environ.get("SCOARY_KEEP", "newest").lower()

# (orthogroup, trait) as the revised manuscript's Table 5 names them.
WANT = [("OG0000157", "Asia"),
        ("OG0000155", "Asia"),
        ("OG0000154", "Europe"),
        ("OG0000159", "Africa"),
        ("OG0000161", "Africa"),
        ("OG0000156", "Africa"),
        ("OG0000152", "Vaccine"),
        # reported in the text but not in Table 5; kept here so their
        # values are printed and can be checked
        ("OG0000151", "Vaccine"),
        ("OG0000004", "Wild")]

STAMP = re.compile(r"_(\d{2}_\d{2}_\d{4}(?:_\d+)?)\.results\.csv$")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")


# Paths matching these are never preferred, whatever their date: a
# superseded or no-tree run must not win a tie on modification time alone.
DEPRECATED = ("superseded", "notree", "no_tree", "_old", "backup", "archive")


def rank(path):
    """2 = preferred, 1 = ordinary, 0 = looks superseded."""
    low = path.lower()
    if any(k in low for k in DEPRECATED):
        return 0
    return 2 if "tree" in os.path.basename(os.path.dirname(low)) else 1


def sortkey(stamp, path):
    """Preference tier first, then DD_MM_YYYY[_HHMMSS], then mtime."""
    r = rank(path)
    if stamp:
        p = stamp.split("_")
        try:
            return (r, 1, p[2], p[1], p[0], p[3] if len(p) > 3 else "")
        except IndexError:
            pass
    return (r, 0, f"{os.path.getmtime(path):020.3f}", "", "", "")


dirs = os.environ.get("SCOARY_DIRS")
dirs = dirs.split(":") if dirs else [f"{ROOT}/11_Scoary",
                                     f"{REV}/22_Scoary_Vaccine"]
files = []
for d in dirs:
    if os.path.isdir(d):
        files += glob.glob(os.path.join(d, "**", "*.results.csv"),
                           recursive=True)
    else:
        say(f"note: {d} is not a directory, skipped")
files = sorted(set(files))
if not files:
    sys.exit("no Scoary *.results.csv found under:\n  " + "\n  ".join(dirs) +
             "\nSet SCOARY_DIRS to the directory holding them.")
say(f"{len(files)} Scoary result files found")

# ---- one file per trait -----------------------------------------------------
by_trait = {}
for f in files:
    base = os.path.basename(f)
    m = STAMP.search(base)
    stamp = m.group(1) if m else None
    trait = base[:m.start()] if m else base.replace(".results.csv", "")
    run = "sensitivity_280" if REV in f else "primary_290"
    by_trait.setdefault((run, trait), []).append((sortkey(stamp, f), stamp, f))

dupes = {k: v for k, v in by_trait.items() if len(v) > 1}
if dupes:
    say(f"\n{len(dupes)} trait(s) have more than one result file:")
    for (run, trait), v in sorted(dupes.items()):
        say(f"  {run} / {trait}")
        for key, stamp, f in sorted(v, reverse=True):
            mark = "KEEP" if (KEEP != "all" and
                              (key, stamp, f) == max(v)) else "    "
            flag = "" if rank(f) else "   <- looks superseded"
            say(f"    {mark} {stamp or '(no timestamp)'}  "
                f"{os.path.relpath(f, ROOT) if f.startswith(ROOT) else f}"
                f"{flag}")
    if KEEP == "all":
        say("\nSCOARY_KEEP=all: using every file, so associations present in "
            "more than one run are counted more than once.")
    else:
        say("\nUsing the newest file per trait. To see the undeduplicated "
            "numbers instead:  SCOARY_KEEP=all python3 scripts/04_scoary_table.py")

used = []
for (run, trait), v in by_trait.items():
    chosen = v if KEEP == "all" else [max(v)]
    for key, stamp, f in chosen:
        used.append(dict(run=run, trait=trait, stamp=stamp, path=f))
pd.DataFrame(used).to_csv(
    os.path.join(D, "02_results", "scoary_files_used.csv"), index=False)
say(f"\n{len(used)} files used, covering "
    f"{len({(u['run'], u['trait']) for u in used})} trait/run combinations")

frames = []
for u in used:
    try:
        t = pd.read_csv(u["path"])
    except Exception as e:                                        # noqa: BLE001
        say(f"  skipped {os.path.basename(u['path'])}: {e}")
        continue
    t["Trait"] = u["trait"]
    t["run"] = u["run"]
    t["source_file"] = (os.path.relpath(u["path"], ROOT)
                        if u["path"].startswith(ROOT) else u["path"])
    frames.append(t)
if not frames:
    sys.exit("no Scoary file could be parsed")
allr = pd.concat(frames, ignore_index=True)
allr.to_csv(os.path.join(D, "02_results", "scoary_all_parsed.csv"), index=False)
say(f"{len(allr)} rows parsed across {allr.Trait.nunique()} traits")


def col(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


c_bonf = col(allr, "Bonferroni_p", "Bonferroni_p_value")
c_emp = col(allr, "Empirical_p", "Empirical_p_value")
c_or = col(allr, "Odds_ratio", "Odds_Ratio")
c_sup = col(allr, "Max_supporting_pairs", "Best_supporting_pairs")
c_opp = col(allr, "Max_opposing_pairs", "Worst_opposing_pairs")
miss = [n for n, c in [("Bonferroni", c_bonf), ("Empirical", c_emp),
                       ("Odds ratio", c_or), ("supporting", c_sup),
                       ("opposing", c_opp)] if c is None]
if miss:
    say(f"WARNING: columns not found: {', '.join(miss)}")
    say("columns present: " + ", ".join(allr.columns))

# a (gene, trait) pair should now appear once; say so if it does not
dup_pairs = allr.duplicated(subset=["Gene", "Trait", "run"]).sum()
if dup_pairs:
    say(f"\nWARNING: {dup_pairs} duplicated (gene, trait) rows remain after "
        f"deduplication. Inspect 02_results/scoary_files_used.csv.")

sig = allr[(allr[c_bonf] < BONF) & (allr[c_emp] < EMP)].copy()
say(f"\npassing both filters (Bonferroni < {BONF}, empirical < {EMP}): "
    f"{len(sig)} gene-trait pairs")
for r, g in sig.groupby("run"):
    say(f"  {r}: {len(g)} pairs, {g.Trait.nunique()} traits, "
        f"{g.Gene.nunique()} orthogroups")
say("manuscript states: 92 pairs, 16 traits, 18 orthogroups "
    "(primary 290-genome run)")
sig.to_csv(os.path.join(D, "03_summary", "scoary_significant.csv"), index=False)

sig["net_support"] = sig[c_sup] - sig[c_opp]
conv = sig[(sig[c_sup] >= MINPAIRS) & (sig.net_support > 0)].copy()
conv = conv.sort_values(["run", "Trait", c_emp])
conv.to_csv(os.path.join(D, "03_summary", "scoary_convergent.csv"),
            index=False)
say(f"\nmeeting the convergence criterion (>= {MINPAIRS} supporting pairs, "
    f"net positive): {len(conv)}")
for _, r in conv.iterrows():
    say(f"  {r['run']:<16} {r['Trait']:<22} {r['Gene']:<12} "
        f"OR={r[c_or]:<10.4g} emp_p={r[c_emp]:<9.4g} "
        f"{int(r[c_sup])}:{int(r[c_opp])}")

# ---- the Table 5 rows -------------------------------------------------------
say("\n" + "=" * 70)
say("TABLE 5 ROWS — looked up by (orthogroup, trait), not by orthogroup")
say("=" * 70)
rows = []
for og, trait in WANT:
    hit = sig[(sig.Gene == og) & (sig.Trait.str.lower() == trait.lower())]
    status = "passes both filters"
    if hit.empty:
        hit = allr[(allr.Gene == og) & (allr.Trait.str.lower() == trait.lower())]
        status = "present but does NOT pass both filters"
    if hit.empty:
        other = sorted(allr.Trait[allr.Gene == og].unique())
        say(f"{og:<12} {trait:<22} NOT FOUND for this trait."
            + (f" Present for: {', '.join(other)}" if other else
               " Orthogroup absent from the output entirely."))
        continue
    r = hit.sort_values(c_emp).iloc[0]
    orv = r[c_or]
    or_tex = ("$\\infty$" if (pd.isna(orv) or orv == float("inf") or orv > 1e6)
              else f"{orv:.2f}")
    conv_ok = (r[c_sup] >= MINPAIRS) and (r[c_sup] - r[c_opp] > 0)
    rows.append(dict(orthogroup=og, trait=trait, run=r["run"],
                     odds_ratio=orv, odds_ratio_tex=or_tex,
                     empirical_p=r[c_emp], bonferroni_p=r[c_bonf],
                     supporting=int(r[c_sup]), opposing=int(r[c_opp]),
                     meets_convergence=conv_ok, status=status))
    say(f"{og:<12} {trait:<22} OR={or_tex:<10} emp_p={r[c_emp]:<9.4g} "
        f"pairs={int(r[c_sup])}:{int(r[c_opp])}  "
        f"{'convergent' if conv_ok else 'NOT convergent'}  [{status}]")

if rows:
    pd.DataFrame(rows).to_csv(
        os.path.join(D, "03_summary", "table5_rows.csv"), index=False)
    say("\nwritten: 03_summary/table5_rows.csv")

say("\nCheck before using any of this:")
say("  - the pair/trait/orthogroup counts above against 92 / 16 / 18")
say("  - that every Table 5 row says 'convergent'")
say("  - 03_summary/scoary_convergent.csv is the authoritative tier-1 list;")
say("    if it disagrees with Table 5 in the manuscript, the table is wrong.")
LOG.close()
