#!/usr/bin/env python3
"""
prepare_itol_datasets.py

Filters and matches iTOL annotation datasets to the exact tip labels of the
final core-genome phylogenetic tree, and builds two additional
professional-grade datasets from Supplementary_File1.xlsx:
  - Viral species colorstrip (LSDV / SPPV / GTPV) -- distinct from host animal
  - Clean leaf labels (Isolate name instead of raw accession)

WHY THIS IS NEEDED
    1) The tree contains only a subset of the original accessions (some were
       excluded during QC). The raw itol_*.txt files still list every
       original accession and must be filtered down, or iTOL will show a
       warning / ignore the extra rows.
    2) Several tree tips are "_Merged" labels representing TWO original SRA
       run accessions (e.g. SRR19090746_SRR19090747_Merged). The raw
       annotation files list each run accession separately, so a direct
       string match against the tree tip fails. This script detects that
       pattern and resolves it automatically.

TWO WAYS TO SUPPLY TIP NAMES (pick one):
    --tree              Use the final Newick/contree tree file directly.
    --tips-from-fasta    Use the concatenated supermatrix FASTA instead --
                         its headers are already the exact future tip
                         labels, so this lets you prepare every dataset
                         WHILE IQ-TREE is still running, no need to wait.

USAGE (run on the cluster):
    # after the tree is done:
    python prepare_itol_datasets.py \
        --tree /path/to/core_genome_ML_tree.contree \
        --supplementary /path/to/Supplementary_File1.xlsx \
        --itol-dir /path/to/folder/with/itol_*.txt \
        --outdir /path/to/output

    # or, while IQ-TREE is still running:
    python prepare_itol_datasets.py \
        --tips-from-fasta /path/to/03_supermatrix/supermatrix.fasta \
        --supplementary /path/to/Supplementary_File1.xlsx \
        --itol-dir /path/to/folder/with/itol_*.txt \
        --outdir /path/to/output
"""
import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import Phylo


def load_tree_tips(tree_path):
    tree = Phylo.read(tree_path, "newick")
    return [t.name for t in tree.get_terminals()]


def load_tips_from_fasta(fasta_path):
    """Alternative tip source: read genome/tip IDs directly from the
    concatenated supermatrix FASTA (03_supermatrix/supermatrix.fasta).
    Useful to prepare datasets in parallel while IQ-TREE is still running,
    since the supermatrix headers are exactly the tree's future tip labels."""
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip())
    return ids


def parse_merged_tip(tip):
    """Return the list of component SRA run accessions if `tip` is a
    '<run1>_<run2>_Merged' tree label, else return [tip] unchanged."""
    m = re.match(r"^(SRR\d+)_(SRR\d+)_Merged$", tip)
    if m:
        return [m.group(1), m.group(2)]
    return [tip]


# ==============================================================================
# HARMONIZED ACADEMIC PALETTE
# One coherent, muted "editorial" color family shared across every ring, so the
# four rings read as a single designed figure instead of four independently
# picked color sets. Old hex (from the original uploaded files) -> new hex.
# Only colorstrip datasets (condition, continent, host) need remapping --
# country_labels is plain black text, viral_species is generated fresh below.
# ==============================================================================
COLOR_REMAP = {
    "itol_condition.txt": {
        "#7d6592": "#6B5876",  # Recombinant -> muted plum
        "#e57b89": "#B8697D",  # Vaccine     -> dusty rose
        "#2f4858": "#3A4750",  # Wild        -> charcoal slate
    },
    "itol_continent.txt": {
        "#d98a01": "#C97B3D",  # Africa -> warm terracotta
        "#4c8c3c": "#3E7C59",  # Asia   -> muted jade
        "#1b6299": "#3B5B8C",  # Europe -> steel blue
    },
    "itol_host.txt": {
        "#1e3900": "#5B4B8A",  # Antidorcas marsupialis -> muted plum (rare host)
        "#3D5A80": "#B98B4E",  # Cattle   -> ochre / sand
        "#4d691a": "#D9A441",  # Giraffe  -> golden amber
        "#81B29A": "#7C9885",  # Goat     -> sage
        "#E07A5F": "#A65B54",  # Sheep    -> clay / terracotta-rose
    },
}

VIRAL_SPECIES_COLORS = {
    "LSDV": "#C1440E",     # burnt vermillion
    "SPPV": "#1D6F8C",     # deep teal-blue
    "GTPV": "#4C7A3D",     # forest olive-green
    "Unknown": "#8C8579",  # warm taupe grey
}


def remap_row_color(src_name, row):
    """Replace the color field (row[0]) of a colorstrip DATA row using
    COLOR_REMAP, if this source file has a remap table."""
    remap = COLOR_REMAP.get(src_name)
    if not remap or not row:
        return row
    new_row = list(row)
    new_row[0] = remap.get(row[0], row[0])
    return new_row


def remap_header_legend_colors(src_name, header_lines):
    """Rewrite the LEGEND_COLORS line (if present) using COLOR_REMAP so the
    on-figure legend matches the recolored data rows."""
    remap = COLOR_REMAP.get(src_name)
    if not remap:
        return header_lines
    new_lines = []
    for line in header_lines:
        if line.startswith("LEGEND_COLORS"):
            parts = line.split("\t")
            parts = [parts[0]] + [remap.get(p, p) for p in parts[1:]]
            line = "\t".join(parts)
        new_lines.append(line)
    return new_lines


def load_itol_dataset(path):
    """Parse an existing iTOL annotation file (any DATASET_* type).
    Returns (header_lines_including_DATA, {leaf_id: [field1, field2, ...]})."""
    header_lines = []
    data = {}
    in_data = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.strip() == "DATA":
                header_lines.append(line)
                in_data = True
                continue
            if not in_data:
                header_lines.append(line)
                continue
            if not line.strip():
                continue
            parts = line.split("\t")
            data[parts[0]] = parts[1:]
    return header_lines, data


def match_tip_to_data(tip, data):
    """Direct match first; fall back to merged-sample component matching.
    Returns the matched row (list of fields) or None."""
    if tip in data:
        return data[tip]
    components = parse_merged_tip(tip)
    if len(components) > 1:
        matches = [data[c] for c in components if c in data]
        if matches:
            first_vals = set(m[0] for m in matches)
            if len(first_vals) > 1:
                print(f"WARNING: merged tip {tip} has inconsistent metadata "
                      f"between components ({first_vals}) -- using the first match",
                      file=sys.stderr)
            return matches[0]
    return None


def write_dataset(header_lines, matched_rows, out_path):
    with open(out_path, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for tip, row in matched_rows.items():
            fh.write(tip + "\t" + "\t".join(row) + "\n")


def build_genbank_lookup(xlsx_path):
    """'Assembeled' sheet: GenBank Accession -> organism name & isolate name."""
    df = pd.read_excel(xlsx_path, sheet_name="Assembeled", header=1)
    df.columns = [str(c).strip() for c in df.columns]
    cols = df.columns.tolist()
    df = df.rename(columns={cols[0]: "Accession", cols[1]: "Organism", cols[2]: "Isolate"})
    lookup = {}
    for _, row in df.iterrows():
        acc = row.get("Accession")
        if pd.isna(acc):
            continue
        isolate = row.get("Isolate")
        lookup[str(acc).strip()] = {
            "organism": "" if pd.isna(row.get("Organism")) else str(row.get("Organism")).strip(),
            "isolate": "" if pd.isna(isolate) else str(isolate).strip(),
        }
    return lookup


def build_sra_lookup(xlsx_path):
    """'Raw_Data' sheet: SRA Run accession -> species group, isolate, description.
    The sheet has three stacked sub-tables (LSDV / SPPV / GTPV), each introduced
    by a lone species-name row -- this function tracks which block it is in."""
    df = pd.read_excel(xlsx_path, sheet_name="Raw_Data", header=1)
    df.columns = [str(c).strip() for c in df.columns]
    run_col = df.columns[1]
    isolate_col = next((c for c in df.columns if c.strip().lower() == "isolate"), None)
    desc_col = next((c for c in df.columns if c.strip().lower() == "description"), None)

    lookup = {}
    current_species = "LSDV"
    for _, row in df.iterrows():
        val = row.get(run_col)
        if pd.isna(val):
            continue
        val = str(val).strip()
        if val in ("SPPV", "GTPV", "LSDV"):
            current_species = val
            continue
        if val == "Run":
            continue
        if re.match(r"^[SED]RR\d+$", val):
            isolate = row.get(isolate_col) if isolate_col else None
            desc = row.get(desc_col) if desc_col else None
            lookup[val] = {
                "species": current_species,
                "isolate": "" if (isolate is None or pd.isna(isolate)) else str(isolate).strip(),
                "description": "" if (desc is None or pd.isna(desc)) else str(desc).strip(),
            }
    return lookup


def classify_species_from_organism(organism_name):
    o = organism_name.lower()
    if "lumpy" in o:
        return "LSDV"
    if "sheep" in o:
        return "SPPV"
    if "goat" in o:
        return "GTPV"
    return "Unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    tip_source = parser.add_mutually_exclusive_group(required=True)
    tip_source.add_argument("--tree", help="Newick/contree tree file (preferred, final run)")
    tip_source.add_argument("--tips-from-fasta",
                             help="Concatenated supermatrix FASTA to source tip IDs from "
                                  "*before* the tree finishes (e.g. 03_supermatrix/supermatrix.fasta). "
                                  "Tip IDs will be identical to the eventual tree, so datasets "
                                  "prepared this way do not need to be regenerated later -- "
                                  "unless you want the belt-and-suspenders double-check.")
    parser.add_argument("--supplementary", required=True)
    parser.add_argument("--itol-dir", required=True,
                         help="Directory containing itol_condition.txt, itol_continent.txt, "
                              "itol_country_labels.txt, itol_host.txt")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    itol_dir = Path(args.itol_dir)

    if args.tree:
        tips = load_tree_tips(args.tree)
        print(f"[INFO] {len(tips)} tips loaded from tree {args.tree}")
    else:
        tips = load_tips_from_fasta(args.tips_from_fasta)
        print(f"[INFO] {len(tips)} tips loaded from supermatrix FASTA {args.tips_from_fasta} "
              f"(tree not yet finished -- these IDs will match the tree exactly)")

    unmatched_report = defaultdict(list)

    # ------------------------------------------------------------------
    # 1) Filter the four existing datasets down to exactly the tree's tips
    # ------------------------------------------------------------------
    existing_files = {
        "itol_condition.txt": "itol_condition_FINAL.txt",
        "itol_continent.txt": "itol_continent_FINAL.txt",
        "itol_country_labels.txt": "itol_country_labels_FINAL.txt",
        "itol_host.txt": "itol_host_FINAL.txt",
    }
    for src_name, out_name in existing_files.items():
        src_path = itol_dir / src_name
        if not src_path.exists():
            print(f"WARNING: {src_path} not found -- skipped", file=sys.stderr)
            continue
        header_lines, data = load_itol_dataset(src_path)
        header_lines = remap_header_legend_colors(src_name, header_lines)
        matched = {}
        for tip in tips:
            row = match_tip_to_data(tip, data)
            if row is not None:
                matched[tip] = remap_row_color(src_name, row)
            else:
                unmatched_report[src_name].append(tip)
        write_dataset(header_lines, matched, outdir / out_name)
        print(f"[DONE] {out_name}: matched {len(matched)}/{len(tips)} tips")

    # ------------------------------------------------------------------
    # 2) Build two bonus datasets directly from Supplementary_File1.xlsx
    # ------------------------------------------------------------------
    genbank_lookup = build_genbank_lookup(args.supplementary)
    sra_lookup = build_sra_lookup(args.supplementary)

    species_colors = VIRAL_SPECIES_COLORS
    species_rows, label_rows = {}, {}

    for tip in tips:
        components = parse_merged_tip(tip)
        species, isolate_label = None, None

        if tip in genbank_lookup:
            info = genbank_lookup[tip]
            species = classify_species_from_organism(info["organism"])
            isolate_label = info["isolate"] or tip

        if species is None:
            for c in components:
                if c in sra_lookup:
                    info = sra_lookup[c]
                    species = info["species"]
                    isolate_label = info["isolate"] or info["description"] or tip
                    break

        if species is None:
            species = "Unknown"
            unmatched_report["viral_species/labels"].append(tip)
        if not isolate_label:
            isolate_label = tip

        isolate_label = isolate_label.replace("\t", " ").strip()
        if len(isolate_label) > 40:
            isolate_label = isolate_label[:37] + "..."

        species_rows[tip] = [species_colors.get(species, "#7f8c8d")]
        label_rows[tip] = [isolate_label]

    species_header = [
        "DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tViral Species",
        "COLOR\t#000000", "STRIP_WIDTH\t30", "MARGIN\t2", "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tViral Species",
        "LEGEND_SHAPES\t1\t1\t1\t1",
        f"LEGEND_COLORS\t{species_colors['LSDV']}\t{species_colors['SPPV']}\t"
        f"{species_colors['GTPV']}\t{species_colors['Unknown']}",
        "LEGEND_LABELS\tLSDV\tSPPV\tGTPV\tUnknown",
        "DATA",
    ]
    write_dataset(species_header, species_rows, outdir / "itol_viral_species_FINAL.txt")
    print(f"[DONE] itol_viral_species_FINAL.txt: {len(species_rows)} tips")

    label_header = ["LABELS", "SEPARATOR TAB", "DATA"]
    write_dataset(label_header, label_rows, outdir / "itol_leaf_labels_FINAL.txt")
    print(f"[DONE] itol_leaf_labels_FINAL.txt: {len(label_rows)} tips")

    # ------------------------------------------------------------------
    # 3) Report anything that could not be matched, for manual review
    # ------------------------------------------------------------------
    report_path = outdir / "unmatched_tips_report.txt"
    with open(report_path, "w") as fh:
        if not unmatched_report:
            fh.write("All tree tips were matched successfully in every dataset.\n")
        for src, tips_list in unmatched_report.items():
            fh.write(f"=== {src}: {len(tips_list)} unmatched tip(s) ===\n")
            for t in tips_list:
                fh.write(f"  {t}\n")
            fh.write("\n")
    total_unmatched = sum(len(v) for v in unmatched_report.values())
    print(f"[DONE] Unmatched-tip report written to {report_path}")
    if total_unmatched:
        print(f"[WARN] {total_unmatched} unmatched entries across all datasets -- review the report.")
    else:
        print("[OK] Every tree tip matched in every dataset.")


if __name__ == "__main__":
    sys.exit(main())
