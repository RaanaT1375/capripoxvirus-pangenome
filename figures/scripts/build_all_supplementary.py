#!/usr/bin/env python3
"""
Build & standardize all 16 Supplementary Excel files (S1-S16)
Unified Q1 Style: Times New Roman, Navy (#002F57) header, zebra rows, auto-width.
"""
import os
import glob
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

OUT_DIR = "revision_2026-11/final_supplementary"
os.makedirs(OUT_DIR, exist_ok=True)

# Unified Q1 styling (matching Figure 1-5 palette)
FONT_TITLE  = Font(name="Times New Roman", size=12, bold=True, color="002F57")
FONT_HEADER = Font(name="Times New Roman", size=11, bold=True, color="FFFFFF")
FONT_DATA   = Font(name="Times New Roman", size=10.5, bold=False, color="111111")
FILL_HEADER = PatternFill(start_color="002F57", end_color="002F57", fill_type="solid")
FILL_ZEBRA  = PatternFill(start_color="F4F6F9", end_color="F4F6F9", fill_type="solid")
FILL_WHITE  = PatternFill(start_color="FFFFFF", end_color="FFFFFF", fill_type="solid")
THIN_SIDE   = Side(border_style="thin", color="D0D7DE")
BORDER_CELL = Border(left=THIN_SIDE, right=THIN_SIDE, top=THIN_SIDE, bottom=THIN_SIDE)

def style_sheet(ws, title_text, header_row=3):
    ws.cell(row=1, column=1, value=title_text)
    ws.cell(row=1, column=1).font = FONT_TITLE
    ws.cell(row=1, column=1).alignment = Alignment(horizontal="left", vertical="center", wrap_text=False)
    ws.row_dimensions[1].height = 24
    ws.row_dimensions[2].height = 8

    max_col = ws.max_column
    max_row = ws.max_row

    # Style header row
    ws.row_dimensions[header_row].height = 26
    for c in range(1, max_col + 1):
        cell = ws.cell(row=header_row, column=c)
        cell.font = FONT_HEADER
        cell.fill = FILL_HEADER
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER_CELL

    # Style data rows
    for r in range(header_row + 1, max_row + 1):
        fill = FILL_ZEBRA if (r - header_row) % 2 == 0 else FILL_WHITE
        for c in range(1, max_col + 1):
            cell = ws.cell(row=r, column=c)
            cell.font = FONT_DATA
            cell.fill = fill
            cell.border = BORDER_CELL
            if isinstance(cell.value, (int, float)):
                cell.alignment = Alignment(horizontal="right", vertical="center")
            else:
                cell.alignment = Alignment(horizontal="left", vertical="center")

    # Auto-fit column widths (skip row 1 title)
    for c in range(1, max_col + 1):
        col_letter = get_column_letter(c)
        max_len = 10
        for r in range(header_row, min(max_row + 1, header_row + 200)):
            val = ws.cell(row=r, column=c).value
            if val is not None:
                max_len = max(max_len, len(str(val)))
        ws.column_dimensions[col_letter].width = min(max_len + 4, 48)

    ws.freeze_panes = f"A{header_row + 1}"

def write_multisheet_xlsx(filename, sheets_dict):
    """
    sheets_dict: {sheet_name: (title_text, df)}
    """
    out_path = os.path.join(OUT_DIR, filename)
    wb = openpyxl.Workbook()
    first = True
    for sh_name, (title_text, df) in sheets_dict.items():
        ws = wb.active if first else wb.create_sheet(title=sh_name[:31])
        ws.title = sh_name[:31]
        first = False
        # Write header at row 3, data from row 4
        for c_idx, col_name in enumerate(df.columns, start=1):
            ws.cell(row=3, column=c_idx, value=str(col_name))
        for r_idx, row in enumerate(df.itertuples(index=False), start=4):
            for c_idx, val in enumerate(row, start=1):
                if pd.isna(val):
                    ws.cell(row=r_idx, column=c_idx, value="")
                else:
                    ws.cell(row=r_idx, column=c_idx, value=val)
        style_sheet(ws, title_text, header_row=3)
    wb.save(out_path)
    print(f"[OK] Saved {out_path} ({len(sheets_dict)} sheet(s))")

def read_table(path):
    sep = "\t" if path.endswith(".tsv") else ","
    return pd.read_csv(path, sep=sep)

# =====================================================================
# PART 1: Fix & Restyle existing S1 - S7 (searches 99_Supplementary/ or .)
# =====================================================================
def find_existing(name):
    for d in ["99_Supplementary", ".", "revision_2026-11"]:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None

# S1: Fix sheet name typo (Assembeled_Data -> Assembled_Data) & A1 titles
p_s1 = find_existing("Supplementary File S1.xlsx")
if p_s1:
    wb = openpyxl.load_workbook(p_s1)
    if "Assembeled_Data" in wb.sheetnames:
        wb["Assembeled_Data"].title = "Assembled_Data"
    for sh in wb.sheetnames:
        ws = wb[sh]
        ws["A1"].value = "Supplementary File S1: Comprehensive metadata, accession numbers, host information, and geographic origins for the 381 compiled Capripoxvirus datasets."
        ws["A1"].font = FONT_TITLE
    out_s1 = os.path.join(OUT_DIR, "Supplementary File S1.xlsx")
    wb.save(out_s1)
    print(f"[OK] Fixed & saved {out_s1}")

# S2 - S7: Read data, apply specific fixes, and restyle uniformly
s2_s7_titles = {
    2: ("Core_Genes_60", "Supplementary File S2: The 60 universal single-copy core orthogroups used for the supermatrix and selection analyses."),
    3: ("Recombination_Blocks", "Supplementary File S3: Gubbins output for the 90 recombination blocks across the 110,321-bp core alignment."),
    4: ("Temporal_Isolates_205", "Supplementary File S4: Per-isolate data underlying the temporal signal analysis (205 dated LSDV isolates)."),
    5: ("Significant_Associations", "Supplementary File S5: The 92 significant gene-trait associations from structure-corrected Scoary pan-GWAS."),
    6: ("Per_Trait_Summary", "Supplementary File S6: Per-trait summary across all 19 evaluated phenotypes in the Scoary pan-GWAS."),
    7: ("Root_Subsampling_MPPA", "Supplementary File S7: Root ancestral states from the balanced continental sub-sampling analysis (MPPA).")
}

for idx, (sh_name, title_txt) in s2_s7_titles.items():
    fname = f"Supplementary File S{idx}.xlsx"
    p_in = find_existing(fname)
    if not p_in:
        print(f"[SKIP] {fname} not found on cluster (upload to 99_Supplementary/ to restyle)")
        continue
    df = pd.read_excel(p_in, header=2)
    if idx == 3 and len(df.columns) > 0:
        df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.strip('"')
    if idx == 7 and "Analysis" in df.columns:
        mask_full = df["Analysis"].astype(str).str.lower() == "full"
        df.loc[mask_full, "States Retained at Root (n)"] = 1
        df.loc[mask_full, "Root Resolution"] = "Asia"
    write_multisheet_xlsx(fname, {sh_name: (title_txt, df)})

# =====================================================================
# PART 2: Build S8 - S16 directly from revision_2026-11/
# =====================================================================
def build_from_dir(filename, folder, title_prefix):
    files = sorted(glob.glob(os.path.join(folder, "**", "*.csv"), recursive=True) +
                   glob.glob(os.path.join(folder, "**", "*.tsv"), recursive=True))
    if not files:
        print(f"[WARN] No CSV/TSV found in {folder} for {filename}")
        return
    sheets = {}
    for f in files:
        sh = os.path.splitext(os.path.basename(f))[0][:31]
        sheets[sh] = (f"{title_prefix} ({sh})", read_table(f))
    write_multisheet_xlsx(filename, sheets)

# S8: Provenance-exclusive orthogroups (19_Provenance_Audit)
build_from_dir("Supplementary File S8.xlsx", "revision_2026-11/19_Provenance_Audit",
               "Supplementary File S8: Provenance-exclusive orthogroups")

# S9: Root inference output (24_Root_Inference)
build_from_dir("Supplementary File S9.xlsx", "revision_2026-11/24_Root_Inference",
               "Supplementary File S9: Non-reversible root inference output")

# S10: DNA-level presence/absence matrix & GWAS (27_DNA_Level_GWAS)
s10_files = [
    ("Annotation_vs_DNA", "revision_2026-11/27_DNA_Level_GWAS/02_results/annotation_vs_dna.csv"),
    ("Old_vs_New_GWAS",   "revision_2026-11/27_DNA_Level_GWAS/03_summary/old_vs_new.csv"),
    ("Scoary_DNA_Sig",    "revision_2026-11/27_DNA_Level_GWAS/03_summary/scoary_dna_significant.csv"),
    ("Scoary_DNA_Conv",   "revision_2026-11/27_DNA_Level_GWAS/03_summary/scoary_dna_convergent.csv"),
    ("DNA_Matrix",        "revision_2026-11/27_DNA_Level_GWAS/03_summary/gene_presence_absence_dna.csv")
]
s10_sheets = {k: (f"Supplementary File S10: Nucleotide-level matrix validation & pan-GWAS ({k})", read_table(v))
              for k, v in s10_files if os.path.exists(v)}
if s10_sheets:
    write_multisheet_xlsx("Supplementary File S10.xlsx", s10_sheets)

# S11: Functional classification & constraint (28_Positive_Results)
s11_files = [
    ("Gene_Constraint_60", "revision_2026-11/28_Positive_Results/02_results/gene_constraint.csv"),
    ("Constraint_by_Class","revision_2026-11/28_Positive_Results/03_summary/constraint_by_class.csv"),
    ("Constraint_Tests",   "revision_2026-11/28_Positive_Results/03_summary/constraint_tests.csv"),
    ("Functional_Rules",   "revision_2026-11/28_Positive_Results/01_inputs/functional_categories.csv")
]
s11_sheets = {k: (f"Supplementary File S11: Functional classification & gene-level constraint ({k})", read_table(v))
              for k, v in s11_files if os.path.exists(v)}
if s11_sheets:
    write_multisheet_xlsx("Supplementary File S11.xlsx", s11_sheets)

# S12: Fixed substitutions & positional clustering (28_Positive_Results)
s12_files = [
    ("Fixed_Substitutions", "revision_2026-11/28_Positive_Results/02_results/fixed_substitutions.csv"),
    ("Clustering_Tests",    "revision_2026-11/28_Positive_Results/02_results/clustering_tests.csv"),
    ("Nameable_Clean",      "revision_2026-11/28_Positive_Results/03_summary/substitutions_clean.csv"),
    ("Flagged_Artefacts",   "revision_2026-11/28_Positive_Results/03_summary/substitutions_flagged.csv")
]
s12_sheets = {k: (f"Supplementary File S12: Fixed substitutions on the ancestral LSDV branch ({k})", read_table(v))
              for k, v in s12_files if os.path.exists(v)}
if s12_sheets:
    write_multisheet_xlsx("Supplementary File S12.xlsx", s12_sheets)

# S13: Recombination breakpoint analysis (28_Positive_Results)
s13_files = [
    ("Breakpoint_Positions", "revision_2026-11/28_Positive_Results/02_results/breakpoint_positions.csv"),
    ("Breakpoint_Density",   "revision_2026-11/28_Positive_Results/02_results/breakpoint_density.csv"),
    ("Recurrent_Intervals",  "revision_2026-11/28_Positive_Results/03_summary/recurrent_intervals.csv"),
    ("Ancestral_Node147",    "revision_2026-11/28_Positive_Results/03_summary/ancestral_breakpoints.csv")
]
s13_sheets = {k: (f"Supplementary File S13: Recombination breakpoint analysis ({k})", read_table(v))
              for k, v in s13_files if os.path.exists(v)}
if s13_sheets:
    write_multisheet_xlsx("Supplementary File S13.xlsx", s13_sheets)

# S14: Lineage-diagnostic core sites (28_Positive_Results)
s14_files = [
    ("Final_Spaced_Panel",  "revision_2026-11/28_Positive_Results/03_summary/panel_final.csv"),
    ("Diagnostic_Panel",    "revision_2026-11/28_Positive_Results/03_summary/diagnostic_panel.csv"),
    ("Recomb_Screened",     "revision_2026-11/28_Positive_Results/02_results/sites_recombination_ok.csv"),
    ("All_Diagnostic_Sites","revision_2026-11/28_Positive_Results/02_results/diagnostic_sites.csv")
]
s14_sheets = {k: (f"Supplementary File S14: Lineage-diagnostic core alignment sites ({k})", read_table(v))
              for k, v in s14_files if os.path.exists(v)}
if s14_sheets:
    write_multisheet_xlsx("Supplementary File S14.xlsx", s14_sheets)

# S15: Ancestral gene-content reconstruction (18_Gene_Content_ASR)
s15_files = [
    ("Compartment_Summary", "revision_2026-11/18_Gene_Content_ASR/03_summary/gene_content_summary.csv"),
    ("Events_Public_237",   "revision_2026-11/18_Gene_Content_ASR/02_results/gene_content_events_genbank237.csv"),
    ("Events_All_290",      "revision_2026-11/18_Gene_Content_ASR/02_results/gene_content_events_all290.csv")
]
s15_sheets = {k: (f"Supplementary File S15: Ancestral gene-content reconstruction ({k})", read_table(v))
              for k, v in s15_files if os.path.exists(v)}
if s15_sheets:
    write_multisheet_xlsx("Supplementary File S15.xlsx", s15_sheets)

# S16: Temporal signal power analysis (29_Temporal_Reconciliation)
s16_files = [
    ("Power_Summary",     "revision_2026-11/29_Temporal_Reconciliation/03_summary/power_summary.csv"),
    ("Temporal_Refit",    "revision_2026-11/29_Temporal_Reconciliation/02_results/temporal_refit.csv"),
    ("Rate_Expectations", "revision_2026-11/29_Temporal_Reconciliation/02_results/rate_expectations.csv")
]
s16_sheets = {k: (f"Supplementary File S16: Temporal signal power analysis ({k})", read_table(v))
              for k, v in s16_files if os.path.exists(v)}
if s16_sheets:
    write_multisheet_xlsx("Supplementary File S16.xlsx", s16_sheets)

print("\nDone! All generated supplementary files are in:", OUT_DIR)
