import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(13, 6.2), dpi=300)
ax.set_xlim(0, 13)
ax.set_ylim(0, 6.2)
ax.axis("off")

# Title banner
ax.text(6.5, 5.82, "Capripoxvirus Pan-Genomic & Population-Genetic Analysis Workflow",
        ha="center", va="center", fontsize=14.5, fontweight="bold", color="#002F57", family="DejaVu Sans")
ax.text(6.5, 5.48, "290 Curated Genomes (240 LSDV  |  34 SPPV  |  16 GTPV)  •  60 Single-Copy Core Genes  •  Release v2.0.1",
        ha="center", va="center", fontsize=9.8, color="#495C7C", family="DejaVu Sans")

def draw_box(x, y, w, h, title, body, header_col="#002F57", fill_col="#F7F9FC"):
    # Main card
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                                  linewidth=1.1, edgecolor=header_col, facecolor=fill_col, zorder=2)
    ax.add_patch(rect)
    # Header bar
    hdr = patches.FancyBboxPatch((x, y + h - 0.42), w, 0.42, boxstyle="round,pad=0.04,rounding_size=0.12",
                                 linewidth=0, facecolor=header_col, zorder=3)
    ax.add_patch(hdr)
    ax.text(x + w/2, y + h - 0.21, title, ha="center", va="center",
            fontsize=9.3, fontweight="bold", color="white", family="DejaVu Sans", zorder=4)
    ax.text(x + w/2, y + (h - 0.42)/2, body, ha="center", va="center",
            fontsize=8.2, color="#1A1A1A", family="DejaVu Sans", linespacing=1.42, zorder=4)

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>,head_width=0.25,head_length=0.45",
                                color="#495C7C", lw=1.4), zorder=5)

# Row 1: Curation -> Pan-genome -> Phylogenomics & Recombination
draw_box(0.35, 3.05, 3.65, 2.05, "1. Data Curation & QC",
         "381 GenBank / ENA datasets\nfastp • Bowtie2 • SPAdes • Prokka\n3-step filtering (46 excluded):\n10 incomplete • 32 contaminated\n4 chimeric amplicon merges\n→ 290 curated genomes",
         header_col="#002F57")

draw_box(4.65, 3.05, 3.65, 2.05, "2. Pan-Genome Architecture",
         "OrthoFinder (206 orthogroups)\nEffectively closed (Heaps γ = 0.068)\n132 strict core • 60 single-copy\nGTR2+FO ancestral gene-content\nProksee macro-synteny\n(terminal diversity enrichment)",
         header_col="#002F57")

draw_box(8.95, 3.05, 3.65, 2.05, "3. Phylogeny, Rooting & Clonal Structure",
         "IQ-TREE partitioned ML (16,215 aa)\nOutgroup-free UNREST rooting\n(LSDV basal in 10/10 replicates)\nParsnp + Gubbins (r/m = 0.13)\nAncestral Node_147 (100 genomes)\nfastBAPS (5 clusters) • PastML",
         header_col="#495C7C")

arrow(4.08, 4.07, 4.57, 4.07)
arrow(8.38, 4.07, 8.87, 4.07)

# Connector down from Row 1 to Row 2
ax.plot([10.77, 10.77, 2.17, 2.17], [3.00, 2.68, 2.68, 2.40], color="#495C7C", lw=1.3, zorder=1)
ax.plot([6.47, 6.47], [2.68, 2.40], color="#495C7C", lw=1.3, zorder=1)
ax.plot([10.77, 10.77], [2.68, 2.40], color="#495C7C", lw=1.3, zorder=1)
arrow(2.17, 2.48, 2.17, 2.38)
arrow(6.47, 2.48, 6.47, 2.38)
arrow(10.77, 2.48, 10.77, 2.38)

# Row 2: Codon Models -> MK Test -> Validation & Diagnostics
draw_box(0.35, 0.25, 3.65, 2.05, "4. Selection & Constraint",
         "HyPhy BUSTED • aBSREL • MEME • RELAX\nPer-site dN/dS unidentifiable in LSDV\n(88–96% sites dS ≈ 0; 48 haplotypes)\nGene-level ω_MG94 structured by\nfunction (Kruskal–Wallis p = 0.0024,\n2.6× transcription vs host-interaction)",
         header_col="#585D26", fill_col="#F9FAF5")

draw_box(4.65, 0.25, 3.65, 2.05, "5. Ancestral Polygenic Adaptation",
         "Polarized McDonald–Kreitman test\nML ASR (PP ≥ 0.90): α = 0.629\n(95% CI 0.26–0.77; CMH p = 0.009)\nRobust to recombination masking\n(α = 0.588) & MAF sweep (1–15%)\n40 clean substitutions in 18 genes",
         header_col="#871746", fill_col="#FDF7F9")

draw_box(8.95, 0.25, 3.65, 2.05, "6. Validation & Diagnostic Markers",
         "Nucleotide-level pan-GWAS (tblastn)\nSeparates annotation dropout &\nassembly provenance from biology\nTemporal signal power analysis\n35 recombination-screened core\nlineage-diagnostic SNP markers",
         header_col="#871746", fill_col="#FDF7F9")

plt.tight_layout(pad=0.3)
plt.savefig("figures/workflow_overview.svg", format="svg", bbox_inches="tight")
plt.savefig("figures/workflow_overview.png", format="png", dpi=300, bbox_inches="tight")
print("[OK] Generated figures/workflow_overview.svg and .png")
