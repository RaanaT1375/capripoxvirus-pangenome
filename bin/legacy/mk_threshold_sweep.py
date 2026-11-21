#!/usr/bin/env python3
"""
Allele-frequency threshold sweep for the polarized McDonald-Kreitman test.

Reuses the exact functions defined in 14_MK_Test/scripts/mk_polarized.py by
exec'ing only its header (imports + constants + function definitions), so the
module's top-level analysis code never runs and NO existing result file is
overwritten.

For each cutoff in {0.00, 0.01, 0.02, 0.05, 0.10, 0.15} the polarized test is
recomputed twice: once with the vaccine-flagged LSDV genomes excluded from the
ingroup (the primary design) and once with them included.

Run from:  /cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis
Outputs:   14_MK_Test/03_summary/mk_threshold_sweep.csv
           14_MK_Test/03_summary/mk_threshold_sweep_pergene.csv
           14_MK_Test/03_summary/mk_threshold_sweep.pdf
"""
import os, sys, glob, collections
import numpy as np
import pandas as pd

BASE = "/cfs/earth/scratch/xpkk/Raana/Pangenome_Analysis"
SRC = os.path.join(BASE, "14_MK_Test/scripts/mk_polarized.py")
OUTDIR = os.path.join(BASE, "14_MK_Test/03_summary")
CUTOFFS = [0.00, 0.01, 0.02, 0.05, 0.10, 0.15]
N_BOOT = 10000
SEED = 20260817

os.chdir(BASE)
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------------------------------------------------
# Load only the header of mk_polarized.py (up to and including the end of
# polarized_gene). The sentinel is the return statement of that function.
# ----------------------------------------------------------------------
with open(SRC) as fh:
    lines = fh.readlines()

end = None
for i, ln in enumerate(lines):
    if ln.strip().startswith("return dict(Pn=Pn"):
        end = i + 1
        break
if end is None:
    sys.exit("ERROR: could not locate the end of polarized_gene() in mk_polarized.py")

ns = {"__name__": "mk_polarized_header"}
exec("".join(lines[:end]), ns)

polarized_gene = ns["polarized_gene"]
cmh = ns["cmh"]
ALN = ns["ALN"]
print(f"Loaded {end} header lines from mk_polarized.py "
      f"(baseline CUT={ns['CUT']}, CONS={ns['CONS']}, COV={ns['COV']})")

# ----------------------------------------------------------------------
# Same group definitions as the original script
# ----------------------------------------------------------------------
sp = pd.read_csv("00_Metadata/species_assignment.csv")
grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV", "SPPV", "GTPV"]}
vacc = set(pd.read_csv("11_Scoary/01_inputs/traits.csv", index_col=0)
             .query("Vaccine == 1").index)

ing_novacc = [g for g in grp["LSDV"] if g not in vacc]
ing_withvacc = list(grp["LSDV"])

print(f"LSDV total            : {len(grp['LSDV'])}")
print(f"Vaccine-flagged (all) : {len(vacc)}")
print(f"Vaccine-flagged, LSDV : {len(grp['LSDV']) - len(ing_novacc)}")
print(f"Ingroup, no vaccine   : {len(ing_novacc)}")
print(f"Ingroup, with vaccine : {len(ing_withvacc)}")
print(f"Outgroups             : SPPV={len(grp['SPPV'])}, GTPV={len(grp['GTPV'])}\n")

files = sorted(glob.glob(f"{ALN}/*_codon.fasta"))
print(f"Alignments found: {len(files)}\n")


def bootstrap_alpha(df, n_boot=N_BOOT, seed=SEED):
    """Resample orthogroups with replacement, recompute alpha each replicate."""
    rng = np.random.default_rng(seed)
    m = df[["Pn", "Ps", "Dn", "Ds"]].to_numpy(float)
    n = len(m)
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        Pn, Ps, Dn, Ds = m[idx].sum(axis=0)
        out[b] = np.nan if (Ps == 0 or Ds == 0 or Dn == 0) else 1 - (Pn / Ps) / (Dn / Ds)
    out = out[~np.isnan(out)]
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if len(out) else (np.nan, np.nan)


rows, pergene = [], []

for include_vacc, ing_ids, label in [
    (False, ing_novacc, "excluded"),
    (True, ing_withvacc, "included"),
]:
    for cut in CUTOFFS:
        ns["CUT"] = cut  # polarized_gene reads CUT from this namespace
        recs = []
        for f in files:
            d = polarized_gene(f, ing_ids, grp["SPPV"], grp["GTPV"])
            if d:
                d["Orthogroup"] = os.path.basename(f).replace("_codon.fasta", "")
                recs.append(d)
        df = pd.DataFrame(recs)
        df["cutoff"] = cut
        df["vaccines"] = label
        pergene.append(df)

        Pn, Ps, Dn, Ds = df[["Pn", "Ps", "Dn", "Ds"]].sum()
        NI = (Pn / Ps) / (Dn / Ds) if Ps and Ds and Dn else np.nan
        alpha = 1 - NI if NI == NI else np.nan

        tabs = [[[round(r.Dn), round(r.Ds)], [round(r.Pn), round(r.Ps)]]
                for r in df.itertuples()]
        try:
            p_cmh, mhor, _ = cmh(tabs)
        except Exception as e:
            p_cmh, mhor = np.nan, np.nan
            print(f"  (CMH failed at cutoff {cut}, vaccines {label}: {e})")

        lo, hi = bootstrap_alpha(df)

        rows.append(dict(cutoff=cut, vaccines_included=include_vacc,
                         n_ingroup=len(ing_ids), n_genes=len(df),
                         Pn=round(Pn, 1), Ps=round(Ps, 1),
                         Dn=round(Dn, 1), Ds=round(Ds, 1),
                         NI=round(NI, 4) if NI == NI else np.nan,
                         alpha=round(alpha, 4) if alpha == alpha else np.nan,
                         alpha_CI_low=round(lo, 4), alpha_CI_high=round(hi, 4),
                         CMH_p=p_cmh, MH_OR=mhor))
        print(f"vaccines {label:8s} | cutoff {cut:.2f} | "
              f"Pn={Pn:7.1f} Ps={Ps:7.1f} Dn={Dn:6.1f} Ds={Ds:6.1f} | "
              f"alpha={alpha:.3f} [{lo:.2f}, {hi:.2f}] | CMH p={p_cmh:.3e}")

res = pd.DataFrame(rows)
res.to_csv(f"{OUTDIR}/mk_threshold_sweep.csv", index=False)
pd.concat(pergene, ignore_index=True).to_csv(
    f"{OUTDIR}/mk_threshold_sweep_pergene.csv", index=False)
print(f"\nWrote {OUTDIR}/mk_threshold_sweep.csv")
print(f"Wrote {OUTDIR}/mk_threshold_sweep_pergene.csv")

# ----------------------------------------------------------------------
# Panel: alpha vs cutoff, one line per vaccine treatment
# ----------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    style = {False: dict(color="#1f4e79", marker="o", label="Vaccine strains excluded"),
             True: dict(color="#b3541e", marker="s", label="Vaccine strains included")}
    for inc, sub in res.groupby("vaccines_included"):
        sub = sub.sort_values("cutoff")
        st = style[bool(inc)]
        ax.plot(sub.cutoff * 100, sub.alpha, lw=1.8, ms=5,
                color=st["color"], marker=st["marker"], label=st["label"])
        ax.fill_between(sub.cutoff * 100, sub.alpha_CI_low, sub.alpha_CI_high,
                        color=st["color"], alpha=0.15, lw=0)

    ax.axhline(0, color="0.4", lw=0.8, ls="--")
    ax.axvline(5, color="0.75", lw=0.8, ls=":")
    ax.set_xlabel("Minor-allele frequency cutoff (%)")
    ax.set_ylabel(r"$\alpha$ (adaptive fixed differences)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/mk_threshold_sweep.pdf")
    fig.savefig(f"{OUTDIR}/mk_threshold_sweep.png", dpi=300)
    print(f"Wrote {OUTDIR}/mk_threshold_sweep.pdf and .png")
except Exception as e:
    print(f"(plot skipped: {e})")
