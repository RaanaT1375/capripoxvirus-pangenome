#!/usr/bin/env python3
"""McDonald-Kreitman test: LSDV polymorphism vs. divergence to SPPV / GTPV."""
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import permutations
from scipy.stats import fisher_exact
import pandas as pd, numpy as np, glob, os, collections, sys

ALN_DIR = "10_Selection_Pressure/01_codon_alignments"
OUT_DIR = "14_MK_Test"
VALID   = set("ACGT")
_tr = {}

def tr(c):
    if c not in _tr:
        _tr[c] = str(Seq(c).translate())
    return _tr[c]

def path_counts(a, b):
    """Nei-Gojobori pathway averaging -> (syn, nonsyn) between two codons."""
    pos = [i for i in range(3) if a[i] != b[i]]
    if not pos:
        return 0.0, 0.0
    npaths = s_tot = n_tot = 0
    for order in set(permutations(pos)):
        cur, s, n, ok = a, 0, 0, True
        for p in order:
            nxt = cur[:p] + b[p] + cur[p+1:]
            if tr(nxt) == "*":
                ok = False; break
            if tr(nxt) == tr(cur): s += 1
            else:                  n += 1
            cur = nxt
        if ok:
            npaths += 1; s_tot += s; n_tot += n
    if npaths == 0:                      # every path passes a stop
        d = len(pos)
        return (d, 0.0) if tr(a) == tr(b) else (0.0, d)
    return s_tot/npaths, n_tot/npaths

def clean(codon):
    return len(codon) == 3 and set(codon) <= VALID and tr(codon) != "*"

def mk_one_gene(path, ingroup, outgroup, min_cov=0.90, out_cons=0.90):
    recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}
    ing = [recs[i] for i in ingroup  if i in recs]
    out = [recs[i] for i in outgroup if i in recs]
    if not ing or not out:
        return None
    L = len(next(iter(recs.values()))) // 3
    Pn = Ps = Dn = Ds = 0.0
    Pn_ns = Ps_ns = 0.0          # no singletons
    n_skip_cov = n_skip_amb = n_poly_div = 0

    for k in range(L):
        sl = slice(3*k, 3*k+3)
        ic = [s[sl] for s in ing]; ic = [c for c in ic if clean(c)]
        oc = [s[sl] for s in out]; oc = [c for c in oc if clean(c)]
        if len(ic) < min_cov*len(ing) or len(oc) < min_cov*len(out):
            n_skip_cov += 1; continue
        icc, occ = collections.Counter(ic), collections.Counter(oc)
        anc, anc_n = occ.most_common(1)[0]
        if anc_n / len(oc) < out_cons:
            n_skip_amb += 1; continue
        major = icc.most_common(1)[0][0]

        if len(icc) > 1:                                   # polymorphic in LSDV
            for v, cnt in icc.items():
                if v == major: continue
                s, n = path_counts(major, v)
                Ps += s; Pn += n
                if cnt > 1:
                    Ps_ns += s; Pn_ns += n
            if major != anc: n_poly_div += 1               # counted as P only
        elif major != anc:                                 # fixed difference
            s, n = path_counts(anc, major)
            Ds += s; Dn += n

    return dict(Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds, Pn_nosingle=Pn_ns, Ps_nosingle=Ps_ns,
                n_codons=L, skip_coverage=n_skip_cov, skip_ambiguous=n_skip_amb,
                poly_and_divergent=n_poly_div, n_in=len(ing), n_out=len(out))

def stats_row(d):
    Pn, Ps, Dn, Ds = d["Pn"], d["Ps"], d["Dn"], d["Ds"]
    tbl = [[round(Dn), round(Ds)], [round(Pn), round(Ps)]]
    p = fisher_exact(tbl)[1] if min(sum(tbl[0]), sum(tbl[1])) > 0 else np.nan
    NI  = (Pn/Ps)/(Dn/Ds) if Ps > 0 and Ds > 0 and Dn > 0 else np.nan
    DoS = (Dn/(Dn+Ds) - Pn/(Pn+Pn*0+Ps)) if (Dn+Ds) > 0 and (Pn+Ps) > 0 else np.nan
    return dict(fisher_p=p, NI=NI, alpha=(1-NI if NI==NI else np.nan), DoS=DoS)

if __name__ == "__main__":
    sp = pd.read_csv("00_Metadata/species_assignment.csv")
    grp = {s: sp.loc[sp.species == s, "Name"].tolist() for s in ["LSDV","SPPV","GTPV"]}
    ann = pd.read_csv("10_Selection_Pressure/04_summary_tables/orthogroup_annotations.csv")

    for og_sp in ["SPPV", "GTPV"]:
        rows = []
        for f in sorted(glob.glob(f"{ALN_DIR}/*_codon.fasta")):
            og = os.path.basename(f).replace("_codon.fasta", "")
            r = mk_one_gene(f, grp["LSDV"], grp[og_sp])
            if r is None: continue
            r.update(Orthogroup=og, outgroup=og_sp); r.update(stats_row(r))
            rows.append(r)
        df = pd.DataFrame(rows)
        ok = df.fisher_p.notna()
        df.loc[ok, "q_BH"] = (df.loc[ok].sort_values("fisher_p").fisher_p.values *
                              ok.sum() / np.arange(1, ok.sum()+1))[
                              np.argsort(np.argsort(df.loc[ok, "fisher_p"].values))]
        df["q_BH"] = df["q_BH"].clip(upper=1.0)
        df = df.merge(ann[["Orthogroup","product"]], on="Orthogroup", how="left")
        df.to_csv(f"{OUT_DIR}/02_results/mk_vs_{og_sp}.csv", index=False)

        tot = df[["Pn","Ps","Dn","Ds"]].sum()
        pooled = stats_row(tot.to_dict())
        print(f"\n{'='*78}\n=== MK: LSDV vs {og_sp}  (n={len(df)} ژن) ===")
        print(f"  مجموع: Pn={tot.Pn:.1f} Ps={tot.Ps:.1f} Dn={tot.Dn:.1f} Ds={tot.Ds:.1f}")
        print(f"  تجمیعی: NI={pooled['NI']:.3f} | alpha={pooled['alpha']:.3f} "
              f"| DoS={pooled['DoS']:+.4f} | Fisher p={pooled['fisher_p']:.3e}")
        print(f"  ژن با Ps=0: {(df.Ps==0).sum()} | با Ds=0: {(df.Ds==0).sum()}")
        print(f"  میانه سایت‌های کنارگذاشته: پوشش {df.skip_coverage.median():.0f}"
              f" | مبهم {df.skip_ambiguous.median():.0f}")
        print(f"  سایت‌های هم‌چندشکل‌هم‌واگرا (محافظه‌کارانه فقط P): "
              f"{df.poly_and_divergent.sum():.0f}")
        sig = df[df.q_BH < 0.05].sort_values("fisher_p")
        print(f"\n  --- ژن‌های معنادار (q<0.05): {len(sig)} ---")
        if len(sig):
            print(sig[["Orthogroup","Pn","Ps","Dn","Ds","NI","DoS","q_BH","product"]]
                  .to_string(index=False, max_colwidth=42))
