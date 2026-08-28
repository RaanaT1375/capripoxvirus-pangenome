import glob, os
from collections import Counter
import numpy as np, pandas as pd
from Bio import AlignIO
from scipy.stats import binomtest

META = "00_Metadata/species_assignment.csv"
ALN  = "10_Selection_Pressure/01_codon_alignments"
OUT  = "14_MK_Test/03_apobec_diagnostic"

meta = pd.read_csv(META)
sp = dict(zip(meta.iloc[:,0].astype(str), meta.iloc[:,1].astype(str)))

def cons(col, min_frac=0.90):
    c = Counter(b for b in col if b in "ACGT")
    if not c: return None
    b, n = c.most_common(1)[0]
    return b if n/sum(c.values()) >= min_frac else None

poly, fixed = Counter(), Counter()          # نوع تغییر
poly_ctx, fixed_ctx = Counter(), Counter()  # تغییر + بافت
bg = Counter()                              # پس‌زمینه: TC و GA در توالی اجماعی LSDV
per_gene = []

for f in sorted(glob.glob(os.path.join(ALN, "OG*_codon.fasta"))):
    og = os.path.basename(f).replace("_codon.fasta","")
    aln = AlignIO.read(f, "fasta")
    ids = [r.id for r in aln]
    arr = np.array([list(str(r.seq).upper()) for r in aln])
    is_l = np.array([sp.get(i,"?")=="LSDV" for i in ids])
    is_o = np.array([sp.get(i,"?") in ("SPPV","GTPV") for i in ids])
    if is_l.sum() < 50 or is_o.sum() < 20:
        print(f"  ! {og}: تخصیص گونه ناقص (LSDV={is_l.sum()} OUT={is_o.sum()}) — رد شد"); continue
    L, O = arr[is_l], arr[is_o]
    n = arr.shape[1]

    # توالی اجماعی LSDV برای تعیین بافت
    lcons = []
    for i in range(n):
        c = Counter(b for b in L[:,i] if b in "ACGT")
        lcons.append(c.most_common(1)[0][0] if c else "-")
    lcons = "".join(lcons)

    # پس‌زمینه: چند درصد از Cها پیش‌ازشان T است؟ چند درصد از Gها پس‌ازشان A؟
    for i in range(n):
        b = lcons[i]
        if b == "C":
            bg["C_total"] += 1
            if i>0 and lcons[i-1]=="T": bg["C_TCctx"] += 1
        elif b == "G":
            bg["G_total"] += 1
            if i<n-1 and lcons[i+1]=="A": bg["G_GActx"] += 1

    gp = gf = gp_ap = gf_ap = 0
    for i in range(n):
        anc = cons(O[:,i])                 # حالت اجدادی از برون‌گروه
        if anc is None: continue
        c = Counter(b for b in L[:,i] if b in "ACGT")
        tot = sum(c.values())
        if tot < 100: continue
        alleles = {b: k/tot for b,k in c.items() if k/tot >= 0.05}
        if not alleles: continue
        prev_, next_ = (lcons[i-1] if i>0 else "-"), (lcons[i+1] if i<n-1 else "-")

        def apo(a, d):
            return (a=="C" and d=="T" and prev_=="T") or (a=="G" and d=="A" and next_=="A")

        if len(alleles) >= 2:              # چندشکل درون LSDV
            for d in alleles:
                if d != anc:
                    poly[f"{anc}>{d}"] += 1; gp += 1
                    if apo(anc,d): poly_ctx[f"{anc}>{d}"] += 1; gp_ap += 1
        else:                              # ثابت
            d = next(iter(alleles))
            if d != anc:
                fixed[f"{anc}>{d}"] += 1; gf += 1
                if apo(anc,d): fixed_ctx[f"{anc}>{d}"] += 1; gf_ap += 1
    per_gene.append(dict(Orthogroup=og, n_poly=gp, n_poly_apobec=gp_ap,
                         n_fixed=gf, n_fixed_apobec=gf_ap))

pd.DataFrame(per_gene).to_csv(f"{OUT}/apobec_per_gene.csv", index=False)

def report(name, cnt, ctx):
    tot = sum(cnt.values())
    ct, ga = cnt.get("C>T",0), cnt.get("G>A",0)
    print(f"\n=== {name} (n={tot}) ===")
    for k,v in sorted(cnt.items(), key=lambda x:-x[1])[:6]:
        print(f"   {k}: {v}  ({100*v/tot:.1f}%)" if tot else "")
    print(f"   C>T + G>A = {ct+ga}  ({100*(ct+ga)/tot:.1f}% از کل)" if tot else "")
    bgC = bg["C_TCctx"]/bg["C_total"] if bg["C_total"] else 0
    bgG = bg["G_GActx"]/bg["G_total"] if bg["G_total"] else 0
    for lab, obs, n_, exp in [("C>T در بافت TC", ctx.get("C>T",0), ct, bgC),
                              ("G>A در بافت GA", ctx.get("G>A",0), ga, bgG)]:
        if n_ >= 10:
            p = binomtest(obs, n_, exp, alternative="greater").pvalue
            print(f"   {lab}: {obs}/{n_} = {obs/n_:.3f} | "
                  f"انتظار از ترکیب باز = {exp:.3f} | "
                  f"غنی‌شدگی = {obs/n_/exp:.2f}× | p = {p:.3g}")
        else:
            print(f"   {lab}: تعداد کم ({n_}) — آزمون‌ناپذیر")

print(f"\nپس‌زمینهٔ ترکیب باز: TC/C = {bg['C_TCctx']}/{bg['C_total']}"
      f" | GA/G = {bg['G_GActx']}/{bg['G_total']}")
report("چندشکلی درون LSDV (P)", poly, poly_ctx)
report("تفاوت ثابت روی شاخهٔ LSDV (D)", fixed, fixed_ctx)
