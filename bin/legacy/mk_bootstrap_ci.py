#!/usr/bin/env python3
import glob, os, json
import numpy as np, pandas as pd

ROOT = "14_MK_Test"
OUT  = "14_MK_Test/03_summary"
NBOOT, SEED = 10000, 20260817
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(SEED)

def alpha_pooled(Pn, Ps, Dn, Ds):
    if Ps <= 0 or Dn <= 0: return np.nan
    return 1.0 - (Pn * Ds) / (Ps * Dn)

def mh_or_rbg(a, b, c, d):
    """a=Dn b=Ds c=Pn d=Ps ؛ CI رابینز-برسلو-گرینلند روی log(OR)"""
    n = a + b + c + d
    ok = n > 0
    a, b, c, d, n = a[ok], b[ok], c[ok], d[ok], n[ok]
    R, S = a * d / n, b * c / n
    Rp, Sp = R.sum(), S.sum()
    if Rp <= 0 or Sp <= 0: return np.nan, np.nan, np.nan, int(ok.sum())
    P, Q = (a + d) / n, (b + c) / n
    var = ((P * R).sum() / (2 * Rp**2)
           + ((P * S).sum() + (Q * R).sum()) / (2 * Rp * Sp)
           + ((Q * S).sum() / (2 * Sp**2)))
    se = np.sqrt(var); lo_r = np.log(Rp / Sp)
    return Rp / Sp, np.exp(lo_r - 1.96 * se), np.exp(lo_r + 1.96 * se), int(ok.sum())

need = {"Pn", "Ps", "Dn", "Ds", "Orthogroup"}
cands = sorted(glob.glob(os.path.join(ROOT, "**", "*.csv"), recursive=True))
print(f"فایل‌های csv یافت‌شده زیر {ROOT}: {len(cands)}")

results = {}
for f in cands:
    base = os.path.basename(f)
    if base.startswith("alpha_"):          # خروجی خودِ این اسکریپت
        continue
    try:
        df = pd.read_csv(f)
    except Exception as e:
        print(f"  ! {f}: خواندن نشد ({e})"); continue
    if not need.issubset(df.columns):
        print(f"  - رد شد (ستون MK ندارد): {f}"); continue

    tag = base.replace(".csv", "")
    g = df[["Pn", "Ps", "Dn", "Ds"]].to_numpy(float)
    g = g[~np.isnan(g).any(axis=1)]
    Pn, Ps, Dn, Ds = g.sum(0); ng = len(g)
    a_obs = alpha_pooled(Pn, Ps, Dn, Ds)

    idx  = rng.integers(0, ng, size=(NBOOT, ng))
    boot = g[idx].sum(axis=1)
    bPn, bPs, bDn, bDs = boot.T
    valid = (bPs > 0) & (bDn > 0)
    a_boot = np.full(NBOOT, np.nan)
    a_boot[valid] = 1.0 - (bPn[valid] * bDs[valid]) / (bPs[valid] * bDn[valid])
    ab = a_boot[~np.isnan(a_boot)]
    lo, hi = np.percentile(ab, [2.5, 97.5])
    p_boot = float((ab <= 0).mean())

    or_mh, or_lo, or_hi, nstr = mh_or_rbg(g[:,2], g[:,3], g[:,0], g[:,1])
    f_ = lambda x: 1 - 1/x if (x is not None and x == x and x > 0) else np.nan
    a_mh, a_mh_lo, a_mh_hi = f_(or_mh), f_(or_lo), f_(or_hi)

    print(f"\n=== {tag} ===")
    print(f"  مسیر: {f}")
    print(f"  ژن‌ها: {ng} | Pn={Pn:.1f} Ps={Ps:.1f} Dn={Dn:.1f} Ds={Ds:.1f}")
    print(f"  ژن‌های بدون تفاوت ثابت (Dn+Ds=0): {(g[:,2]+g[:,3] == 0).sum()}")
    print(f"  alpha (pooled)       = {a_obs:.3f}")
    print(f"  بوت‌استرپ ژنی 95% CI  = [{lo:.3f}, {hi:.3f}]"
          f"  (میانه {np.median(ab):.3f}, سوگیری {np.mean(ab)-a_obs:+.3f})")
    print(f"  سهم تکرارهای alpha<=0 = {p_boot:.4f} | نامعتبر: {NBOOT-len(ab)}")
    print(f"  alpha_MH             = {a_mh:.3f}  RBG 95% CI = [{a_mh_lo:.3f}, {a_mh_hi:.3f}]")
    print(f"  OR_MH                = {or_mh:.3f}  CI = [{or_lo:.3f}, {or_hi:.3f}] (لایه: {nstr})")

    results[tag] = dict(source=f, n_genes=int(ng), Pn=Pn, Ps=Ps, Dn=Dn, Ds=Ds,
                        alpha=a_obs, boot_lo=lo, boot_hi=hi,
                        boot_median=float(np.median(ab)), boot_p_le0=p_boot,
                        n_invalid=int(NBOOT-len(ab)),
                        alpha_MH=a_mh, MH_lo=a_mh_lo, MH_hi=a_mh_hi,
                        OR_MH=or_mh, OR_lo=or_lo, OR_hi=or_hi)
    np.save(f"{OUT}/{tag}_alpha_bootstrap.npy", ab)

if not results:
    print("\n⚠ هیچ فایلی با ستون‌های Pn/Ps/Dn/Ds/Orthogroup پیدا نشد.")
else:
    pd.DataFrame(results).T.to_csv(f"{OUT}/alpha_confidence_intervals.csv")
    json.dump(results, open(f"{OUT}/alpha_confidence_intervals.json","w"), indent=2)
    print(f"\n✓ {len(results)} تحلیل ذخیره شد در {OUT}/alpha_confidence_intervals.csv"
          f" (seed={SEED}, B={NBOOT})")
