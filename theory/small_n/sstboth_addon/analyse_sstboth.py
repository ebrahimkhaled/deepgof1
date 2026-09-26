"""analyse_sstboth.py -- Studies 1 and 2 re-analysed with the Stukel arm in the SstBoth form, by the SAME rules as
02_analyse.py and 04_analyse_trim.py (PREDECLARATION.md, PREDECLARATION_TRIM.md), side by side with the stored
joint-form results. Pairing: every regenerated dataset must reproduce the stored HL and joint-Stukel p-values.

Run:  python analyse_sstboth.py   -> RESULTS_SSTBOTH.json
"""
import os, glob, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); SN = os.path.dirname(HERE)
A = 0.05
def load(d):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, "*.csv")))], ignore_index=True)
def rate(x):
    x = np.asarray(x, bool); r = x.mean(); return r, np.sqrt(max(r * (1 - r), 1e-4) / len(x))

out = {"pairing": {}, "undefined": {}, "S1": {}, "S2": {}}
# ------------------------------------------------------------------ Study 1
O = load(os.path.join(SN, "rows")); N = load(os.path.join(HERE, "rows_s1"))
M = O.merge(N, on=["cell", "rep"], suffixes=("", "_new"))
assert len(M) == len(O) == len(N) == 11600, (len(M), len(O), len(N))
pair = dict(rows=len(M),
            HL_maxdiff=float((M.chk_p_HL - M.p_HL).abs().max()),
            Stukel_joint_full_maxdiff=float((M.chk_p_Stj - M.p_Stukel.fillna(1)).abs().max()),
            Stukel_joint_full_na_stored=int(M.p_Stukel.isna().sum()),
            Stukel_joint_loo_maxdiff=float((M.chk_loo_Stj - M.loo_Stukel).abs().max()),
            corrupted_equal=bool((M.corrupted.fillna(0) == M.corrupted_new.fillna(0)).all()),
            events_equal=bool((M.events == M.events_new).all()))
out["pairing"]["S1"] = pair; print("Study 1 pairing:", pair)
M["pS_j"] = M.p_Stukel; M["pS_m"] = M.p_Stm.fillna(1.0)       # undefined SstBoth: p = 1 (non-rejection)
M["lS_j"] = M.loo_Stukel; M["lS_m"] = M.loo_Stm
out["undefined"]["S1"] = dict(full_share=float(M.full_undef.mean()),
                              full_by_n={int(n): float(g.full_undef.mean()) for n, g in M.groupby("n")},
                              loo_any_share=float((M.loo_undef > 0).mean()),
                              loo_any_by_n={int(n): float((g.loo_undef > 0).mean()) for n, g in M.groupby("n")},
                              loo_deletions_share=float(M.loo_undef.sum() / M.n.sum()))
print("Study 1 SstBoth undefined:", out["undefined"]["S1"])
null = M[M.C == 0]; alt = M[M.C > 0]
q1 = {}; q2 = {}; pw = {}; ph = {}; dg = {}
print("\nQ1 size (joint -> SstBoth), * = above .075")
for n in (30, 50):
    clean = null[(null.n == n) & (null.type == "none")]
    for ty in ("none", "decimal", "sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        rj, sj = rate(s.pS_j <= A); rm, sm = rate(s.pS_m <= A)
        q1["%d|%s" % (n, ty)] = dict(joint=rj, sstboth=rm, se=sm, misbehaves_joint=bool(rj > .075), misbehaves_sstboth=bool(rm > .075))
        # Q2: LOO size vs clean base size (+2 SE), each form against its own clean size
        res = {}
        for f, pc, lc in (("joint", "pS_j", "lS_j"), ("sstboth", "pS_m", "lS_m")):
            c0, c0se = rate(clean[pc] <= A); r, se = rate(s[lc] <= A)
            lim = c0 + 2 * np.sqrt(c0se ** 2 + se ** 2)
            res[f] = dict(loo_size=r, clean_base=c0, limit=lim, holds=bool(r <= lim))
        q2["%d|%s" % (n, ty)] = res
        print("  n=%d %-8s size %.3f%s -> %.3f%s | LOO %.3f (clean %.3f %s) -> %.3f (clean %.3f %s)" %
              (n, ty, rj, "*" if rj > .075 else " ", rm, "*" if rm > .075 else " ",
               res["joint"]["loo_size"], res["joint"]["clean_base"], "ok" if res["joint"]["holds"] else "FAIL",
               res["sstboth"]["loo_size"], res["sstboth"]["clean_base"], "ok" if res["sstboth"]["holds"] else "FAIL"))
print("\npower base / LOO (joint | SstBoth)")
for n in (30, 50):
    for C in (0.5, 1.0, 1.5):
        for ty in ("none", "decimal"):
            s = alt[(alt.n == n) & (alt.C == C) & (alt.type == ty)]
            d = dict(joint=(rate(s.pS_j <= A)[0], rate(s.lS_j <= A)[0]), sstboth=(rate(s.pS_m <= A)[0], rate(s.lS_m <= A)[0]))
            pw["%d|%.1f|%s" % (n, C, ty)] = d
            print("  n=%d C=%.1f %-8s %.3f/%.3f | %.3f/%.3f" % (n, C, ty, *d["joint"], *d["sstboth"]))
print("\npost hoc calibrated alpha' (clean LOO size .05)")
for n in (30, 50):
    clean = null[(null.n == n) & (null.type == "none")]
    for f, pc, lc in (("joint", "pS_j", "lS_j"), ("sstboth", "pS_m", "lS_m")):
        ap = float(np.quantile(clean[lc], A))
        sz = {ty: rate(null[(null.n == n) & (null.type == ty)][lc] <= ap)[0] for ty in ("decimal", "sign", "outcome")}
        base = {ty: rate(null[(null.n == n) & (null.type == ty)][pc] <= A)[0] for ty in ("decimal", "sign", "outcome")}
        pwr = {C: (rate(alt[(alt.n == n) & (alt.C == C) & (alt.type == "none")][pc] <= A)[0],
                   rate(alt[(alt.n == n) & (alt.C == C) & (alt.type == "none")][lc] <= ap)[0]) for C in (0.5, 1.0, 1.5)}
        ph["%d|%s" % (n, f)] = dict(alpha_prime=ap, size=sz, base=base, power_base_vs_cal={str(k): v for k, v in pwr.items()})
        print("  n=%d %-7s alpha' %.3f  decimal size %.3f (base %.3f)" % (n, f, ap, sz["decimal"], base["decimal"]))
print("\ndiagnostic: most-rescuing deletion = corrupted record, among base rejections (SstBoth)")
for n in (30, 50):
    for ty in ("decimal", "sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        rj = s[s.pS_j <= A]; rm = s[s.pS_m <= A]
        hj = float((rj.amax_Stukel == rj.corrupted).mean()) if len(rj) else float("nan")
        hm = float((rm.amax_Stm == rm.corrupted).mean()) if len(rm) else float("nan")
        dg["%d|%s" % (n, ty)] = dict(joint=(len(rj), hj), sstboth=(len(rm), hm))
        print("  n=%d %-8s joint %3d rej, %.0f%% | SstBoth %3d rej, %.0f%%" % (n, ty, len(rj), 100 * hj, len(rm), 100 * hm))
out["S1"] = dict(Q1=q1, Q2=q2, power=pw, posthoc=ph, diagnostic=dg)

# ------------------------------------------------------------------ Study 2
O = load(os.path.join(SN, "rows_trim")); N = load(os.path.join(HERE, "rows_s2"))
M = O.merge(N, on=["cell", "rep"], suffixes=("", "_new"))
assert len(M) == len(O) == len(N) == 11600, (len(M), len(O), len(N))
pair = dict(rows=len(M), HL_maxdiff=float((M.chk_p_HL - M.p_HL).abs().max()),
            HL_trim_maxdiff=float((M.chk_t_HL - M.t_HL).abs().max()),
            Stukel_joint_full_maxdiff=float((M.chk_p_Stj.fillna(1) - M.p_Stukel).abs().max()),
            Stukel_joint_trim_maxdiff=float((M.chk_t_Stj.fillna(1) - M.t_Stukel).abs().max()),
            trimmed_equal=bool((M.trimmed == M.trimmed_new).all()), events_equal=bool((M.events == M.events_new).all()))
out["pairing"]["S2"] = pair; print("\nStudy 2 pairing:", pair)
M["pS_j"] = M.p_Stukel; M["tS_j"] = M.t_Stukel
M["pS_m"] = M.p_Stm.fillna(1.0); M["tS_m"] = M.t_Stm.fillna(1.0)
out["undefined"]["S2"] = dict(full_share=float(M.full_undef.mean()), trimmed_share=float(M.t_undef.mean()),
                              full_by_n={int(n): float(g.full_undef.mean()) for n, g in M.groupby("n")})
print("Study 2 SstBoth undefined:", out["undefined"]["S2"])
null = M[M.C == 0]; c1 = {}; c2 = {}; c3 = {}; c4 = {}
print("\n(1) protection: decimal slip, base -> trimmed, clean base (joint | SstBoth)")
for n in (30, 50):
    cl = null[(null.n == n) & (null.type == "none")]; de = null[(null.n == n) & (null.type == "decimal")]
    res = {}
    for f, pc, tc in (("joint", "pS_j", "tS_j"), ("sstboth", "pS_m", "tS_m")):
        c0, s0 = rate(cl[pc] <= A); b, _ = rate(de[pc] <= A); r, s = rate(de[tc] <= A)
        res[f] = dict(clean_base=c0, decimal_base=b, decimal_trim=r, se=s, holds=bool(r <= c0 + 2 * np.sqrt(s0 ** 2 + s ** 2)),
                      in_band=bool(.035 <= r <= .065))
    c1[str(n)] = res
    print("  n=%d joint %.3f -> %.3f (clean %.3f) %s | SstBoth %.3f -> %.3f (clean %.3f) %s" %
          (n, res["joint"]["decimal_base"], res["joint"]["decimal_trim"], res["joint"]["clean_base"], "ok" if res["joint"]["holds"] else "FAIL",
           res["sstboth"]["decimal_base"], res["sstboth"]["decimal_trim"], res["sstboth"]["clean_base"], "ok" if res["sstboth"]["holds"] else "FAIL"))
print("\n(2) no harm when clean (base/trimmed; ! = beyond 2 paired SE)")
for n in (30, 50):
    for C in (0.0, 0.5, 1.0, 1.5):
        s = M[(M.n == n) & (M.C == C) & (M.type == "none")]; res = {}
        for f, pc, tc in (("joint", "pS_j", "tS_j"), ("sstboth", "pS_m", "tS_m")):
            b, _ = rate(s[pc] <= A); r, _ = rate(s[tc] <= A)
            se = np.sqrt(max(np.mean((s[tc] <= A) != (s[pc] <= A)), 1e-4) / len(s))
            res[f] = dict(base=b, trim=r, paired_se=se, holds=bool(abs(r - b) <= 2 * max(se, 1e-3)))
        c2["%d|%.1f" % (n, C)] = res
        print("  n=%d C=%.1f joint %.3f/%.3f%s | SstBoth %.3f/%.3f%s" % (n, C, res["joint"]["base"], res["joint"]["trim"], "" if res["joint"]["holds"] else "!",
                                                                     res["sstboth"]["base"], res["sstboth"]["trim"], "" if res["sstboth"]["holds"] else "!"))
print("\n(3) mechanism: size among null decimal datasets whose slip was NOT deleted")
for n in (30, 50):
    m0 = null[(null.n == n) & (null.type == "decimal")]; mm = m0[m0.trimmed != m0.corrupted]
    c3[str(n)] = dict(n_missed=len(mm), joint=float((mm.tS_j <= A).mean()), sstboth=float((mm.tS_m <= A).mean()))
    print("  n=%d  %d missed: joint %.3f  SstBoth %.3f" % (n, len(mm), c3[str(n)]["joint"], c3[str(n)]["sstboth"]))
print("\n(4) sign / outcome (base -> trimmed)")
for n in (30, 50):
    for ty in ("sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        c4["%d|%s" % (n, ty)] = dict(joint=[float((s.pS_j <= A).mean()), float((s.tS_j <= A).mean())],
                                     sstboth=[float((s.pS_m <= A).mean()), float((s.tS_m <= A).mean())])
        print("  n=%d %-8s joint %.3f->%.3f | SstBoth %.3f->%.3f" % (n, ty, *c4["%d|%s" % (n, ty)]["joint"], *c4["%d|%s" % (n, ty)]["sstboth"]))
out["S2"] = dict(check1=c1, check2=c2, check3=c3, check4=c4)
json.dump(out, open(os.path.join(HERE, "RESULTS_SSTBOTH.json"), "w"), indent=1)
print("\nwrote RESULTS_SSTBOTH.json")
