"""make_tables.py -- every generated number of the paper's tables and text, from the deposited files.

Writes tables/tables_generated.tex (the LaTeX macros the manuscript \\input's) and
results/liu/TABLE2_SCORES.csv (the per-block size and matched power behind Table 3 and Figure 5).
Run from the archive root (or set DEEPGOF_ARCHIVE):  python tables/make_tables.py

  * The Liu benchmark (Table 3) is recomputed from the per-replicate p-values: each test is evaluated
    on the datasets where it returns a p-value (complete-case rule, uniform over every test); size =
    share of those nulls at p <= .05, matched power = randomized_matched_power(block nulls, cell
    alternatives), both per block (setting x n x corpus), then averaged. Under the original
    convention (NA counts as no rejection; columns *_na1 of TABLE2_SCORES.csv) the arms computed in
    the first run reproduce results/liu/AUGSCORE_INTERIM.csv to 1e-16 (asserted below).
  * DeepGOF-1 is bolded in Table 3 (label and numbers).
  * Arms: Liu et al. (2024)'s tests (HL, mHL, Osius-Rojek, Stukel = their SstBoth, Stute-Zhu, the
    projection test), the HL variants (Pigeon-Heyse, HL equal-width), GiViTI, BAGofT and DeepGOF-1.
    mHL and GiViTI come from theory/addon_arms/liu_rows (joined on kind, rep, C; HL re-checked on every
    row); the projection test is proj_liu_raw from theory/addon_arms/proj_rows_liu (Liu et al.'s
    statistic on the raw covariates). The first run's projection column (proj_superseded in
    augment_pvalues.csv, a different weight) is used only to check that the add-on rows are the same
    datasets: the recomputed first weight (proj_old) must equal it wherever it was recomputed.
    BAGofT is read from bagoft/confirm/.
  * An arm is "inflated" when its mean size exceeds SIZE_CAP = .075.
  * The FDIC table reports the published partition tests one by one.
"""
import glob
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import stdtr

PAP = os.path.abspath(os.environ.get("DEEPGOF_ARCHIVE", "."))
sys.path.insert(0, os.path.join(PAP, "benchmark", "common"))
from matched_power import randomized_matched_power, ALPHA  # noqa: E402

OUT = os.path.join(PAP, "tables", "tables_generated.tex")
RES = os.path.join(PAP, "results", "liu")
ADD = os.path.join(PAP, "theory", "addon_arms")
CONF = os.path.join(PAP, "bagoft", "confirm")
CS = [0.25, 0.5, 0.75, 1.0]

frag = []


def emit(name, body):
    frag.append("%% ---- %s ----\n\\newcommand{\\%s}{%%\n%s}\n" % (name, name, body))


# ================= Table 2: the Liu benchmark, from the per-replicate p-values ================
STORED = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Stukel", "Osius-Rojek"]
ARMS = STORED + ["Stute-Zhu", "DEEP", "mHL", "GiViTI", "proj"]
rows = []
for d in sorted(glob.glob(os.path.join(RES, "S*_n*_s*"))):
    tag = os.path.basename(d)
    S, N = tag.split("_")[0], int(tag.split("_")[1][1:])
    b = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    a = pd.read_csv(os.path.join(d, "augment_pvalues.csv")).rename(
        columns={"SZ": "Stute-Zhu", "proj_superseded": "proj_stored"})
    p = pd.read_csv(os.path.join(d, "deep_pvalues.csv"))
    for f in (b, a, p):
        f["_C"] = f["C"].fillna(-1.0)
    m = b.merge(a[["kind", "rep", "_C", "Stute-Zhu", "proj_stored"]], on=["kind", "rep", "_C"])
    m = m.merge(p[["kind", "rep", "_C", "p_deep"]], on=["kind", "rep", "_C"]).rename(
        columns={"p_deep": "DEEP"})
    ad = pd.read_csv(os.path.join(ADD, "liu_rows", tag + ".csv"))
    ad["_C"] = ad["C"].fillna(-1.0)
    m = m.merge(ad[["kind", "rep", "_C", "mHL", "GiViTI", "HL_check"]],
                on=["kind", "rep", "_C"], how="left")
    pr = pd.concat([pd.read_csv(f) for f in
                    glob.glob(os.path.join(ADD, "proj_rows_liu", tag, "chunk_*.csv"))])
    pr = pr[pr.status == "ok"].copy()
    pr["_C"] = pr["C"].fillna(-1.0)
    m = m.merge(pr[["kind", "rep", "_C", "proj_liu_raw", "proj_old"]],
                on=["kind", "rep", "_C"], how="left").rename(columns={"proj_liu_raw": "proj"})
    m = m[m.kind != 1]
    # pairing checks: the add-on rows are the stored datasets
    c = m[m.HL_check.notna()]
    assert np.allclose(c.HL_check, c.HL, rtol=0, atol=1e-12), "HL pairing failed in " + tag
    c = m[m.proj_old.notna()]
    assert np.allclose(c.proj_old, c.proj_stored, rtol=0, atol=1e-12), "proj pairing " + tag
    nul = m[m.kind == 0]
    for arm in ARMS:
        pn = nul[arm].to_numpy(float)
        pnc = pn[np.isfinite(pn)]
        # COMPLETE-CASE RULE, uniform over every test: each test is evaluated on the datasets where
        # it returns a p-value (nulls and alternatives alike). The *_na1 columns keep the original
        # convention (no p-value = no rejection) for the reproduction check against
        # AUGSCORE_INTERIM.csv and for the one sentence of the paper that reports it.
        size = float(np.mean(pnc <= ALPHA))
        size_na1 = float(np.mean(np.nan_to_num(pn, nan=1.0) <= ALPHA))
        for cc in CS:
            pa = m[(m.kind == 2) & np.isclose(m["_C"], cc)][arm].to_numpy(float)
            pac = pa[np.isfinite(pa)]
            rows.append(dict(block=tag, setting=S, n=N, arm=arm, C=cc, size=size,
                             matched=randomized_matched_power(pnc, pac)[0],
                             raw=float(np.mean(pac <= ALPHA)),
                             size_na1=size_na1,
                             matched_na1=randomized_matched_power(pn, pa)[0]))
R = pd.DataFrame(rows)
assert R.block.nunique() == 24

# the stored arms must reproduce the deposited scores exactly (under the original convention)
old = pd.read_csv(os.path.join(RES, "AUGSCORE_INTERIM.csv"))
j = old.merge(R[R.arm != "proj"], on=["block", "arm", "C"], suffixes=("_o", "_n"))  # proj is the corrected arm
assert len(j) == 24 * 4 * 7, len(j)          # HL, HLw, PH, Stukel, OR, SZ, DEEP
assert (j.size_o - j.size_na1).abs().max() < 1e-12 and (j.matched_o - j.matched_na1).abs().max() < 1e-12
# what the complete-case rule changes, test by test (printed for the record)
_cc = R.groupby("arm")[["size", "size_na1", "matched", "matched_na1"]].mean()
print("complete-case vs no-p-value-as-no-rejection (overall):\n", _cc.round(4).to_string())

FAMILY = {"proj": "covariate-space", "Stukel": "link", "Stute-Zhu": "covariate-space",
          "DEEP": "partition (learned)", "HL": "partition", "Pigeon-Heyse": "partition",
          "HL-equalwidth": "partition", "Osius-Rojek": "global", "mHL": "partition",
          "GiViTI": "calibration curve", "BAGofT": "partition (adaptive)"}
NAME = {"proj": "projection", "DEEP": "DeepGOF-1",
        "Stukel": "Stukel"}   # the SstBoth form of Liu et al. (2024) is named in the paper's text

size = R.groupby("arm")["size"].mean()
pw_n = R.pivot_table(index="arm", columns="n", values="matched", aggfunc="mean")
overall = R.groupby("arm")["matched"].mean()
cellsize = R.drop_duplicates(["block", "arm"]).groupby(["arm", "setting", "n"])["size"].mean()
cellmax = cellsize.groupby("arm").max()

# BAGofT joins from its own confirmatory run (bagoft/confirm/)
crow = pd.read_csv(os.path.join(CONF, "CONFIRM_ROW.csv")).iloc[0]
size["BAGofT"] = crow["size"]
overall["BAGofT"] = crow["overall"]
for _n in (50, 100, 500):
    pw_n.loc["BAGofT", _n] = crow["pw%d" % _n]
sz = pd.read_csv(os.path.join(CONF, "SIZE_BY_N.csv"))
_sz_sets = sz[(sz.arm == "BAGofT") & (sz.setting != "all") & (sz.setting != "trend")]
cellmax["BAGofT"] = float(_sz_sets["size"].max())
overall = overall.sort_values(ascending=False)

SIZE_CAP = 0.075
inflated = [a for a in overall.index if size[a] > SIZE_CAP]
valid = [a for a in overall.index if a not in inflated]
assert "DEEP" in valid and "mHL" in valid and inflated == ["Osius-Rojek"]


def liu_row(arm):
    nm = NAME.get(arm, arm)
    if arm in inflated:
        nm += " (inflated)"
    # DeepGOF-1 is bolded: its label and every number
    b = (lambda s: "\\textbf{%s}" % s) if arm == "DEEP" else (lambda s: s)
    return "%s & %s & %s & %s & %s & %s & %s \\\\" % (
        b(nm), b("%.3f" % size[arm]), FAMILY[arm], b("%.3f" % pw_n.loc[arm, 50]),
        b("%.3f" % pw_n.loc[arm, 100]), b("%.3f" % pw_n.loc[arm, 500]), b("%.3f" % overall[arm]))


emit("tableliu", "\n".join([liu_row(a) for a in valid] + ["\\midrule"]
                           + [liu_row(a) for a in inflated]) + "\n")

ov_valid = overall[valid]
_rankvalid = int((ov_valid > overall["DEEP"]).sum() + 1)
emit("rankdeepvalid", str(_rankvalid))
emit("narmsvalid", str(len(valid)))
emit("narmsinfl", str(len(inflated)))
CARD = ("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
        "fifteen sixteen seventeen eighteen nineteen twenty").split()
ORD = ("zeroth first second third fourth fifth sixth seventh eighth ninth tenth eleventh twelfth "
       "thirteenth fourteenth fifteenth sixteenth seventeenth eighteenth nineteenth "
       "twentieth").split()
emit("narmsword", CARD[len(overall)])
emit("narmsvalidword", CARD[len(valid)])
emit("narmsinflword", CARD[len(inflated)])
emit("rankdeepvalidword", ORD[_rankvalid])
_rk = {}
for _n, _nm in ((50, "fifty"), (100, "hundred"), (500, "fivehundred")):
    _rk[_n] = int((pw_n.loc[valid, _n] > pw_n.loc["DEEP", _n]).sum() + 1)
    emit("rankdeepn" + _nm, str(_rk[_n]))
    emit("rankdeepn" + _nm + "word", ORD[_rk[_n]])
assert _rk[50] == _rk[100], "n=50 and n=100 ranks differ: rewrite the 'n<=100' sentence"
emit("powdeep", "%.3f" % overall["DEEP"])
emit("rankdeep", str(int((overall > overall["DEEP"]).sum() + 1)))
emit("narms", str(len(overall)))
emit("sizecap", "%.3f" % SIZE_CAP)

# the calibrated arms ahead of DeepGOF-1 overall, in order (the text names them)
ahead = [a for a in valid if overall[a] > overall["DEEP"]]
print("calibrated tests ahead of DeepGOF-1 overall:", ahead)
assert ahead == ["proj"], ahead       # under the complete-case rule only the projection test
# the one sentence on the original convention (no p-value = no rejection): Stukel's overall and rank
_ov1 = R.groupby("arm")["matched_na1"].mean()
_ov1["BAGofT"] = crow["overall"]
_sz1 = R.groupby("arm")["size_na1"].mean()
_sz1["BAGofT"] = crow["size"]
_valid1 = [a for a in _ov1.index if _sz1[a] <= SIZE_CAP]
emit("powstukeldeclared", "%.3f" % _ov1["Stukel"])
emit("rankstukeldeclaredword", ORD[int((_ov1[_valid1] > _ov1["Stukel"]).sum() + 1)])
emit("powproj", "%.3f" % overall["proj"])
emit("powstukel", "%.3f" % overall["Stukel"])
emit("sizeproj", "%.3f" % size["proj"])
emit("sizestukel", "%.3f" % size["Stukel"])
emit("sizegiviti", "%.3f" % size["GiViTI"])
emit("powgiviti", "%.3f" % overall["GiViTI"])
emit("sizemhl", "%.3f" % size["mHL"])
emit("sizemhlmax", "%.3f" % cellmax["mHL"])
_mc = cellsize.loc["mHL"]
assert _mc.idxmax() == ("S1", 500)
emit("sizemhlsfour", "%.3f" % _mc.loc[("S4", 500)])
emit("sizegivitimax", "%.3f" % cellmax["GiViTI"])

# family margins (DeepGOF-1 minus each rival, per n)
dg = pw_n.loc["DEEP"]


def marg(r):
    return "%.2f/%.2f/%.2f" % tuple(dg[n] - pw_n.loc[r, n] for n in (50, 100, 500))


emit("margHL", marg("HL"))
emit("margHLW", marg("HL-equalwidth"))
emit("margSZ", marg("Stute-Zhu"))
emit("margGiViTI", marg("GiViTI"))
emit("gaptoproj", "%.2f/%.2f/%.2f" % tuple(pw_n.loc["proj", n] - dg[n] for n in (50, 100, 500)))
emit("gapstukelfivehundred", "%.2f" % (pw_n.loc["Stukel", 500] - dg[500]))

# ---- Monte Carlo error of the overall column: paired bootstrap over nulls and alternatives ----
# (resample each block's nulls and each cell's alternatives with replacement, the same indices for
# every arm; 200 replicates). This is the estimator behind \deepseboot in the earlier versions.
BOOT_ARMS = ["DEEP", "Stukel", "proj", "GiViTI"]
blocks = {}
for d in sorted(glob.glob(os.path.join(RES, "S*_n*_s*"))):
    tag = os.path.basename(d)
    b = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    p = pd.read_csv(os.path.join(d, "deep_pvalues.csv"))
    for f in (b, p):
        f["_C"] = f["C"].fillna(-1.0)
    m = b.merge(p[["kind", "rep", "_C", "p_deep"]], on=["kind", "rep", "_C"]).rename(
        columns={"p_deep": "DEEP"})
    a = pd.read_csv(os.path.join(d, "augment_pvalues.csv"))
    a["_C"] = a["C"].fillna(-1.0)
    m = m.merge(a[["kind", "rep", "_C"]], on=["kind", "rep", "_C"])
    ad = pd.read_csv(os.path.join(ADD, "liu_rows", tag + ".csv"))
    ad["_C"] = ad["C"].fillna(-1.0)
    m = m.merge(ad[["kind", "rep", "_C", "GiViTI"]], on=["kind", "rep", "_C"], how="left")
    pr = pd.concat([pd.read_csv(f) for f in
                    glob.glob(os.path.join(ADD, "proj_rows_liu", tag, "chunk_*.csv"))])
    pr = pr[pr.status == "ok"].copy()
    pr["_C"] = pr["C"].fillna(-1.0)
    m = m.merge(pr[["kind", "rep", "_C", "proj_liu_raw"]], on=["kind", "rep", "_C"],
                how="left").rename(columns={"proj_liu_raw": "proj"})
    blocks[tag] = (m[m.kind == 0][BOOT_ARMS].to_numpy(float),
                   [m[(m.kind == 2) & np.isclose(m._C, cc)][BOOT_ARMS].to_numpy(float) for cc in CS])
rng = np.random.default_rng(20260925)
NB = 200
boot = np.zeros((NB, len(BOOT_ARMS)))
for r in range(NB):
    acc = np.zeros(len(BOOT_ARMS))
    k = 0
    for tag, (N0, AL) in blocks.items():
        i0 = rng.integers(0, len(N0), len(N0))
        for A in AL:
            ia = rng.integers(0, len(A), len(A))
            for jj in range(len(BOOT_ARMS)):
                _n0, _a = N0[i0, jj], A[ia, jj]          # complete-case rule, as in the table
                acc[jj] += randomized_matched_power(_n0[np.isfinite(_n0)], _a[np.isfinite(_a)])[0]
            k += 1
    boot[r] = acc / k
bse = dict(zip(BOOT_ARMS, boot.std(axis=0, ddof=1)))
d_st = boot[:, 0] - boot[:, 1]
emit("diffstukel", "%.3f" % (overall["Stukel"] - overall["DEEP"]))
emit("sediffstukel", "%.3f" % d_st.std(ddof=1))
emit("zdiffstukel", "%.1f" % ((overall["Stukel"] - overall["DEEP"]) / d_st.std(ddof=1)))
d_gv = boot[:, 0] - boot[:, 3]
emit("diffgiviti", "%.3f" % (overall["DEEP"] - overall["GiViTI"]))
emit("sediffgiviti", "%.3f" % d_gv.std(ddof=1))
print("paired bootstrap SEs:", {k: round(v, 4) for k, v in bse.items()},
      "| DEEP-Stukel SE %.4f | DEEP-GiViTI SE %.4f" % (d_st.std(ddof=1), d_gv.std(ddof=1)))

# ---- level stability: per-block SD, and paired t on per-block distance from nominal ---------
bb = R.drop_duplicates(["block", "arm"]).pivot(index="block", columns="arm", values="size")
sd = bb.std()
emit("sdsizedeep", "%.3f" % sd["DEEP"])
emit("sdsizeproj", "%.3f" % sd["proj"])
emit("sdsizestukel", "%.3f" % sd["Stukel"])
calib_blocks = [a for a in valid if a in bb.columns and a != "DEEP"]
assert sd["DEEP"] == sd[[a for a in bb.columns]].min()
tt = {}
for a in calib_blocks:
    dd = (bb["DEEP"] - .05).abs() - (bb[a] - .05).abs()
    t = dd.mean() / (dd.std(ddof=1) / np.sqrt(len(dd)))
    tt[a] = (t, 2 * stdtr(len(dd) - 1, -abs(t)), dd.mean())
tmin = max(tt.values(), key=lambda v: v[0])       # the weakest (least negative) t
tmax = min(tt.values(), key=lambda v: v[0])
assert all(v[0] < 0 for v in tt.values())
emit("tstatweak", "%.1f" % tmin[0])
emit("tstatstrong", "%.1f" % tmax[0])
emit("pstatweak", ("$%.3f$" % tmin[1]).replace("0.", "."))
_weakarm = max(tt, key=lambda k: tt[k][0])
assert _weakarm == "mHL", _weakarm       # the text names the arm with the weakest t
emit("tstatproj", "%.1f" % tt["proj"][0])
emit("distgainproj", "%.3f" % -tt["proj"][2])
print("paired t (DEEP vs calibrated arms):", {a: round(v[0], 2) for a, v in tt.items()})

# ---- Setting 2: rise in matched-level power from null to the strongest departure --------------
s2 = R[(R.setting == "S2") & (R.C == 1.0)].groupby("arm")["matched"].mean() - ALPHA
for arm, mac in (("proj", "riseproj"), ("Stute-Zhu", "risestutezhu"), ("HL", "risehl"),
                 ("mHL", "risemhl"), ("Stukel", "risestukel"), ("Osius-Rojek", "riseosius"),
                 ("GiViTI", "risegiviti"), ("Pigeon-Heyse", "riseph"),
                 ("HL-equalwidth", "risehlw"), ("DEEP", "risedeep")):
    emit(mac, "%+.2f" % s2[arm])
# BAGofT's Setting 2 rise is read from the confirmatory recomputation, on its own raw scale
cf = pd.read_csv(os.path.join(CONF, "REVIEW_RECOMPUTED.csv")).set_index("name")["value"]
emit("risebagoft", "%+.2f" % float(cf["rise_BAGofT"]))
print("Setting 2 rises:", s2.round(3).to_dict())

# ---- one corpus against two: how much it moves the overall column ----------------------------
R["corpus"] = R.block.str.split("_s").str[1]
one = R[R.corpus == "101000"].groupby("arm")["matched"].mean()
emit("corpusgap", "%.3f" % float((one - R.groupby("arm")["matched"].mean()).abs().max()))

# ---- why the Osius-Rojek level fails: the design, not the sample size -------------------------
orow = R[R.arm == "Osius-Rojek"].drop_duplicates(subset=["block"])[["block", "setting", "n", "size"]]
good = orow[orow.setting.isin(["S1", "S2"])]["size"]
emit("osiussizegood", "$%.3f$--$%.3f$" % (good.min(), good.max()))
emit("osiussizebadthree", "%.3f" % orow[(orow.setting == "S3") & (orow.n == 500)]["size"].mean())
emit("osiussizebadfour", "%.3f" % orow[(orow.setting == "S4") & (orow.n == 500)]["size"].mean())
ov = pd.read_csv(os.path.join(RES, "OSIUS_VARIANCE.csv"))
emit("osiusvargood", "$%.0f$--$%.0f$" % (ov[ov.setting <= 2].pct_var_left.min(),
                                         ov[ov.setting <= 2].pct_var_left.max()))
emit("osiusvarthree", "%.1f" % ov[ov.setting == 3].pct_var_left.median())
emit("osiusvarfour", "%.1f" % ov[ov.setting == 4].pct_var_left.median())
emit("sizedeep", "%.4f" % size["DEEP"])
emit("sizeosius", "%.3f" % size["Osius-Rojek"])

# the per-block scores, kept for Figure 8 and for checking
R.to_csv(os.path.join(RES, "TABLE2_SCORES.csv"), index=False)

# ================= Table: the confirmatory BAGofT comparison (unchanged) =======================
cp = pd.read_csv(os.path.join(CONF, "CONFIRM_PAIRED.csv")).set_index("group")
PRIMARY = "all 48 power cells (O2, primary)"
BGROWS = [(PRIMARY, "all forty-eight power cells"),
          ("n = 50", "$n=50$"), ("n = 100", "$n=100$"), ("n = 500", "$n=500$"),
          ("Setting 1", "Setting 1"), ("Setting 2", "Setting 2"),
          ("Setting 3", "Setting 3"), ("Setting 4", "Setting 4"),
          ("C = 0.25", "$C=0.25$"), ("C = 0.5", "$C=0.5$"),
          ("C = 0.75", "$C=0.75$"), ("C = 1", "$C=1$")]


def pfmt(p):
    if p >= 1e-4:
        return ("$%.4f$" % p).replace("0.", ".")
    e = int(np.floor(np.log10(p)))
    return "$%.1f\\times10^{%d}$" % (p / 10.0 ** e, e)


rows = []
for key, label in BGROWS:
    r = cp.loc[key]
    p = r["p_exact"] if key == PRIMARY else r["p_holm"]
    rows.append("%s & %d & %.3f & %.3f & %d:%d & %s \\\\" % (
        label, int(r["pairs"]), r["bagoft_rate"], r["deepgof_rate"],
        int(r["only_bagoft"]), int(r["only_deepgof"]), pfmt(p)))
    if key == PRIMARY:
        rows.append("\\midrule")
emit("tablebagoft", "\n".join(rows) + "\n")

pr_ = cp.loc[PRIMARY]
emit("bagpow", "%.3f" % crow["overall"])
emit("bagsize", "%.3f" % crow["size"])
emit("bagpowfifty", "%.3f" % crow["pw50"])
emit("bagpairs", "{:,}".format(int(pr_["pairs"])))
emit("bagse", "%.3f" % crow["se_overall"])
emit("bagdiscord", "%d:%d" % (int(pr_["only_bagoft"]), int(pr_["only_deepgof"])))
emit("bagmcnemar", pfmt(pr_["p_exact"]))
emit("bagptwo", pfmt(cp.loc["Setting 2", "p_holm"]))
emit("bagpnfifty", pfmt(cp.loc["n = 50", "p_holm"]))
szall = sz[sz.setting == "all"]
emit("bagnullsize", "%.3f" % crow["size"])
emit("deepnullsize", "%.3f" % float(np.average(szall[szall.arm == "DeepGOF-1"]["size"],
                                               weights=szall[szall.arm == "DeepGOF-1"]["nulls"])))
emit("bagsdsize", "%.3f" % float(_sz_sets["size"].std()))
emit("bagsizelo", "%.3f" % float(_sz_sets["size"].min()))
emit("bagsizehi", "%.3f" % float(_sz_sets["size"].max()))
tr = sz[(sz.arm == "BAGofT") & (sz.setting == "trend")].iloc[0]
for _n, _nm in ((50, "fifty"), (100, "hundred"), (500, "fivehundred")):
    emit("bagsizen" + _nm, "%.3f" % float(szall[(szall.arm == "BAGofT") & (szall.n == _n)]["size"].iloc[0]))
emit("bagslope", "%.2f" % tr["slope_logn"])
emit("bagslopeci", "$%.2f$--$%.2f$" % (tr["lo"], tr["hi"]))
emit("bagslopep", ("$%.3f$" % tr["p_trend"]).replace("0.", "."))
trd = sz[(sz.arm == "DeepGOF-1") & (sz.setting == "trend")].iloc[0]
emit("deepslopep", ("$%.2f$" % trd["p_trend"]).replace("0.", "."))


def cfnum(k, fmt="%.3f"):
    return fmt % float(cf[k])


emit("bagmatchedbag", cfnum("matched_power_bagoft"))
emit("bagmatcheddeep", cfnum("matched_power_deep"))
emit("bagmatchedgap", cfnum("matched_gap"))
emit("bagmatcheddiscord", "%d:%d" % (int(cf["matched_only_bagoft"]), int(cf["matched_only_deep"])))
emit("bagmatchedp", pfmt(float(cf["matched_mcnemar_p"])))
emit("bagrawgap", cfnum("raw_gap"))
emit("bagclusters", "{:,}".format(int(cf["n_clusters"])))
emit("bagmcnemarcl", pfmt(float(cf["cluster_p"])))
emit("bagseboot", cfnum("bag_se_bootstrap", "%.3f"))
emit("deepseboot", cfnum("deep_se_bootstrap", "%.3f"))
emit("bagsetimes", cfnum("se_ratio", "%.1f"))
_pkg = pd.read_csv(os.path.join(PAP, "bagoft", "explore", "EXPLORATORY_FULL_TABLE.csv"))
emit("bagsecspkglo", "%.0f" % _pkg.groupby("n").secs.median().min())
emit("bagsecspkghi", "{:,.0f}".format(_pkg.groupby("n").secs.median().max()))
emit("bagsdsame", cfnum("bag_sd_12", "%.3f"))
emit("deepsdsame", cfnum("deep_sd_12", "%.3f"))
emit("noisefloor", cfnum("noise_floor_200", "%.3f"))
emit("bagchisq", cfnum("bag_chisq", "%.1f"))
emit("bagchisqp", ("$%.2f$" % float(cf["bag_chisq_p"])).replace("0.", "."))
emit("bagstrictfifty", cfnum("strict_size_50"))
emit("bagstricthundred", cfnum("strict_size_100"))
emit("bagstrictfivehundred", cfnum("strict_size_500"))
emit("bagslopepstrict", ("$%.3f$" % float(cf["strict_slope_p"])).replace("0.", "."))
emit("gridnullstotal", "{:,}".format(int(cf["grid_nulls_total"])))
emit("gridnullsheld", "{:,}".format(int(cf["grid_nulls_held"])))
emit("gridse", cfnum("grid_cell_se", "%.3f"))
ct = pd.read_csv(os.path.join(CONF, "COST.csv")).set_index("n")["median_secs"]
emit("bagsecslo", "%.0f" % ct.min())
emit("bagsecshi", "%.0f" % ct.max())

# ================= Table: certificate margins (verified in benchmark/theory/theory_verify_numbers.R) ===
cert = [
    ("step / threshold", 2.200, 1.075),
    ("curvature stripe (omitted $x^{2}$)", 1.729, 1.644),
    ("local bump", 1.506, 1.493),
    ("diffuse calibration drift", 1.464, 1.463),
    ("interaction saddle", 1.188, 1.141),
]
rows = ["%s & $%.3f$ & $%.3f$ \\\\" % r for r in cert]
rows.append("\\midrule")
rows.append("most negative direction found (checkerboard) & \\multicolumn{2}{c}{$-0.281$} \\\\")
emit("tablecert", "\n".join(rows) + "\n")

# ================= Table: FDIC application, the published partition tests one by one ==========
FD = os.path.join(PAP, "results", "fdic")
ob = pd.read_csv(os.path.join(FD, "observed_pvalues.csv"))
dgo = pd.read_csv(os.path.join(FD, "deepgof_observed.csv"))
dgs = pd.read_csv(os.path.join(FD, "deepgof_size.csv"))
ob = ob.merge(dgo[["state", "p_deepgof"]], on="state").merge(
    dgs[["state", "size"]].rename(columns={"size": "dg_size"}), on="state")
PART = ["HL", "Pigeon-Heyse", "HL-equalwidth"]
fsize = {}
for s in ob.state:
    P = pd.read_csv(os.path.join(FD, "boot_pvalues_%s.csv" % s))
    fsize[s] = {t: float((P[t].fillna(1.0) <= ALPHA).mean()) for t in PART}   # NA = no rejection
rows = []
for _, r in ob.sort_values("n", ascending=False).iterrows():
    f = fsize[r.state]
    rows.append("%s & %d & %d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\" % (
        r.state, r.n, r.events, f["HL"], f["Pigeon-Heyse"], f["HL-equalwidth"],
        min(r.HL, r["Pigeon-Heyse"]), r["HL-equalwidth"], r.p_deepgof, r.dg_size))
emit("tablefdic", "\n".join(rows) + "\n")
allsz = pd.DataFrame(fsize).T
emit("fdicpartmax", "%.3f" % allsz.to_numpy().max())
emit("fdicpartga", "$%.3f$--$%.3f$" % (allsz.loc["GA"].min(), allsz.loc["GA"].max()))
emit("fdicdglo", "%.4f" % ob.dg_size.min())
emit("fdicdghi", "%.4f" % ob.dg_size.max())
emit("fdicninside", str(int(((allsz >= .035) & (allsz <= .065)).to_numpy().sum())))
emit("fdicncells", str(int(allsz.size)))

# ================= grid mean of the deployed level (Section 5) =================================
v = pd.read_csv(os.path.join(PAP, "results", "s3", "VERDICT.csv"))
vd = v[v["arm"] == "DEEP"]
assert len(vd) == 60
emit("sizegrid", "%.4f" % vd["raw_size"].mean())

io.open(OUT, "w", encoding="utf-8", newline="\n").write(
    "%% GENERATED by tables/make_tables.py -- do not edit by hand.\n" + "\n".join(frag))
print("wrote", OUT)
print("valid:", [(a, round(overall[a], 4), round(size[a], 4)) for a in valid])
print("inflated:", [(a, round(overall[a], 4), round(size[a], 4), round(cellmax[a], 3)) for a in inflated])
print("rank DEEP among calibrated %d of %d; by n %s" % (_rankvalid, len(valid), _rk))
print("FDIC partition sizes:\n", allsz.round(4))
