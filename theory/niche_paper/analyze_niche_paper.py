"""analyze_niche_paper.py -- the niche simulation study as reported in the audit manuscript.

Reads only:
  theory/niche_paper/NICHE_ROWS.csv             one row per dataset (16,000): design keys, DeepGOF-1 (DG),
                                                Hosmer-Lemeshow (HL), mHL, GiViTI, Stukel_marginal (Liu et al.'s
                                                SstBoth, the paper's Stukel), proj_liu_raw (Liu et al.'s statistic,
                                                raw covariates). The table joins the stored rows of the run with
                                                the add-on rows of theory/addon_arms (niche_rows/, proj_rows_niche/)
                                                on (scen, rep); the join checked the event count and HL on every
                                                row, and the recomputed first projection weight against the stored
                                                one wherever it was recomputed.
  theory/addon_arms/niche_rows_bagoft/*.csv     BAGofT p (p.value), n = 200 only, reps 1-100

Estimands (the logic of theory/niche/analyze_niche.R):
  * NA -> p = 1 (a failed test cannot reject), for nulls and alternatives alike.
  * size = raw rejection at .05 on the 2,000 nulls at each n, binomial SE.
  * matched power: each arm randomized to exactly .05 on its OWN nulls at the same n (matched_phi,
    a line-for-line port of analyze_niche.R, itself a port of probit_common.randomized_matched_power);
    per-replicate phi in {0, gamma, 1}.
  * paired difference DG - arm = mean(phi_DG - phi_arm), paired SE = sd(.)/sqrt(R).
  * decision rule (stated in the paper): in a cell, DG exceeds EVERY other arm by more than 2 paired SE;
    overall, in no cell is any arm above DG by more than 2 paired SE. Arms in the rule: HL, mHL, Stukel,
    GiViTI, projection. BAGofT (100 reps, n = 200 only) is reported separately, not in the rule.
Writes theory/niche_paper/{niche_power.csv, niche_diff.csv, niche_size.csv, niche_bagoft.csv, niche_rule.csv,
niche_na_share.csv, RESULT.md} and niche_generated.tex (macros and table body the manuscript \\input's).
Run:  python theory/niche_paper/analyze_niche_paper.py   (any working directory)
"""
import glob
import io
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ADD = os.path.join(os.path.dirname(HERE), "addon_arms")
TEX = os.path.join(HERE, "niche_generated.tex")
ALPHA = 0.05


def read_dir(d, pat="*.csv"):
    # keep_default_na=False: the family "null" must not be read as a missing value
    return pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", "", "NaN"])
                      for f in sorted(glob.glob(os.path.join(d, pat)))], ignore_index=True)


D = pd.read_csv(os.path.join(HERE, "NICHE_ROWS.csv"), keep_default_na=False, na_values=["NA", "", "NaN"])
D = D[["scen", "fam", "C", "n", "rep", "events", "HL", "DG", "Stukel_marginal", "mHL", "GiViTI", "proj_liu_raw"]]
assert len(D) == 16000 and not D.duplicated(["scen", "rep"]).any()
D = D.rename(columns={"Stukel_marginal": "Stukel", "proj_liu_raw": "proj"})
ARMS = ["DG", "proj", "GiViTI", "Stukel", "HL", "mHL"]
RIVALS = ARMS[1:]
na_share = D.groupby("scen")[ARMS].apply(lambda g: g.isna().mean())
for a in ARMS:
    D[a] = D[a].fillna(1.0)


def matched_phi(pn, pa, alpha=ALPHA):
    t = np.unique(pn)
    F = np.searchsorted(np.sort(pn), t, side="right") / len(pn)
    ok = np.where(F <= alpha)[0]
    if len(ok) == 0:
        g = alpha / F[0] if F[0] > 0 else 0.0
        return g * (pa <= t[0]), t[0], g
    i = ok[-1]
    t0, F0 = t[i], F[i]
    if i + 1 >= len(t):
        return (pa <= t0).astype(float), t0, 0.0
    t1, F1 = t[i + 1], F[i + 1]
    g = 0.0 if F1 <= F0 else (alpha - F0) / (F1 - F0)
    return (pa <= t0).astype(float) + g * ((pa > t0) & (pa <= t1)), t0, g


# ---- sizes ----
SZ = []
for n in (200, 500):
    Z = D[(D.fam == "null") & (D.n == n)]
    assert len(Z) == 2000
    for a in ARMS:
        r = float((Z[a] <= ALPHA).mean())
        SZ.append(dict(n=n, arm=a, R=len(Z), size=r, se=np.sqrt(r * (1 - r) / len(Z))))
SZ = pd.DataFrame(SZ)

# ---- matched power and paired differences ----
cells = (D[D.fam != "null"][["fam", "C", "n"]].drop_duplicates()
         .sort_values(["fam", "n", "C"]).reset_index(drop=True))
assert len(cells) == 12
PW, DF = [], []
for _, cc in cells.iterrows():
    X = D[(D.fam == cc.fam) & (D.C == cc.C) & (D.n == cc.n)]
    N0 = D[(D.fam == "null") & (D.n == cc.n)]
    assert len(X) == 1000
    phi = {a: matched_phi(N0[a].to_numpy(float), X[a].to_numpy(float))[0] for a in ARMS}
    for a in ARMS:
        PW.append(dict(fam=cc.fam, C=cc.C, n=cc.n, arm=a, R=len(X), matched=phi[a].mean(),
                       se=phi[a].std(ddof=1) / np.sqrt(len(X)),
                       raw=float((X[a] <= ALPHA).mean())))
    for a in RIVALS:
        d = phi["DG"] - phi[a]
        se = d.std(ddof=1) / np.sqrt(len(X))
        DF.append(dict(fam=cc.fam, C=cc.C, n=cc.n, arm=a, d=d.mean(), se=se,
                       z=(d.mean() / se) if se > 0 else 0.0))
PW = pd.DataFrame(PW)
DF = pd.DataFrame(DF)

# ---- the decision rule ----
rule = []
for _, cc in cells.iterrows():
    Y = DF[(DF.fam == cc.fam) & (DF.C == cc.C) & (DF.n == cc.n)]
    closest = Y.loc[Y.d.idxmin()]
    rule.append(dict(fam=cc.fam, C=cc.C, n=cc.n, lead_all=bool((Y.z > 2).all()),
                     any_above=bool((Y.z < -2).any()), n_rivals_lead=int((Y.z > 2).sum()),
                     closest=closest.arm, closest_d=closest.d, closest_se=closest.se, closest_z=closest.z))
rule = pd.DataFrame(rule)
n_lead = int(rule.lead_all.sum())
n_above = int(rule.any_above.sum())
holds = n_lead >= 1 and n_above == 0
print(rule.to_string())
print("cells where DG leads every arm by > 2 paired SE: %d of 12; cells with an arm > 2 SE above DG: %d"
      % (n_lead, n_above))

# ---- BAGofT, n = 200, 100 reps per cell, reported separately ----
B = read_dir(os.path.join(ADD, "niche_rows_bagoft"), "*_r*.csv")
assert not B.duplicated(["scen", "rep"]).any()
B = B.merge(D[["scen", "rep", "events", "DG"]], on=["scen", "rep"], suffixes=("", "_st"))
assert (B.events == B.events_st).all()
B["p"] = B["p"].fillna(1.0)
BN = B[B.scen == "null_n200"]
DN = D[(D.fam == "null") & (D.n == 200)]
BG = []
for _, cc in cells[cells.n == 200].iterrows():
    X = B[(B.fam == cc.fam) & (B.C == cc.C)]
    phb, t0, g = matched_phi(BN.p.to_numpy(float), X.p.to_numpy(float))
    phd = matched_phi(DN.DG.to_numpy(float), X.DG.to_numpy(float))[0]   # DG at its own 2,000-null threshold
    d = phd - phb
    BG.append(dict(fam=cc.fam, C=cc.C, n=200, R=len(X), bag_matched=phb.mean(),
                   bag_se=phb.std(ddof=1) / np.sqrt(len(X)), bag_raw=float((X.p <= ALPHA).mean()),
                   dg_same=phd.mean(), d=d.mean(), d_se=d.std(ddof=1) / np.sqrt(len(X)), t0=t0, gamma=g))
BG = pd.DataFrame(BG)
bag_null = float((BN.p <= ALPHA).mean())
bag_null_se = np.sqrt(bag_null * (1 - bag_null) / len(BN))
print(BG.round(3).to_string())
print("BAGofT null size %.3f (%.3f) on %d nulls" % (bag_null, bag_null_se, len(BN)))

os.makedirs(HERE, exist_ok=True)
PW.to_csv(os.path.join(HERE, "niche_power.csv"), index=False)
DF.to_csv(os.path.join(HERE, "niche_diff.csv"), index=False)
SZ.to_csv(os.path.join(HERE, "niche_size.csv"), index=False)
BG.to_csv(os.path.join(HERE, "niche_bagoft.csv"), index=False)
rule.to_csv(os.path.join(HERE, "niche_rule.csv"), index=False)
na_share.to_csv(os.path.join(HERE, "niche_na_share.csv"))

# ---- the LaTeX fragment ----
frag = []


def emit(name, body):
    frag.append("%% ---- %s ----\n\\newcommand{\\%s}{%%\n%s}\n" % (name, name, body))


def f3(x):
    return ("%.3f" % x).replace("0.", ".", 1) if x >= 0 else ("-%.3f" % -x).replace("0.", ".", 1)


FAM = {"inter": "interaction", "quad": "quadratic"}
rows = []
for _, cc in cells.iterrows():
    z = PW[(PW.fam == cc.fam) & (PW.C == cc.C) & (PW.n == cc.n)].set_index("arm")
    rl = rule[(rule.fam == cc.fam) & (rule.C == cc.C) & (rule.n == cc.n)].iloc[0]
    bg = BG[(BG.fam == cc.fam) & (BG.C == cc.C)]
    bag = f3(bg.bag_matched.iloc[0]) if cc.n == 200 else "--"
    mark = "$^{*}$" if rl.lead_all else ""
    # DeepGOF-1 is bolded (its header in the manuscript, every number here)
    rows.append("%s & %.1f & %d & \\textbf{%s}%s & %s & %s & %s & %s & %s & %s \\\\" % (
        FAM[cc.fam], cc.C, cc.n, f3(z.loc["DG", "matched"]), mark, f3(z.loc["proj", "matched"]),
        f3(z.loc["GiViTI", "matched"]), f3(z.loc["Stukel", "matched"]), f3(z.loc["HL", "matched"]),
        f3(z.loc["mHL", "matched"]), bag))
rows.append("\\midrule")
for n in (200, 500):
    z = SZ[SZ.n == n].set_index("arm")
    bag = f3(bag_null) if n == 200 else "--"
    rows.append("null size & 0 & %d & \\textbf{%s} & %s & %s & %s & %s & %s & %s \\\\" % (
        n, f3(z.loc["DG", "size"]), f3(z.loc["proj", "size"]), f3(z.loc["GiViTI", "size"]),
        f3(z.loc["Stukel", "size"]), f3(z.loc["HL", "size"]), f3(z.loc["mHL", "size"]), bag))
emit("tableniche", "\n".join(rows) + "\n")

CARD = "zero one two three four five six seven eight nine ten eleven twelve".split()
emit("nichenlead", str(n_lead))
emit("nichenleadword", CARD[n_lead])
emit("nichenabove", str(n_above))
# DG minus each arm: range of paired differences, and how many cells DG leads each arm by > 2 SE
for a, m in (("proj", "proj"), ("GiViTI", "giviti"), ("Stukel", "stukel"), ("HL", "hl"), ("mHL", "mhl")):
    y = DF[DF.arm == a]
    emit("nichelead" + m, str(int((y.z > 2).sum())))
    emit("nichedmin" + m, f3(y.d.min()))
    emit("nichedmax" + m, f3(y.d.max()))
# the closest arm overall (smallest mean paired difference)
closest_overall = DF.groupby("arm").d.mean().idxmin()
emit("nicheclosest", {"proj": "the projection test", "GiViTI": "GiViTI", "Stukel": "Stukel's test",
                      "HL": "the Hosmer--Lemeshow test", "mHL": "mHL"}[closest_overall])
# the cells that fail the rule, with their closest arm
fails = rule[~rule.lead_all]
emit("nichenfail", str(len(fails)))
desc = []
for _, r in fails.iterrows():
    nm = {"proj": "the projection test", "GiViTI": "GiViTI", "Stukel": "Stukel's test",
          "HL": "Hosmer--Lemeshow", "mHL": "mHL"}[r.closest]
    desc.append("%s, $C=%.1f$, $n=%d$ (%s, $%s\\pm%s$)" % (
        FAM[r.fam], r.C, r.n, nm, ("+" if r.closest_d >= 0 else "") + f3(r.closest_d), f3(r.closest_se)))
emit("nichefailcells", "; ".join(desc) if desc else "none")
# sizes
sz = SZ.set_index(["arm", "n"])
for a, m in (("DG", "dg"), ("proj", "proj"), ("GiViTI", "giviti"), ("Stukel", "stukel"), ("HL", "hl"),
             ("mHL", "mhl")):
    emit("nichesize" + m, "%s/%s" % (f3(sz.loc[(a, 200), "size"]), f3(sz.loc[(a, 500), "size"])))
emit("nichesizese", f3(np.sqrt(.05 * .95 / 2000)))
# DG power range
dgp = PW[PW.arm == "DG"].matched
emit("nichedgmin", f3(dgp.min()))
emit("nichedgmax", f3(dgp.max()))
# BAGofT
emit("nichebagnull", f3(bag_null))
emit("nichebagnullse", f3(bag_null_se))
emit("nichebagsemax", f3(BG.bag_se.max()))
emit("nichebagmin", f3(BG.bag_matched.min()))
emit("nichebagmax", f3(BG.bag_matched.max()))
emit("nichebagnlead", str(int((BG.d / BG.d_se > 2).sum())))
emit("nichebagnabove", str(int((BG.d / BG.d_se < -2).sum())))
emit("nichebagdmin", f3(BG.d.min()))
emit("nichebagdmax", f3(BG.d.max()))
# highest point estimate of matched power among ALL arms in the cell (BAGofT included at n = 200)
hi = []
for _, cc in cells.iterrows():
    m = PW[(PW.fam == cc.fam) & (PW.C == cc.C) & (PW.n == cc.n)].set_index("arm").matched
    others = list(m.drop("DG"))
    bg = BG[(BG.fam == cc.fam) & (BG.C == cc.C)]
    if cc.n == 200:
        others.append(float(bg.bag_matched.iloc[0]))
    mx = max(others)
    hi.append(dict(fam=cc.fam, C=cc.C, n=cc.n,
                   state="highest" if m["DG"] > mx + 1e-12 else ("tie" if abs(m["DG"] - mx) < 1e-12 else "behind")))
hi = pd.DataFrame(hi)
emit("nichenhighest", str(int((hi.state == "highest").sum())))
emit("nichenhighestword", CARD[int((hi.state == "highest").sum())])
emit("nichentie", str(int((hi.state == "tie").sum())))
emit("nichenbehind", str(int((hi.state == "behind").sum())))
emit("nichenquadhighest", str(int(((hi.fam == "quad") & (hi.state == "highest")).sum())))
emit("nichenprojhigher", str(int((DF[DF.arm == "proj"].d > 0).sum())))
emit("nichentests", CARD[len(ARMS) + 1])          # the six arms of the rule plus BAGofT
print(hi.to_string())
io.open(TEX, "w", encoding="utf-8", newline="\n").write(
    "%% GENERATED by theory/niche_paper/analyze_niche_paper.py -- do not edit by hand.\n" + "\n".join(frag))

# ---- RESULT.md ----
L = ["# Niche study, audit-paper arm set (analyze_niche_paper.py)", "",
     "Decision rule: DG exceeds every other arm (proj, GiViTI, Stukel, HL, mHL) by > 2 paired SE in a cell,",
     "and in no cell is any arm above DG by > 2 paired SE.", "",
     "cells satisfying the lead: **%d of 12**; cells with an arm > 2 SE above DG: **%d**; rule holds overall: %s"
     % (n_lead, n_above, holds), "", "## Matched power (MC SE)", ""]
L.append("| fam | C | n | " + " | ".join(ARMS) + " |")
L.append("|---|---|---|" + "---|" * len(ARMS))
for _, cc in cells.iterrows():
    z = PW[(PW.fam == cc.fam) & (PW.C == cc.C) & (PW.n == cc.n)].set_index("arm")
    L.append("| %s | %.1f | %d | %s |" % (cc.fam, cc.C, cc.n, " | ".join(
        "%.3f (%.3f)" % (z.loc[a, "matched"], z.loc[a, "se"]) for a in ARMS)))
L += ["", "## DG minus arm, matched (paired SE) [z]", "", "| fam | C | n | " + " | ".join(RIVALS) + " |",
      "|---|---|---|" + "---|" * len(RIVALS)]
for _, cc in cells.iterrows():
    z = DF[(DF.fam == cc.fam) & (DF.C == cc.C) & (DF.n == cc.n)].set_index("arm")
    L.append("| %s | %.1f | %d | %s |" % (cc.fam, cc.C, cc.n, " | ".join(
        "%+.3f (%.3f) [%+.1f]" % (z.loc[a, "d"], z.loc[a, "se"], z.loc[a, "z"]) for a in RIVALS)))
L += ["", "## Null size (binomial SE)", "", SZ.round(4).to_markdown(index=False), "",
      "## BAGofT (n = 200, 100 reps per cell; matched on its own 100 nulls)", "",
      BG.round(4).to_markdown(index=False), "", "BAGofT null size %.3f (%.3f), %d nulls" % (bag_null, bag_null_se, len(BN)),
      "", "## Rule per cell", "", rule.round(4).to_markdown(index=False), "",
      "## NA share per arm and scenario (NA -> p = 1)", "", na_share.round(4).to_markdown()]
io.open(os.path.join(HERE, "RESULT.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("wrote", TEX)
