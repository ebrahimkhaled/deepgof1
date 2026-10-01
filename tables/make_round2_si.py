"""make_round2_si.py -- the five supporting studies of the second referee round, written from their result files
(nothing hand-typed). Writes esm/round2_si.tex (Supporting Information Sections 21-25, input after the FDIC
section) and ../round2_macros.tex (the numbers the main text quotes). Re-run after any of the studies changes.

  ablation       theory/ablation/ABLATION_POWER.csv, ABLATION_SIZE.csv
  localization   theory/localization/LOCALIZATION.csv
  external       theory/external/EXTERNAL_LEVEL.csv, EXTERNAL_POWER.csv
  plasmode       medical_example/PLASMODE2_LEVEL.csv, PLASMODE2_POWER.csv; timing screen_SUPPORT/*_M0.csv
  corruption     theory/corruption/CORRUPTION_FALSE_ALARMS_MT.csv (+ rows_rivals/ for projection and BAGofT)
"""
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
TH, ME = os.path.join(ROOT, "theory"), os.path.join(ROOT, "medical_example")
f3 = lambda v: "--" if v is None or (isinstance(v, float) and not np.isfinite(v)) else ("%.3f" % v).replace("0.", ".", 1) if v < 1 else "%.2f" % v
f2 = lambda v: "--" if v is None or (isinstance(v, float) and not np.isfinite(v)) else ("%.2f" % v).replace("0.", ".", 1) if v < 1 else "1.00"
pc = lambda v: "%d\\%%" % round(100 * v)
M = {}                                                     # macros for the main text
out = []

# ---------------------------------------------------------------- 21 ablation
A = pd.read_csv(os.path.join(TH, "ablation", "ABLATION_POWER.csv"))
Z = pd.read_csv(os.path.join(TH, "ablation", "ABLATION_SIZE.csv"))
ST = ["axes", "allpairs", "combined", "ss", "max"]
def verdict(d, se):
    return np.where(d > 2 * se, "ahead", np.where(d < -2 * se, "behind", "level"))
A["vs_ss"] = verdict(A.net_minus_ss, A.se_ss); A["vs_max"] = verdict(A.net_minus_max, A.se_max)
DES = [("liu", "Liu et al.\\ Settings 2--4"), ("niche", "two active among ten"), ("ushape", "U-shape, five covariates")]
rows = []
for k, lab in DES:
    a = A[A.design == k]; z = Z[Z.design == k]
    cnt = lambda col: "/".join(str(int((a[col] == v).sum())) for v in ("ahead", "behind", "level"))
    rows.append("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s--%s \\\\" % (
        lab, len(a), *[f3(a[s].mean()) for s in ST], cnt("vs_ss"), cnt("vs_max"),
        f3(z[ST].min().min()), f3(z[ST].max().max())))
    M["ab" + k.capitalize() + "Net"] = f3(a["axes"].mean()); M["ab" + k.capitalize() + "SS"] = f3(a["ss"].mean())
    M["ab" + k.capitalize() + "Max"] = f3(a["max"].mean()); M["ab" + k.capitalize() + "All"] = f3(a["allpairs"].mean())
    M["ab" + k.capitalize() + "Comb"] = f3(a["combined"].mean())
    for v in ("ahead", "behind", "level"):
        M["ab" + k.capitalize() + "SS" + v.capitalize()] = str(int((a.vs_ss == v).sum()))
        M["ab" + k.capitalize() + "Max" + v.capitalize()] = str(int((a.vs_max == v).sum()))
    M["ab" + k.capitalize() + "Cells"] = str(len(a))
M["abNicheGain"] = f3(A[A.design == "niche"]["axes"].mean() - A[A.design == "niche"].ss.mean())
nz = A[A.design == "niche"].copy()
cell_rows = []
for _, r in pd.concat([A[A.design == "niche"], A[A.design == "ushape"]]).iterrows():
    dep = {"inter": "interaction", "quad": "quadratic"}.get(r.fam, r.fam)
    cell_rows.append("%s & %s & %g & %d & %s & %s & %s & %s & %s & $%+.3f$ (%s) & %s \\\\" % (
        "ten" if r.design == "niche" else "five", dep, r.C, r.n, f3(r["axes"]), f3(r.allpairs), f3(r.combined),
        f3(r.ss), f3(r["max"]), r.net_minus_ss, f3(r.se_ss), f3(r.axes_hit)))
out.append(r"""\section{The network against simple statistics of the same map}\label{sec:ablationS}

Does the trained network add anything to the residual map, or would a plain statistic of the same $36$ cells do as
well? On every dataset one set of $B=199$ parametric-bootstrap refits gives five p-values: the shipped network on the
map of the axis rule (the default test), the largest network score over every pair of covariates (the all-pairs
reading), the smaller of these two p-values calibrated exactly over the $B+1$ draws (the combined reading), the sum of
the squared cells of the axis-rule map (a $\chi^{2}$ statistic over a partition of two covariates) and the largest
absolute cell of the same map. Each statistic is placed at exactly level $.05$ on its own nulls of the same design and
$n$ by the randomized threshold of the main text, and the paired difference in matched rejection, network minus each
simple statistic, is taken over the same datasets. Designs: Liu et al.'s Settings 2--4 at $n=50,100,500$ (nulls $500$,
alternatives $250$ per strength $C\in\{.25,.5,.75,1\}$); the design of Section~7 of the main text, ten standard normal
covariates with an omitted interaction or quadratic on the first two ($C\in\{.3,.5,.8\}$, $n=200,500$; nulls $600$,
alternatives $150$); and a pure U-shape in one of five covariates ($C\in\{.3,.6\}$, $n=200,500$; nulls $1{,}000$,
alternatives $400$).

\begin{table}[h]
\centering
\caption{Mean matched-level power per design (each cell weighted equally), and the number of cells in which the
network (axis rule) is ahead of, behind, or level with the sum of squares and the largest cell by more than two paired
standard errors. Last column: range of raw null size over the five statistics and all sample sizes.}
\label{tab:ablation}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}lrccccccccc@{}}
\toprule
design & cells & network & all pairs & combined & $\Sigma$ squares & max cell & vs $\Sigma$ sq. & vs max & size \\
 & & (axis rule) & & & & & a/b/l & a/b/l & \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[h]
\centering
\caption{Per cell, the designs with many covariates: matched-level power of the five statistics, the paired difference
network minus sum of squares (standard error), and the share of datasets in which the axis rule selected a pair
containing the misfit covariate(s).}
\label{tab:ablationcells}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}llrrcccccrc@{}}
\toprule
covariates & departure & $C$ & $n$ & network & all pairs & combined & $\Sigma$ sq. & max & net $-$ $\Sigma$ sq. & axes hit \\
\midrule
""" + "\n".join(cell_rows) + r"""
\bottomrule
\end{tabular}
\end{table}

On the benchmark of Liu et al.\ the network and the $\chi^{2}$ of the same map are level in most cells: the learned
score adds nothing there, and the power comes from the map and its calibration. With two active covariates among ten,
the setting of Section~7 of the main text, the network is ahead of the sum of squares in half of the cells and behind
in none, because the misfit is concentrated in a few cells and a sum of squares spreads its power over all $36$; the
largest-cell statistic, which reads only one cell, is weaker still. In the U-shape design the axis rule misses the
misfit covariate and every statistic of its map stays near the level; the all-pairs reading recovers it.
""")

# ---------------------------------------------------------------- 22 localization
L = pd.read_csv(os.path.join(TH, "localization", "LOCALIZATION.csv"), keep_default_na=False)
LAB = {"null": "none", "ushape": "U-shape in $x_1$", "thresh": "threshold in $x_2$", "inter": "interaction $x_1x_3$"}
TOOL = {"DG axis rule": "DeepGOF-1, axis rule", "DG all-pairs": "DeepGOF-1, all pairs", "Lin-Wei-Ying": "Lin--Wei--Ying",
        "spline screen": "spline screen"}
lr = []
for (fam, C), g in L.groupby(["family", "C"], sort=False):
    for _, r in g.iterrows():
        lr.append("%s & %s & %s & %d & %s & %s & %s & %s \\\\" % (LAB[fam], "--" if fam == "null" else "%g" % C, TOOL[r.tool],
                                                              r.n, f3(r.rej), "--" if fam == "null" else pc(r.right),
                                                              "--" if fam == "null" else pc(r.located), f3(r["false"])))
    lr.append("\\addlinespace")
g = lambda fam, C, tool: L[(L.family == fam) & (np.isclose(L.C, C)) & (L.tool == tool)].iloc[0]
r1 = g("inter", 1.0, "DG all-pairs"); M["locInterRej"] = f2(r1.rej); M["locInterRight"] = pc(r1.right)
M["locInterLocated"] = pc(r1.located)
M["locInterOtherRej"] = "%s and %s" % (f2(g("inter", 1.0, "Lin-Wei-Ying").rej), f2(g("inter", 1.0, "spline screen").rej))
M["locInterOtherRightMax"] = pc(max(g("inter", 1.0, "Lin-Wei-Ying").right, g("inter", 1.0, "spline screen").right,
                                  g("inter", .6, "Lin-Wei-Ying").right, g("inter", .6, "spline screen").right))
M["locNullFalseAll"] = f3(g("null", 0.0, "DG all-pairs")["false"]); M["locNullFalseAxes"] = f3(g("null", 0.0, "DG axis rule")["false"])
M["locUSplineLocated"] = pc(g("ushape", .4, "spline screen").located); M["locUAllLocated"] = pc(g("ushape", .4, "DG all-pairs").located)
M["locThreshSpline"] = pc(g("thresh", 4.0, "spline screen").located); M["locThreshAll"] = pc(g("thresh", 4.0, "DG all-pairs").located)
out.append(r"""\section{Locating the misfit}\label{sec:locS}

Five standard normal covariates, $n=500$, $\eta_{0}=-0.5+0.5x_{1}+0.8x_{2}+0.6x_{3}+0.4x_{4}+0.3x_{5}$, and a model fitted
linear in all five. Truths: $\eta_{0}+C(x_{1}^{2}-1)$ (U-shape), $\eta_{0}+C\max(x_{2}-1,0)$ (threshold),
$\eta_{0}+Cx_{1}x_{3}$ (interaction), and $\eta_{0}$ ($500$ null datasets; $300$ per alternative cell). One set of
$B=199$ refits per dataset serves every tool. DeepGOF-1: its p-value, the pair it displays, and the cells beyond the
simultaneous $95\%$ threshold of that map (the $95\%$ quantile of its largest absolute cell over the bootstrap
replicates), read only when the test rejects. The true region is both extreme sixths of $x_{1}$ (U-shape), the top two
sixths of $x_{2}$ (threshold) and the four corner blocks of the $(x_{1},x_{3})$ map (interaction). Lin--Wei--Ying: for
each covariate the supremum of the cumulative residual over its ranks, calibrated by the same refits, and the spline
screen: for each covariate the likelihood-ratio test of a natural spline with three degrees of freedom in place of the
linear term; both flag covariates at the Bonferroni level $.05/5$ and reject when they flag any.

\begin{table}[h]
\centering
\caption{Locating the misfit ($n=500$). Rejects: rejection rate at $.05$. Right: the displayed pair contains the true
covariate(s) (DeepGOF-1), or the flagged set is exactly the true covariate(s) (the other two). Located: the tool rejects
and every flag lies in the truth (for DeepGOF-1, at least one flagged cell and every flagged cell inside the true
region). False: the tool flags something outside the truth; under the null, any flag. Monte Carlo standard error at most
$.029$.}
\label{tab:loc}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}llllcccc@{}}
\toprule
truth & $C$ & tool & datasets & rejects & right & located & false \\
\midrule
""" + "\n".join(lr[:-1]) + r"""
\bottomrule
\end{tabular}
\end{table}

A missed interaction of two covariates is invisible to tools that look at one covariate at a time, and the all-pairs
map both detects it and names the pair. Curvature or a threshold in one covariate is a one-covariate problem, and the
spline screen locates it best; the axis rule, which ranks covariates by their linear effect, locates neither the
U-shape nor the interaction.
""")

# ---------------------------------------------------------------- 23 external validation
EL = pd.read_csv(os.path.join(TH, "external", "EXTERNAL_LEVEL.csv"), index_col=0)
EP = pd.read_csv(os.path.join(TH, "external", "EXTERNAL_POWER.csv"))
KEEP = [("DG_axes", "DG-1 axes"), ("DG_allpairs", "DG-1 all"), ("DG_combined", "DG-1 comb."), ("proj", "proj."),
        ("pearson", "Pearson"), ("HL", "HL"), ("stukel", "Stukel"), ("giviti", "GiViTI")]
ELr = ["%d & %s \\\\" % (int(n), " & ".join(f3(EL.loc[t, n]) for t, _ in KEEP)) for n in EL.columns]
FAM = {"large": "overall rate", "slope": "slope", "ushape": "U-shape", "thresh": "threshold", "inter": "interaction"}
EPr = []
for (fam, C), g2 in EP.groupby(["fam", "C"], sort=False):
    v = g2.set_index("test").matched
    best = max(v[t] for t, _ in KEEP)
    EPr.append("%s & %g & %s \\\\" % (FAM[fam], C, " & ".join(("\\textbf{%s}" if v[t] == best else "%s") % f2(v[t]) for t, _ in KEEP)))
gp = lambda fam, C, t: EP[(EP.fam == fam) & np.isclose(EP.C, C) & (EP.test == t)].matched.iloc[0]
others = [t for t, _ in KEEP if not t.startswith("DG")]
for fam, C, key in [("thresh", 4, "Thresh"), ("ushape", .8, "U"), ("ushape", .4, "Uweak"), ("inter", 1, "Inter")]:
    dg = max(gp(fam, C, t) for t in ("DG_axes", "DG_allpairs", "DG_combined"))
    M["ext" + key + "DGbest"] = f2(dg); M["ext" + key + "DG"] = f2(gp(fam, C, "DG_combined"))
    M["ext" + key + "Giv"] = f2(gp(fam, C, "giviti"))
    M["ext" + key + "Best"] = f2(max(gp(fam, C, t) for t in others))
M["extSlopeGiv"] = f2(gp("slope", .6, "giviti")); M["extSlopeDG"] = f2(gp("slope", .6, "DG_combined"))
M["extLargeProj"] = f2(gp("large", .4, "proj")); M["extLargeGiv"] = f2(gp("large", .4, "giviti"))
M["extLargeDG"] = f2(gp("large", .4, "DG_combined"))
M["extInterWeakDG"] = f2(gp("inter", .6, "DG_combined")); M["extInterWeakBest"] = f2(max(gp("inter", .6, t) for t in others))
dgl = EL.loc[["DG_axes", "DG_allpairs", "DG_combined"]]
M["extLevelDGlo"] = f3(dgl.min().min()); M["extLevelDGhi"] = f3(dgl.max().max())
oth = EL.loc[[t for t in others]]
M["extLevelOthHi"] = f3(oth.max().max())
out.append(r"""\section{External validation of a published risk model}\label{sec:extS}

A published model with fixed coefficients, $\eta_{0}=-0.5+0.5x_{1}+0.8x_{2}+0.6x_{3}+0.4x_{4}+0.3x_{5}$, is checked on new
patients with five independent standard normal covariates. The null is $y_{i}\sim\mathrm{Bernoulli}(p_{i})$ with
$p_{i}=\mathrm{logit}^{-1}(\eta_{0i})$ given; nothing is refitted. Truths in the validation population: $\eta_{0}+C$
(overall rate, $C\in\{.2,.4\}$), $C\eta_{0}$ (an overfitted model with calibration slope $C\in\{.8,.6\}$),
$\eta_{0}+C(x_{1}^{2}-1)$ (U-shape, $C\in\{.4,.8\}$), $\eta_{0}+C\max(x_{2}-1,0)$ (threshold, $C\in\{2,4\}$) and
$\eta_{0}+Cx_{1}x_{3}$ (interaction, $C\in\{.6,1\}$). Null datasets: $500$ at $n=250$ and $500$, $300$ at $n=1{,}000$;
alternatives: $300$ per cell at $n=500$.

With $p$ fixed, a statistic calibrated by drawing $y^{*}\sim\mathrm{Bernoulli}(p)$ and recomputing it is an exact Monte
Carlo test at any $n$ \citep{besagclifford1989}. DeepGOF-1 (its map over the covariates, standardized by the given $p$;
the axis rule uses the published coefficients), the projection statistic and Pearson's statistic (the statistic of the
Osius--Rojek test) are calibrated in this way from $B=199$ shared draws. The others use their usual external
references: Hosmer--Lemeshow with ten groups against $\chi^{2}_{10}$; Stukel's two directions added to the offset
$\eta_{0}$, likelihood ratio against $\chi^{2}_{2}$; and the GiViTI calibration belt in its external form, which reduces
to the test of calibration intercept and slope when the selected calibration curve is linear.

\begin{table}[h]
\centering
\caption{External validation: size at nominal $.05$ (standard error about $.010$ at $500$ datasets and $.013$ at
$300$).}
\label{tab:extlevel}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}r""" + "c" * len(KEEP) + r"""@{}}
\toprule
$n$ & """ + " & ".join(l for _, l in KEEP) + r""" \\
\midrule
""" + "\n".join(ELr) + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[h]
\centering
\caption{External validation, $n=500$: matched-level power (each test placed at exactly $.05$ on its own $500$ nulls;
$300$ datasets per cell, standard error at most $.029$). Bold: the highest in the row.}
\label{tab:extpower}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}lr""" + "c" * len(KEEP) + r"""@{}}
\toprule
miscalibration & $C$ & """ + " & ".join(l for _, l in KEEP) + r""" \\
\midrule
""" + "\n".join(EPr) + r"""
\bottomrule
\end{tabular}
\end{table}

Miscalibration that is the same for every patient, a wrong overall rate or an overfitted slope, is found best by the
tests built along the predicted risk. Miscalibration that differs between patients with the same predicted risk is
invisible to them, and DeepGOF-1, which reads the residuals over the covariates, has the most power on every such cell
but the weaker interaction, where every test is weak.
""")

# ---------------------------------------------------------------- 24 SUPPORT plasmode and timing
PL = pd.read_csv(os.path.join(ME, "PLASMODE2_LEVEL.csv")); PP = pd.read_csv(os.path.join(ME, "PLASMODE2_POWER.csv"))
PT = [("DG_axes", "DG-1 axes"), ("DG_all", "DG-1 all"), ("DG_comb", "DG-1 comb."), ("proj", "proj."), ("StukelJ", "Stukel J"),
      ("StukelS", "SstBoth"), ("HL", "HL"), ("GiViTI", "GiViTI")]
lv = PL.pivot(index="n", columns="test", values="size"); reps = PL.groupby("n").reps.max()
PLr = ["%s & %d & %s \\\\" % ("{:,}".format(n).replace(",", "{,}"), reps[n], " & ".join(f3(lv.loc[n, t]) if t in lv and np.isfinite(lv.loc[n, t]) else "--" for t, _ in PT)) for n in lv.index]
pw = PP.pivot(index="n", columns="test", values="matched"); praw = PP.pivot(index="n", columns="test", values="raw")
PPr = ["%s & %s \\\\" % ("{:,}".format(n).replace(",", "{,}"), " & ".join("%s (%s)" % (f2(pw.loc[n, t]), f2(praw.loc[n, t])) for t, _ in PT)) for n in pw.index]
full = lv.loc[lv.index.max()]
M["plFullDGlo"] = f3(min(full["DG_axes"], full["DG_all"], full["DG_comb"])); M["plFullDGhi"] = f3(max(full["DG_axes"], full["DG_all"], full["DG_comb"]))
M["plFullReps"] = str(int(reps[lv.index.max()]))
sub = lv.loc[[n for n in lv.index if n < lv.index.max()]]
M["plDGlo"] = f3(sub[["DG_axes", "DG_all", "DG_comb"]].min().min()); M["plDGhi"] = f3(sub[["DG_axes", "DG_all", "DG_comb"]].max().max())
M["plProjPowK"] = f2(pw.loc[1000, "proj"]); M["plDGPowK"] = f2(pw.loc[1000, ["DG_axes", "DG_all", "DG_comb"]].max())
M["plClassPowK"] = f2(pw.loc[1000, ["StukelJ", "StukelS", "HL", "GiViTI"]].max())
TD = pd.read_csv(os.path.join(ME, "screen_SUPPORT", "deepgof_timing_M0.csv")); TP = pd.read_csv(os.path.join(ME, "screen_SUPPORT", "projection_M0.csv"))
ps = float(TP.seconds.iloc[0]); da = float(TD[TD.reading == "axes"].seconds.iloc[0]); dl = float(TD[TD.reading == "allpairs"].seconds.iloc[0])
dc = float(TD[TD.reading == "combined"].seconds.iloc[0])
M["tmProjHours"] = "%d~h~%d~min" % (ps // 3600, round((ps % 3600) / 60)); M["tmProjSecs"] = "{:,}".format(round(ps)).replace(",", "{,}")
M["tmProjGB"] = "%.1f" % (float(TP.peak_MB.iloc[0]) / 1000); M["tmDGaxes"] = "%d" % round(da); M["tmDGall"] = "%d" % round(dl)
M["tmRatio"] = "{:,}".format(int(round(ps / da, -2))).replace(",", "{,}"); M["tmRatioAll"] = "%d" % int(round(ps / dl, -1))
M["tmProjP"] = f3(float(TP.p.iloc[0]))
# split-sample external validation on SUPPORT (support_external_split.R)
XS = pd.read_csv(os.path.join(ME, "screen_SUPPORT", "external_split.csv"))
nspl = XS.split.nunique(); M["xsSplits"] = str(nspl); M["xsNval"] = "{:,}".format(int(XS.n_val.iloc[0])).replace(",", "{,}")
XT = [("DG_combined", "DG-1 comb."), ("DG_axes", "DG-1 axes"), ("DG_allpairs", "DG-1 all"), ("GiViTI", "GiViTI"),
      ("HL", "HL"), ("Stukel", "Stukel")]
xsr = []
for m in ("M0", "M2"):
    z = XS[XS.model == m]
    rej = {t: int((z[t] <= .05).sum()) for t, _ in XT}
    xsr.append("%s & %s & %s & %s & %s \\\\" % (m, " & ".join("%d" % rej[t] for t, _ in XT),
                                               " & ".join(f3(z[t].median()) for t in ("DG_combined", "GiViTI")),
                                               f2(z.citl.median()), f2(z.slope.median())))
    for t, nm in [("DG_combined", "DG"), ("GiViTI", "Giv"), ("HL", "HL"), ("Stukel", "Stuk")]:
        M["xs" + {"M0": "Mzero", "M2": "Mtwo"}[m] + nm] = str(rej[t])
    M["xs" + {"M0": "Mzero", "M2": "Mtwo"}[m] + "DGmed"] = f3(z.DG_combined.median()); M["xs" + {"M0": "Mzero", "M2": "Mtwo"}[m] + "Givmed"] = f2(z.GiViTI.median())
    M["xs" + {"M0": "Mzero", "M2": "Mtwo"}[m] + "Slope"] = f2(z.slope.median()); M["xs" + {"M0": "Mzero", "M2": "Mtwo"}[m] + "Citl"] = f2(z.citl.median())
M["xsMzeroResp"] = pc(XS[XS.model == "M0"].top_pair.str.contains("resp").mean())
out.append(r"""\section{The SUPPORT design: level, power, computing time and a split-sample validation}\label{sec:plasmodeS}

Each replicate draws $n$ of the $8{,}873$ complete-case patients without replacement, keeps their real covariates, and
simulates in-hospital death from a true model fitted on all of them: the null is M0, linear in every covariate, and the
alternative is M1, with a natural spline in mean arterial pressure (a U-shaped risk). The analyst fits M0. DeepGOF-1's
three readings share one set of $B=199$ refits; the projection test ($B=199$), Stukel's joint and SstBoth tests, the
Hosmer--Lemeshow test and GiViTI were run up to $n=1{,}000$, and at the full $n=8{,}873$ DeepGOF-1 alone (every patient,
outcomes redrawn).

\begin{table}[h]
\centering
\caption{SUPPORT design: size at nominal $.05$ (standard error about $.013$ at $300$ datasets, $.015$ at $200$ and
$.022$ at $100$).}
\label{tab:plasmodelevel}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}rr""" + "c" * len(PT) + r"""@{}}
\toprule
$n$ & datasets & """ + " & ".join(l for _, l in PT) + r""" \\
\midrule
""" + "\n".join(PLr) + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[h]
\centering
\caption{SUPPORT design, U-shaped risk of blood pressure: matched-level power (raw power in brackets); $200$ datasets
at $n=250$ and $500$, $150$ at $1{,}000$.}
\label{tab:plasmodepower}
\small\setlength{\tabcolsep}{2pt}
\begin{tabular}{@{}r""" + "c" * len(PT) + r"""@{}}
\toprule
$n$ & """ + " & ".join(l for _, l in PT) + r""" \\
\midrule
""" + "\n".join(PPr) + r"""
\bottomrule
\end{tabular}
\end{table}

\emph{Computing time.} On M0 with all $8{,}873$ patients and $B=199$, on one core of the same machine under the same
load, DeepGOF-1 took """ + M["tmDGaxes"] + r"""~s (all-pairs and combined readings """ + M["tmDGall"] + r"""~s), and
\texttt{projection.gof()} took """ + M["tmProjSecs"] + r"""~s (""" + M["tmProjHours"] + r"""), with a peak memory of
""" + M["tmProjGB"] + r"""~GB; both returned $p=""" + M["tmProjP"] + r"""$.

\emph{Split-sample external validation.} The public release has no hospital or date identifier, so the split is
random: half of the $8{,}873$ patients develop a model, which is frozen and checked on the other half ($""" + M["xsNval"] + r"""$
patients), repeated over """ + M["xsSplits"] + r""" random splits. The frozen models are M0 and M2 of the main text.
DeepGOF-1 is \texttt{deepgof1.external()} ($B=199$) on the ten covariates (sex coded $0/1$); GiViTI is the calibration
test in its external form; Hosmer--Lemeshow uses ten groups against $\chi^{2}_{10}$; Stukel's two terms are added to
the frozen linear predictor.

\begin{table}[h]
\centering
\caption{SUPPORT, split-sample external validation: number of the """ + M["xsSplits"] + r""" splits in which each test
rejects at $.05$; median p-values of DeepGOF-1 (combined reading) and GiViTI; median calibration intercept (slope fixed
at one) and calibration slope on the validation half.}
\label{tab:extsplit}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}l""" + "c" * len(XT) + r"""cccc@{}}
\toprule
model & """ + " & ".join(l for _, l in XT) + r""" & med.\ p DG-1 & med.\ p GiViTI & intercept & slope \\
\midrule
""" + "\n".join(xsr) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# ---------------------------------------------------------------- 25 corrupted records at n = 1,000
CF = pd.read_csv(os.path.join(TH, "corruption", "CORRUPTION_FALSE_ALARMS_MT.csv"))
RV = glob.glob(os.path.join(TH, "corruption", "rows_rivals", "*.csv"))
rv = pd.concat([pd.read_csv(f) for f in RV]) if RV else pd.DataFrame(columns=["test", "file", "p", "match"])
CL = [("logit_clean_n1000", "1{,}000", "none"), ("logit_C1_r001_n1000", "1{,}000", "$\\times4$, 1 record"),
      ("logit_C1_r002_n1000", "1{,}000", "$\\times4$, 2 records"), ("logit_C1_r005_n1000", "1{,}000", "$\\times4$, 5 records"),
      ("logit_C1_r010_n1000", "1{,}000", "$\\times4$, 10 records"),
      ("clean", "500", "none"), ("x4", "500", "$\\times4$, 3 records"), ("x8", "500", "$\\times8$, 3 records")]
cr = []
for key, n, lab in CL:
    r = CF[CF.cell == key].iloc[0]
    def rival(t, col):
        if n == "500":
            return f3(r[col]), int(r.reps)
        z = rv[(rv.test == t) & (rv.file == key)]
        return (f3((z.p.dropna() <= .05).mean()), len(z)) if len(z) else ("--", 0)
    pj, npj = rival("proj", "proj"); bg, nbg = rival("BAGofT", "BAGofT")
    wn = lambda v, k: "--" if k == 0 else "%s (%d)" % (v, k)
    cr.append("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\" % (
        n, lab, "{:,}".format(int(r.reps)).replace(",", "{,}"), f3(r.DeepGOF1), f3(r["HL.G10"]), f3(r["Stk.joint"]),
        f3(r.GiViTI), wn(pj, npj), wn(bg, nbg)))
    if key == "logit_C1_r010_n1000":
        cr.append("\\addlinespace")
g10 = CF[CF.cell == "logit_C1_r010_n1000"].iloc[0]; g5 = CF[CF.cell == "logit_C1_r005_n1000"].iloc[0]
M["coDGten"] = f3(g10.DeepGOF1); M["coStukten"] = f3(g10["Stk.joint"]); M["coGivten"] = f3(g10.GiViTI); M["coHLten"] = f3(g10["HL.G10"])
M["coDGfive"] = f3(g5.DeepGOF1); M["coStukfive"] = f3(g5["Stk.joint"]); M["coGivfive"] = f3(g5.GiViTI)
M["coDGclean"] = f3(CF[CF.cell == "logit_clean_n1000"].iloc[0].DeepGOF1)
x8 = CF[CF.cell == "x8"].iloc[0]; M["coDGxeight"] = f2(x8.DeepGOF1); M["coProjxeight"] = f2(x8.proj); M["coStukxeight"] = f2(x8["Stk.joint"])
for key, t, nm in [("logit_C1_r010_n1000", "proj", "coProjten"), ("logit_C1_r005_n1000", "proj", "coProjfive"),
                   ("logit_C1_r010_n1000", "BAGofT", "coBagten"), ("logit_clean_n1000", "proj", "coProjclean"),
                   ("logit_clean_n1000", "BAGofT", "coBagclean")]:
    z = rv[(rv.test == t) & (rv.file == key)]
    M[nm] = f3((z.p.dropna() <= .05).mean()) if len(z) else "--"
pr = rv[rv.test == "proj"].groupby("file").p.apply(lambda s: (s <= .05).mean())
M["coProjlo"] = f3(pr.min()); M["coProjhi"] = f3(pr.max())
M["coProjN"] = str(int(rv[rv.test == "proj"].groupby("file").size().min()))
M["coBagN"] = str(int(rv[rv.test == "BAGofT"].groupby("file").size().min()))
# the paired comparison with DeepGOF-1 on the same datasets: the largest discordant count either way
dgr = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(TH, "corruption", "rows_edge_mt", "9_*.csv"))])
disc = []
for (t, fl), z in rv.groupby(["test", "file"]):
    m = z.merge(dgr[dgr.file == fl][["rep", "DG"]], on="rep")
    disc += [int(((m.p <= .05) & (m.DG > .05)).sum()), int(((m.p > .05) & (m.DG <= .05)).sum())]
M["coDiscMax"] = str(max(disc))
out.append(r"""\section{A few corrupted records at $n=500$ and $1{,}000$}\label{sec:corrupt1000S}

The design of the corrupted-record studies of the directed calibration tests, so that the stored p-values of the other
tests are on the same datasets: $x\sim U(-3,3)$, $d\sim\mathrm{Bernoulli}(.5)$,
$y\sim\mathrm{Bernoulli}\{\mathrm{logit}^{-1}(.6x+.5d)\}$, fitted as $y\sim x+d$, the correct model. After the outcome is
drawn, $k$ records chosen at random have $x$ multiplied by $4$ or $8$; the outcome is right and the covariate is wrong,
so every rejection is a false alarm caused by the corruption. Each dataset was rebuilt from its stored seed and checked
against the stored number of events and fitted slope (no mismatch), and DeepGOF-1 ($B=199$; with one pair of covariates
its readings coincide) was run on it. The projection test ($B=199$) and BAGofT (package defaults) at $n=1{,}000$ were run
on the first $200$ and $100$ datasets of their cells (number in brackets).

\begin{table}[h]
\centering
\caption{False-alarm rate at $.05$ with corrupted covariate values. Standard error about $.007$ at $1{,}000$ datasets,
$.015$ at $200$ and $.022$ at $100$. HL: Hosmer--Lemeshow, ten groups; Stukel: joint score test; GiViTI: internal form.}
\label{tab:corrupt1000}
\small\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}llrcccccc@{}}
\toprule
$n$ & corruption & datasets & DeepGOF-1 & HL & Stukel & GiViTI & projection & BAGofT \\
\midrule
""" + "\n".join(cr) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

# the main-text version of the corrupted-record table: n = 1,000 only, the tests of the paper
main = [r for r in cr if r.startswith("1{,}000")]
M["tablecorruptmain"] = "\n".join(main)
with open(os.path.join(HERE, "..", "corrupt1000_main.tex"), "w", encoding="utf-8") as fh:
    fh.write("% generated by esm/make_round2_si.py -- do not edit by hand\n" + "\n".join(main) + "\n")

with open(os.path.join(HERE, "round2_si.tex"), "w", encoding="utf-8") as fh:
    fh.write("% generated by make_round2_si.py -- do not edit by hand\n" + "\n".join(out))
with open(os.path.join(HERE, "..", "round2_macros.tex"), "w", encoding="utf-8") as fh:
    fh.write("% generated by esm/make_round2_si.py -- do not edit by hand\n")
    for k, v in M.items():
        fh.write("\\newcommand{\\%s}{%s}\n" % (k, v))
print("wrote round2_si.tex and round2_macros.tex;", len(M), "macros")
for k, v in M.items():
    print("  %-22s %s" % (k, v))
