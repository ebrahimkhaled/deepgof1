"""make_si_tables.py -- the Supporting Information section on very small samples and wrong links, written from
LEVEL.csv, POWER.csv and POWER_COVARIATE_MISFIT.csv (never hand-typed). Output: smalln_link_si.tex, copied to
manuscript/sim/esm/ (the SI \\input's it)."""
import pandas as pd, pathlib, shutil
H = pathlib.Path(__file__).resolve().parent
L = pd.read_csv(H / "LEVEL.csv"); P = pd.read_csv(H / "POWER.csv"); C = pd.read_csv(H / "POWER_COVARIATE_MISFIT.csv")
T = [("DG", "DG-1"), ("proj", "proj."), ("StukelJ", "Stukel J"), ("StukelS", "SstBoth"), ("GiViTI", "GiViTI"), ("HL", "HL")]
f3 = lambda v: ("%.3f" % float(v)).lstrip("0") if float(v) >= 0 else "-" + ("%.3f" % -float(v)).lstrip("0")
num = lambda s: float(str(s).split()[0]); na = lambda s: float(str(s).split("[")[1].rstrip("]"))
LINK = {"probit": "probit", "cloglog": "cloglog", "stukel_heavy": "Stukel, heavier tails", "stukel_asym": "Stukel, asymmetric"}

o = [r"""\section{Very small samples and wrong links}\label{sec:smallnlink}

The design has two standard normal covariates and linear predictor $\eta=0.25+1.5x_{1}+x_{2}$, about half
events; the fitted model is the logistic regression in $x_{1}$ and $x_{2}$, fitted by a plain \texttt{glm()}.
Null datasets: $2{,}000$ at $n=20$ and $30$ and $1{,}000$ at $n=50$ and $100$. Wrong links: the outcome is
drawn with probability $F(\eta)$ for $F$ the probit, the complementary log-log, and Stukel's generalized
logistic with heavier tails ($\alpha_{1}=\alpha_{2}=-1$) and with asymmetric tails ($\alpha_{1}=0.5$,
$\alpha_{2}=-1$), at $n=20$, $30$, $50$ and $100$. Misfit in the covariates: $\eta+C\{x_{1}^{2}-1\}$ (omitted
quadratic) or $\eta+Cx_{1}x_{2}$ (omitted interaction), $C\in\{0.5,1\}$, at $n=20$, $30$ and $50$. Every
alternative cell has $1{,}000$ datasets. The tests are DeepGOF-1 (shipped network, $B=199$; with two
covariates both readings coincide), the projection test ($B=199$), Stukel's joint score test and its SstBoth
form, GiViTI through \texttt{givitiR}, and the Hosmer--Lemeshow test. Power is at matched level $.05$: each
test is placed at exactly $.05$ on its own null datasets at the same $n$ by the randomized threshold of the
main text, and is evaluated on the datasets where it returns a p-value.

\begin{table}[h]
\centering
\caption{Size at nominal $.05$ on the null datasets; in brackets, the share of datasets with no p-value
where it is not zero. Separated: share of fits that did not converge or had a coefficient above $8$ in
absolute value. Band: $[.035,.065]$.}
\label{tab:smalln_size}
\small
\begin{tabular}{@{}rrcccccc@{}}
\toprule
$n$ & separated & DG-1 & proj. & Stukel J & SstBoth & GiViTI & HL \\
\midrule"""]
for r in L.itertuples():
    cells = []
    for k, _ in T:
        s = getattr(r, k)
        cells.append(f3(num(s)) + (" [%.2f]" % na(s) if na(s) > 0 else ""))
    o.append("%d & %s & %s \\\\" % (r.n, f3(r.separated), " & ".join(cells)))
o.append(r"""\bottomrule
\end{tabular}
\end{table}

\begin{table}[h]
\centering
\caption{Matched-level power against misfit in the covariates. Monte Carlo standard error at most $.016$.}
\label{tab:smalln_power}
\small
\begin{tabular}{@{}lrrcccccc@{}}
\toprule
misfit & $C$ & $n$ & DG-1 & proj. & Stukel J & SstBoth & GiViTI & HL \\
\midrule""")
for r in C.sort_values(["fam", "C", "n"], ascending=[False, True, True]).itertuples():
    o.append("%s & %.1f & %d & %s \\\\" % ("quadratic" if r.fam == "quad" else "interaction", r.C, r.n,
                                          " & ".join(f3(getattr(r, k)) for k, _ in T)))
o.append(r"""\bottomrule
\end{tabular}
\end{table}

\begin{table}[h]
\centering
\caption{Matched-level power when the true link is not the logit. Monte Carlo standard error at most $.015$.}
\label{tab:smalln_link}
\small
\begin{tabular}{@{}lrcccccc@{}}
\toprule
true link & $n$ & DG-1 & proj. & Stukel J & SstBoth & GiViTI & HL \\
\midrule""")
for lk in LINK:
    for r in P[P.link == lk].sort_values("n").itertuples():
        o.append("%s & %d & %s \\\\" % (LINK[lk], r.n, " & ".join(f3(getattr(r, k)) for k, _ in T)))
o.append(r"""\bottomrule
\end{tabular}
\end{table}

DeepGOF-1 is the only test whose size lies inside the band at every $n$
(Table~\ref{tab:smalln_size}). Against misfit in the covariates it has more power than the Hosmer--Lemeshow
test in ten of the twelve cells and less than the projection test in all twelve
(Table~\ref{tab:smalln_power}); with only two covariates, both carrying signal, this is the order the local
power comparison of the main text predicts. Against a wrong link DeepGOF-1 and the projection test stay close
to their level (Table~\ref{tab:smalln_link}): a wrong link bends the risk along the linear predictor and
leaves little pattern over the covariate ranks, while Stukel's test and GiViTI look along the linear
predictor and are the tests to use; even they reach appreciable power only for the asymmetric link at
$n\ge50$. The scripts and every per-dataset p-value are in \texttt{theory/smalln\_link/} of the deposit.
""")
out = H / "smalln_link_si.tex"
out.write_text("\n".join(o), encoding="utf-8", newline="\n")
dst = H.parents[1] / "manuscript" / "sim" / "esm" / "smalln_link_si.tex"
shutil.copy2(out, dst)
print("wrote", out, "and", dst)
