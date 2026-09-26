"""05_certify.py -- assemble the smoothing certificate from the net (02), the tail (03) and the truth (04).

FIXED test.  At net point t_j the Clopper-Pearson lower bound q_lo,j at level 1 - ETA_NET/J gives
Lambda_j = Phi^{-1}(q_lo,j).  By the Gaussian-shift lemma, Phi^{-1} q(t) >= Lambda_j - |t - t_j| for every t,
so q(t) >= beta on the union of the intervals [t_j - r_j, t_j + r_j], r_j = Lambda_j - Phi^{-1}(beta) >= 0.
With the tail intervals of 03 (q >= 1 - EPS_TAIL) and Lemma 2 beyond t_global, the certified t is the left
end of the connected component of that union which contains +infinity.  Exact, no grid.

B199 test.  Power = E[Bin(9; 199, W)] >= sum_l dg_l q_{c_l}(t) (layer cake, THEOREM_SMOOTH.md Lemma 4), with
each q_{c_l}(t) >= Phi(max_j [Lambda_{j,l} - |t - t_j|]) and CP level 1 - ETA_NET/(J m).  The bound is
evaluated on a grid of spacing DT and, being (0.4)-Lipschitz in t, is reduced by 0.2 DT between grid points.

Run: python -u 05_certify.py            (writes SMOOTH_CERT.json and prints the table)
"""
import os, sys, json
import numpy as np
from scipy.stats import norm
from smooth_core import (ORDER, T_BAR, LEVELS, LAYER_W, ETA_NET, ETA_CRIT, EPS_TAIL, HERE, RHO,
                         T_EMP_PAPER, cp_lower, read_jsonl)

# the paper's usage-aware CROWN certificate (Online Resource 1, Table 'Certified against true sample sizes')
N_CROWN = {"step_x1": 3789, "diffuse": 4311, "stripe_x1": 5881, "saddle": 8044, "random": 13918}
T_CROWN = {"step_x1": 30.78, "diffuse": 32.83, "stripe_x1": 38.34, "saddle": 44.84, "random": 58.99}
N_TRUE_PAPER = {"step_x1": 72, "diffuse": 103, "stripe_x1": 122, "saddle": 204, "random": 729}
DT = 0.001

def n_of_t(t):
    return (t / RHO) ** 2

def component_left_end(intervals, tail_from, tail_closed):
    """Left end of the connected component containing +inf of  U intervals  U  [tail_from, inf)."""
    if not tail_closed:
        return None
    ivs = sorted([iv for iv in intervals if iv[1] >= iv[0]], key=lambda x: x[0])
    left = tail_from
    changed = True
    while changed:
        changed = False
        for a, b in ivs:
            if a < left and b >= left - 1e-12:
                left = a; changed = True
    return left

def certify_direction(d, eta_net):
    net = sorted(read_jsonl(os.path.join(HERE, "net_%s.jsonl" % d)), key=lambda r: r["t"])
    tail = read_jsonl(os.path.join(HERE, "tail_%s.jsonl" % d))
    tsum = json.load(open(os.path.join(HERE, "tail_%s_summary.json" % d)))
    tail_closed = bool(tsum["closed"])
    t_start = tsum["t_start"]
    J = len(net); m = len(LEVELS)
    ts = np.array([r["t"] for r in net]); Ms = np.array([r["M"] for r in net])
    out = {}
    # ---- FIXED ----
    a = eta_net / J
    lam = np.array([norm.ppf(cp_lower(r["n_fixed"], r["M"], a)) for r in net])
    for beta_ in (.95, .50):
        zb = norm.ppf(beta_)
        ivs = [(t - (l - zb), t + (l - zb)) for t, l in zip(ts, lam) if l >= zb]
        assert norm.ppf(1 - EPS_TAIL) >= zb
        tc = component_left_end(ivs, t_start, tail_closed)
        out[("fixed", beta_)] = tc
    # ---- B199 ----
    a2 = eta_net / (J * m)
    lam2 = np.array([[norm.ppf(cp_lower(r["n_levels"][l], r["M"], a2)) for l in range(m)] for r in net])   # J x m
    grid = np.arange(0.0, t_start + DT / 2, DT)
    # L_l(t) = max_j (lam2[j,l] - |t - t_j|), computed in blocks
    P = np.zeros(len(grid))
    for l in range(m):
        Ll = np.full(len(grid), -np.inf)
        for j in range(J):
            Ll = np.maximum(Ll, lam2[j, l] - np.abs(grid - ts[j]))
        P += LAYER_W[l] * norm.cdf(Ll)
    P_tail = (1 - EPS_TAIL) * LAYER_W.sum()
    for beta_ in (.95, .50):
        ok = P >= beta_ + 0.2 * DT
        if not tail_closed or P_tail < beta_:
            out[("b199", beta_)] = None; continue
        # last grid index that fails; certified from the next one
        bad = np.where(~ok)[0]
        if len(bad) == 0: tc = 0.0
        elif bad.max() + 1 >= len(grid): tc = t_start          # only the deterministic tail certifies
        else: tc = grid[bad.max() + 1]
        out[("b199", beta_)] = float(tc)
    diag = dict(J=J, m=m, net_t_min=float(ts.min()), net_t_max=float(ts.max()), tail_closed=tail_closed,
                t_global=tsum["t_global"], n_tail_intervals=tsum["n_intervals"],
                lambda_fixed_max=float(lam.max()), draws_total=int(Ms.sum()))
    return out, diag

if __name__ == "__main__":
    rows = []; diags = {}
    truth = {}
    for d in ORDER:
        for r in read_jsonl(os.path.join(HERE, "truth_%s.jsonl" % d)):
            truth[(d, r["version"], r["target"])] = r
    for d in ORDER:
        res, diag = certify_direction(d, ETA_NET)
        res01, _ = certify_direction(d, 0.01)
        diags[d] = diag
        for (ver, beta_), tc in res.items():
            tr = truth.get((d, ver, beta_))
            row = dict(direction=d, version=ver, target=beta_, t_cert=tc, n_cert=n_of_t(tc) if tc is not None else None,
                       t_cert_eta01=res01[(ver, beta_)],
                       t_true=tr["t"] if tr else None, se_t_true=tr["se_t"] if tr else None,
                       n_true=n_of_t(tr["t"]) if tr else None)
            if tc is not None and tr:
                row["ratio_n"] = (tc / tr["t"]) ** 2
                row["ratio_n_ci95"] = [(tc / (tr["t"] + 2 * tr["se_t"])) ** 2, (tc / (tr["t"] - 2 * tr["se_t"])) ** 2]
            if ver == "fixed" and beta_ == .95:
                row.update(n_crown=N_CROWN[d], t_crown=T_CROWN[d], n_true_paper=N_TRUE_PAPER[d], t_emp_paper=T_EMP_PAPER[d],
                           ratio_crown_paper=N_CROWN[d] / N_TRUE_PAPER[d],
                           ratio_n_vs_paper_truth=(n_of_t(tc) / N_TRUE_PAPER[d]) if tc is not None else None)
            rows.append(row)
    json.dump(dict(eta_net=ETA_NET, eta_crit=ETA_CRIT, eps_tail=EPS_TAIL, t_bar=T_BAR, levels=LEVELS,
                   layer_weights=LAYER_W.tolist(), rows=rows, diagnostics=diags),
              open(os.path.join(HERE, "SMOOTH_CERT.json"), "w"), indent=1)
    print("\nFIXED test (S > %.4f), target .95   [eta = %.3f net + %.3f critical value]" % (T_BAR, ETA_NET, ETA_CRIT))
    print("%-10s %8s %8s | %8s %8s %8s | %8s %8s %8s %8s" % ("direction", "t_true", "SE", "t_CROWN", "t_smooth", "t(eta.01)", "n_true", "n_CROWN", "n_smooth", "ratio"))
    for r in rows:
        if r["version"] == "fixed" and r["target"] == .95:
            print("%-10s %8.3f %8.4f | %8.2f %8.3f %8.3f | %8.0f %8.0f %8.0f %8.2f  [95%% CI %.2f-%.2f]  (paper: n_true %d, CROWN ratio %.1f)"
                  % (r["direction"], r["t_true"], r["se_t_true"], r["t_crown"], r["t_cert"], r["t_cert_eta01"], r["n_true"], r["n_crown"],
                     r["n_cert"], r["ratio_n"], r["ratio_n_ci95"][0], r["ratio_n_ci95"][1], r["n_true_paper"], r["ratio_crown_paper"]))
    for ver, beta_ in [("fixed", .50), ("b199", .95), ("b199", .50)]:
        print("\n%s test, target %.2f" % (ver.upper(), beta_))
        for r in rows:
            if r["version"] == ver and r["target"] == beta_:
                if r["t_cert"] is None:
                    print("%-10s  NOT CERTIFIED (tail open or layer cake below target)" % r["direction"]); continue
                print("%-10s t_true %8.3f (SE %.4f)  t_smooth %8.3f  n_true %8.1f  n_smooth %8.1f  ratio %.2f [%.2f-%.2f]"
                      % (r["direction"], r["t_true"], r["se_t_true"], r["t_cert"], r["n_true"], r["n_cert"], r["ratio_n"],
                         r["ratio_n_ci95"][0], r["ratio_n_ci95"][1]))
    print("\ndiagnostics:")
    for d, g in diags.items():
        print("  %-10s J=%d net points (%.1f M draws) on [%.2f, %.2f]; tail %s to t_global %.0f with %d CROWN intervals; max Lambda %.2f"
              % (d, g["J"], g["draws_total"] / 1e6, g["net_t_min"], g["net_t_max"], "closed" if g["tail_closed"] else "OPEN",
                 g["t_global"], g["n_tail_intervals"], g["lambda_fixed_max"]))
