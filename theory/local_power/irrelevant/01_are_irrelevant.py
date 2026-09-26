"""01_are_irrelevant.py -- local sample-size ratio n*_arm / n*_DG at power one half as a function of the
number k of irrelevant covariates, in the niche design of theory/niche/run_niche.R.

For each (variant, k, departure): the limit power of every arm on the same drift h = tau v g, tau_50 by
bisection, and the ratio (tau_50^arm / tau_50^DG)^2 (above 1: DeepGOF-1 needs fewer observations).
Arms: DG (deployed, B = 199), HL, Stukel, and proj with the benchmark's FIRST weight (proj_kernel, a length
factor |x_i - x_l||x_j - x_l| on standardised covariates). That weight is superseded: the projection column of the
paper is 08_proj_liu_raw.py (Liu et al.'s angle-only weight on the raw covariates, PROJ_LIU_RAW.json).
Variants: 'study' (c_irr = 0.15, the niche design) and 'pure' (c_irr = 0: the irrelevant covariates are
absent from the null predictor, the case of the proposition).

Resumable: every finished record is appended to ARE_RAW.jsonl and skipped on rerun; several processes
may run on disjoint k lists (python -u 01_are_irrelevant.py 0 4  runs k = 0 and 4).
Assemble ARE_IRRELEVANT.json from the jsonl with:  python 01_are_irrelevant.py assemble
"""
import os, sys, json, time
import numpy as np
from ir_core import *

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "ARE_RAW.jsonl")
KS = [0, 1, 2, 4, 6, 8]
VARIANTS = {"study": dict(cirr=C_IRR, seeds=[1, 2], M_proj=1000),
            "pure": dict(cirr=0.0, seeds=[1], M_proj=1000)}

def load_done():
    done = {}
    if os.path.exists(RAW):
        for line in open(RAW):
            try:
                r = json.loads(line)
            except ValueError:
                continue                                        # a truncated last line (power loss) is dropped
            done[(r["variant"], r["k"], r["dep"])] = r
    return done

def assemble(done):
    out = {"design": "x ~ N(0, I_{2+k}); eta0 = -0.5 + x1 + 0.8 x2 + c_irr sum_{j>=3} x_j; local departure tau g(x1,x2)/sqrt(n)",
           "reading": "ratio = (tau50_arm / tau50_DG)^2 = n*_arm / n*_DG at limit power one half, alpha = .05; "
                      "DG and proj at the deployed B = 199",
           "k": KS, "variants": {}}
    for vname in VARIANTS:
        out["variants"][vname] = {dep: {str(k): done[(vname, k, dep)] for k in KS if (vname, k, dep) in done} for dep in DEPS}
    json.dump(out, open(os.path.join(HERE, "ARE_IRRELEVANT.json"), "w"), indent=1)
    print("wrote ARE_IRRELEVANT.json with %d records" % len(done))

args = sys.argv[1:]
if args and args[0] == "assemble":
    assemble(load_done()); sys.exit(0)
my_vars = [a[2:] for a in args if a.startswith("v=")] or list(VARIANTS)      # v=pure restricts the variant
my_ks = [int(a) for a in args if not a.startswith("v=")] or KS
VARIANTS = {vn: VARIANTS[vn] for vn in my_vars}
done = load_done()
print("bench %.1e | %d records already done | this process: variants %s, k in %s" % (check_bench(), len(done), my_vars, my_ks))

def log(rec):
    with open(RAW, "a") as fh:
        fh.write(json.dumps(rec) + "\n")

for vname, V in VARIANTS.items():
    for k in my_ks:
        if all((vname, k, dep) in done for dep in DEPS):
            continue
        t0 = time.time()
        des = GridDesign(k, cirr=V["cirr"])
        dgp = DGPower(des); hlp = HLPower(des)
        size_dg = dgp.power(np.zeros(K * K))
        print("[%s k=%d] design + DG null law %.0f s; DG size %.4f" % (vname, k, time.time() - t0, size_dg)); sys.stdout.flush()
        # quadratic-form tests at reference samples: kernels once per seed, drift directions cached, Z dropped
        qfs = {"proj": []}
        for sd in V["seeds"]:
            xp, Xp, pi0p, vp = ref_sample(k, V["M_proj"], 2000 * sd + k, cirr=V["cirr"])
            for arm, q, xx, vv in (("proj", QuadForm(proj_kernel(xp), Xp, vp, seed=sd + 100), xp, vp),):
                for dep in DEPS:
                    q.prep(vv * gfun(dep, xx[:, 0], xx[:, 1]), dep)
                q.drop(); qfs[arm].append((sd, q))
            print("[%s k=%d] seed %d kernels %.0f s" % (vname, k, sd, time.time() - t0)); sys.stdout.flush()
        for dep in DEPS:
            if (vname, k, dep) in done:
                continue
            a, ahl, zh, b = des.drift(h_local(dep))
            d1 = des.delta_dg(a, b); dhl = des.drift_hl(ahl, b); ncp1 = des.ncp_stukel(zh, b)
            rec = dict(variant=vname, cirr=V["cirr"], k=k, p=2 + k, dep=dep,
                       delta_norm=float(np.linalg.norm(d1)), stukel_ncp1=ncp1,
                       hl_drift_norm=float(np.linalg.norm(dhl)), size={"DG": size_dg}, tau50={}, arms={})
            rec["tau50"]["DG"] = tau_for(lambda t: dgp.power(t * d1))
            rec["tau50"]["HL"] = tau_for(lambda t: hlp.power(t, dhl))
            rec["tau50"]["Stukel"] = tau_for(lambda t: power_stukel(t, ncp1))
            rec["size"]["HL"] = hlp.power(0, dhl); rec["size"]["Stukel"] = ALPHA
            for arm in ("proj",):
                per = []
                for sd, q in qfs[arm]:
                    pw = (lambda t, q=q: q.power_bootstrap(t, dep))
                    per.append(dict(seed=sd, tau50=tau_for(pw), size=pw(0.0), signal=q.signal(dep),
                                    noise=q.noise(), snr=q.signal(dep) / q.noise(), df=q.df,
                                    max_share=q.max_share, n_neg=q.n_neg, n_keep=q.n_keep,
                                    cross_sd=q.cross_sd(dep), sum_lam=q.sum_lam, sum_lam2=q.sum_lam2))
                rec["arms"][arm] = per
                rec["tau50"][arm] = float(np.mean([r["tau50"] for r in per]))
                rec["size"][arm] = float(np.mean([r["size"] for r in per]))
            t = rec["tau50"]["DG"]
            rec["ratio_vs_DG"] = {arm: (rec["tau50"][arm] / t) ** 2 for arm in rec["tau50"] if arm != "DG"}
            rec["power_at_tauDG"] = {"HL": hlp.power(t, dhl), "Stukel": power_stukel(t, ncp1)}
            for arm in ("proj",):
                rec["power_at_tauDG"][arm] = float(np.mean([q.power_bootstrap(t, dep) for sd, q in qfs[arm]]))
            rec["secs"] = time.time() - t0
            log(rec); done[(vname, k, dep)] = rec
            r = rec["ratio_vs_DG"]
            print("  %-5s tau50 DG %6.2f | n*/n*_DG: HL %5.2f St %5.2f proj %5.2f  (%.0f s)"
                  % (dep, t, r["HL"], r["Stukel"], r["proj"], rec["secs"]))
            sys.stdout.flush()
assemble(load_done())
