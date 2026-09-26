"""03_trim_miss.py -- the only non-trivial term of Theorem 9: P{rho(X) != r}, the chance that the
covariate-only rule fails to delete the corrupted record. It depends on the covariate law and the
corruption mechanism ONLY (no outcomes, no fit), so it is computed exactly up to Monte Carlo error by
simulating covariates alone. Sample median and MAD, as the rule uses (the population-quantity figure in
the first draft of Corollary 9a was about 30x too small; adversarial review 2026-09-23)."""
import numpy as np

rng = np.random.default_rng(20260923)
R, CUT = 400000, 4.5
def o_of(X):
    med = np.median(X, axis=1, keepdims=True)
    mad = 1.4826 * np.median(np.abs(X - med), axis=1, keepdims=True)
    return (np.abs(X - med) / mad).max(2)

for n in (30, 50):
    X = rng.normal(size=(R, n, 2))
    oc = o_of(X)
    touch = (oc.max(1) > CUT).mean()
    # decimal slip on x1 of a uniformly chosen record (record 0 w.l.o.g.; the law is exchangeable)
    Xs = X.copy(); Xs[:, 0, 0] *= 10
    o = o_of(Xs)
    trim = np.where(o.max(1) > CUT, o.argmax(1), -1)
    miss = (trim != 0).mean()
    big = np.abs(X[:, 0, 0]) > 0.45                         # slips that move the value past ~4.5 robust sd
    print("n=%d  clean data touched: %.4f | decimal slip missed: %.4f  (given |x|>.45: %.4f ; given |x|<=.45: %.4f)"
          % (n, touch, miss, (trim[big] != 0).mean(), (trim[~big] != 0).mean()))
