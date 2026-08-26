"""probit_common.py -- the matched-level device and the NA policy, in ONE place. (2026-08-21)

Both the per-block scorer and the cross-seed verdict import from here, so the two can never drift
apart on the comparison rule -- which is the part a referee will check hardest.

TWO THINGS WERE LEARNED THE HARD WAY IN THE FIRST BLOCK; both are encoded here.

1. NaN p-values are REAL at these sample sizes. At n=50, 1.71% of datasets get an NA from at least
   one of the eleven tests (Pigeon-Heyse 1.71%, HL-equalwidth .39%, DEF.poly3 .39%, DEF.poly2 .16%,
   DEF.stukel .16%). This is not a nuisance to be dropped -- it is a measurement of small-n
   fragility, and gate S1 was pre-declared to report exactly it.

   POLICY: an NA is a test that FAILED, and a test that failed cannot reject. NaN -> p = 1.0,
   applied identically to nulls and alternatives so no arm gains or loses from the convention.
   The rate itself is carried into the summary as `na_rate` and reported, never silently dropped.

2. R's `mean()` propagates NA, so the Cauchy combination inherited every one of those failures
   (1.71% NA) -- and an NaN threshold silently poisoned the matched-level device: np.unique sorts
   NaN LAST, `(p <= nan)` is False everywhere so its empirical F reads 0.0, that 0.0 passes an
   `F <= alpha` filter, and taking the last passing index picked the NaN as the threshold. Every
   comparison against it was then False and the arm scored a power of exactly 0.0000.

   POLICY: the CCT is recomputed here from the eleven deposited p-values with `nanmean`, which
   combines the components that DID work. NA rate falls to 0. The naive R column is kept as
   `cct_naive` for the record, never as a comparator.
"""
import numpy as np

ALPHA = 0.05


def clean(p):
    """NaN -> 1.0 (a failed test cannot reject). Returns (cleaned, na_rate)."""
    p = np.asarray(p, dtype=float)
    bad = ~np.isfinite(p)
    out = p.copy()
    out[bad] = 1.0
    return out, float(bad.mean())


def narobust_cct(P):
    """Cauchy-combine an (ndat x 11) block of p-values, skipping components that failed."""
    P = np.clip(np.asarray(P, dtype=float), 1e-12, 1 - 1e-12)
    with np.errstate(invalid="ignore"):
        t = np.tan((0.5 - P) * np.pi)
    allbad = ~np.isfinite(t).any(axis=1)
    m = np.full(len(P), np.nan)
    ok = ~allbad
    m[ok] = np.nanmean(np.where(np.isfinite(t), t, np.nan)[ok], axis=1)
    c = 0.5 - np.arctan(m) / np.pi
    c[allbad] = 1.0                      # every component failed: the combination cannot reject
    return c


def randomized_matched_power(p_null, p_alt, alpha=ALPHA):
    """Place an arm at EXACTLY `alpha` on the null sample, then read its power.

    Discrete p-values (the deployed bootstrap lives on a 1/200 grid) cannot generally hit alpha
    exactly, and rounding down would quietly hand an advantage to whichever arm happens to land
    closer. The randomized test fixes that: reject outright below t0, and with probability gamma on
    the boundary atom, chosen so the null rejection rate is alpha to the digit. Applied identically
    to every arm, so the power comparison is fair by construction.

    Returns (power, F0, t0, gamma, na_rate).
    """
    pn, na_n = clean(p_null)
    pa, na_a = clean(p_alt)
    t = np.unique(pn)
    assert np.isfinite(t).all(), "non-finite threshold survived clean() -- refusing to score"
    F = np.searchsorted(np.sort(pn), t, side="right") / len(pn)
    ok = np.where(F <= alpha)[0]
    na = float(max(na_n, na_a))
    if len(ok) == 0:                       # even the smallest atom over-rejects
        g = alpha / F[0] if F[0] > 0 else 0.0
        return float(g * (pa <= t[0]).mean()), 0.0, float(t[0]), float(g), na
    i = ok[-1]
    t0, F0 = t[i], F[i]
    if i + 1 >= len(t):                    # nothing above it: no randomization possible
        return float((pa <= t0).mean()), float(F0), float(t0), 0.0, na
    t1, F1 = t[i + 1], F[i + 1]
    g = 0.0 if F1 <= F0 else (alpha - F0) / (F1 - F0)
    base = (pa <= t0).mean()
    atom = ((pa > t0) & (pa <= t1)).mean()
    return float(base + g * atom), float(F0), float(t0), float(g), na
