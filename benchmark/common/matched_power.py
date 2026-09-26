"""matched_power.py -- the matched-level device and the NA policy, in ONE place.

Every scorer and the table generator import from here, so they can never drift apart on the
comparison rule. The two functions are copied verbatim from the module the deposited numbers were
computed with.

NA POLICY: an NA p-value is a test that FAILED, and a test that failed cannot reject. NaN -> p = 1.0,
applied identically to nulls and alternatives so no arm gains or loses from the convention. The rate
itself is carried into the summaries as `na_rate`.
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
