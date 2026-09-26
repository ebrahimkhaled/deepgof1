# W4 — Finite-sample power: what is certified, what is exact, and why the gap is where it is

**Decision (the author, 2026-09-23): option A, the bracket.**

## The two results
1. **A sound finite-sample power guarantee (probabilistic CROWN, `07_pcrown.py` → `PCERT.json`).**
   For map noise $e\sim N(0,I)$: every intermediate bound is a linear function of $e$ minus 4.47 of its
   Gaussian sd (union bound over the 3,200 ReLU units, total .025); the final linear lower function's
   Gaussian tail takes the other .025. So power ≥ .95 is **proved** at the certified signal $t$. Checked on
   200,000 draws per direction: bounds were violated 0.35% of the time against a budget of 2.5%, and the
   Monte Carlo power at the certified point is 1.000.

   | direction | guaranteed n* | true n* (simulation) | ratio |
   |---|---|---|---|
   | step | 4,113 | 72 | 57× |
   | diffuse | 4,726 | 103 | 46× |
   | stripe | 6,420 | 122 | 53× |
   | saddle | 8,960 | 204 | 44× |
   | random | 15,609 | 729 | 21× |
   (signal ‖μ‖ = 0.5 per √n)

2. **The exact local-power formula (W1, Theorem L)**, computed from the frozen weights before any data,
   and matched by simulation up to n = 32,000 (`local_power/CONVERGENCE.json`).

Together these bracket the truth. The first is a guarantee. The second is the precise large-sample value.

## Why no tighter guarantee is available with present tools (measured, not asserted)
| route | best it can do | evidence |
|---|---|---|
| norm-product constants | 10⁶–10⁷× | `01_constants.py` |
| CROWN over the 95% noise ball | 90–150× | `03_certify.py` |
| α-CROWN (final layer) | about 2% better | `04_alpha_crown.py` |
| **exact** worst case over the ball | 4.4–6.7× floor | `05_ball_limit.py`, by attack |
| probabilistic CROWN (ours) | **21–57×** | `07_pcrown.py` |
| α-CROWN with every intermediate bound optimised | about 1.5× gain | `09_alpha_full.py` |
| CROWN with **perfect** intermediate bounds | ≈10× (stripe) | `10_oracle_bounds.py` (oracle diagnostic) |
| Lipschitz + Gaussian concentration, blocks | ≥ 25.9 → ≈50× **for any block certificate** | `11_lipsdp.py`: head 1.823 (sampled 1.786), block A sampled 4.04 |

The limit is structural. At the signal strengths that matter, more than 2,000 of the 3,200 ReLU units
are undecided under noise. Any method that replaces an undecided unit by a linear relaxation loses a
little per unit, and the losses add up. Beating it needs case-splitting over those units, which grows
exponentially.

## How to present it (defensive framing)
"A finite-sample power guarantee exists and is computable from the frozen weights; for the canonical
departures it certifies 95% power at 21–57 times the sample size simulation shows is needed. The exact
large-sample power is given by Theorem L. The gap between the two is the cost of certifying a
3,200-unit piecewise-linear function under noise, and we show it is not an artefact of our bounding: even
perfect intermediate bounds leave a factor of ten." Supplement, not headline.

## Update 2026-09-23: branch-and-bound done (`bab/`, README.md, CERT_UA_*.json)
* **Usage-aware charging** (the probability budget pays only for relaxation statements actually used; a
  sound refinement, validated on 200,000 draws per direction) tightens every certificate. Certified t:
  stripe 40.06 → 38.00, step 32.07 → 30.78, saddle 47.33 → 44.84, diffuse 34.37 → 32.83,
  random 62.47 → 58.99. **The paper's range becomes 19–53× (was 21–57×).** Use these numbers.
* **Branch-and-bound cannot reach 10× (measured):**
  * Input splitting needs at least 24 of the 36 directions pinned.
  * Pre-activation splitting (BaB-prob style) makes the bound worse, because the loose deep-layer forms
    overlap across the split.
  * Even an oracle for the band-failure probability reaches only t ≈ 35.
  * The limit is relaxation error compounding through 10 layers under bands about 4.2 sd wide.
* Survey (2026-09-23): no published neural-network test proves any finite-sample power guarantee.
  Minwoo Kim et al. (Stat. Comput. 2025) prove no theorems. Finite-sample power results exist only as
  rate theorems for kernel tests (MMDAgg, KSDAgg, Kim–Balakrishnan–Wasserman 2022).

## Update 2026-09-24: code audit of the usage-aware charge (proof review, Proposition S.12)
Question from the proof review: does the code compute the "used" set U as the least set containing the
final-form units and closed under "appears in the CROWN back-substitution form of a member", with each unit
charged the strongest requirement over all forms it enters?

Finding (read `bab/bab_core.py::cell_pcrown_ua`, called by `bab/07_certify.py`, `bab/03_bab.py`,
`bab/08_validate_strat.py`):
* Each unit's band (l, u), and hence its relaxation, is computed once and is the same in every form.
* The code charges a CLOSED SUPERSET of the least closed set, not the least set itself:
  - every stable-active ReLU unit is charged Phi(-a_lo/s_lo) and every stable-inactive unit
    Phi(-a_up/s_up), whether or not the unit appears in any form that matters;
  - an unstable unit is charged 2 Phi(-z) iff it meets a negative coefficient in the back-substitution of
    ANY computed form (`used[s] |= un & (Ln != 0).any(0)`), i.e. the lower and upper forms of every unit in
    every layer plus the final form, not only the forms of members of U;
  - an unstable unit met only with nonnegative coefficients is charged 0 (its lower line 0 or x is valid
    for every x).
* Max-over-forms charging: yes. The only form-dependent requirement is for unstable units (negative
  coefficient needs l <= x <= u; nonnegative needs nothing), and `any` over all forms charges the two-sided
  requirement if any form needs it. Stable units have the same requirement in every form.
* The charged set is closed (every form of every unit is scanned) and contains U, so the charge is >= the
  charge of U: the certificate is SOUND, possibly slightly conservative relative to the least closed set.
* Consequence: no code change, no rerun. The reported numbers (certified t: step 30.78, diffuse 32.83,
  stripe 38.34 one region / 38.00 with 33 regions, saddle 44.84, random 58.99; ratios 19-53x) stand.
* The ESM proof of Proposition S.12 now defines U as the least closed set with max-over-forms charging and
  states that the implementation charges a closed superset of it.
