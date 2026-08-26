# 03 — C1 DEFECTS + the deployment gap (adversarial review, 2026-08-17)

Source: a 10-agent review (4 independent JASA lenses → adjudicated synthesis → 4 figure
concepts → figure judge). **Every finding below was re-verified by hand against the source
before being written down.** Line numbers are `code/c1_confirmatory.R` unless stated.

C1's verdict is not withdrawn — the pre-declared criteria were met as written. What is withdrawn
is the *interpretation*: the comparison was not size-matched, the run is not instrumented, and
the evidence certifies a procedure a user cannot actually run.

---

## D1 — The comparison is not size-matched *(verified)*

```r
:155  cut  <- quantile(ense(Scal), .95)      # ensemble: independent cell EV_CAL
:156  size <- mean(ense(Snul) > cut)         # realized .059 / .040 / .065
:161  ccut <- quantile(ccts(EV_NUL), .95)    # CCT: the SAME cell its power is scored on
```

The ensemble is judged at its *realized* size; CCT-full is pinned to exactly .050 on the cell it
is then scored against. Both rules must be calibrated on `EV_CAL` and both realized sizes reported
on `EV_NUL`.

| n | reported gain | ensemble realized size | size-matched estimate |
|---|---|---|---|
| 200 | +.118 | .059 | ≈ **+.090** |
| 500 | +.041 | .040 | ≈ **+.057** ⬆ |
| 1000 | +.057 | .065 | ≈ **+.042** |

**Not all downside** — at n=500 the ensemble ran *conservative* (.040 vs .050), so the correction
moves that gain up. Fixing this makes the n=500 result stronger, not weaker.

## D2 — Nothing per-replicate was stored *(verified)*

`:172` is the only writer: 12 rounded numbers per n. The C1 corpus, the 20 fitted nets, the scaler
attributes, the feature matrices and every per-replicate score were discarded. This violates the
standing project rule that every simulation deposits a `*_pvalues*.csv`, and it means every
referee question — standard errors, size-matching, a sign test, an ROC — costs an ~18-hour rerun
instead of a re-score.

## D3 — The evaluation null is 2-covariate only *(verified)*

Training null: `:107` `four <- (i %% 10) < 3` → 70/30 two/four-covariate.
Evaluation null: `../learned-def/R/scenarios.R` `gen_data("null")` → the 4-covariate branch fires
only for `ALT_MULTI`, so the null cell is 2-covariate throughout. `multivar_quad` and
`multivar_int` are therefore scored against a threshold calibrated on the wrong design.
Fix: add a 4-covariate null cell (~72 min at n=1000) and calibrate per design, or label the
headline as a 10-family 2-covariate result.

## D4 — ⛔ RETRACTED: "the undetectable, detected"

| corpus | probit power, n=1000 |
|---|---|
| exploratory (`nnet_stability_n1000.csv`) | .147 |
| kind-holdout (`..._holdout.csv`) | **.080** |
| C1 confirmatory | .170 |

Individual-test ceiling ≈ .078–.081. In one of three trained ensembles probit sits *on* the
ceiling. Worse, `01_C1_PREDECLARATION.md` §"What C1 does NOT test" states plainly: **"No
per-family claims from C1."** The claim broke the project's own frozen rule and appeared in the
README, MANIFEST, memory, tracker and the money figure. Withdrawn from all of them.

**What survives, and it is better as mechanism than as a boast:** a departure prior containing
index-direction terms (`indexq`/`indexs` in `make_family()` — functions of η₀ alone, i.e. exactly
the probit direction) buys detection of link misspecification that no battery member achieves.
Magnitude is corpus-variable (.08–.17 at n=1000). State the mechanism, report the range.

## D5 — 🔴 THE BIGGEST RISK: the evidence certifies a test nobody can run

Every size number in the project — .059/.040/.065 — comes from `evcell("cal","null",2000)`, drawn
from the **true** null generator with known β = (0.3, 0.8, −0.5) and a known covariate law. A user
has one dataset and must calibrate by parametric bootstrap from β̂.

The only bootstrap-calibration measurement anywhere in the family,
`../learned-def/results/bootstrap_calibration_table.csv`: **size .083 at α = .05** — on 120
replicates, on a reduced 4-feature rule. Underpowered, but pointing the wrong way.

If the shipped 33-feature ensemble bootstraps at .08, then (a) the power gains largely evaporate
(an 8% test against a 5% test), and (b) the intended abstract sentence — *"validity does not rest
on the learner; size is exact by construction"* — is **false as shipped**, in the same register
that already cost this family two fit-rejections.

**Run the deployed-calibration study before writing a word.** If size is inflated there, that is
the paper's most important finding and it leads the paper.

## D6 — Two smaller corrections

- **"the strongest classical rule" is not defensible at n=200.** `results/smalln_curve.csv` and
  `results/oracle_room.csv` show gamma-tilt (.3448) and ridge (.3378) above CCT-full (.3075).
  The gamma-tilt edge is itself noise (sign test 7–1–4, p = .549), but the honest n=200 gain
  against best-of-three is ≈ **+.07**, not +.118. Say it before a referee does.
  Always call CCT-full "the pre-declared comparator", never "the strongest".
- **"one forward pass" is wrong** (`02_C3_LITERATURE.md` §A). Under bootstrap calibration the
  use-time cost is B refits of the whole battery. Fix before it reaches an abstract.
- **Good news from the ablation check:** `diverse_training_n1000.csv` gives ridge-ALL .7620 and a
  single net .7628 vs CCT .7459 and the 20-net ensemble .7955 — ridge captures ~a third of the
  gain, the ensemble adds the other two-thirds. The "deep" in the title is earned; the ablation
  is mandatory but winnable.
