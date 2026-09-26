# One Stukel test everywhere: SstBoth (Liu et al. 2024) — result, 2026-09-26

Referee report `manuscript/sc_audit/REFEREE_FINAL.md`, item B1. The paper cites one Stukel test
(`stukel1988, liu2024`) and Table 2 and the niche already use SstBoth. Three places still used the joint
2-df score form: `tab:are`, `fig:irrelevant`, and the two corrupted-record studies (`tab:corrupt`,
`tab:trim`, ESM `tab:studyone`, `tab:trimS`). All of them are recomputed here with SstBoth, which is
`ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))`: the sum of the squares of the two marginal Stukel
score z-statistics (`statmod::glm.scoretest`), referred to chi-square(2). **The manuscript was not edited.**

## 0. Headline

* **Pairing is exact.** Every one of the 23,200 small-n datasets was regenerated from its seed and reproduced
  the stored HL and joint-Stukel p-values exactly (max difference 0.0; LOO maximum, trimmed index and trimmed
  p-values too). The 18,000 local-power simulation datasets reproduced the stored joint p to 1.3e-14.
* **The SstBoth limit is not chi-square(2).** The two marginal scores are correlated after the fit
  (rho = −0.613 in the tab:are design; −0.613 at every k in the niche design). So under H0 the statistic
  tends to 1.613 U1² + 0.387 U2², and the deployed chi-square(2) rule has **limit size .064**, not .05.
  Simulation agrees: .071/.063/.066 at n = 2,000/8,000/32,000 (SE .008); niche .075 ± .006 (n = 200),
  .068 ± .006 (n = 500). At n ≤ 50 the clean size is .061/.077 (n = 30) and .087/.062 (n = 50) in the
  two studies.
* **tab:are, Stukel column:** DeepGOF-1 now needs fewer observations than Stukel on **two** of five
  departures (bump, threshold), not three. SstBoth is more efficient on the interaction as well.
* **fig:irrelevant:** SstBoth's ratio is **0.75–0.98 (interaction) and 0.90–1.22 (quadratic)**, against
  1.06–1.39 and 1.28–1.74 for the joint form. It is **below one** for the interaction at every k and for the
  quadratic up to k = 2. The caption's "Hosmer–Lemeshow and Stukel stay above one at every k" is no longer
  true.
* **Study 1 (tab:corrupt):** the decimal slip still breaks Stukel at both n (.104*, .146*). Under SstBoth,
  **sign (n = 30, .079) and outcome (n = 30, .086) also exceed .075**, and so does the **clean** cell at
  n = 50 (.087). This reflects SstBoth's own over-size, not the corruption: each is within about 2 SE of
  that cell's clean size. The LOO theorem check still holds in all cells.
* **Study 2 (tab:trim): the pre-declared check (1) FAILS for SstBoth at n = 50.** Trimmed decimal-slip size
  .089 against a clean size of .062, where the limit is .0856 (2.3 SE). It passes at n = 30 (.086 vs clean
  .077). Check (2) (no harm when clean) passes in all 8 Stukel cells. So "the rule passed its pre-declared
  check in every cell" is no longer true with SstBoth as the Stukel arm.

## 1. Pre-declarations and whether a verdict changes

* `theory/small_n/PREDECLARATION.md` **names the form**: "Stukel (joint score, `ebrahim.gof` 2.8.0)".
  `PREDECLARATION_TRIM.md` names no form. It says "Same four tests" as the first study, so it inherits the
  joint form. **Reporting SstBoth in both studies is therefore a post hoc deviation from both
  pre-declarations.** It is computed on exactly the same stored datasets, by the same pre-declared
  rules, and must be reported as such. Suggested wording is in §5.
* Pre-declared verdicts about Stukel:

| check | joint (pre-declared form) | SstBoth | verdict changes? |
|---|---|---|---|
| Q1 misbehaves (>.075) under decimal slip | yes, both n (.079, .142) | yes, both n (.104, .146) | no |
| Q1 under sign / outcome | no (max .071) | **yes at n = 30** (.079, .086); n = 50 no (.062, .071) | **yes** (but clean n = 30 is .061, n = 50 is .087: this is the test's own over-size) |
| Q2 LOO size ≤ clean + 2 SE | 8/8 cells hold | 8/8 hold | no |
| Trim check (1), protection | passes n = 30 and n = 50 | passes n = 30, **FAILS n = 50** (.089 > .0856) | **yes** |
| Trim check (2), no harm when clean | 8/8 | 8/8 | no |
| Trim check (3), mild slips missed (size) | .077 / .124 | .097 / .129 | numbers only |
| Trim check (4), sign/outcome | harmless | harmless (n = 30 and n = 50 sign at about .08 before and after: over-size, unchanged by the rule) | no |

  The Study 2 failure is marginal. It is 2.3 SE above the clean size, and the clean reference in that
  study (.062) is lower than in Study 1 (.087) at the same n. Pooled over both studies (2,000 clean nulls)
  the clean size is .0745, and .089 would lie within 2 SE of it. But the pre-declared rule compares with
  the same study's clean size, and by that rule it fails.

## 2. SstBoth undefined (NA): share and handling

As in the package, SstBoth is undefined when one Stukel column is identically zero, which happens when no
fitted risk lies on one side of .5. It then counts as non-rejection (p = 1), the Table 2 convention. For
the LOO maximum, an undefined deletion gives p = 1, the original `all_p` rule `loo[is.na(loo)] <- 1`.

| where | share undefined |
|---|---|
| tab:are / fig:irrelevant limits | 0 (both sides of .5 have positive mass) |
| local-power simulation (CONV_RAW, 18,000 datasets, n ≥ 2,000) | 0 |
| niche simulation (addon rows, k = 8, n = 200/500) | 0 (`na` column of `niche_addon_rates.csv`) |
| Study 1, full data | 2.5 % (3.6 % at n = 30, 1.5 % at n = 50) |
| Study 1, at least one LOO deletion undefined (so LOO p = 1) | 7.0 % of datasets (10.0 % n = 30, 4.0 % n = 50); 2.3 % of all deletions |
| Study 2, full / trimmed data | 2.3 % / 2.3 % (3.0 % n = 30, 1.6 % n = 50) |

The joint form never had this problem: it drops the empty direction and uses 1 df.

## 3. (a) tab:are — local sample-size ratio (τ50^Stukel/τ50^DG)²

Derivation (script `01_tab_are.py`, core `sst_core.py`). Let u = n^{-1/2} z'(y − π̂) for the two Stukel
columns. Under π_n = π0 + τh/√n, u → N(τ s(h), Ω), with Ω and s(h) exactly as in ESM §8.4, because the
joint and marginal forms share the score. SstBoth is ‖D^{-1/2}u‖², D = diag Ω, so the limit is
‖W‖², W ~ N(τ D^{-1/2}s(h), R), R = D^{-1/2} Ω D^{-1/2}. This is a generalised noncentral chi-square,
Σ λ_i (U_i + c_i)², with λ = 1 ∓ |ρ| and c = Λ^{-1/2} Qᵀ μ. P(‖W‖² > χ²_{2,.95}) is computed by
one-dimensional quadrature (checked against 2·10⁶ Monte Carlo draws: .20284 vs .20292). The joint form was
recomputed alongside and reproduces FORMULA.json's τ50 to every digit. DeepGOF-1's τ50 is unchanged
(FORMULA.json).

| departure | DG τ50 | joint τ50 | joint ratio (printed) | **SstBoth τ50** | **SstBoth ratio** | limit power at DG τ50, joint → SstBoth |
|---|---|---|---|---|---|---|
| quadratic | 8.76 | 7.70 | 0.77 | 6.48 | **0.55** | .615 → .743 |
| interaction | 11.53 | 12.45 | 1.17 | 10.46 | **0.82** | .439 → .576 |
| cubic | 12.91 | 7.84 | 0.37 | 10.62 | **0.68** | .917 → .704 |
| bump | 13.28 | 22.16 | 2.79 | 21.46 | **2.61** | .205 → .235 |
| threshold | 16.70 | 78.44 | 22.1 | 66.06 | **15.6** | .067 → .093 |

Ω = [[.01369, −.01347], [−.01347, .03523]], ρ = −.613, limit size .0643.

**Simulation check** (`02_sim_check_tabare.py`). SstBoth was computed on the stored datasets of
`03_convergence.py` (`CONV_RAW.jsonl`, seeded; `SIM_RAW.json` of `02_simulate.py` stores no seeds and
cannot be regenerated). Rejection rate at DG τ50, with z = (sim − limit)/SE:

| cell | limit | n = 2,000 | n = 8,000 | n = 32,000 | joint form, limit / sim at 2k, 8k, 32k |
|---|---|---|---|---|---|
| null | .064 | .071 (+0.8) | .063 (−0.2) | .066 (+0.2) | .050 / .051 .059 .059 |
| quadratic | .743 | .783 (+3.0) | .793 (+3.9) | .761 (+1.3) | .615 / .634 .667 .622 |
| interaction | .576 | .597 (+1.4) | .579 (+0.2) | .609 (+2.2) | .439 / .436 .405 .451 |
| cubic | .704 | .180 (−43) | .469 (−15) | .615 (−5.8) | .917 / .378 .774 .875 |
| bump | .235 | .219 (−1.3) | .267 (+2.3) | .250 (+1.1) | .205 / .252 .250 .237 |
| threshold | .093 | .225 (+10) | .140 (+4.2) | .118 (+2.4) | .067 / .134 .096 .076 |

The null limit (.064) is confirmed. For the cubic and threshold, **both** Stukel forms are far from
their limits at n ≤ 32,000. The cubic is badly over-predicted (the joint form too: .378 at n = 2,000
against .917), because Stukel's added covariates are quadratic in η and the cubic's τ50 is not yet
local. The threshold is under-predicted. This was already true of the printed joint column. The ESM
sentence "The rival formulas ... match simulation under the null and show a similar finite-n bias at
n = 2,000" understates it for Stukel on the cubic and the threshold (see §6).

### Ready-to-paste rows (main `tab:are` and ESM `tab:areS`: only the Stukel column changes)

```
quadratic, $x_{1}^{2}-1$ & 1.62 & 0.55 & 0.62 (.01) \\
interaction, $x_{1}x_{2}$ & 2.52 & 0.82 & 0.73 (.02) \\
cubic, $x_{1}^{3}-3x_{1}$ & 1.01 & 0.68 & 0.26 (.01) \\
local bump & 3.36 & 2.61 & 0.87 (.03) \\
threshold, $1\{x_{2}>1\}$ & 43.7 & 15.6 & 0.67 (.01) \\
```
ESM `tab:areS` (short labels): `quadratic & 1.62 & 0.55 & 0.62 (.01)`, `interaction & 2.52 & 0.82 & 0.73 (.02)`,
`cubic & 1.01 & 0.68 & 0.26 (.01)`, `bump & 3.36 & 2.61 & 0.87 (.03)`, `threshold & 43.7 & 15.6 & 0.67 (.01)`.

Suggested caption addition (tab:are): "Stukel: the SstBoth form of \citet{liu2024}. Its two score
statistics are correlated ($\rho=-.61$), so referred to $\chi^{2}_{2}$ its limit size is $.064$, not
$.05$."

## 4. (b) fig:irrelevant — Stukel curve with SstBoth

Script `03_irrelevant.py` (ir_core.GridDesign unchanged; DG τ50 from ARE_IRRELEVANT.json; the joint ratios
are reproduced exactly). Study variant (c_irr = .15):

| k | 0 | 1 | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|---|
| interaction, joint (printed) | 1.063 | 1.098 | 1.138 | 1.214 | 1.308 | 1.394 |
| **interaction, SstBoth** | **0.749** | **0.773** | **0.801** | **0.855** | **0.921** | **0.981** |
| quadratic, joint (printed) | 1.281 | 1.332 | 1.392 | 1.494 | 1.625 | 1.738 |
| **quadratic, SstBoth** | **0.903** | **0.939** | **0.980** | **1.052** | **1.144** | **1.223** |

Pure variant (c_irr = 0): SstBoth 0.749 (k = 0) and 0.735 (k ≥ 1) for the interaction; 0.903 and 0.888
for the quadratic. It is flat in k, as Proposition S.10a(i) says. The k = 0 vs k ≥ 1 step comes from
DeepGOF-1's τ50, as in the joint column.

**Figure:** `theory/stukel_sstboth/fig_are_irrelevant.pdf` (+ `.png`). It is made by
`fig_are_irrelevant_sstboth.py`, a copy of the original script with the same style (Arial 8 pt, 84 × 58 mm,
greyscale line style + marker, no tight bbox); only the Stukel series is read from IRR_SSTBOTH.json. The
original in `theory/local_power/irrelevant/` is untouched. Copy the new PDF over
`manuscript/sc_audit/figures/fig_are_irrelevant.pdf` when the text is updated.

**Predicted vs simulated at k = 8** (niche rows, `Stukel_marginal`, 1,000 reps/cell; exact finite-n drift):

| cell | SstBoth pred | SstBoth sim (SE) | z | joint pred / sim / z (printed basis) |
|---|---|---|---|---|
| inter C.3 n200 | .121 | .110 (.010) | −1.1 | .085 / .081 / −0.5 |
| inter C.3 n500 | .207 | .208 (.013) | +0.1 | .143 / .152 / +0.8 |
| inter C.5 n200 | .232 | .178 (.012) | −4.5 | .167 / .134 / −3.1 |
| inter C.5 n500 | .467 | .411 (.016) | −3.6 | .364 / .309 / −3.8 |
| inter C.8 n200 | .565 | .360 (.015) | −13.5 | .492 / .283 / −14.7 |
| inter C.8 n500 | .916 | .767 (.013) | −11.1 | .888 / .672 / −14.5 |
| quad C.3 n200 | .171 | .167 (.012) | −0.4 | .123 / .124 / +0.1 |
| quad C.3 n500 | .331 | .292 (.014) | −2.7 | .248 / .227 / −1.6 |
| quad C.5 n200 | .455 | .301 (.015) | −10.6 | .385 / .241 / −10.6 |
| quad C.5 n500 | .824 | .647 (.015) | −11.7 | .779 / .569 / −13.4 |
| quad C.8 n200 | .917 | .478 (.016) | −27.8 | .907 / .439 / −29.8 |
| quad C.8 n500 | 1.000 | .895 (.010) | −10.8 | 1.000 / .855 / −13.0 |
| null n200 | .064 | .075 (.006) | +1.8 | .05 / .061 |
| null n500 | .064 | .068 (.006) | +0.6 | .05 / .057 |

The same conclusion holds: the limit over-predicts by **3.6 to 28** SE in every cell with C ≥ 0.5 (joint:
3.1 to 30). So the ratio in the figure is still favourable to Stukel.

## 5. (c) Corrupted-record studies — new numbers and ready-to-paste tables

Scripts: `theory/small_n/sstboth_addon/recompute_sstboth.R` (4 workers, 928 resumable 25-replicate
blocks in `rows_s1/`, `rows_s2/`, 5 min), `analyse_sstboth.py` → `RESULTS_SSTBOTH.json`,
`ANALYSE_SSTBOTH.log`.

### tab:corrupt (main) and tab:studyone (ESM): Stukel column only

| n | corruption | HL | Stukel joint (printed) | **Stukel SstBoth** | DG-1 |
|---|---|---|---|---|---|
| 30 | none | .046 | .048 | **.061** | .056 |
| 30 | decimal | .049 | .079* | **.104*** | .086* |
| 30 | sign | .033 | .063 | **.079*** | .057 |
| 30 | outcome | .034 | .071 | **.086*** | .038 |
| 50 | none | .049 | .061 | **.087*** | .055 |
| 50 | decimal | .041 | .142* | **.146*** | .118* |
| 50 | sign | .045 | .051 | **.062** | .046 |
| 50 | outcome | .041 | .054 | **.071** | .056 |

```
30 & none & .046 & .061 & .056 \\
30 & decimal slip & .049 & .104* & .086* \\
30 & sign error & .033 & .079* & .057 \\
30 & flipped outcome & .034 & .086* & .038 \\
50 & none & .049 & .087* & .055 \\
50 & decimal slip & .041 & .146* & .118* \\
50 & sign error & .045 & .062 & .046 \\
50 & flipped outcome & .041 & .071 & .056 \\
```
(ESM `tab:studyone` uses the labels decimal / sign / outcome; same numbers.) The caption's "Only the
decimal slip, a leverage point, breaks tests; sign errors and flipped outcomes do not" must be qualified.
Suggested: "... breaks tests beyond their clean size; Stukel's test (SstBoth) is above $.075$ already
on clean data at $n=50$, and its sign and outcome cells at $n=30$ are within two standard errors of its
clean size." One option is to put the asterisk only on the contaminated cells (the pre-declared
criterion) and mark .087 with a dagger ("clean size above .075").

Other Study 1 numbers (ESM):
* LOO size, Stukel: .003–.013 (joint .002–.012). All 8 Stukel cells within 2 SE of the clean size. The
  "24 of 24 cells" statement stands. (The "between .000 and .027" range in the text is wider than any
  HL/Stukel/DG cell in either form: the max is .012 for the joint form and .013 for SstBoth. .027 seems to
  come from the removed fourth arm. Check.)
* Price at n = 50, C = 1.5, clean: Stukel **.48 → .24** (joint .47 → .24).
* Post hoc α′: Stukel α′ = .26 at n = 50, decimal-slip size **.061 (from .146)** (joint .059 from .142).
* Diagnostic: at n = 50, decimal, the most-rescuing deletion is the corrupted record in **63 %** of
  Stukel's base rejections (joint 61 %, DG-1 57 %) → range "57–63 %" (if the printed 57–65 % came from the
  removed arm). Under sign/outcome at n = 50 it is at most **18 %** for SstBoth (joint 13 %, so the printed "at most 16 %"
  presumably included the removed arm) → "at most 18 %". At n = 30 it reaches 24 % (outcome; joint 28 %).

### tab:trim (main) and tab:trimS (ESM): Stukel rows

| n | test | slip, no rule | slip, with rule | clean | check (1) |
|---|---|---|---|---|---|
| 30 | Stukel joint (printed) | .091 | .066 | .053 | pass |
| 30 | **Stukel SstBoth** | **.113** | **.086** | **.077** | pass (limit .1015) |
| 50 | Stukel joint (printed) | .145 | .073 | .064 | pass |
| 50 | **Stukel SstBoth** | **.161** | **.089** | **.062** | **FAIL** (limit .0856; +2.3 SE) |

```
30 & Stukel & .113 & .086 & .077 \\
30 & DeepGOF-1 & .088 & .047 & .055 \\
50 & Stukel & .161 & .089 & .062 \\
50 & DeepGOF-1 & .111 & .055 & .055 \\
```
Other Study 2 numbers:
* Check (2): SstBoth base/trimmed, n = 30: .077/.077, .103/.100, .277/.277, .353/.353; n = 50: .062/.063,
  .203/.203, .430/.437, .490/.487. All within 2 paired SE. "24 of 24 cells" stands.
* Check (3), mild slips missed: SstBoth **.097 / .129** (n = 30/50; joint .077/.124).
* Check (4): SstBoth sign .080→.078 (n = 30), .080→.080 (n = 50); outcome .085→.084, .063→.063.

## 6. Every manuscript sentence whose Stukel number or statement changes

Line numbers are those of `manuscript/sc_audit/deepgof1_sc_audit.tex` and `esm/ESM_audit.tex` as read on
2026-09-26. They may shift under the concurrent edit.

**Main text**
1. l.95–96 (abstract) "one corrupted covariate value among fifty patients breaks this test and Stukel's": still
   true (.146 at n = 50). No change needed.
2. l.715–735 `tab:are` caption + Stukel column: new column (§3). Add the SstBoth definition and its limit
   size .064.
3. l.748–750 "it needs fewer than Stukel's link test on the interaction, the bump and the threshold, and
   *more* on the quadratic and the cubic, up to 2.7 times as many observations on the cubic" →
   "it needs fewer than Stukel's test on the bump ($2.61$) and the threshold ($15.6$), and *more* on the
   quadratic, the interaction and the cubic, up to $1.8$ times as many observations on the quadratic
   ($0.55$)."
4. l.769–770 "The same is true of the Hosmer–Lemeshow and Stukel tests, whose ... added terms are built from
   the fitted linear predictor": still true for SstBoth (flat in the pure variant). No change.
5. l.794–797 "Stukel's test needs $1.06$–$1.39$ (interaction) and $1.28$–$1.74$ (quadratic) times as many
   observations as DeepGOF-1, and the Hosmer–Lemeshow test 2.3–3.7 times; both drift up slowly with k" →
   "Stukel's test needs $0.75$–$0.98$ (interaction) and $0.90$–$1.22$ (quadratic) times as many observations
   as DeepGOF-1: fewer for the interaction at every k and for the quadratic up to $k=2$. The Hosmer–Lemeshow
   test needs 2.3–3.7 times as many; both ratios drift up slowly with k ..." The next sentence, "DeepGOF-1
   therefore needs fewer observations than the projection test only once about two to four covariates ...",
   stays, but no longer implies that DeepGOF-1 beats every rival there.
6. l.817–820 "Stukel's limit fails first: it over-predicts the simulated power by 3 to 30 standard errors in
   every cell with C ≥ 0.5" → "by 4 to 28 standard errors". The rest stands ("favourable to Stukel"). Add
   "its limit size is .064; the simulated null size is .075 ± .006 (n = 200) and .068 ± .006 (n = 500)" if
   the null sizes are quoted.
7. l.826–840 `fig:irrelevant` caption: "Hosmer–Lemeshow and Stukel stay above one at every k" →
   "Hosmer–Lemeshow stays above one at every k; Stukel's test (SstBoth) is below one for the interaction at
   every k and for the quadratic up to k = 2". Swap in the new PDF.
8. l.1272–1283 `tab:corrupt`: new Stukel column (§5). Caption sentence on sign/outcome must be qualified.
9. l.1286–1290 "A slipped decimal pushes Stukel's test and DeepGOF-1 above the pre-declared .075 at both
   sample sizes ... Sign errors and flipped outcomes break nothing by that criterion, although at n = 30 one
   cell comes close (.071)" → the decimal part stands. Replace the last clause: "Sign errors and flipped
   outcomes do not raise any test above its own clean size; Stukel's test in the SstBoth form is above .075
   in those cells at n = 30 (.079, .086) and on clean data at n = 50 (.087), because its chi-square(2)
   reference is too small (limit size .064; Section 3.2)." Add the one-line deviation note: "The
   pre-declaration named the joint score form of Stukel's test; both studies are reported here with the
   SstBoth form of Table 2, computed on the same stored datasets (the joint form gives .079/.142 under the
   slip)."
10. l.1328–1333 `tab:trim`: new Stukel rows (§5).
11. l.1338–1342 "The rule passed its pre-declared check in every cell: each broken test returned to within
    two standard errors of its own clean size ... Stukel's test improves but stays above that band (.066 at
    n = 30 and .073 at n = 50 with the rule, against clean sizes of .053 and .064)" → "The rule passed its
    pre-declared check for DeepGOF-1 at both sample sizes and for Stukel's test at n = 30. For Stukel's test
    at n = 50 it failed narrowly: .089 with the rule against a clean size of .062 (2.3 standard errors; its
    clean size in the first study was .087). Stukel's test improves but stays above the band (.086 and .089
    with the rule, against clean sizes of .077 and .062)."
12. l.1343–1346 "left DeepGOF-1 at its nominal size but not Stukel's test (.124 at n = 50)" → "(.129 at n = 50)".
13. l.1363–1366 "for Stukel's test it is not": stands.
14. l.1484–1485 (Discussion) "fewer than Stukel's link test on three of five" → "on two of five".
15. l.1517 / l.1519 "Stukel's test improves but stays above the band": stands.
16. l.1320–1322 (first repair) "size between .000 and .027": see the §5 note. The HL/Stukel/DG range is .000–.013.
    "at n = 50 and C = 1.5 DeepGOF-1 falls from .55 to .26": unaffected.

**ESM**
17. l.1965–1981 §8.4 (`sec:rivals`): replace "a score test with 2 degrees of freedom" and the noncentral
    chi-square limit with:
    "*Stukel's test* in the SstBoth form of \citet{liu2024} \citep{stukel1988}. Its two added covariates
    are [z(x) as now]. Each is tested alone: $Z_{j}=u_{j}/\Omega_{jj}^{1/2}$ with
    $u=n^{-1/2}\sum_i z(x_i)\{y_i-\hat\pi_i\}$, and the statistic $Z_{1}^{2}+Z_{2}^{2}$ is referred to
    $\chi^{2}_{2}$. [Ω and s(h) as now.] With $D=\mathrm{diag}(\Omega)$, the statistic tends in law to
    $\|W\|^{2}$, $W\sim N(\tau D^{-1/2}s(h),\,D^{-1/2}\Omega D^{-1/2})$: a sum
    $(1+|\rho|)(U_{1}+c_{1})^{2}+(1-|\rho|)(U_{2}+c_{2})^{2}$ of two independent noncentral terms, where
    $\rho$ is the correlation of the two scores after the fit, evaluated by one-dimensional quadrature.
    In the design of \S\ref{sec:localnum}, $\rho=-.61$, so the limit size of the deployed rule is
    $.064$, not $.05$; simulation gives $.071$, $.063$ and $.066$ at $n=2{,}000$, $8{,}000$ and $32{,}000$.
    The statistic is undefined when no fitted risk lies on one side of $.5$, which has probability tending
    to zero here, and is then counted as a non-rejection."
18. l.2041–2049 `tab:areS`: new Stukel column. l.2053–2054 "it needs fewer than Stukel's test on the
    interaction, the bump and the threshold, and more on the quadratic and the cubic" → "on the bump and the
    threshold, and more on the quadratic, the interaction and the cubic".
19. l.2030–2032 "The rival formulas ... match simulation under the null and show a similar finite-n bias at
    n = 2,000" → for Stukel this is inaccurate in both forms. Suggested: "The Hosmer–Lemeshow formula
    matches simulation; Stukel's matches it under the null and for the quadratic, interaction and bump
    (within 3.9 SE), but not for the cubic (simulated .18, .47, .62 at n = 2,000, 8,000, 32,000 against a
    limit of .70) or the threshold (.23, .14, .12 against .09): its added covariates are quadratic in η,
    and these amplitudes are not yet local."
20. l.1755 and l.1807 (invariance of the Stukel limit in k): hold for SstBoth (its Ω and s(h) are the joint
    form's). The pure-variant ratio is flat for k ≥ 1. No change.
21. l.3564 Study 1 "Tests: Hosmer–Lemeshow, Stukel": add "(the pre-declaration named the joint score form;
    reported here in the SstBoth form of the main benchmark, computed on the same datasets)".
22. l.3572–3590 `tab:studyone`: new Stukel column (§5).
23. l.3596–3600 "breaks Stukel's test and DeepGOF-1 ... Sign errors and mis-keyed outcomes break no test" →
    qualify as in item 9.
24. l.3602–3604 LOO: "between .000 and .027" (see §5 note); "24 of 24" stands.
25. l.3606–3608 "$.47$ to $.24$ for Stukel" → "$.48$ to $.24$".
26. l.3618–3619 "$.059$ for Stukel (from $.142$)" → "$.061$ for Stukel (from $.146$)".
27. l.3623–3626 diagnostic "57–65 % ... at most 16 %" (n = 50) → with SstBoth "57–63 %" (DG-1 57 %, Stukel 63 %)
    and "at most 18 %" under sign/outcome (Stukel outcome 18 %).
28. l.3648–3656 `tab:trimS`: new Stukel rows. Caption "Every cell passes" → "Every cell passes except
    Stukel's test at n = 50 (.089 against a limit of .086)".
29. l.3663–3666 (2): the example numbers are for DG and stand. "24 of 24" stands.
30. l.3668–3673 (3): "$.077/.124$ for Stukel" → "$.097/.129$".
31. l.3674–3677 "Stukel's test is sensitive even to mild covariate errors": stands.
32. l.3681–3686 Verdict: "Stukel's test passes the pre-declared check but stays above the nominal band (.066
    and .073 with the rule, and .124 on the mild slips ...)" → "Stukel's test fails the pre-declared check
    narrowly at n = 50 and stays above the nominal band (.086 and .089 with the rule, and .129 on the mild
    slips the rule misses at n = 50); its clean size is itself .06–.09 at these n."

**Not affected:** Table 2, the niche table and every niche sentence (already SstBoth); every HL, DG-1 and
projection number; the abstract; the cover letter's "Stukel" mentions (REFEREE_FINAL minor 28 is separate).

## 7. Files

`theory/stukel_sstboth/`: `sst_core.py` (limit law + finite-sample statistic), `01_tab_are.py` →
`TAB_ARE_SSTBOTH.json`, `TAB_ARE.log`; `02_sim_check_tabare.py` → `SIM_CHECK_TABARE.json`, `.log`;
`03_irrelevant.py` → `IRR_RAW.jsonl`, `IRR_SSTBOTH.json`, `PRED_VS_SIM_SSTBOTH.json`, `IRRELEVANT.log`;
`fig_are_irrelevant_sstboth.py` → `fig_are_irrelevant.pdf/.png`.
`theory/small_n/sstboth_addon/`: `recompute_sstboth.R`, `rows_s1/` (464), `rows_s2/` (464),
`RUN_SSTBOTH.log`, `analyse_sstboth.py` → `RESULTS_SSTBOTH.json`, `ANALYSE_SSTBOTH.log`.
Stored rows, RESULTS*.json, the pre-declarations and their hashes were only read, never written.
