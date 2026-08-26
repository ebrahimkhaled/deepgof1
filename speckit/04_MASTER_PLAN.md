# 04 — MASTER PLAN: how this research gets done well, and where it lands

**Written 2026-08-17**, after gate C1 passed, gate C3 cleared the literature, and a 10-agent
adversarial review withdrew part of C1's interpretation (`03_C1_DEFECTS.md`).

This document is the executor's plan. It supersedes `00_VISION.md` wherever they disagree. Its
job is to remove every decision so the work can be run without re-litigating anything.

---

## §0. The one-paragraph state of the project

A goodness-of-fit test for logistic regression, **pretrained once offline** on 100 randomly
generated misspecification families. At use time it reads 33 numbers from the user's fitted model
(11 classical GOF z-scores + risk-ordered grouped residuals + the same residuals sorted + hl +
max|·|), passes them through a frozen 20-network ensemble, and turns the score into a p-value by
Monte Carlo rank. Under a pre-declared confirmatory protocol it beat the pre-declared comparator
at all three sample sizes. **But** the comparison was not size-matched, nothing per-replicate was
stored, and — decisively — **every size number was computed against a null reference drawn from
the true model with known β, which no analyst has.** The programme's next act is to measure the
test as a user would actually run it. Everything else waits on that.

## §1. The thesis, in the form it will be defended

> Existing learned goodness-of-fit tests spend the analyst's sample to learn. This one spent a
> simulator's.

Three claims follow, each falsifiable, each with the experiment that decides it:

| # | Claim | Decided by |
|---|---|---|
| T1 | Amortized omnibus GOF of a *fitted parametric regression* under a randomized departure prior is new | §2 gate G-LIT (C3 + the two missed literatures) |
| T2 | The user pays **no sample** for the learning, so the advantage is largest where a sample split hurts most (small n) | C-BOOT at n = 200, then the n-ladder |
| T3 | Validity does not rest on the learner | **C-BOOT** — and if it fails, T3 is deleted and the paper changes |

**Never claimed** (each already exists): learned GOF testing; Monte Carlo calibration
(Barnard 1963; NPLM); neural networks as GOF alternatives; amortized inference in general.

## §2. The gate ladder — ordered, with pre-declared criteria AND fallbacks

Each gate is declared *before* it runs. A gate that fails is reported as it fell; the fallback is
already written down so no result gets renegotiated after the fact. This is the discipline that
cost the LEGofT paper a 78% shrinkage when it was skipped once.

### G-BOOT — deployed calibration ⟵ RUNNING, and it gates everything

**Question.** Does the test hold its level under per-dataset parametric bootstrap from β̂?
**Design.** `code/c_boot_deployed.R`. n = 200, B = 999, 1200 null replicates (600 two-covariate,
600 four-covariate), 200 replicates on each of the 12 families. Exact rank p-value
p = (1 + #{S_b ≥ S_obs})/(B + 1). **Both rules are deployed identically from the same B refits**,
so the power comparison is size-matched by construction and defect D1 cannot recur.
**Criteria.** (a) realized size ∈ [.035, .065] with the Clopper–Pearson interval covering .05;
(b) the null p-value distribution passes a KS test against U(0,1) at α = .01;
(c) the deployed mean-12 gain over the comparator is positive.
**Fallbacks, declared now.**
- Size inflated (> .065, CI excluding .05) → **that becomes the paper's lead finding.** The paper
  is re-framed as *"pretrained GOF statistics are powerful but bootstrap calibration is the
  binding constraint"*, the abstract's validity sentence is deleted, and the target is JCGS.
- Size fine, gain gone → the paper is a **negative result about deployment** and still publishable
  at JCGS/CSDA, but the amortization pitch is retired.
- Size fine, gain intact → proceed, and this becomes Figure 3, the paper's survival figure.

**Cost note (measured, not estimated).** One bootstrap step (refit + full battery) costs 111 ms at
n = 200, 656 ms at n = 500, 4.80 s at n = 1000. B = 999 therefore costs 111 s / 11 min / 80 min per
replicate. **The deployed study is only affordable at n = 200 at full B** — which is exactly the
size the paper should lead with. At n = 500 use B = 299 and 600 reps; at n = 1000, B = 199 and 300
reps, reported as supporting evidence with their wider intervals shown.

### G-RERUN — the instrumented confirmatory run

Re-run C1 with three fixes: both rules calibrated on the same held-out cell; **every per-replicate
score, the 20 fitted nets, and the scaler saved**; a four-covariate null cell added. ~20 h wall on
12 workers. Nothing new is learned — but every later referee question becomes a re-score instead of
a re-run, and the standing project rule on depositing per-replicate p-values is honoured.
**Criterion:** the size-matched gain remains ≥ +.03 at n = 200. **Fallback:** if it does not, the
headline is whatever the size-matched number is; no re-runs with new seeds.

### G-TRANSFER (C8) — does the advantage survive outside the training design?

**This is the gate that decides the paper's ceiling**, and the one a referee will reach for.
**Design.** p = 3, 5, 10 covariates; correlated X (AR(1), ρ = 0.5); one categorical covariate;
event rates ~5% and ~15%; plus the two empirical designs (Pima, LBW). All through the deployed
calibration. **Criterion:** positive mean gain in at least 5 of the 7 new design classes, with
size in range in all 7. **Fallback, in writing before the run:** if the gain lives only inside the
training design class, the amortization framing is retired, the contribution becomes *a recipe for
building a design-specific pretrained test*, and the target is JCGS with no further debate.
⚠️ Note the tension: the leakage autopsy proved these features **can** read design fingerprints
(`lev_sd` |coef| 462; +.509 on 4-cov vs +.005 on 2-cov). Broadening the design prior may dilute
exactly the signal that is winning. Expect it; do not be surprised by it.

### G-THEORY — two propositions, and only two

- **Proposition 1** (days, no compute). The randomized library is a prior Π over departures; the
  ensemble estimates dP_Π/dP₀ on the 33-feature map; therefore the calibrated rule approximates
  the **Π-average-power-optimal test** (Wald 1943; King 1988; Andrews & Ploberger 1994). One page.
  Its value is rhetorical and large: it converts every honest loss (cloglog −.141,
  multivar_int −.091) from an apology into a **corollary** — a Π-optimal test is *supposed* to
  trade power on directions the prior weights lightly.
- **Proposition 2** (1–2 weeks). Pivotality / transportability: conditions under which the
  bootstrap null distribution of the score is asymptotically free of β. This one object closes
  the calibration hole *and* the design-transfer hole.
- **Explicitly skipped:** any excess-risk-to-power regret bound. At these sample sizes it is
  vacuous and no referee will verify it. Say so in one sentence rather than pretending.

### G-LIT — discharge the literature debts

Fetch and read in full the six sources still marked **S** in `02_C3_LITERATURE.md`, plus the two
literatures that sweep missed: the **specification-testing tradition** (Ledwina 1994; Kallenberg &
Ledwina 1997; Aerts–Claeskens–Hart 1999; Bierens & Ploberger 1997; Horowitz & Spokoiny 2001;
Janssen 2000; Kim–Ramdas–Singh–Wasserman 2021) and **LF2I**, which already amortizes frequentist
testing *including the critical value*. LF2I is the single most dangerous omission — read it first
and position against it explicitly.

### G-SOFTWARE — `deepgof()` in `ebrahim.gof`

Frozen ensemble weights shipped with the package; the deployed bootstrap as the default p-value
route; a `reference =` argument for users who *can* simulate their design. Ship before submission —
software appearing between submission and first report reads as strength.

## §3. Venue: the decision, and the rule that decides it

**First target: JCGS.** A computational construction with an honest deployment study is JCGS house
style, not a deficiency. Realistic acceptance ≈ 40%.

| Rank | Venue | Why | When |
|---|---|---|---|
| 1 | **JCGS** | contribution is a computational construction; deployment studies are welcome | after G-BOOT, G-RERUN, G-TRANSFER, G-SOFTWARE — ~3–4 months |
| 2 | **JASA T&M** | the reach | only if G-BOOT **and** G-TRANSFER both come back clean, **plus** both propositions — ~6–9 months, ≈20% |
| 3 | **CSDA** | simulation-friendly, Q1-adjacent, forgiving on theory | if JCGS bounces |
| 4 | **Statistics and Computing** | genuinely fits *this* paper — but they desk-rejected EDGE on fit, so ranked below CSDA on history | fallback |
| — | ~~Biometrika~~, ~~JRSS-B~~ | higher theory bar than JASA for an empirically-led paper. **Off the ladder.** | never |
| — | ~~Statistics in Medicine~~ | LEGofT owns that lane | never |

**The decision rule, fixed now so it is not re-argued under optimism:**

> Make the JASA-vs-JCGS call **at a single checkpoint, after G-BOOT and G-TRANSFER report**.
> If both are clean and both propositions are written → JASA T&M. If either scopes the claim →
> JCGS, immediately, without further debate.

JASA today is 2–3%: T&M expects theoretical content, and a simulation-led paper is desk-rejected on
remit, not on quality. That is a **fit** problem — and this family has already lost months to two
fit-rejections (Statistics & Computing, then ADAC). Do not repeat it out of ambition.

## §4. The manuscript

**Lead with n = 200, not n = 1000.** It is the only sample size whose win survives a sign test
(10–2, p = .019); it keeps a large size-matched gain; it is where the no-sample-splitting mechanism
argument is strongest (a rival that splits the sample pays most there); and it reproduces in 6.4
minutes for a reproducibility reviewer. n = 1000 is a ceiling test in which half the families pin
at 1.000.

Sections: Introduction (CARS) → The pretrained test (the departure library; the 33-feature map; the
ensemble; deployed calibration) → Theory (Props 1–2) → Simulation (**size before power**, always) →
Design transfer → Two applications → Discussion with limitations in our own words.

**Standing writing rules for this paper.** Always "the pre-declared comparator", never "the
strongest classical rule" — at n = 200 gamma-tilt (.3448) and ridge (.3378) both exceed CCT-full
(.3075). Report at least one setting where the method is not best (there are three; keep them
visible). State the honest two-component story for generalization: generic residual-curve reading
(survives kind-holdout: step keeps 88%) **plus** prior coverage. And never write "one forward pass"
— under bootstrap calibration the use-time cost is B refits of the whole battery.

## §5. The figure program

JASA and JCGS both impose **no graphics conservatism** — colour is free online, print colour is
$400/figure, figures sit in-body at ≥300 dpi. The binding constraint is the page budget, not taste.

**The rule that separates distinctive from amateur:** a figure is *one coordinate system, or a
lattice of identical ones*. Bold **marks** are always safe — arrows, dumbbells, 100 small
multiples, shaded acceptance regions, dot strips of 20 initializations. Bold **page furniture** is
what reads as a conference poster: tinted banner bands, numbered stations, spine arrows linking
panels, editorial sentences typeset inside the plotting region. *Test: cover the axes and tick
labels; if what remains still looks like a slide, it is a poster.* One new visual grammar per paper
is affordable; two is not. And every figure should contain **the object that would have falsified
its own claim** — the pre-declared threshold, the size window, the ceiling. For a paper whose whole
method is pre-declaration, drawing the acceptance region is the most persuasive thing on the page.

| # | Figure | Carries | Status |
|---|---|---|---|
| 1 | **The departure library** — 10×10 lattice of the 100 training families, ordered by how much of each departure is visible in the risk ordering | the prior is arbitrary shapes, not rival models | data exists |
| 2 | **The residual-curve gallery** — 4×4 lattice, 50 bins, three single replicates in grey + the null band + the 1000-replicate mean | *what the net actually reads*; **nobody in the GOF literature draws this** | curves exist at n=1000 only |
| 3 | **Size under deployment** — every pre-declared condition, deployed vs idealized, Clopper–Pearson intervals, plus the null p-value ECDF against the diagonal | the paper's survival claim | ⟵ G-BOOT produces it |
| 4 | **The confirmatory ledger** — per-family arrows, both rules size-matched, the pre-declared window drawn | the headline, honestly | needs G-RERUN numbers |
| 5 | **Design transfer** — gain against distance from the training prior, crossing point marked | pre-empts the question that would otherwise sink the paper | ⟵ G-TRANSFER |
| 6 | **The ladder** — every ingredient alone failed, with the measured ceilings and the CNN loss | the ablation, kept on the page with its failures | data exists |

**Kill `fig_money.png` as it stands** — two lines and six bars, every value of which sits in a CSV.
A referee writes "Figure 1 is a table." **Do not rebuild any figure on the current numbers**: the
instrumented re-run moves every headline value, and three of the four proposed designs hard-code
the old ones.

## §5b. Compute: the battery was 4.4x slower than it needed to be (fixed 2026-08-19)

Profiling the deployed run found **98% of every bootstrap step is the GOF battery** (142.9 ms of
145.2 ms; the glm refit is 2.3 ms). Two causes, neither arithmetic:

1. `run.all.gof(include_slow = FALSE)` computes **21 tests**; the feature vector uses **9**.
   Twelve results were computed and discarded on each of ~3.6 million evaluations.
2. `eval_cheap()` also calls `run.all.gof(tests = "le-Cessie", include_slow = TRUE)` at 22.8 ms —
   and **le-Cessie is not in the 33-feature vector**. It is the LEGofT paper's test, not this one's.

`code/fast_battery.R` restricts the battery to the 9 tests actually used. **Verified on 80 datasets
across both designs and two sample sizes: 0 mismatches, max |difference| = 0.000e+00** — the same
functions on the same inputs, so this is not an approximation. **155.4 ms → 33.5 ms per feature row;
the full bootstrap step 157.7 ms → 35.8 ms (4.4x).** It is additive: the shared `eval_cheap()` that
LEGofT and the benchmark paper depend on is untouched.

**Consequence for the plan.** G-RERUN and G-TRANSFER cost about a quarter of what §2 assumed —
n = 500 at B = 999 becomes affordable, and G-TRANSFER's seven design classes go from prohibitive to
routine. Re-time before scheduling. Always run `verify_fast_battery()` after any package update.

**A separate, exact speedup for le-Cessie itself** (`code/lecessie_algebraic.R`, not needed here but
valuable to `ebrahim.gof`): the package forms the N×N hat matrix and computes
`sum(diag(R.tmp %*% R.tmp))` as a full N³ product. Two identities remove that — H is rank p, so
(I−H)′R(I−H) expands to O(N²p); and since R.cor is symmetric, tr(R.tmp²) = μ′(R.cor∘R.cor)μ in
O(N²). **4.4x at n = 200, max relative difference 2.9e-15 over 75 datasets**, and the gap widens
with n because it removes an O(N³) term. Worth porting into the package.

## §6. Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| Deployed size inflated | **fatal to the current framing** | G-BOOT running; fallback declared in §2 |
| Advantage confined to the training design | **fatal to amortization** | G-TRANSFER with a written fallback |
| Desk-reject on remit at JASA | high | JCGS first; checkpoint rule in §3 |
| A rival ships first (Chen & Lei is 4 months old) | medium | re-check for released code near submission; the separating axis is structural, not incremental |
| Self-similarity flag at the desk screen (3 sibling manuscripts) | medium | declare all three in the cover letter — self-declaring turns a flag into routine disclosure |
| Over-claiming, again | **the historical failure mode** | the claims ledger in §1; "size before power"; keep the losses visible |

## §7. Kill criteria — when to stop

Stop and write this up as a methods-and-failures paper if **both**: G-BOOT shows deployed size
outside [.03, .08] *and* a corrected calibration removes the mean gain at n = 200. In that case the
honest contribution is the failure catalogue (§ manifest phases 0–3) plus the recipe plus the
measured ceilings — real, publishable, and considerably more useful to the field than a fourth
paper claiming a small power win.
