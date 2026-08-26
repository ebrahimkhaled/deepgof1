# 29 — THE JASA MANUSCRIPT PLAN (Theory & Methods)

**Decided by the author 2026-08-25: target = JASA T&M.** This document is the executor-ready plan:
every structural decision is made here so the writing session only writes. Venue mechanics come
from the `jasa-submission` skill (template, `\anon` toggle, agsm.bst, ACC scaffold — all bundled
there and NOT restated here).

**Honest odds, stated once:** JASA T&M is a reach — my estimate 15–25% (vs ~60% at S&C). What
makes it defensible now: BAGofT (Zhang, Ding & Yang 2023, JASA) proves the venue takes
adaptive/learned GOF tests; our theory package is smaller than theirs in count but different in
kind (exact finite-sample validity; a **computable** consistency certificate — theirs has rate
conditions no user can check); and JASA T&M "strongly encourages" a real-data illustration, which
we have (FDIC). The cascade on rejection is §10. The rejection risk is priced in; proceed.

---

## 1. The claim ladder (what the paper claims, in this order, never reordered)

1. **Mechanism** — a goodness-of-fit statistic can be *learned once, offline* and deployed frozen;
   the analyst's own parametric bootstrap supplies the level, so **validity never depends on the
   network being right** (Lemma S.2 idealized; Prop S.3 in general — our p-value is Beran's
   first-stage prepivoted root by Monte Carlo).
2. **Theory** — consistency is *certified, not assumed*: the recession function of the frozen
   network is computable in one zero-bias forward pass (Lemma S.1), power → 1 whenever the
   certificate is positive (Thm S.4), and the blind cone is measured and characterized as
   sub-resolution oscillation (Cor S.5). The representation limit is stated, not hidden (Cor S.6).
3. **Evidence** — exact-size behaviour where classical tests are miscalibrated (58/60 pre-declared
   cells; external benchmark size .0506, the closest to nominal of 13 arms); **beats every
   partition-based test at every n** on the reproduced Liu et al. (2024) benchmark
   (+.13/+.15/+.06 vs HL at n=50/100/500); honest global rank **7th of 13** — a small-sample
   instrument, stated as such.
4. **Reproducibility** — the versioned DeepGOF-1 release: training code, corpora, frozen weights,
   base-R deployment verified to 9.3e-09 against PyTorch, per-replicate p-values for every table.

**Banned claims (standing, from the record):** "first neural/pretrained GOF test" in any wording
(Wilson & Engel 1990 → Sigut 2006 → Simić 2021 → Kim et al. 2025); any power claim not at matched
level; "exact by construction" without the pivotality condition; the pre-refinement Cor-4 wording
("undetectable by ANY map-based test") — use Cor S.6's tightness/no-consistency form only.

## 2. Title — numbered options (author picks)

1. **"DeepGOF-1: A Pretrained Convolutional Goodness-of-Fit Test for Logistic Regression with
   Certified Validity"** — name+colon pattern (+49% CNCI in the family's title study), versioned
   identity the author wants, property-led. **Recommended.**
2. "A Goodness-of-Fit Test That Is Trained Once: Exact Calibration for a Learned Statistic" —
   mechanism-led, no version name.
3. "Pretrained Goodness-of-Fit Testing for Logistic Regression: Validity Without Trusting the
   Network" — contrast-led.

Keywords (3–6, none in the title): *model misspecification; parametric bootstrap; Monte Carlo
test; residual diagnostics; consistency; calibration*. (Adjust to avoid overlap with the chosen
title.)

## 3. Abstract skeleton (≤200 words, unstructured, no citations)

- Sentences 1–2 (field problem): classical GOF tests for logistic regression are miscalibrated
  exactly where they are most needed — small samples, rare events.
- Sentences 3–5 (contribution): a convolutional statistic trained once on simulated departures,
  reading misfit as a picture over covariate ranks; calibration by the analyst's own bootstrap, so
  the level does not depend on the network; the prepivoting identity.
- Sentences 6–7 (theory): consistency certified by a computable functional of the frozen weights;
  the blind cone measured.
- Sentences 8–9 (evidence, with numbers): size within the nominal band in 58/60 pre-declared
  cells; dominates every partition-based test at every sample size on a published benchmark;
  honest overall placement + small-sample framing.
- Sentence 10: versioned release enabling successor tests (DeepGOF-2).

## 4. Section architecture with page budget (≤35 pp double-spaced; target 31–33)

| § | content | pp | figures/tables in-body |
|---|---|---|---|
| 1 | Introduction: CARS — miscalibration territory → learned-test lineage gap (level never controllable: Simić's own words) → this paper occupies it. Contributions enumerated (4). | 3.5 | — |
| 2 | The test: residual map, axis rule, network, deployment-matched training, bootstrap rank p-value. Practitioner code box (3 lines of R). | 4 | **Fig 1** (architecture, two-lane), **Fig 2** (map gallery = fig3_maps) |
| 3 | Theory: A1–A6 compressed; Lemma S.1–S.2, Prop S.3, Thm S.4, Cor S.5–S.6 **statements only** + the certificate table (margins 1.08–2.20 vs blind cone −0.28). Proofs → supplement. | 4 | **Table 1** (certificate margins) |
| 4 | Simulation I — size: the pre-declared grid; 58/60; the CCT/Osius size failures shown, not hidden. | 2.5 | **Fig 3** (= fig6_s3grid right panel or both) |
| 5 | Simulation II — the published benchmark: Liu et al. settings verbatim; beats all partition tests at every n; proj strongest overall (their remedy works); le-Cessie +.954 completeness observation (neutral wording per `academic-writing` §3.2 — Härdle/Liang could referee); rank 7/13 stated. | 4.5 | **Fig 4** (= fig8_liu), **Table 2** (power by family × n) |
| 6 | Real data: FDIC bank-failure application (phase 22 deposit). | 2 | **Table 3** (FDIC results) |
| 7 | DeepGOF-1 as a versioned object: release contents, base-R deployment (9.3e-09), how DeepGOF-2 improves it (retraining levers). | 1.5 | — |
| 8 | Discussion: small-sample instrument; representation limit; what the certificate buys practitioners; limitations ledger. | 2 | — |
| — | References (~45 entries) | 3 | — |

Figure budget: 4 figures + 3 tables. fig5_probit, fig7_scenarios, fig1_breakdown, fig2, fig4
move to the **supplement** (they support, not carry). Optional new mini-figure: the blind-cone
checkerboard vs a natural direction (cheap, striking; decide at draft time).

## 5. The supplement (separate anonymized file, per skill)

S1 = theory_supplement.tex content (post-workflow fixes), renumbered to the skill's skeleton.
S2 = supplementary simulations (scenario lattice fig7, probit/cloglog fig5, S3 full grid).
S3 = training protocol details + corpus recipe (the manifest's "recipe" section, prose form).
S4 = reproducibility pointers. Keep the supplement anonymous-buildable too.

## 6. Double-anonymity traps (this family's specific ones)

- The six excluded own-family arms (EF/DEF/A.*): anonymous build says *"six additional tests under
  development by the authors, excluded from all comparisons by pre-declaration; reference blinded
  for review."* Never name EF/EDGE in the anonymous PDF.
- The EDGE arXiv preprint exists and is citable — cite in the with-authors build only.
- Release repo name (`ebrahim.gof`) → anonymous build refers to "the released package (name
  blinded)". ACC materials anonymized "to the extent possible" at revision stage.
- Acknowledgments, funding, ORCID: with-authors build only (the `\anon` toggle handles this).
- FDIC data: public API, no identity leak. Kasem/Alexandria affiliations: strip from prose.

## 7. Referee anticipation (write these paragraphs before they ask)

| likely referee angle | pre-emption in draft |
|---|---|
| BAGofT authors: "we did adaptive GOF with 4 theorems" | §3 contrast: their validity asymptotic w/ rate conditions; ours exact-by-rank + computable certificate; cite generously |
| "power is only 7th of 13" | §5 leads with the family win + small-n margins; §8 the small-sample-instrument sentence; matched-level methodology emphasized (Simić's numbers weren't) |
| "consistency is conditional" | the condition is *checkable* — that is the feature, not the bug; blind cone measured at <7.4e-5 |
| "CNN is a black box" | Lemma S.1: at signal strength the network IS its zero-bias linearization — the certificate is interpretable |
| "why K=6, why these axes" | gate S4 pre-declared resolution sweep; axis rule stated + Cor S.6 covers its cost |
| Härdle/Liang (if crossed to them): the le-Cessie observation | neutral completeness framing; their proj confirmed strongest — we validate their remedy |

## 8. Writing order (per `academic-writing` §6.2 — backwards)

1. Tables + figure captions (each caption a declarative claim).
2. §5 → §4 → §6 (results prose from the deposits; every number traced to a CSV).
3. §3 (statements from the verified supplement).
4. §2 (mechanism).
5. §8 discussion → §1 introduction → abstract → title.
6. Passes: claim-consistency (one term per concept: "map", "certificate", "level"), hedge
   calibration, anonymity grep (`Ebrahim|EDGE|EF test|ebrahim.gof|Alexandria|Kasem`), then the
   full `academic-writing` §8 revision ladder.

## 9. Gates before submission (the "100%" list)

- [ ] **G1** Adversarial workflow findings fixed in `theory_supplement.tex` (REAL ones), recompiled
- [ ] **G2** DeepGOF-1 release tree built per speckit/27 (also = the ACC zip at revision stage)
- [ ] **G3** Manuscript drafted per §4, both `\anon` builds compile, ≤35 pp verified
- [ ] **G4** Every number in the paper regenerated from a deposited CSV (no orphan numbers)
- [ ] **G5** Anonymity grep clean on the anonymous PDF + supplement
- [ ] **G6** Multi-agent editor pass (separate from the proof workflow) on the full draft
- [ ] **G7** Cover letter: significance ¶, no-proceedings declaration, gen-AI disclosure per policy
- [ ] **G8** Similarity check (family siblings overlap risk: EDGE/LEGofT prose) — reword, don't cite-pad
- [ ] Author-side: CRediT roles, ORCID, corresponding author fixed

## 10. Cascade (pre-decided, so a rejection costs days not weeks)

JASA reject → **Statistics and Computing** (`springer-stat-comput` skill): sn-jnl `iicol`,
author-year, proofs return to an Online Resource, abstract 150–250, cover letter names Liu et al.
and Kim et al. as the journal's own precedents. The §4–§7 content is venue-independent by design;
only the shell and the abstract re-target. Estimated conversion: 1–2 sessions.

## 11. Timeline estimate

| step | effort |
|---|---|
| G1 proof fixes | 0.5–1 session (workflow-dependent) |
| G2 release tree | 1 session |
| Draft §§ per §8 | 3–4 sessions |
| Editor pass + gates G4–G8 | 1–2 sessions |
| **Total to submission-ready** | **~6–8 sessions** |
