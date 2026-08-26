<!-- 14_FINANCE_APPLICATION.md - produced by the 4-agent finance workflow (wf_820b3c84-b21),
2026-08-21. VERIFICATION STATUS: the FDIC BankFind endpoints were re-verified BY HAND (148 failures
in FAILYR:2009; 335 GA banks at REPDTE:20081231 - both live, keyless). The EBA Staff Paper N.24
4/2026 grade-size figures and the BCBS WP14 page references are AGENT-verified only: OPUS MUST OPEN
AND CONFIRM EACH before citation. The ECB-Jeffreys nuance in section 1 is the referee-proofing
line - never claim regulators mandate HL. -->

# Finance Application Section — Spec (Executor-Ready)

## 1. Problem Statement (the referee-proof version)

Write these five sentences, with these citations, in this order:

1. Supervisory backtesting of probability-of-default (PD) models is performed at the rating-grade level on small samples: in the EBA's COREP C08.05 backtesting dataset (2,250 rating-grade observations, 45 EU banks, corporate SME, 2017–2024), the 5th percentile of obligors per grade is 34 and the 25th percentile is 260, with 0–2 defaults at those percentiles (EBA Staff Paper Series N.24 – 4/2026, Casellina, Chionsini, Kopp & Riabi, Table 3, p. 11).
2. The model-level calibration test documented by the Basel Committee (BCBS Working Paper No. 14, *Studies on the Validation of Internal Rating Systems*, May 2005, p. 52) and embedded in bank practice and standard software (MATLAB Risk Management Toolbox `risk.validation.hosmerLemeshowTest`, which cites BCBS WP14; the R `PDtoolkit` package; IFRS 9 validation frameworks such as Moldoveanu et al. 2022, which list Hosmer–Lemeshow first among 12-month-PD calibration tests) is the Hosmer–Lemeshow chi-square.
3. Its small-sample unreliability is acknowledged qualitatively by the very sources that document it: BCBS WP14 warns that "if the p_i are very small the rate of convergence to the chi²-distribution may be low" (p. 52); Tasche (2006, arXiv physics/0606071, p. 24) states "there may be an issue of bad approximation for rating grades with small numbers of borrowers"; Nehrebecka (BIS IFC Bulletin, 2016) documents slow chi-square convergence at low PDs on real Polish master scales.
4. No published source quantifies the true size of the Hosmer–Lemeshow family at the sample sizes and event rates where it is actually applied; the low-default-portfolio validation problem is a formally recognized open problem (Pluto & Tasche 2005/2006; BCBS Newsletter No. 6, September 2005; Döhler 2010, arXiv:1006.4968: "default data is usually sparse and sample sizes are often small").
5. We measure that quantity — true size .061–.0845 at n = 50–200 (up to 69% above nominal, worst at 15% event rates) and .1008 (double nominal) on a compound realistic design — and supply a pretrained test that is exactly sized in this regime with higher power (+.090 at n = 100 deployed; +.071, 7 SE, on the compound design), an advantage that shrinks as n grows and vanishes by n ≈ 400, which is theory-consistent and maps the tool precisely onto the lower quartile of real EU rating grades.

**Mandatory wording rule (referee-proofing):** never write "regulators mandate Hosmer–Lemeshow." The ECB's 2019 *Instructions for reporting the validation results of internal models* (Section 2.5.3.1) prescribe the one-sided Jeffreys test at grade level. HL's role is the Basel-documented, industry-standard **model-level (multi-grade)** calibration test. The paper **complements** the Jeffreys test (one-sided PD-underestimation check); it does not replace it (two-sided model-adequacy check). State this distinction explicitly in one paragraph.

## 2. Dataset and Pre-Declared Cells

### Primary: FDIC crisis-era bank-failure early-warning logit

- **Source (free, keyless, verified live):** FDIC BankFind Suite API. Covariates from `https://api.fdic.gov/banks/financials?filters=REPDTE:20081231` (8,395 banks at 2008-12-31; fields ASSET, EQ, NIMY, NPERFV, ROA and hundreds more). Outcomes from `https://api.fdic.gov/banks/failures` (148 failures in 2009 alone). Join key: CERT.
- **Pre-declared cell rule (fixed before any outcome is examined; write it verbatim in the paper):**
  - Cell 1 (primary): all FDIC-insured banks headquartered in Georgia (`PSTALP:GA`) filing a 2008-12-31 call report; outcome = appears in `/failures` with FAILDATE between 2009-01-01 and 2011-12-31. Georgia is externally motivated as the publicly known epicenter of crisis-era community-bank failures — the choice pre-dates and is independent of any test result.
  - Replication cells: every US state with ≥ 15 failures in 2009–2011, same rule, all reported.
  - Model (pre-declared, 5 covariates, CAMELS-proxy standard in the Cole & White 2012 lineage): failure ~ equity/assets + nonperforming/assets + ROA + brokered-deposit share + CRE concentration.
- **Why it survives attack:** the small n is institutional, not manufactured — a state supervisor fitting this exact logit faces ~200–350 banks and 20–60 events. Everything reproduces from two public API calls; ship the download script.

### Backup: UCI Statlog German Credit (n = 1,000, CC BY 4.0; CRAN mirror `scorecard::germancredit`, MIT)

- **Pre-declared cell rule:** segment by loan PURPOSE — the segmentation every scorecard validator actually performs — keep ALL purpose categories with n ≥ 50 (new car ~234, radio/TV ~280, furniture ~181, used car ~103; ~30% bad rate → 30–85 events per cell), report every retained cell. Pre-declared 5-covariate scorecard logit per cell.
- **Role:** face-validity anchor — HL on a scorecard logit is literally what validation manuals prescribe.
- **Action item before writing:** recount the 700/300 good/bad split from the raw file (corroborated but not parsed by the researcher).

### Optional multiplicity engine: UCI Taiwan credit-card default (n = 30,000, 22.12% defaults, CC BY 4.0)

Full pre-declared cross of sex × education(1–3) × age-decade, keep every cell with 50 ≤ n ≤ 300 and ≥ 10 events, report ALL cells. Turns the size distortion into a distribution over dozens of real cells — structurally immune to cherry-picking. Recount the 6,636 events from the raw CSV first. Skip sovereign-default data (no covariates) and give Altman's 66-firm sample (`MixGHD::bankruptcy`, p = 2) at most a one-paragraph cameo.

## 3. What the Analysis Shows — Both Branches

**If the thesis holds:** on each pre-declared cell, (a) the parametric-bootstrap Monte Carlo (null model refit on the real design matrix) shows the classical battery's rejection rate at nominal .05 lands in the .06–.10 band measured in simulation — the real-covariate geometry reproduces the size distortion; (b) the learned test and the classical battery disagree on adequacy calls in a predictable direction (classical rejects more often, consistent with its measured liberality); (c) power at matched level on perturbed-null injections (see §6) favors the learned test. Headline sentence: "on real rating-segment-sized cells, the classical battery's excess rejections are of exactly the magnitude our simulations predict."

**If it does not hold:** report it. Specifically: if the real-design bootstrap sizes come out near nominal (the textbook-design distortion fails to transfer), the paper reports the transfer failure per cell and reframes the contribution as simulation-quantified + conditions-for-distortion; if the tests agree everywhere, report agreement and let the simulation section carry the size claim alone. No cell is dropped, no rule is changed post hoc — the pre-declared rules above are the commitment device, and the all-cells reporting is the proof.

## 4. Venue Call

| Rank | Venue | Why | Cost (Egypt) |
|---|---|---|---|
| 1 | **EJOR** (submit here first) | Q1, IF 6.0, ABDC A*; explicit "Stochastics and Statistics" T&M track; the strongest credit-scoring tradition anywhere (Baesens/Lessmann/Crook/Calabrese school; Lessmann et al. 2015 is the field's standard reference, 216 RePEc citations); the LDP-validation framing defuses its desk-reject rule against pure statistics papers | Zero — no Elsevier–Egypt APC deal exists; subscription route + arXiv green OA |
| 2 | IJF | Calibration of probability forecasts is its core identity (IF ~5.7–8); mandatory code+data since July 2023 rewards the honest-boundary reporting; choose it over EJOR only if the final application is genuinely out-of-time PD forecasting | Zero (same route) |
| 3 | ESWA | IF ~9.4, ML-friendly fallback if both reject on scope | Zero on subscription route (gold APC $3,490 only if chosen) |

Eliminated: Journal of Credit Risk (Q4, IF ~0.27), JBF (empirical-economics culture), Quantitative Finance (IF ~1.9), biostat route (Biometrical Journal requires life-science motivation; the finance application actively disqualifies). Hedge if free gold OA becomes a priority: Annals of Operations Research (Springer — 100% APC covered by the STDF/EKB agreement, verified), borderline on strict JCR Q1.

**EJOR hard constraints for the executor:** ~30 pages all-inclusive (proofs + extra simulations to online supplement), structured abstract 50–250 words, ≤ 5 keywords with the first from the official EJOR list, data availability statement.

## 5. Related-Work Additions (finance block)

1. **Blöchlinger & Leippold (2011), Management Science 57(3), 487–505** — the closest predecessor and MUST be engaged head-on: they also beat HL in credit-default validation, but their angle is dependence-robustness and level/shape decomposition; they do not quantify small-sample size. Position: their edge is dependence, ours is exact small-n validity + power. Benchmark against their test in at least one scenario.
2. **BCBS WP14 (2005)** — regulatory documentation of HL + the convergence warning + Fig. 9 (correlation destroys binomial-test size).
3. **ECB (2019) Instructions, §2.5.3.1** — the Jeffreys grade-level prescription; cite to get the one-sided/two-sided positioning right.
4. **EBA Staff Paper N.24 – 4/2026 (Casellina et al.)** — the grade-size distribution (Table 3) that justifies n = 50–200 line-by-line; also their note that the binomial test "is known to inflate the type I error rate."
5. **Tasche (2006), arXiv physics/0606071** — the standard validation survey; the "bad approximation for rating grades with small numbers of borrowers" quote; the PIT/conditional argument that makes independence defensible in the IFRS 9 setting.
6. **LDP canon:** Pluto & Tasche (2005/2006, Springer ch. DOI 10.1007/3-540-33087-9_5); BCBS Newsletter No. 6 (Sep 2005); Döhler (2010, arXiv:1006.4968).
7. **Nehrebecka (2016, BIS IFC Bulletin)** — central-bank practice: HL "the most popular test of the entire rating system," plus its small-count pathologies on real master scales.
8. **Moldoveanu et al. (2022)** — IFRS 9 practice evidence only (minor venue; cite as evidence of practice, not authority).
9. **Benchmark battery in the application:** HL, binomial, Jeffreys, Spiegelhalter, Blöchlinger–Leippold — matching the tests named across BCBS/ECB/EBA/Tasche.

## 6. Honest Risks and Pre-Committed Answers

1. **"Real data has no ground truth — disagreement between tests proves nothing."** Pre-committed answer, three layers: (a) the size claim never relies on truth — it comes from the parametric bootstrap under the null *on the real design matrix*, which is a statement about the tests, not about the world; (b) power on real data comes from **pre-declared perturbed-null injection**: fit the null logit on a cell, simulate outcomes from a fixed, published perturbation family (omitted quadratic / omitted interaction at fixed amplitudes) over the real covariates, and measure rejection rates — semi-synthetic, fully reproducible, truth known by construction; (c) frame every real-cell disagreement descriptively ("the classical battery rejects k of m cells; the size-corrected test rejects j; the gap matches the measured liberality"), never as "the classical test is wrong about this bank."
2. **Default correlation.** BCBS WP14 and EBA (2026) both show independence-based tests over-reject under asset/serial correlation, and Blöchlinger–Leippold market dependence-robustness. Pre-commit: the test is trained and valid for the conditional/PIT regime (independence given the cycle — Tasche's argument, the IFRS 9 setting); run one correlated-defaults robustness scenario; state openly that heavy asset-correlation regimes are out of scope and a separate, known problem (cite WP14 Fig. 9, EBA 2026).
3. **"You replace a regulator-prescribed test."** No — the prescribed grade-level test (Jeffreys) is one-sided underestimation; ours is two-sided model adequacy at the model/multi-grade level where HL lives. One explicit paragraph, §1 wording rule.
4. **Cherry-picking accusation.** Fully pre-declared cells, externally motivated choices (Georgia = known crisis epicenter; purpose segmentation = standard scorecard practice), all-cells reporting, two-API-call reproducibility, code + data shipped.
5. **"Advantage vanishes at n = 400 — so what?"** Own it: the EBA size distribution shows the lower quartile of real EU rating grades sits at ≤ 260 obligors; the tool targets exactly that quartile and says so. The vanishing boundary is measured, theory-consistent (classical tests are asymptotically near-optimal), and stated in the abstract.
6. **Pretrained-on-simulation suspicion ("learned test overfits its training distribution").** The deployed procedure recalibrates per-dataset via parametric bootstrap, so validity never depends on the training distribution — only power does; the real-data cells are the out-of-distribution check, and branch two of §3 is the pre-committed fallback if it fails.