# 27 — DeepGOF-1: THE VERSIONED RELEASE (spec, 2026-08-23)

Author's decision (2026-08-23): ship the test as a **numbered first model** with the training code,
the training corpus generator, the benchmark harness and every deposited p-value, so that other
researchers can retrain it and publish **DeepGOF-2**. This document is the release plan.

## 1. Why version it — the strategic argument, stated once

The measured position is: correct size on an external benchmark, detects the off-index departure the
partition family and Stukel miss, **5th of 11 on power**. Presented as a finished method that is a
mediocre result. Presented as **DeepGOF-1**, it is a documented starting point with a published
route to improvement — which is how every deep-learning architecture has ever been released.

⭐ **And the harness is the bigger asset.** We can run 11 published tests plus `proj` and Stute-Zhu
across Liu et al.'s four settings, 9 scenarios x 5 severities, and a 60-cell design grid, with every
per-replicate p-value deposited. **That is a reusable evaluation platform for ANY new GOF test.**
Releasing it makes this paper the reference implementation and the leaderboard. It also reframes
"le-Cessie beats us" correctly: we found that, with tooling we are giving away.

## 2. Naming (fixed)

| context | form |
|---|---|
| prose, title, abstract | **DeepGOF-1** |
| R function | `deepgof1()` (with `deepgof()` aliased to the current default version) |
| shipped weights object | `deepgof1_weights` |
| future models | DeepGOF-2, -3 … same interface, new weights |

## 3. ⚠️ Title — the ONE disagreement with the author's draft wording

The author proposed *"DeepGOF1: a trial for making Deep Learned Goodness of Fit test"*.
**"A trial for making" must not survive into the title.** `academic-writing` section 2.6 bans exactly
this register ("A study of", "An investigation into", "An attempt at"), and the Statistics and
Computing desk screen reads title + abstract + conclusion + references only (`springer-stat-comput`
section 4). A tentative title reads as unconfident work before a referee sees a single number.

**The `-1` already carries the humility.** Recommended:

> **DeepGOF-1: a pretrained goodness-of-fit test for logistic regression**

Front-loaded searchable terms, 10 words, no filler. The honest limits go in the abstract's closing
sentence and the results table, where they are EVIDENCE rather than hedging.

## 4. Release contents

```
deepgof1/
  README.md                  <- section 5; the front door
  R/                         deepgof1(), the base-R forward pass, the map builder
  inst/extdata/              deepgof1_weights + scaler  (18,273 numbers, no Python at deploy)
  training/
    01_corpus.R              generate the training corpus (the ONE file a v2 author edits)
    02_train.py              PyTorch training; writes weights + scaler
    03_export_to_R.py        the exporter, with the row/column-major note
    VERIFY.R                 base-R forward pass vs PyTorch, must print < 1e-6
  benchmark/
    liu/                     Liu et al. (2024) settings 1-4, all 13 comparators
    scenarios/               9 scenarios x 5 severities x 3 n
    grid/                    the 60-cell design grid
    run_all.sh               regenerates every number in the paper
  results/                   every per-replicate p-value, as deposited
  speckit/                   the frozen pre-declarations for every gate
```

**Licence:** MIT for code, CC-BY for the deposited results. Archived at Zenodo with a concept DOI
(always resolves to the newest version) plus a version DOI for the DeepGOF-1 release exactly.

## 5. README must answer four questions, in this order

1. **How do I use the test?** Three lines of R, and the ~5 s runtime. Nothing about training.
2. **How was it trained?** The corpus, the architecture, the deployment-matched null class, the one
   command that reproduces the shipped weights.
3. **How do I make DeepGOF-2?** Edit `training/01_corpus.R`, rerun `02_train.py`, run `VERIFY.R`,
   then `benchmark/run_all.sh` — and **compare against the same frozen criteria we did**. The
   pre-declarations are shipped so a successor cannot move the goalposts.
4. **What are the known weaknesses?** Section 6, verbatim, with numbers. A release that hides its
   failure modes wastes the successor's time.

## 6. FUTURE WORK — five directions, each with a measurement behind it

Not speculation. Every one is an inference from a deposited result, and each names the gate.

1. **Representation is the bottleneck, not the learner.** p = 3 designs lose **3-4x** what p = 2 does
   (gate S3: D2 -.146, D3 -.132 vs D4 -.036), purely because the 2-D grid DISCARDS a covariate while
   every classical test uses all of them. A permutation-invariant set network over
   (x_i, pi-hat_i, r_i) removes the projection entirely and is the single highest-value change.
2. **Training-corpus variance dominates design comparisons.** One cell moved **.0219 -> .1786**
   across three corpora (gate EM) while an analytic rival moved .005. Corpus-bagging over M corpora
   shrinks that like sqrt(M) and would make small design effects measurable at all.
3. **Selection over the classical battery is unexploited.** The Cauchy combination loses **~.057** to
   its own best member, measured four times. Gate SEL failed (-.0354) for a DIAGNOSED reason: the
   selector was trained to predict which arm fires most EXTREMELY, not which has most POWER at
   matched level. The gap is still there and the objective is fixable.
4. **Conditioning on the sufficient statistic is an OPEN question, not a closed one.** Gate CF was
   ruled an INVALID measurement: its own baseline failed to reproduce phase 12, and the diagnosis
   showed the AUC depends on classifier training size and had not converged (.509 -> .607 from 1.5k
   to 30k maps). Re-run at N ~ 30,000 before drawing any conclusion.
5. **The departure prior is a first-class object and should be published, tuned and versioned.**
   The same test scores **+.090 under one departure library and -.071 under another** (gate S3 vs
   gate S2). A DeepGOF-2 could improve power by changing ONLY the prior, with no architecture change
   at all.

⛔ **Forbidden in this section:** implying any of these is expected to work. Direction 3 already
failed once and direction 4's probe was invalid — say so.

## 7. What the paper must state about versioning

- The abstract names it DeepGOF-1 and says a first model.
- The results table shows the full ranking including the four tests that beat it.
- A short "Reproducibility and extension" subsection points at the Zenodo DOI and says, plainly,
  that the pre-declared criteria ship with the code so successors are judged the same way.
