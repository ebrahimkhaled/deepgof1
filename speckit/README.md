# The frozen pre-declarations — read this first

These thirty files are the **goalposts as they stood before each result was known**. They
ship unedited, including the parts that turned out to be wrong, because that is the only way
a later version can be judged by the same standard this one was.

**Nothing in this folder has been revised after the fact.** Where a pre-declaration was later
contradicted by the evidence, the contradiction is recorded here rather than fixed there.

## Claims made in these files that the paper does NOT make

Two phrases survive in the pre-declarations that the finished paper retired. If you are
checking the paper against these files, these are the differences you will find, and the
paper is the corrected side of each.

**1. "The first pretrained goodness-of-fit test."**
Appears in `26_LIU_BENCHMARK.md` (positioning paragraph). **Retired: it is false.** At least
four earlier tests are pretrained in the same sense, closest among them Simić (2021) for
normality and Kim, Genton, Huser and Castruccio (*Statistics and Computing*, 2025). The paper
claims exact finite-sample validity for a composite fitted-regression null instead, which is
the part those precedents do not have.

**2. "Size .0506, the closest to nominal of 13 arms."**
Appears in `29_JASA_MANUSCRIPT_PLAN.md` (headline plan). **Retired: it was Monte Carlo
noise.** The margin over the next arm was about one fiftieth of its own standard error. Three
arms are tied on the mean level; they differ in *dispersion*, and that is what the paper
claims: a block-to-block standard deviation of .004 against .012 (le-Cessie) and .011 (proj),
a paired difference of roughly four standard errors.

## Pre-declared criteria this version FAILED

Kept here deliberately — a successor that quietly drops them is not being judged the same way.

- **Gate SEL** (`23_SEL_PREDECLARATION.md`) — selection over the classical battery. Failed at
  −.0354. Diagnosed: the selector was trained to predict which arm *fires most extremely*,
  not which has the most matched-level power. The objective is fixable; the opportunity it
  was aimed at is still open.
- **Gate S3** (`21_S3_PREDECLARATION.md`) — deficit against the best published test per cell.
  Passed on p = 2 designs (+.009, +.003, parity) and **failed on p = 3** (−.060, −.082). This
  is the two-axis representation limit, and it is formalized in the paper as Corollary S.6
  rather than left for a referee to find.
- **The near-fibre probe** (`24_NEARFIBER_PROBE.md`) — ruled an *invalid measurement*, not a
  negative result. Its diagnostic AUC had not converged (.509 → .607 between 1,500 and 30,000
  maps), so it measures nothing. Re-run at N ≈ 30,000 before drawing any conclusion from it.

## What a successor must not do

Edit these files. Add `NN_*.md` for a new gate, and record the verdict against the old ones
as they are written here — including the three failures above.
