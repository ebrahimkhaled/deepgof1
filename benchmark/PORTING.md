# Porting the benchmark scripts

The scripts in `liu/`, `scenarios/`, and `grid/` are shipped **verbatim** — they are the exact
files that produced every deposited p-value in `../results/`, and provenance beats cosmetic
portability. They reference two machine-specific locations that a fresh machine must remap.

## The two constants to change

1. **`harness/config.R`, line `JASA$base <- ...`** — point it at the parent of this
   `harness/` directory (the harness files are copied here: `config.R`, `harness_core.R`,
   `methods_registry.R`, `scenarios_registry.R`, `scenarios.R`, `methods.R`). The scripts
   `source()` these via `JASA$study` / `JASA$base` paths; after remapping, the copies in
   `harness/` resolve.
2. **`harness/config.R`, `JASA$pkg_tree <- ...`** — the harness sources the `ebrahim.gof`
   R package from a git working tree. Clone `github.com/ebrahimkhaled/ebrahim.gof`
   (version ≥ 2.4.1 — earlier versions carry a le-Cessie transpose bug) and point this at
   the clone; or adapt `jasa_init()` to `library(ebrahim.gof)` once ≥ 2.4.1 is on CRAN.

Each top-level script also carries absolute `ROOT`/`PAP` constants in its first lines; set them
to your clone of this release. Everything else is environment-variable driven (`AU_*`, `DG_*`,
`JASA_SMOKE=1` for smoke tests) and documented in each script's header.

## What you can and cannot cheaply reproduce

- `theory/` certificate scripts: **seconds**. `training/VERIFY.R`: seconds.
- One LIU block at n = 50–100: minutes to ~1 h on 12 cores.
- The full deposit (24 LIU blocks incl. n = 500, the 60-cell grid, the scenario lattice):
  **multiple CPU-days** — the n = 500 augmentation alone ran 4.5–8 h per block. The deposited
  per-replicate p-values in `../results/` exist precisely so nobody has to rerun this to check
  a number: every table regenerates from the CSVs alone via the `*_score.py` / `*_verdict.py`
  scripts, which are pure post-processing.

## Determinism

Every dataset in every block is generated from a deterministic seed recorded in the script
(`set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i)` and analogues), so any single
replicate can be re-derived in isolation and compared against its row in the deposited CSV.
