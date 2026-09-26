# calibration_gap -- progress log

One line per completed step, newest at the bottom.

- 2026-09-16 06:20 -- Read the five scripts. Smoke-test outputs (M=2000, D=4: `per_dataset.csv`,
  `cell_summary.csv`, `run_meta.rds`, `scores_n*_rep1.csv.gz`) moved to `smoke_M2000/` so they are
  not mixed with the real run. `verify.rds` kept in place (it is the verification, not a run output).
- 2026-09-16 06:22 -- `02_pilot.R` rewritten as CHECKPOINTED and RESUMABLE. Protocol unchanged:
  design D1, seeds `901000+i` for the dataset and `SEED0+i+3e6/6e6+1e7*cellcode` for the two laws,
  tau = .30, gamma = .01, alpha = .05, the same separation rule (Firth refit when any |beta_j| > 8
  or IRLS did not converge; an outright failure scores +Inf), the same functionals. What changed is
  only bookkeeping: (a) the moment a dataset finishes its row is written atomically to
  `rows/n<N>_rep<i>.csv` (temp file + rename) and appended to `per_dataset.csv` under a lock
  directory (`.pdlock`, `dir.create` is atomic on NTFS), header written if the file is absent;
  (b) on start, every (n, rep, M) already present in `rows/` or in `per_dataset.csv` is skipped;
  (c) `per_dataset.csv`, `cell_summary.csv` and `run_meta.rds` are rebuilt from the row files at the
  end of every invocation; (d) the `scores_n<N>_rep1.csv.gz` dump is written by the worker that does
  rep 1, only if the file does not already exist; (e) new env vars `CG_CHUNK` (network batch size,
  default 2500 -- halves peak worker memory versus the old fixed 5000) and `CG_OUTDIR` (output
  directory, for dry runs); (f) a `progress.log` heartbeat every 4 chunks, because on a Windows PSOCK
  cluster the workers' output never reaches the master's stdout.
- 2026-09-16 06:23 -- `04_report.R` patched to match: it still prefers `rows/`, but falls back to
  `per_dataset.csv`, filters to the run's M, de-duplicates on (n, rep); the hard-coded "M = 100,000
  / D = 40 / 12 workers" narrative in section 3 and the "protocol D = 40" line now come from
  `run_meta.rds`, and the missing `mt$NW` that would have errored is guarded.
- 2026-09-16 06:26 -- TIMING PROBE under the real current load (48 R processes from the two other
  jobs, 24 logical cores, CPU pinned at 100%): `CG_M=4000 CG_D=8 CG_W=8 CG_CHUNK=2500` into a
  scratch directory. One wave of 8 datasets took **1.08 min at n = 50** and **1.33 min at n = 200**
  (8000 scores per dataset, i.e. 8.1 and 10.0 ms per score -- 1.3 to 1.8x slower per worker than the
  4-worker smoke test, as expected with 8 workers added to a saturated machine).
- 2026-09-16 06:29 -- RESUME VERIFIED: one row file deleted and its line removed from
  `per_dataset.csv`; the rerun reported "resume: 15 dataset rows already on disk", recomputed that
  one dataset only, left the n = 200 cell untouched ("already complete"), and rebuilt
  `per_dataset.csv` back to 16 rows.
- 2026-09-16 06:30 -- **RUN SIZE FIXED BEFORE LAUNCH: M = 100,000, D = 20, CG_W = 8, CG_CHUNK = 2500,
  cells n in {50, 200}.** Arithmetic, straight from the probe. D = 20 on 8 workers is 3 waves per
  cell (8 + 8 + 4).
  * At the target **M = 200,000** (x50 the probe): 54.0 min per wave at n = 50 and 66.7 min at
    n = 200, so 3 x (54.0 + 66.7) = **362 min = 6.0 hours** -- this EXCEEDS the ~5 hour budget.
  * The pre-declared fallback is therefore used: **M = 100,000** (x25 the probe): 27.0 min per wave
    at n = 50 and 33.3 min at n = 200, so 3 x (27.0 + 33.3) = **181 min = 3.0 hours**.
  * DKW term at M = 100,000: 2*sqrt(log(4/0.01)/(2M)) = **0.01095** (the .011 of the protocol), against
    .0077 that M = 200,000 would have bought.
  * D is NOT lowered; it stays at the protocol's 20 per cell, 40 datasets in all.
- 2026-09-16 06:29:33 -- REAL RUN LAUNCHED DETACHED, master pid 48572 (see `RUN_INSTRUCTIONS.md` for
  the exact command and for how to resume it). All 8 workers were computing by 06:30:00.
- 2026-09-16 06:39 -- Throughput on the live run, measured over the first 8 minutes from
  `progress.log`: **1670 scores/s aggregate** across the 8 workers at n = 50 (per-worker spread 4.2
  to 8.6 ms per score). That is about twice the rate the probe implied, because a probe wave lasting
  only one minute is dominated by cluster start-up and by each worker loading `brglm2`, the weights
  and the scaler. On this better number the n = 50 cell is ~50 min and the n = 200 cell ~65 min, so
  **the run should finish about 08:30**, not 09:30.
- 2026-09-16 06:40 -- Reconsidered and **kept M = 100,000**. On the corrected throughput the target
  M = 200,000 would have been about 3.8 h, i.e. inside the 5 h budget; but (a) the protocol requires
  M and D to be fixed *before* launch and they were, from the best measurement available at that
  moment, and (b) conditions are getting worse, not better -- free memory is down to 1.4 GB of 32 GB,
  there are now 50 R worker processes on 24 cores, and a THIRD unrelated job (`analyse_selftest.R`)
  started at 06:15, so doubling my own footprint's residency time would put the other jobs at risk
  for a 0.0033 improvement in the DKW term. Recorded here so the author knows the option exists: a
  rerun at M = 200,000 (DKW .0077) is feasible in under 4 h on a quieter machine, and would require a
  fresh `rows/` directory because rows are matched on (n, rep, M).
- 2026-09-16 06:48 -- FIRST DATASETS LANDED. The first wave of 8 finished between 06:42 and 06:46,
  13.0 to 16.4 min per dataset at n = 50; `rows/` and `per_dataset.csv` both gained the rows the
  moment each finished, with six near-simultaneous appends and no corruption of the file. Expected
  finish now ~08:15-08:30 (n = 50 cell ~45 min, n = 200 cell ~60 min, 3 waves each).
- 2026-09-16 06:50 -- Resume-key matching re-verified against the LIVE checkpoint files: R stores M
  in the csv as `1e+05`, and `sprintf("%d", ...)` recovers `100000` from it, so a finished dataset
  (50, 7) matches the master's key and an unfinished one (50, 2) does not. A resume will therefore
  skip exactly what is done.
- 2026-09-16 06:50 -- `RUN_INSTRUCTIONS.md` written: resume command, how to tell whether the run is
  alive, the expected finish, the `04_report.R` command, and the warning not to kill the other jobs'
  R processes.
