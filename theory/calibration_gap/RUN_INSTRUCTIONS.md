# Calibration-gap pilot -- running, resuming, reporting

Folder: `PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/calibration_gap`
R: `"/c/Program Files/R/R-4.4.1/bin/Rscript.exe"` (R 4.4.1)

## The run that is in flight

Launched 2026-09-16 06:29:33, detached, from Git Bash in that folder:

```bash
cd "/path/to/project/paper_deepgof/theory/calibration_gap"
(CG_M=100000 CG_D=20 CG_W=8 CG_N=50,200 CG_CHUNK=2500 \
  nohup "/c/Program Files/R/R-4.4.1/bin/Rscript.exe" 02_pilot.R > RUN.log 2>&1 &)
```

| | |
|---|---|
| M (scores per law per dataset) | 100,000 |
| D (datasets per cell) | 20 |
| cells | n = 50, n = 200 |
| workers | 8 (the machine is also running ~48 R processes from two other jobs) |
| network batch | 2500 maps |
| DKW at gamma = .01 | 0.01095 |
| master process | `Rscript.exe` pid 48572 (children: 8 `slaveRSOCK` workers) |
| **expected finish** | **about 08:30 on 2026-09-16**, and in any case before 09:30 |

Two finish estimates, because the second is better. The pre-launch one, from a probe of 8-worker
waves at M = 4,000, was 3.0 h (finish ~09:30). Once the real run was up, its own `progress.log` gave
**1670 scores/s aggregate** at n = 50, roughly twice what the probe implied -- a one-minute probe wave
is dominated by cluster start-up and by each worker loading `brglm2`, the weights and the scaler. At
that rate the n = 50 cell is ~50 min and the n = 200 cell ~65 min, so **finish ~08:30**.

Why M = 100,000 and not the target 200,000: at launch, the only measurement available was the probe,
on which D = 20 on 8 workers is 3 waves per cell and M = 200,000 came to 3 x (54.0 + 66.7) = 362 min
= 6.0 h, over the ~5 h budget, so the pre-declared fallback M = 100,000 (3.0 h) was fixed and
recorded before launch. The corrected throughput would have allowed M = 200,000 in about 3.8 h; it
was not adopted, because M and D are required to be fixed before launch and because free memory on
this machine had fallen to 1.4 GB of 32 GB with a third unrelated job starting at 06:15. **D was not
lowered** -- it stays at the protocol's 20 per cell.

## Is it alive, and how far along?

```bash
cd "/path/to/project/paper_deepgof/theory/calibration_gap"
tail -5 RUN.log                  # master: header, then one line per finished cell
tail -20 progress.log            # workers: START / heartbeat every 10,000 scores / DONE, per dataset
ls rows | wc -l                  # finished datasets (40 when the run is complete)
wc -l per_dataset.csv            # header + one line per finished dataset
```

`RUN.log` only grows when a whole cell finishes -- on a Windows PSOCK cluster the workers' `cat()`
never reaches the master's stdout. `progress.log` is the live one; it gains a line roughly every
1.5-2.5 min per worker.

## If the run dies -- resume

Re-run the **identical** command. Every dataset whose (n, rep, M) row is already in `rows/` or in
`per_dataset.csv` is skipped, so only what is missing is recomputed; nothing is recomputed twice and
no completed work is lost.

```bash
cd "/path/to/project/paper_deepgof/theory/calibration_gap"
(CG_M=100000 CG_D=20 CG_W=8 CG_N=50,200 CG_CHUNK=2500 \
  nohup "/c/Program Files/R/R-4.4.1/bin/Rscript.exe" 02_pilot.R >> RUN.log 2>&1 &)
```

The first lines it prints tell you what it decided:

```
resume: 27 dataset rows already on disk
  cell n= 50: 0 of 20 datasets to compute
  cell n=200: 13 of 20 datasets to compute
```

Two things to know before resuming:

* If a stale lock directory `.pdlock` is left behind by a hard kill, delete it first
  (`rm -rf .pdlock`); otherwise the first append waits two minutes and then proceeds anyway.
* Do not change `CG_M` on a resume. Rows are matched on (n, rep, **M**), so a different M makes
  every dataset look unfinished and the final rebuild keeps only rows at the new M.

To stop the run: kill the master `Rscript.exe` and its 8 children. **Do not kill the other R
processes on this machine** -- `run_M_rivals.R` and `T30_keycells_NB499.R` are two unrelated jobs.

## When it finishes -- the report

`02_pilot.R` already rebuilds `per_dataset.csv`, `cell_summary.csv` and `run_meta.rds` at the end of
every invocation. `PILOT_RESULTS.md` is rendered by:

```bash
cd "/path/to/project/paper_deepgof/theory/calibration_gap"
"/c/Program Files/R/R-4.4.1/bin/Rscript.exe" 04_report.R
```

It is safe to run at any moment, including mid-run: it reports on whatever rows exist and says how
many datasets that was.

## Files

| file | what |
|---|---|
| `00_verify_and_bench.R` | forward-pass and fast-path verification, timing (already passed, `verify.rds`) |
| `01_batched_net.R` | batched re-expression of the shipped forward pass |
| `02_pilot.R` | the pilot; checkpointed, resumable |
| `03_profile.R` | per-component timing |
| `04_report.R` | renders `PILOT_RESULTS.md` |
| `rows/n<N>_rep<i>.csv` | per-dataset checkpoint, written the moment a dataset finishes |
| `per_dataset.csv` | all rows; appended live, rebuilt from `rows/` at the end |
| `cell_summary.csv`, `run_meta.rds` | rebuilt at the end of every invocation |
| `scores_n<N>_rep1.csv.gz` | the M scores of each law for rep 1 of each cell |
| `RUN.log`, `progress.log` | master log, worker heartbeat |
| `smoke_M2000/` | the earlier M = 2,000, D = 4 smoke test, kept out of the way |
