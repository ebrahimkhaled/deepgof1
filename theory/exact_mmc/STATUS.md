# exact_mmc -- STATUS log

One-day feasibility study: can DeepGOF-1 be made exactly valid at finite n by replacing the
plug-in parametric bootstrap with Dufour's confidence-set maximised Monte Carlo (CMMC)?

Folder owned by this study only. The calibration_gap pilot is running concurrently in the
sibling folder with 12 workers; nothing here touches it. This study uses <= 6 workers.

---

- 2026-09-16 -- step 0. Read AE_FINITE_SAMPLE_BOUND.md section 5(d), rdeploy/deepgof.R +
  deepgof_forward.R, calibration_gap/00_verify_and_bench.R + 01_batched_net.R + 02_pilot.R,
  code/universal_design.R, code/s3_grid.R, code/s3_score.py. Confirmed the deployed p-value is
  p = (1 + #{b : T*_b >= T_obs}) / (B+1), reject for LARGE T. Folder created.

- 2026-09-16 07:30 -- step 1. `10_core.R` written (CMMC machinery), `11_bench.R` run. ALL CHECKS
  PASS: batched net vs PyTorch test_output.csv 4.6e-11; local null generator RNG-identical to
  universal_design.R's gen_design_null (max|diff| 0), captured theta0 reproduces p0 to 2.2e-16;
  fast fit/map/score vs shipped .deepgof_map + deepgof_score over 240 comparisons, max|diff|
  2.3e-15, zero Firth-flag mismatches; pB() from the CRN matrix reproduces the shipped-path
  p-value exactly (0.320 vs 0.320, count 63); early abandonment leaves the count unchanged;
  theta0 in C_gamma on 400/400 null datasets (nominal >= .99).

- 2026-09-16 07:40 -- step 2. `12_geom.R`: geometry of the Bernstein set at n=100, 40 null
  datasets, 200 random directions each. FINDINGS THAT SHAPE THE STUDY:
    * radial extent of C_gamma, in Wald-SE units: 3.8 (tightest direction) to 34.7 (loosest),
      median over directions 5.4;
    * C_gamma is UNBOUNDED on 3 of 40 datasets (directions along which the score saturates
      inside the Bernstein thresholds). Dufour needs a compact set -- so the search is now
      truncated at RCAP = 200 Wald SEs and the truncation flagged;
    * an 11^3 grid over the BOUNDING BOX of C_gamma puts a median of 6.5 points inside the set
      (min 0, max 150). The box grid is a very poor maximiser; this is a result, not a bug.

- 2026-09-16 07:56 -- step 3. Profiling: one pB evaluation (199 refits + network) is ~0.65 s at
  theta-hat, but ~3 s at parameter values where the replicates are near-separated and the
  deployed policy falls back to brglmFit (14 ms/fit vs glm.fit's 2 ms). Two smoke datasets ran
  past 13 min. A per-dataset search BUDGET of 480 s was added, with a `trunc` flag: on a flagged
  dataset p_exact is a strict lower bound and exactness is void there.

- 2026-09-16 07:58 -- step 4. MAIN RUN LAUNCHED. `13_run.R`, 6 workers, 300 reps/arm queued
  interleaved (null i, alt i), 205-minute queuing deadline, checkpointed to rows/*.csv.

