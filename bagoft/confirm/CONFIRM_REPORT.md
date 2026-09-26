# BAGofT confirmatory arm -- report

Pre-declared analysis (PREDECLARATION.md), `analyze_confirm.py`. Implementation `bagoft_fast()` (identical to BAGofT 1.0.0, VERIFY_FAST.log).

Datasets: 2400 nulls, 4800 alternatives, 0 dropped as degenerate, 0 errors (0.00%).

## O1 -- the Table 2 row

| size | n=50 | n=100 | n=500 | overall | MC SE (overall) | inflated? |
|---|---|---|---|---|---|---|
| 0.048 | 0.373 | 0.590 | 0.838 | **0.600** | 0.0042 | False |

## O2 -- primary paired test: **DeepGOF-1 AHEAD**

| group | pairs | BAGofT | DeepGOF-1 | only BAGofT | only DeepGOF-1 | exact p | Holm p |
|---|---|---|---|---|---|---|---|
| all 48 power cells (O2, primary) | 4800 | 0.595 | 0.635 | 389 | 581 | 7.641e-10 | -- |
| n = 50 | 1600 | 0.344 | 0.460 | 105 | 291 | 2.579e-21 | 1.805e-20 |
| n = 100 | 1600 | 0.598 | 0.619 | 142 | 176 | 0.06407 | 0.1922 |
| n = 500 | 1600 | 0.843 | 0.825 | 142 | 114 | 0.09131 | 0.1922 |
| Setting 1 | 1200 | 0.759 | 0.840 | 9 | 106 | 3.693e-22 | 2.955e-21 |
| Setting 2 | 1200 | 0.756 | 0.578 | 311 | 97 | 3.1e-27 | 3.1e-26 |
| Setting 3 | 1200 | 0.227 | 0.335 | 45 | 175 | 2.699e-19 | 1.619e-18 |
| Setting 4 | 1200 | 0.637 | 0.786 | 24 | 203 | 1.682e-36 | 1.85e-35 |
| C = 0.25 | 1200 | 0.284 | 0.413 | 37 | 192 | 1.977e-26 | 1.779e-25 |
| C = 0.5 | 1200 | 0.577 | 0.648 | 78 | 164 | 3.373e-08 | 1.686e-07 |
| C = 0.75 | 1200 | 0.741 | 0.728 | 130 | 114 | 0.3369 | 0.3369 |
| C = 1 | 1200 | 0.777 | 0.749 | 144 | 111 | 0.04487 | 0.1795 |

## O4 -- sizes on this run's own nulls (replicates 13-212)

BAGofT 0.048, DeepGOF-1 0.052. O2 compares raw decisions at p <= .05; the more liberal test is favoured by at most roughly the difference.

