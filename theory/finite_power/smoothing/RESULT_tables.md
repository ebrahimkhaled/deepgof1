Certified values rounded up: t_cert to the printed digit, n* smoothing = ceil(4 t_cert^2). Ratios from unrounded values.

| direction | true t* (SE) | true n* | n* CROWN (paper) | ratio CROWN | t_cert smoothing | n* smoothing | ratio smoothing [2 SE] |
|---|---|---|---|---|---|---|---|
| step | 4.274 (0.006) | 73 | 3,789 | 51.9 | 4.327 | 75 | **1.02** [1.02, 1.03] |
| diffuse | 5.105 (0.006) | 104 | 4,311 | 41.4 | 5.154 | 107 | **1.02** [1.01, 1.02] |
| stripe | 5.547 (0.006) | 123 | 5,881 | 47.8 | 5.606 | 126 | **1.02** [1.02, 1.03] |
| saddle | 7.156 (0.007) | 205 | 8,044 | 39.3 | 7.213 | 209 | **1.02** [1.01, 1.02] |
| random | 13.529 (0.013) | 732 | 13,918 | 19.0 | 13.701 | 751 | **1.03** [1.02, 1.03] |

Paper's own numbers for reference (threshold = MC point estimate 2.4616, 20,000 draws): true n* = step 72, diffuse 103, stripe 122, saddle 204, random 729; CROWN ratios 52.6, 41.9, 48.2, 39.4, 19.1.

**Fixed threshold, power 1/2**

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 2.312 (0.003) | 21.4 | 2.342 | 22 | **1.03** [1.02, 1.03] |
| diffuse | 3.048 (0.004) | 37.2 | 3.088 | 39 | **1.03** [1.02, 1.03] |
| stripe | 3.341 (0.004) | 44.6 | 3.389 | 46 | **1.03** [1.02, 1.03] |
| saddle | 4.689 (0.005) | 88.0 | 4.774 | 92 | **1.04** [1.03, 1.04] |
| random | 8.839 (0.008) | 312.5 | 9.001 | 325 | **1.04** [1.03, 1.04] |

**Bootstrap rule B = 199, power .95**

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 4.302 (0.005) | 74.0 | 4.429 | 79 | **1.06** [1.05, 1.06] |
| diffuse | 5.145 (0.006) | 105.9 | 5.268 | 112 | **1.05** [1.04, 1.05] |
| stripe | 5.594 (0.006) | 125.1 | 5.730 | 132 | **1.05** [1.04, 1.05] |
| saddle | 7.232 (0.007) | 209.2 | 7.432 | 221 | **1.06** [1.05, 1.06] |
| random | 13.657 (0.013) | 746.0 | 14.083 | 794 | **1.06** [1.06, 1.07] |

**Bootstrap rule B = 199, power 1/2**

| direction | true t* (SE) | true n* | t_cert | n* smoothing | ratio [2 SE] |
|---|---|---|---|---|---|
| step | 2.326 (0.004) | 21.6 | 2.438 | 24 | **1.10** [1.09, 1.11] |
| diffuse | 3.066 (0.004) | 37.6 | 3.206 | 42 | **1.09** [1.09, 1.10] |
| stripe | 3.356 (0.004) | 45.1 | 3.507 | 50 | **1.09** [1.09, 1.10] |
| saddle | 4.719 (0.005) | 89.1 | 4.914 | 97 | **1.08** [1.08, 1.09] |
| random | 8.886 (0.009) | 315.9 | 9.305 | 347 | **1.10** [1.09, 1.10] |

Sensitivity to the Monte Carlo budget (fixed threshold, .95): certified t at eta = .001 vs .01: step 4.327 / 4.323, diffuse 5.154 / 5.150, stripe 5.606 / 5.602, saddle 7.213 / 7.210, random 13.701 / 13.697.

| direction | net points J | draws (M) | net range | tail intervals | t_global (rounded up) | max Lambda |
|---|---|---|---|---|---|---|
| step | 260 | 6.3 | [1.27, 250] | 19 | 6,944 | 3.28 |
| diffuse | 271 | 7.4 | [1.52, 249] | 21 | 10,430 | 3.44 |
| stripe | 278 | 8.0 | [1.66, 250] | 22 | 8,835 | 3.42 |
| saddle | 300 | 10.0 | [2.14, 250] | 31 | 12,862 | 3.66 |
| random | 387 | 18.2 | [4.05, 249] | 50 | 20,908 | 3.65 |

Slack of the fixed .95 certificate over the truth: step 9.5 SE / 1.23 % in t / 2.48 % in n, diffuse 8.3 SE / 0.95 % in t / 1.91 % in n, stripe 9.7 SE / 1.06 % in t / 2.13 % in n, saddle 8.8 SE / 0.80 % in t / 1.60 % in n, random 13.1 SE / 1.27 % in t / 2.55 % in n. Relaxation ratios against the matched truth: 51.9, 41.4, 47.8, 39.3, 19.0 (range 19--52).

Tail chains: c_max = 5.8321; smallest CROWN lower bound per direction: step 247.940, diffuse 97.419, stripe 39.515, saddle 5.953, random 10.658.
