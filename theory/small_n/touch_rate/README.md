Touch rate and miss probability of the covariate-only deletion rule (Online Resource 1, Section 16).

touch_sim.py simulates covariates only (20,000 replicates per cell; laws normal, uniform, t4,
lognormal, Bernoulli, and the three-covariate design of Liu et al.'s Setting 2) and writes cells/.

touch_analyse.py reads cells/ and writes TOUCH_RATE.json: the touch rate on clean data, the miss
probability of a decimal slip (overall, gross |x| > .45, mild), three analytic approximations for
normal covariates (the naive union bound, the same with the finite-sample MAD factor, and an
integral over a lognormal law for the scaled MAD with mean (n - .8)/n and variance 1.361/n, within
15% of the simulation in every normal cell and within 5% from n = 100), and the adjusted
thresholds (a Bonferroni threshold, which touches 39% of clean datasets at d = 10, n = 30, and the
empirical 95th percentile of max o_i, which holds the touch rate at .04-.05).

power_check_d10.R is the paired false-fire check reported in Section 16 (DeepGOF-1 only): the
per-replicate rows are in rows_power/, the run log in power_check.log, and the summary in the
"power_cost" entry of TOUCH_RATE.json.
