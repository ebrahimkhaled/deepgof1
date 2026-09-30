SUPPORT application (Section 9 of the manuscript): code, data and results.

Data: data/support2.csv, the public SUPPORT2 release (UCI Machine Learning Repository,
https://archive.ics.uci.edu/dataset/880/support2), 9,105 patients; 8,873 complete cases are used.

Requirements: R 4.4.1 with ebrahim.gof 2.9.0 or later (CRAN), splines, givitiR; Python 3 with numpy
and matplotlib for the figure. Version 2.9.0 is needed: earlier versions return p = 1 for any model
with a spline term, because their bootstrap refit failed on such terms. Run the scripts in order from
this folder:
  01_M0.R               M0 and the shipped DeepGOF-1 test (p = .005, axes meanbp & num.co)
  02_repairs.R          the fixed repair list (a)-(d); repair (b), ns(meanbp, 3), is chosen
  03_fixed_axes.R       M0 and M1 read on the axes meanbp x num.co fixed in advance (p = .005, .90)
  04_axis_robustness.R  M0 on every pair of the four top-ranked covariates (all p = .005)
  05_classical_tests.R  the classical tests on M0 and M1 (Table 8)
  06-09                 the 28-pair reading of M0, M1, M2 and M3 and the ladder rungs M2, M3
  10_clinical_impact.R  likelihood-ratio tests, optimism-corrected AUC and calibration slope,
                        predicted risks for patient profiles
  11_export_maps.R      the three maps of the figure (M0, M1, M2 on blood pressure x respiratory
                        rate) and their simultaneous threshold; 13_figure.py draws Figure 6
  12_allpairs.R         every model under both readings of deepgof1(); the all-pairs reading gives
                        p = .005, .005, .020, .21 for M0-M3 (Table 8)
Every random step has a fixed seed. results/screen_SUPPORT/ holds the outputs as produced for the
manuscript; rerunning the scripts overwrites them with identical values.
