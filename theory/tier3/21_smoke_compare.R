# 21_smoke_compare.R -- smoke rows vs the timing pilot (tier3_timing/pilot_AC_w22.csv) on the shared data sets
T3 <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3"
R <- do.call(rbind, lapply(list.files(file.path(T3, "rows_AB_smoke"), full.names = TRUE), read.csv))
P <- read.csv(file.path(T3, "../tier3_timing/pilot_AC_w22.csv"))
kR <- with(R, ifelse(design == "liu", sprintf("liu_%d_%d_%d_%.2f_%d", S, n, kind, C, rep), sprintf("niche_%d_%d_%d", n, ci, rep)))
kP <- with(P, ifelse(design == "liu", sprintf("liu_%d_%d_%d_%.2f_%d", S, n, kind, C, rep), sprintf("niche_%d_%d_%d", n, ci, rep)))
m <- match(kR, kP); cat("smoke rows:", nrow(R), "| shared with pilot:", sum(!is.na(m)), "\n")
L <- !is.na(m) & R$design == "liu"
for (v in c("p_ss", "p_score", "p_score_asym")) cat(sprintf("liu %-13s identical to pilot on %d of %d (max |diff| %.2g)\n", v,
  sum(R[[v]][L] == P[[v]][m[L]], na.rm = TRUE), sum(L), max(abs(R[[v]][L] - P[[v]][m[L]]), na.rm = TRUE)))
cat(sprintf("liu p_dg == pilot p_dg on %d of %d\n", sum(R$p_dg[L] == P$p_dg[m[L]]), sum(L)))
G <- !is.na(m)
cat(sprintf("GRP with the pilot seed equals pilot p_grp to 1e-12 relative on %d of %d (NA in both: %d)\n",
  sum(abs(R$p_grp_pilotseed[G] - P$p_grp[m[G]]) <= 1e-12 * pmax(1e-300, abs(P$p_grp[m[G]])), na.rm = TRUE), sum(G), sum(is.na(R$p_grp_pilotseed[G]) & is.na(P$p_grp[m[G]]))))
cat("events vs pilot:", all(R$events[G] == P$events[m[G]]), "\n")
print(table(R$design, R$score_df)); cat("refit failures:", sum(R$n_fail), "| score NA in bootstrap:", sum(R$sc_na_boot), "\n")
print(R[is.na(R$p_grp), c("design", "S", "n", "kind", "C", "rep", "events", "grp_msg")])
cat("\nmean seconds by design/n:\n"); print(aggregate(cbind(t_combo, t_nets, t_grp) ~ design + n, R, mean))
cat("\nnetwork p-values (first rows):\n"); print(head(R[, c("scen", "kind", "C", "p_dg", grep("^p_(Control|NoInter|NoQuad)", names(R), value = TRUE))], 6), digits = 3)
cat("\nGRP mismatches (ours with pilot seed vs pilot):\n")
X <- data.frame(key = kR[G], ours = R$p_grp_pilotseed[G], ours_newseed = R$p_grp[G], pilot = P$p_grp[m[G]])
print(head(X[is.na(X$ours) | is.na(X$pilot) | X$ours != X$pilot, ], 25), digits = 4)
