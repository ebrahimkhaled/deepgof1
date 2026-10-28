# 31_C_smoke_check.R -- quick agreement checks on rows_C_smoke (or rows_C with argument "full")
T3 <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3"
d <- if (length(commandArgs(TRUE)) && commandArgs(TRUE)[1] == "full") "rows_C" else "rows_C_smoke"
R <- do.call(rbind, lapply(list.files(file.path(T3, d), full.names = TRUE), read.csv))
ok <- R$status == "ok"; print(table(R$status)); print(table(S = R$S, fallback = R$fallback))
b <- ok & R$fallback == 0
cat("SstBoth vs stored, both terms defined: max |diff|", max(abs(R$p_sst_both[b] - R$p_stukel_stored[b]), na.rm = TRUE),
    "| stored NA there:", sum(is.na(R$p_stukel_stored[b])), "\n")
cat("SstBoth vs LogisticDx: max |diff|", max(abs(R$p_sst_both[b] - R$p_sst_ldx[b]), na.rm = TRUE), "\n")
cat("stored Stukel NA wherever fallback > 0:", all(is.na(R$p_stukel_stored[ok & R$fallback > 0])), "\n")
cat("Osius-Rojek hand vs stored:", max(abs(R$p_or_hand[ok] - R$p_or_stored[ok])), "| LogisticDx vs stored:",
    max(abs(R$p_or_ldx[ok] - R$p_or_stored[ok]), na.rm = TRUE), "| LogisticDx NA:", sum(is.na(R$p_or_ldx[ok])), "\n")
print(table(firth = R$firth))
big <- ok & abs(R$p_or_hand - R$p_or_stored) > 1e-10
cat("rows with Osius-Rojek |hand - stored| > 1e-10:", sum(big), "\n"); print(R[big, c("block", "kind", "C", "rep", "firth", "p_or_hand", "p_or_ldx", "p_or_stored")])
cat("max |hand - stored| on non-Firth fits:", max(abs(R$p_or_hand - R$p_or_stored)[ok & R$firth == 0]), "\n")
cat("min_pq of the discrepant rows:", R$min_pq[big], "\n")
