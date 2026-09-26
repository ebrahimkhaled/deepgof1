# =============================================================================
# analyze_bagoft_niche.R -- rejection rates of the BAGofT arm on the niche, n = 200.
#
# 1. Pairing: every row of niche_rows_bagoft/ is joined to theory/niche/rows/ on (scen, rep),
#    and the event count and (fam, C, n) must agree on every row; the script stops otherwise.
# 2. The p-value used is p (= BAGofT p.value, splits combined by the mean), the arm's headline.
#    p_min (= p.value3) is reported as a secondary column. bag_pmin_stat is the statistic, never
#    a p-value, and is not used.
# 3. Raw rejection: p <= .05, with the binomial SE; a failed run (NA) cannot reject.
# 4. Matched level: each alternative cell is placed at exactly .05 on BAGofT's OWN null rows at
#    n = 200 (null_n200, the same reps), by the randomized rule of theory/niche/analyze_niche.R
#    (matched_phi, copied below); NA counts as p = 1 there, as in that script. With 100 null
#    rows the null-estimated threshold itself is noisy: its value and the randomization weight
#    are printed in the table.
# Writes niche_bagoft_rates.csv.  Run: Rscript analyze_bagoft_niche.R
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
NRO  <- "/path/to/project/paper_deepgof/theory/niche/rows"
ALPHA <- 0.05

fs <- list.files(Sys.getenv("AB_ROWS", file.path(HERE, "niche_rows_bagoft")), pattern = "_r[0-9]+[.]csv$", full.names = TRUE)
B  <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
if (anyDuplicated(B[, c("scen", "rep")])) stop("duplicate (scen, rep) rows")
S  <- do.call(rbind, lapply(list.files(NRO, pattern = "_n200_r[0-9]+[.]csv$", full.names = TRUE),
                            utils::read.csv, stringsAsFactors = FALSE))[, c("scen", "fam", "C", "n", "rep", "events")]
M  <- merge(S, B, by = c("scen", "rep"), suffixes = c("_st", ""))
if (nrow(M) != nrow(B)) stop(sprintf("%d BAGofT rows but %d joined to the stored rows", nrow(B), nrow(M)))
bad <- M$events_st != M$events | M$fam_st != M$fam | M$C_st != M$C | M$n_st != M$n
if (any(bad)) { print(utils::head(M[bad, ])); stop("pairing FAILED") }
cat(sprintf("joined %d BAGofT rows to the stored niche rows; events identical on every row\n", nrow(M)))

matched_phi <- function(pn, pa, alpha = ALPHA) {           # from theory/niche/analyze_niche.R
  t <- sort(unique(pn)); F <- findInterval(t, sort(pn)) / length(pn)
  ok <- which(F <= alpha)
  if (!length(ok)) { g <- if (F[1] > 0) alpha / F[1] else 0; return(list(phi = g * (pa <= t[1]), t0 = t[1], g = g)) }
  i <- max(ok); t0 <- t[i]; F0 <- F[i]
  if (i + 1 > length(t)) return(list(phi = as.numeric(pa <= t0), t0 = t0, g = 0))
  t1 <- t[i + 1]; F1 <- F[i + 1]; g <- if (F1 <= F0) 0 else (alpha - F0) / (F1 - F0)
  list(phi = as.numeric(pa <= t0) + g * as.numeric(pa > t0 & pa <= t1), t0 = t0, g = g)
}
raw <- function(p) { r <- mean(!is.na(p) & p <= ALPHA); c(r, sqrt(r * (1 - r) / length(p))) }
na1 <- function(p) ifelse(is.na(p), 1, p)

N0 <- M[M$scen == "null_n200", ]
out <- list()
for (sc in sort(unique(M$scen))) {
  X <- M[M$scen == sc, ]
  rp <- raw(X$p); rm <- raw(X$p_min)
  row <- data.frame(scen = sc, fam = X$fam[1], C = X$C[1], n = X$n[1], R = nrow(X),
                    na_share = mean(is.na(X$p)), reps = sprintf("%d-%d", min(X$rep), max(X$rep)),
                    raw_p = rp[1], raw_p_se = rp[2], raw_pmin = rm[1], raw_pmin_se = rm[2])
  if (nrow(N0)) {
    mp <- matched_phi(na1(N0$p), na1(X$p)); mm <- matched_phi(na1(N0$p_min), na1(X$p_min))
    row$matched_p <- mean(mp$phi); row$matched_p_se <- stats::sd(mp$phi) / sqrt(nrow(X))
    row$matched_p_t0 <- mp$t0; row$matched_p_gamma <- mp$g
    row$matched_pmin <- mean(mm$phi); row$matched_pmin_se <- stats::sd(mm$phi) / sqrt(nrow(X))
    row$null_R <- nrow(N0)
  }
  row$mean_time_s <- mean(X$time_s)
  out[[sc]] <- row
}
O <- do.call(rbind, out); rownames(O) <- NULL
utils::write.csv(O, Sys.getenv("AB_OUT", file.path(HERE, "niche_bagoft_rates.csv")), row.names = FALSE)
print(O, digits = 3)
cat("\nraw = p <= .05 with binomial SE; matched = randomized to exactly .05 on BAGofT's own null_n200 rows",
    "(SE conditional on the null threshold, which rests on", nrow(N0), "null rows)\n")
