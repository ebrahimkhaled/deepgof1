## fdic_download.R -- the finance application's data, from two public FDIC API calls.
##
## The data rule was fixed before any test was run (the declaration is available from the author).
## SHIPPED WITH THE PAPER: everything
## below reproduces from two keyless GET requests against the FDIC BankFind Suite API, so a referee
## can rebuild the analysis dataset without credentials, a login, or a data-use agreement.
##
##   call 1  /banks/financials?filters=REPDTE:20081231   -> 8,395 banks' 2008-12-31 call reports
##   call 2  /banks/failures                             -> every FDIC-resolved failure, with FAILDATE
##   join    CERT
##
## OUTCOME (pre-declared): failed = appears in /failures with FAILDATE in 2009-01-01..2011-12-31.
## COVARIATES (pre-declared, 5, CAMELS-proxy in the Cole & White 2012 lineage) -- every one an
## explicit ratio of raw call-report items, so nothing depends on an FDIC-derived convenience field:
##   eq_ta   = EQ / ASSET                       capital
##   npa_ta  = (P3ASSET + P9ASSET) / ASSET      asset quality (past-due 30-89 + 90+)
##   roa     = ROA                              earnings (already a ratio)
##   bro_dep = BRO / DEP                        funding fragility (brokered-deposit share)
##   cre_ta  = (LNRECONS + LNRENRES) / ASSET    concentration (construction + nonres CRE)
##
## Nothing here looks at an outcome. The cell rule, the model, and the covariates were all fixed
## before any of this was run.

## run from the archive root (or set DEEPGOF_ARCHIVE); writes application/fdic/data/
ARC <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
OUT <- file.path(ARC, "application", "fdic", "data"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
API <- "https://api.fdic.gov/banks"

fetch_all <- function(path, filters, fields, page = 5000L) {
  rows <- list(); off <- 0L
  repeat {
    url <- sprintf("%s/%s?filters=%s&fields=%s&limit=%d&offset=%d&format=json",
                   API, path, utils::URLencode(filters, reserved = TRUE), fields, page, off)
    txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
    if (is.null(txt)) stop("FDIC API unreachable for ", path)
    j <- jsonlite::fromJSON(paste(txt, collapse = ""), simplifyDataFrame = TRUE)
    d <- j$data$data
    if (is.null(d) || !nrow(d)) break
    rows[[length(rows) + 1L]] <- d
    off <- off + nrow(d)
    if (off >= j$meta$total) break
  }
  do.call(rbind, rows)
}

cat("call 1: 2008-12-31 call reports ...\n"); flush.console()
FIN <- fetch_all("financials", "REPDTE:20081231",
                 "CERT,STALP,ASSET,EQ,ROA,BRO,DEP,P3ASSET,P9ASSET,LNRECONS,LNRENRES")
cat(sprintf("  %d banks\n", nrow(FIN)))

cat("call 2: failures ...\n"); flush.console()
FAIL <- fetch_all("failures", "FAILDATE:[2009-01-01 TO 2011-12-31]", "CERT,FAILDATE,NAME,PSTALP")
cat(sprintf("  %d failures 2009-2011\n", nrow(FAIL)))

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
d <- data.frame(
  cert    = num(FIN$CERT),
  state   = as.character(FIN$STALP),
  asset   = num(FIN$ASSET),
  eq_ta   = num(FIN$EQ)      / num(FIN$ASSET),
  npa_ta  = (num(FIN$P3ASSET) + num(FIN$P9ASSET)) / num(FIN$ASSET),
  roa     = num(FIN$ROA),
  bro_dep = num(FIN$BRO)     / num(FIN$DEP),
  cre_ta  = (num(FIN$LNRECONS) + num(FIN$LNRENRES)) / num(FIN$ASSET),
  stringsAsFactors = FALSE)
d$failed <- as.integer(d$cert %in% num(FAIL$CERT))

keep <- is.finite(d$eq_ta) & is.finite(d$npa_ta) & is.finite(d$roa) &
        is.finite(d$bro_dep) & is.finite(d$cre_ta) & d$asset > 0
cat(sprintf("complete-covariate banks: %d of %d (%.2f%% dropped for missing/zero-asset)\n",
            sum(keep), nrow(d), 100*(1-mean(keep))))
d <- d[keep, ]

write.csv(d, file.path(OUT, "fdic_2008_banks.csv"), row.names = FALSE)
write.csv(FAIL, file.path(OUT, "fdic_failures_2009_2011.csv"), row.names = FALSE)
st <- as.data.frame(table(d$state[d$failed == 1]), stringsAsFactors = FALSE)
names(st) <- c("state", "failures")
st$banks <- as.integer(table(d$state)[st$state])
st <- st[order(-st$failures), ]
write.csv(st, file.path(OUT, "fdic_state_counts.csv"), row.names = FALSE)

cat("\nstates with >= 15 failures (the pre-declared replication cells):\n")
print(st[st$failures >= 15, ], row.names = FALSE)
cat(sprintf("\nPRIMARY CELL (Georgia): %d banks, %d failures (%.1f%% event rate)\n",
            sum(d$state == "GA"), sum(d$state == "GA" & d$failed == 1),
            100*mean(d$failed[d$state == "GA"])))
