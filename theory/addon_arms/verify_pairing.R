# =============================================================================
# verify_pairing.R -- are the add-on rows computed on the SAME datasets as the stored rows, and
# which Stukel form do the stored rows hold?
#
# For a spread of stored blocks it regenerates datasets with the add-on scripts' OWN generator
# code (the liu_defs / niche_defs expressions parsed out of liu_addon.R and niche_addon.R, so the
# code checked is the code that runs) and compares, dataset by dataset:
#   HL        recomputed  vs the stored HL column      -> must be identical (pairing)
#   Stukel    joint form and marginal form vs the stored Stukel column -> tells which form it is
# It also checks, on a few datasets, that the persistent-session GiViTI call returns the same
# p-value as run.all.gof(tests = "GiViTI"), which is ebrahim.gof's own gof_giviti().
# Output: verify_pairing.log (console) and verify_pairing.csv (every compared dataset).
# Run:    Rscript verify_pairing.R
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
RES  <- "/path/to/project/paper_deepgof/release/deepgof1/results/liu"
NRO  <- "/path/to/project/paper_deepgof/theory/niche/rows"
source(file.path(HERE, "addon_common.R"))

grab <- function(file, names) {                     # evaluate only the named top-level assignments
  ex <- parse(file)
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                    ((is.name(e[[2]]) && as.character(e[[2]]) %in% names) ||
                     (is.call(e[[2]]) && identical(e[[2]][[2]], as.name("SC")) && "SC" %in% names)))
    eval(e, globalenv())
}
grab(file.path(HERE, "liu_addon.R"), "liu_defs");   LE <- new.env(); eval(liu_defs, LE)
grab(file.path(HERE, "niche_addon.R"), c("CH", "SEEDB", "NNULL", "NALT", "SC", "niche_defs")); NE <- new.env(); eval(niche_defs, NE)

same <- function(a, b) is.na(a) & is.na(b) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(b)))
out <- list()

## ---- Liu grid: 6 blocks, every setting, every n, both corpora; nulls and all four C ----
LB <- c("S1_n50_s101000", "S2_n100_s101000", "S3_n500_s101000", "S4_n50_s102000", "S2_n500_s102000",
        "S3_n100_s102000")
for (b in LB) {
  m  <- regmatches(b, regexec("^S([0-9])_n([0-9]+)_s([0-9]+)$", b))[[1]]
  SET <- as.integer(m[2]); NOBS <- as.integer(m[3]); SEED0 <- as.integer(m[4])
  st <- utils::read.csv(file.path(RES, b, "blk_pvalues.csv"), check.names = FALSE)
  st$ci <- ifelse(st$kind == 2, match(st$C, c(0, 0.25, 0.5, 0.75, 1.0)), 0L)
  pick <- rbind(st[st$kind == 0 & st$rep <= 15, ], st[st$kind == 2 & st$rep <= 5, ])
  for (j in seq_len(nrow(pick))) {
    k <- pick[j, ]
    a <- LE$liu_one(SET, NOBS, SEED0, k$kind, k$rep, if (k$kind == 2) k$C else NA_real_, k$ci, dogv = FALSE)
    out[[length(out) + 1L]] <- data.frame(bench = "liu", block = b, kind = k$kind, rep = k$rep, C = k$C,
      HL_stored = k$HL, HL_check = a$HL_check, Stukel_stored = k$Stukel,
      Stukel_joint = a$Stukel_joint, Stukel_marginal = a$Stukel_marginal)
  }
}

## ---- niche: 5 stored row files spread over the scenarios ----
NF <- c("null_n200_r0001.csv", "null_n500_r1981.csv", "inter_C0.3_n200_r0501.csv", "quad_C0.8_n500_r0981.csv",
        "inter_C0.8_n500_r0021.csv")
for (f in NF) {
  st <- utils::read.csv(file.path(NRO, f))
  sc <- SC[SC$scen == st$scen[1], ]
  for (j in seq_len(nrow(st))) {
    s <- SEEDB + sc$id * 100000L + st$rep[j]
    set.seed(s); dat <- NE$gen(sc); fit <- NE$fit_pol(dat)
    a <- arms_for_fit(fit, do_giviti = FALSE)
    out[[length(out) + 1L]] <- data.frame(bench = "niche", block = f, kind = NA, rep = st$rep[j], C = st$C[j],
      HL_stored = st$HL[j], HL_check = a[["HL_check"]], Stukel_stored = st$Stukel[j],
      Stukel_joint = a[["Stukel_joint"]], Stukel_marginal = a[["Stukel_marginal"]])
  }
}
V <- do.call(rbind, out)
utils::write.csv(V, file.path(HERE, "verify_pairing.csv"), row.names = FALSE)

for (bn in c("liu", "niche")) {
  v <- V[V$bench == bn, ]
  cat(sprintf("\n[%s] %d datasets from %d stored blocks\n", bn, nrow(v), length(unique(v$block))))
  cat(sprintf("  HL identical to stored            : %d / %d   (max |diff| %.3g)\n", sum(same(v$HL_check, v$HL_stored)),
              nrow(v), max(abs(v$HL_check - v$HL_stored), na.rm = TRUE)))
  cat(sprintf("  Stukel stored == joint form       : %d / %d   (max |diff| %.3g)\n", sum(same(v$Stukel_joint, v$Stukel_stored)),
              nrow(v), max(abs(v$Stukel_joint - v$Stukel_stored), na.rm = TRUE)))
  cat(sprintf("  Stukel stored == marginal form    : %d / %d   (max |diff| %.3g)\n", sum(same(v$Stukel_marginal, v$Stukel_stored)),
              nrow(v), max(abs(v$Stukel_marginal - v$Stukel_stored), na.rm = TRUE)))
}

## ---- GiViTI: persistent session vs run.all.gof's own gof_giviti, 6 datasets ----
cat("\n[GiViTI] session call vs run.all.gof(tests = 'GiViTI')\n")
for (k in 1:3) for (bn in c("liu", "niche")) {
  fit <- if (bn == "liu") { set.seed(101000 + 2L*500000L + k); d <- LE$gen_liu(100L, 2L, 0); LE$fit_pol(LE$FRM(d), d)$fit }
         else { sc <- SC[SC$scen == "quad_C0.8_n200", ]; set.seed(SEEDB + sc$id * 100000L + k); NE$fit_pol(NE$gen(sc)) }
  a  <- arms_for_fit(fit, do_giviti = TRUE)
  rr <- suppressMessages(run.all.gof(fit, tests = "GiViTI", include_slow = TRUE, install = "no"))
  cat(sprintf("  %-5s rep %d: session %.15g | run.all.gof %.15g | identical %s | %.2f s\n", bn, k, a[["GiViTI"]],
              rr$p_value[1], identical(unname(a[["GiViTI"]]), unname(as.numeric(rr$p_value[1]))), a[["giviti_secs"]]))
}
giviti_close()
cat("\nverify_pairing DONE", format(Sys.time()), "\n")
