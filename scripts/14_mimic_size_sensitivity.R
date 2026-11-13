suppressMessages({library(lavaan); library(MASS)})
args <- commandArgs(TRUE); x <- read.csv(args[1]); outf <- args[2]; set.seed(20261113)
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
fm <- c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr","aic","bic")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
sink(outf)
## 1. loadings in logit units from pretrained models (var(g) = 1)
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); li <- which(pe$op == "=~"); lam <- setNames(pe$est[li], pe$rhs[li])[vl]
Vlam <- vcov(fp)[paste0("g=~", vl), paste0("g=~", vl)]
cat("Pretrained unstandardised loadings (logit units, var(g)=1)\n"); print(round(lam, 3))
## 2. pairs
idx <- match(x$base, x$model); ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base
sets <- list(all = ok, strict = ok & isF(x$type) & x$merged == "False" & isP(x$type[idx]))
sets$strict[is.na(sets$strict)] <- FALSE
for (s in names(sets)) {
  C <- x[sets[[s]], ]; Bm <- x[idx[sets[[s]]], ]
  G <- as.matrix(C[, vl]) - as.matrix(Bm[, vl]); colnames(G) <- vl
  L <- aggregate(G, by = list(base = C$base), FUN = mean); M <- as.matrix(L[, vl]); nL <- nrow(M)
  mbar <- colMeans(M); V <- cov(M) / nL
  gls <- function(lm) { Vi <- solve(V); mu <- as.numeric(t(lm) %*% Vi %*% mbar / (t(lm) %*% Vi %*% lm)); r <- mbar - lm * mu
    c(mu = mu, Q = as.numeric(t(r) %*% Vi %*% r)) }
  g0 <- gls(lam)
  draws <- t(replicate(2000, gls(mvrnorm(1, lam, Vlam))))
  cat("\n=== Pair set:", s, "| pairs:", nrow(G), "| lineages:", nL, "===\n")
  cat("Child types:\n"); print(table(C$type)); cat("Base types:\n"); print(table(Bm$type))
  cat("Lineage-clustered mean gain (logit units) and SE\n"); print(round(rbind(gain = mbar, se = sqrt(diag(V))), 3))
  cat("Gain in pretrained-SD units:\n"); print(round(mbar / sapply(P[, vl], sd), 3))
  cat(sprintf("GLS proportional-to-loadings model (pure g gain): mu = %.3f, Q = %.2f, df = 5, p = %.3g\n", g0["mu"], g0["Q"], pchisq(g0["Q"], 5, lower.tail = FALSE)))
  cat("With loading uncertainty (2000 draws): Q 2.5/50/97.5% =", round(quantile(draws[, "Q"], c(.025, .5, .975)), 2),
      "| p at those Q =", signif(pchisq(quantile(draws[, "Q"], c(.975, .5, .025)), 5, lower.tail = FALSE), 3), "\n")
  ## 3. MIMIC stacked (Ritchie et al. 2015 Models A/B/C); bases appear once
  S <- rbind(transform(C[, c("model","base",vl)], ft = 1), transform(unique(Bm[, c("model","base",vl)]), ft = 0))
  S$base[S$ft == 0] <- S$model[S$ft == 0]
  mods <- list(A = paste(m4, "\ng ~ ft"),
               B = paste(m4, "\ng ~ ft\nifeval_lg + musr_lg + mmlupro_lg ~ ft"),
               C = paste(m4, "\n", paste(vl, collapse = " + "), " ~ ft"))
  fits <- lapply(mods, function(mm) cfa(mm, data = S, estimator = "MLR", cluster = "base", std.lv = TRUE))
  cat("\nMIMIC on stacked children + bases (n =", nrow(S), ", clusters =", length(unique(S$base)), ")\n")
  print(round(t(sapply(fits, fitMeasures, fm)), 3))
  print(lavTestLRT(fits$C, fits$B, fits$A, method = "satorra.bentler.2001"))
  for (k in names(fits)) { p2 <- parameterEstimates(fits[[k]]); p2 <- p2[p2$op == "~" & p2$rhs == "ft", c("lhs","est","se","ci.lower","ci.upper","pvalue")]; cat("\nModel", k, "paths from ft\n"); print(p2, digits = 3) }
}
## 4. cohort: size control
P2 <- P[!is.na(P$upload_year) & P$params > 0, ]
P2$c24 <- as.numeric(P2$upload_year == 2024); P2$c25 <- as.numeric(P2$upload_year >= 2025); P2$lp <- log(P2$params)
cat("\n=== Cohort with log-parameter control (n =", nrow(P2), ") ===\n")
cm <- list(A = paste(m4, "\ng ~ c24 + c25 + lp"),
           C = paste(m4, "\ng ~ lp\n", paste(vl, collapse = " + "), " ~ c24 + c25"),
           A_full = paste(m4, "\ng ~ c24 + c25\n", paste(vl, collapse = " + "), " ~ lp"),
           C_full = paste(m4, "\n", paste(vl, collapse = " + "), " ~ c24 + c25 + lp"))
cf <- lapply(cm, function(mm) cfa(mm, data = P2, estimator = "MLR", cluster = "org", std.lv = TRUE))
print(round(t(sapply(cf, fitMeasures, fm)), 3))
cat("\nA vs C (size -> g only)\n"); print(lavTestLRT(cf$C, cf$A, method = "satorra.bentler.2001"))
cat("\nA_full vs C_full (size -> every indicator)\n"); print(lavTestLRT(cf$C_full, cf$A_full, method = "satorra.bentler.2001"))
p3 <- parameterEstimates(cf$C_full); print(p3[p3$op == "~" & p3$rhs %in% c("c24","c25"), c("lhs","rhs","est","se","pvalue")], digits = 3)
## 5. invariance sequences: size band and sensitivity exclusions
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
seqfit <- function(dd, lab) { tab <- t(sapply(names(eqs), function(k) { f <- try(cfa(m4, data = dd, group = "grp", estimator = "MLR", cluster = "cl", group.equal = eqs[[k]]), silent = TRUE)
    if (inherits(f, "try-error") || !lavInspect(f, "converged")) return(rep(NA, 5)); round(fitMeasures(f, fm[1:5]), 3) }))
  cat("\n---", lab, "| n:", paste(names(table(dd$grp)), table(dd$grp), collapse = "; "), "\n"); print(tab)
  cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), " Delta RMSEA:", round(diff(tab[, "rmsea.robust"]), 3), "\n") }
mkrole <- function(d) { d <- d[isP(d$type) | isF(d$type), ]; d$grp <- factor(ifelse(isP(d$type), "pretrained", "finetuned"), levels = c("pretrained","finetuned")); d$cl <- d$base; d }
mkcoh <- function(d) { d <- d[isP(d$type) & !is.na(d$upload_year), ]; d$grp <- factor(ifelse(d$upload_year <= 2023, "c2022_23", ifelse(d$upload_year == 2024, "c2024", "c2025"))); d$cl <- sub("/.*", "", d$model); d }
dup <- duplicated(round(x[, v], 4))
cat("\n=== Sensitivity: near-duplicate score vectors (identical to 4 dp):", sum(dup), "; clipped:", sum(x$any_clip == 1), "; merged flag:", sum(x$merged == "True"), "===\n")
xs <- list(baseline = x, no_clip = x[x$any_clip == 0, ], no_merged_flag = x[x$merged == "False", ], no_dup = x[!dup, ],
           all_three = x[x$any_clip == 0 & x$merged == "False" & !dup, ])
for (nm in names(xs)) { seqfit(mkrole(xs[[nm]]), paste("role |", nm)); seqfit(mkcoh(xs[[nm]]), paste("cohort |", nm)) }
band <- mkcoh(x); band <- band[band$params >= 0.5 & band$params <= 15, ]; seqfit(band, "cohort | size band 0.5-15B")
## 6. variance ratios (range-restriction check)
cat("\n=== Observed SD by group and ratio (cluster bootstrap 95% CI) ===\n")
vr <- function(d, ref, foc) { cls <- unique(d$cl)
  st <- function(dd) sapply(vl, function(j) sd(dd[dd$grp == foc, j]) / sd(dd[dd$grp == ref, j]))
  b <- replicate(1000, { pick <- sample(cls, length(cls), replace = TRUE); st(d[unlist(lapply(pick, function(z) which(d$cl == z))), ]) })
  rbind(sd_ref = sapply(vl, function(j) sd(d[d$grp == ref, j])), sd_foc = sapply(vl, function(j) sd(d[d$grp == foc, j])), ratio = st(d),
        lo = apply(b, 1, quantile, .025, na.rm = TRUE), hi = apply(b, 1, quantile, .975, na.rm = TRUE)) }
cat("\nfinetuned / pretrained\n"); print(round(vr(mkrole(x), "pretrained", "finetuned"), 3))
cat("\nc2024 / c2022_23\n"); print(round(vr(mkcoh(x), "c2022_23", "c2024"), 3))
cat("\nc2025 / c2022_23\n"); print(round(vr(mkcoh(x), "c2022_23", "c2025"), 3))
sink()
