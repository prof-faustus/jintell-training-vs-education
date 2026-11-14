suppressMessages({library(lavaan); library(MASS)})
# Review tests of conditions that the main analyses would otherwise take for granted:
# (A) validity of the declared model types; (B) evaluation format (chat template) as a confound;
# (C) homogeneity of descendant-base pairs; (D) independence of lineages within organisations;
# (E) dependence of the g-change estimate on the linking constraint (Robitzsch & Luedtke, 2023);
# (F) per-benchmark gain intervals; (G) coverage of the base-selection explanation; (H) precision mismatch in pairs.
set.seed(20261114)
x <- read.csv("derived/hf_models_v2.csv"); e <- read.csv("derived/hf_extra_v1.csv"); x <- merge(x, e, by = "model", sort = FALSE)
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
Vlam <- vcov(fp)[paste0("g=~", vl), paste0("g=~", vl)]
Theta <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Theta)
wfun <- function(l) as.numeric(solve(t(l) %*% Ti %*% l) %*% t(l) %*% Ti)
sdP <- sapply(P[, vl], sd)
idx <- match(x$base, x$model)
ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl]); colnames(G) <- vl

# Estimators of the change in the common factor from a vector of mean gains m (logits), loadings l (pretrained g SD units)
estm <- function(M, l) { mbar <- colMeans(M); n <- nrow(M); V <- cov(M) / n; Vi <- solve(V)
  gls <- as.numeric(t(l) %*% Vi %*% mbar / (t(l) %*% Vi %*% l)); r <- mbar - l * gls
  c(bartlett = sum(wfun(l) * mbar), gls = gls, uls = sum(l * mbar) / sum(l^2), median = median(mbar / l),
    Q = as.numeric(t(r) %*% Vi %*% r), n = n) }
lmeans <- function(G, cl) as.matrix(aggregate(as.data.frame(G), by = list(cl = cl), FUN = mean)[, vl])
block <- function(lab, G, cl, R = 2000, bootcl = cl) {
  M <- lmeans(G, cl); if (nrow(M) < 12) { cat(sprintf("\n[%s] pairs %d, lineages %d: too few lineages\n", lab, nrow(G), nrow(M))); return(invisible(NULL)) }
  e0 <- estm(M, lam); u <- unique(bootcl)
  b <- t(replicate(R, { s <- sample(u, length(u), replace = TRUE)
    ii <- unlist(lapply(s, function(z) which(bootcl == z))); cc <- unlist(lapply(seq_along(s), function(j) paste(j, cl[bootcl == s[j]])))
    tryCatch(estm(lmeans(G[ii, , drop = FALSE], cc), mvrnorm(1, lam, Vlam))[1:4], error = function(z) rep(NA, 4)) }))
  q <- apply(b, 2, quantile, c(.025, .975), na.rm = TRUE)
  cat(sprintf("\n[%s] pairs %d, lineages %d, bootstrap units %d (usable %d of %d)\n", lab, nrow(G), nrow(M), length(u), sum(complete.cases(b)), R))
  cat(sprintf("  Q(5) = %.1f (chi2 p = %.2g)\n", e0["Q"], pchisq(e0["Q"], 5, lower.tail = FALSE)))
  c0 <- M[, 1] / lam[1] - as.numeric(M[, -1] %*% lam[-1]) / sum(lam[-1]^2); tt <- t.test(c0)
  cat(sprintf("  IFEval-implied minus rest-implied factor change: %+.3f, t(%d) = %.2f, p = %.2g\n", mean(c0), nrow(M) - 1, tt$statistic, tt$p.value))
  for (k in c("bartlett","gls","uls","median")) cat(sprintf("  %-9s %+.3f  95%% CI [%+.3f, %+.3f]\n", k, e0[k], q[1, k], q[2, k]))
  ci <- sapply(seq_along(vl), function(j) { t0 <- t.test(M[, j] / sdP[j]); c(t0$estimate, t0$conf.int) })
  cat("  benchmark gains (pretrained SD, t 95% CI across lineages):", paste(sprintf("%s %+.2f [%+.2f, %+.2f]", v, ci[1, ], ci[2, ], ci[3, ]), collapse = "; "), "\n")
  invisible(e0) }
fitseq <- function(dd, cl) { eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
  fits <- lapply(eqs, function(q) cfa(m4, data = dd, group = "grp", estimator = "MLR", cluster = cl, group.equal = q))
  tab <- t(sapply(fits, function(f) round(fitMeasures(f, c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr")), 3)))
  print(tab); cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), "\n")
  print(lavTestLRT(fits$configural, fits$metric, fits$scalar, fits$strict, method = "satorra.bentler.2001")) }

sink("out/r_review_tests_v1.txt")
cat("=== (A) Declared model type against the chat-template flag (all 4496 models) ===\n")
print(table(type = x$type, chat_template = x$chat_template))
cat(sprintf("Share with chat template: chat models %.3f; domain fine-tunes %.3f; pretrained %.3f; continuously pretrained %.3f\n",
  mean(x$chat_template[x$type == "chat models (RLHF, DPO, IFT, ...)"] == "True"), mean(x$chat_template[x$type == "fine-tuned on domain-specific datasets"] == "True"),
  mean(x$chat_template[x$type == "pretrained"] == "True"), mean(x$chat_template[x$type == "continuously pretrained"] == "True")))
cat("Pretrained models without upload date (excluded from the cohort contrast):", sum(isP(x$type) & is.na(x$upload_year)), "\n")

cat("\n=== (B) Evaluation format: models evaluated without a chat template only ===\n")
A <- x[(isP(x$type) | isF(x$type)) & x$chat_template == "False", ]; A$grp <- factor(ifelse(isP(A$type), "pretrained", "finetuned"), levels = c("pretrained","finetuned"))
cat("Training-role contrast, no chat template | n:", paste(table(A$grp), collapse = "/"), "| lineage clusters:", length(unique(A$base)), "\n"); fitseq(A, "base")
write.csv(x[(isP(x$type) | isF(x$type)) & x$chat_template == "False", names(read.csv("derived/hf_models_v2.csv", nrows = 1))], "derived/hf_models_v2_notemplate.csv", row.names = FALSE)
nt <- C$chat_template == "False" & B$chat_template == "False"
block("Strict pairs, descendant and base both without chat template", G[nt, ], C$base[nt])
block("Strict pairs, descendant with chat template", G[C$chat_template == "True", ], C$base[C$chat_template == "True"])

cat("\nDepth slopes (lineage-bootstrap 95% CI) by evaluation format of the descendant, pretrained SD per generation step\n")
for (tp in c("True","False")) { s0 <- C$chat_template == tp; d0 <- data.frame(G[s0, ] / rep(sdP, each = sum(s0)), gen = C$gen[s0], cl = C$base[s0])
  sl <- function(dd) sapply(vl, function(j) coef(lm(dd[[j]] ~ dd$gen))[2]); b0 <- sl(d0); u <- unique(d0$cl)
  bb <- t(replicate(1000, { s <- sample(u, length(u), replace = TRUE); dd <- do.call(rbind, lapply(s, function(z) d0[d0$cl == z, ])); tryCatch(sl(dd), error = function(z) rep(NA, 6)) }))
  q <- apply(bb, 2, quantile, c(.025, .975), na.rm = TRUE)
  cat(sprintf("  chat template %s (pairs %d, lineages %d, depths %s):", tp, sum(s0), length(u), paste(range(d0$gen), collapse = "-")), paste(sprintf("%s %+.2f [%+.2f, %+.2f]", v, b0, q[1, ], q[2, ]), collapse = "; "), "\n") }

cat("\n=== (C) Homogeneity of descendant-base pairs ===\n")
same_arch <- C$arch == B$arch; same_size <- abs(C$params / B$params - 1) <= .05; same_prec <- C$precision == B$precision; orig <- C$weight_type == "Original"
cat(sprintf("Of %d strict pairs: same architecture %d; parameter count within 5%% %d; same evaluation precision %d; full weights (not adapter) %d; all four %d\n",
  nrow(C), sum(same_arch), sum(same_size, na.rm = TRUE), sum(same_prec), sum(orig), sum(same_arch & same_size & same_prec & orig, na.rm = TRUE)))
h <- same_arch & same_size & same_prec & orig; h[is.na(h)] <- FALSE
block("Strict pairs with same architecture, size within 5%, same precision, full weights", G[h, ], C$base[h])

cat("\n=== (D) Organisation-level dependence ===\n")
Ar <- x[isP(x$type) | isF(x$type), ]; Ar$grp <- factor(ifelse(isP(Ar$type), "pretrained", "finetuned"), levels = c("pretrained","finetuned")); Ar$borg <- sub("/.*", "", Ar$base)
cat("Training-role contrast clustered by organisation of the base model | clusters:", length(unique(Ar$borg)), "\n"); fitseq(Ar, "borg")
block("Strict pairs, bootstrap resampling base-model organisations", G, C$base, bootcl = sub("/.*", "", C$base))

Ml <- lmeans(G, C$base); lin <- sort(unique(C$base))   # aggregate() orders lineages as sort() does
Mo <- as.matrix(aggregate(as.data.frame(Ml), by = list(o = sub("/.*", "", lin)), FUN = mean)[, vl])
qo <- estm(Mo, lam)["Q"]; Mc <- sweep(Mo, 2, colMeans(Mo))
qnull <- replicate(5000, { l <- mvrnorm(1, lam, Vlam); s <- sample(nrow(Mc), nrow(Mc), replace = TRUE); tryCatch(estm(Mc[s, ], l)["Q"], error = function(z) NA) })
cat(sprintf("Q computed on organisation means (%d organisations): Q(5) = %.1f; empirical pure-g null (organisation bootstrap, 5000 replicates): 95th %.1f, 99th %.1f, max %.1f, p = %.4f\n",
  nrow(Mo), qo, quantile(qnull, .95, na.rm = TRUE), quantile(qnull, .99, na.rm = TRUE), max(qnull, na.rm = TRUE), (1 + sum(qnull >= qo, na.rm = TRUE)) / (1 + sum(!is.na(qnull)))))

# One-df contrast fixed by the hypothesis (IFEval is the benchmark that instruction tuning targets): under a pure g gain the factor
# change implied by IFEval (gain / loading) equals the change implied by the other five (ULS on the five). Organisations are the units.
con <- function(Mx, l) Mx[, 1] / l[1] - as.numeric(Mx[, -1] %*% l[-1]) / sum(l[-1]^2)
for (lev in c("lineage","organisation")) { Mx <- if (lev == "lineage") Ml else Mo; c0 <- con(Mx, lam)
  bb <- replicate(5000, { s <- sample(nrow(Mx), nrow(Mx), replace = TRUE); mean(con(Mx[s, , drop = FALSE], mvrnorm(1, lam, Vlam))) })
  tt <- t.test(c0); cat(sprintf("IFEval-implied minus rest-implied factor change, %s units (n = %d): mean %+.3f; t(%d) = %.2f, p = %.2g; bootstrap 95%% CI [%+.3f, %+.3f]; units with a positive contrast %d of %d\n",
    lev, nrow(Mx), mean(c0), nrow(Mx) - 1, tt$statistic, tt$p.value, quantile(bb, .025), quantile(bb, .975), sum(c0 > 0), nrow(Mx))) }

cat("\n=== (E) Linking-constraint sensitivity (strict pairs) ===\n")
e0 <- block("Strict pairs, all", G, C$base)
M <- lmeans(G, C$base); mbar <- colMeans(M)
jk <- t(sapply(seq_along(vl), function(i) { l <- lam[-i]; m <- mbar[-i]; Tm <- Theta[-i, -i]; Tim <- solve(Tm)
  w <- as.numeric(solve(t(l) %*% Tim %*% l) %*% t(l) %*% Tim); c(bartlett = sum(w * m), uls = sum(l * m) / sum(l^2), median = median(m / l)) }))
rownames(jk) <- paste("without", v); print(round(jk, 3))
cat(sprintf("Leave-one-benchmark-out range: Bartlett [%.3f, %.3f]; ULS [%.3f, %.3f]; median [%.3f, %.3f]\n",
  min(jk[, 1]), max(jk[, 1]), min(jk[, 2]), max(jk[, 2]), min(jk[, 3]), max(jk[, 3])))
cat("Implied factor change per benchmark (mean gain / loading):", paste(sprintf("%s %+.3f", v, mbar / lam), collapse = "; "), "\n")

cat("\n=== (F) Per-benchmark within-lineage gains, pretrained SD units, t-based 95% CI across 48 lineages ===\n")
for (j in seq_along(vl)) { t0 <- t.test(M[, j] / sdP[j]); cat(sprintf("%-8s %+.3f [%+.3f, %+.3f]  p = %.2g\n", v[j], t0$estimate, t0$conf.int[1], t0$conf.int[2], t0$p.value)) }

cat("\n=== (G) Coverage of the base-selection explanation ===\n")
fs <- lavPredict(fp)[, 1]; names(fs) <- P$model; z <- (fs - mean(fs)) / sd(fs)
Fall <- x[isF(x$type), ]; hasb <- Fall$base %in% names(z)
cat(sprintf("Fine-tuned models (domain + chat): %d; with a pretrained base on the leaderboard: %d (%.1f%%)\n", nrow(Fall), sum(hasb), 100 * mean(hasb)))
w <- table(Fall$base[hasb]); zb <- z[names(w)]
cat(sprintf("Descendant-weighted mean factor score of those bases: %+.3f pretrained SD (unweighted mean of the %d bases %+.3f)\n", sum(w * zb) / sum(w), length(w), mean(zb)))

cat("\n=== (H) Evaluation precision of pair members ===\n"); print(table(descendant = C$precision, base = B$precision))
sink()
