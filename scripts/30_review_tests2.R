suppressMessages({library(lavaan); library(MASS)})
# Second set of review tests: (T) evaluation time of descendant and base as a confound of within-lineage gains;
# (E) sensitivity to the clipping constant of the logit scale; (L) linearity of each benchmark in the common factor,
# including whether curvature combined with different latent levels could produce the training-role non-invariance.
set.seed(20261116)
x <- read.csv("derived/hf_models_v2.csv"); nz <- read.csv("derived/hf_norm_v1.csv"); stopifnot(identical(nz$model, x$model))
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
lg <- function(s, eps) { p <- pmin(pmax(s / 100, eps), 1 - eps); log(p / (1 - p)) }
prep <- function(eps) { y <- x; for (b in v) y[[paste0(b, "_lg")]] <- lg(nz[[paste0(b, "_norm")]], eps); y }   # normalised scores, as in 02_prep_logit.py
cat_ <- function(...) cat(sprintf(...))
cat("Check: recomputed logits at eps = .005 equal the analysed logits:", isTRUE(all.equal(as.matrix(prep(.005)[, vl]), as.matrix(x[, vl]), tolerance = 1e-9)), "\n")

core <- function(y) {   # loadings, pairs and the within-lineage estimators for a data set y
  P <- y[isP(y$type), ]; P$org <- sub("/.*", "", P$model)
  fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
  pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
  Th <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti)
  idx <- match(y$base, y$model)
  ok <- y$gen >= 1 & !is.na(idx) & y$model != y$base & isF(y$type) & y$merged == "False" & isP(y$type[idx]); ok[is.na(ok)] <- FALSE
  C <- y[ok, ]; B <- y[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
  list(P = P, fp = fp, lam = lam, w = w, C = C, B = B, G = G, sdP = sapply(P[, vl], sd)) }
lmean <- function(G, cl) as.matrix(aggregate(as.data.frame(G), by = list(cl = cl), FUN = mean)[, -1])
contrast <- function(M, lam) M[, 1] / lam[1] - as.numeric(M[, -1] %*% lam[-1]) / sum(lam[-1]^2)
summ <- function(k, sel, lab) { G <- k$G[sel, , drop = FALSE]; cl <- k$C$base[sel]; M <- lmean(G, cl)
  bart <- as.numeric(M %*% k$w); u <- unique(cl)
  bb <- replicate(2000, { s <- sample(seq_len(nrow(M)), nrow(M), replace = TRUE); mean(bart[s]) })
  cc <- contrast(M, k$lam); tt <- t.test(cc)
  cat_("  [%s] pairs %d, lineages %d: Bartlett factor change %+.3f [%+.3f, %+.3f]; IFEval gain %+.2f SD; IFEval contrast %+.2f, t(%d) = %.2f, p = %.2g\n",
    lab, nrow(G), nrow(M), mean(bart), quantile(bb, .025), quantile(bb, .975), mean(M[, 1]) / k$sdP[1], mean(cc), nrow(M) - 1, tt$statistic, tt$p.value) }
k0 <- core(x)

sink("out/t_review_tests2_v1.txt")
cat("=== (T) Evaluation time of descendant and base ===\n")
gap <- as.numeric(as.Date(substr(k0$C$submit, 1, 10)) - as.Date(substr(k0$B$submit, 1, 10))) / 365.25
cat_("Gap (descendant minus base submission date, years): median %.2f, IQR %.2f to %.2f, range %.2f to %.2f; descendant evaluated later in %.1f%% of pairs\n",
  median(gap), quantile(gap, .25), quantile(gap, .75), min(gap), max(gap), 100 * mean(gap > 0))
pp <- data.frame(k0$G / rep(k0$sdP, each = nrow(k0$G)), bart = as.numeric(k0$G %*% k0$w), con = contrast(k0$G, k0$lam), gap = gap, cl = k0$C$base)
u <- unique(pp$cl); outc <- c(vl, "bart", "con")
sl <- function(d) sapply(outc, function(o) coef(lm(d[[o]] ~ d$gap))[2])
b0 <- sl(pp); bb <- t(replicate(2000, { s <- sample(u, length(u), replace = TRUE); d <- do.call(rbind, lapply(s, function(z) pp[pp$cl == z, ])); sl(d) }))
q <- apply(bb, 2, quantile, c(.025, .975))
cat("Slope of each pair gain on the evaluation gap (per year), lineage-bootstrap 95% CI:\n")
for (j in seq_along(outc)) cat_("  %-10s %+.3f [%+.3f, %+.3f]\n", c(v, "Bartlett", "IFEval contrast")[j], b0[j], q[1, j], q[2, j])
summ(k0, abs(gap) <= 90 / 365.25, "pairs evaluated within 90 days of each other")
summ(k0, abs(gap) <= 30 / 365.25, "pairs evaluated within 30 days of each other")

cat("\n=== (E) Clipping constant of the logit scale ===\n")
inv <- function(d, cl) { fits <- lapply(list("", "loadings", c("loadings","intercepts")), function(q) cfa(m4, data = d, group = "grp", estimator = "MLR", cluster = cl, group.equal = q))
  lr <- lavTestLRT(fits[[1]], fits[[2]], fits[[3]], method = "satorra.bentler.2001"); cf <- sapply(fits, fitMeasures, "cfi.robust")
  sprintf("metric dchi2 %.1f, dCFI %+.3f; scalar dchi2 %.1f, dCFI %+.3f", lr[2, "Chisq diff"], cf[2] - cf[1], lr[3, "Chisq diff"], cf[3] - cf[2]) }
for (eps in c(.001, .005, .01, .02)) { y <- prep(eps)
  A <- y[isP(y$type) | isF(y$type), ]; A$grp <- factor(ifelse(isP(A$type), "pretrained", "finetuned"), levels = c("pretrained","finetuned"))
  Pc <- y[isP(y$type) & !is.na(y$upload_year), ]; Pc$grp <- cut(Pc$upload_year, c(2021, 2023, 2024, 2025), labels = c("c2022_23","c2024","c2025")); Pc$org <- sub("/.*", "", Pc$model)
  cat_("eps = %.3f\n  training role: %s\n  cohorts:       %s\n", eps, inv(A, "base"), inv(Pc, "org"))
  summ(core(y), rep(TRUE, nrow(k0$C)), sprintf("within lineages, eps %.3f", eps)) }

cat("\n=== (L) Linearity of each benchmark in the common factor ===\n")
# Rest score: loading-weighted mean of the other five standardised benchmarks (pretrained scaling). Quadratic and cubic terms test curvature;
# local slopes at the 10th and 90th percentiles of the rest score show how much the effective loading changes over the observed range.
A <- x[isP(x$type) | isF(x$type), ]; A$grp <- ifelse(isP(A$type), "pretrained", "finetuned")
P <- A[A$grp == "pretrained", ]; mu <- colMeans(P[, vl]); sdv <- sapply(P[, vl], sd); Z <- sweep(sweep(as.matrix(A[, vl]), 2, mu), 2, sdv, "/")
lin <- function(d, i) { r <- as.numeric(d$Z[, -i] %*% k0$lam[-i]) / sum(k0$lam[-i]); y <- d$Z[, i]
  f1 <- lm(y ~ r); f3 <- lm(y ~ r + I(r^2) + I(r^3)); qs <- quantile(r, c(.1, .9)); cf <- coef(f3)
  slope <- function(z) cf[2] + 2 * cf[3] * z + 3 * cf[4] * z^2
  c(dR2 = summary(f3)$r.squared - summary(f1)$r.squared, s10 = slope(qs[1]), s90 = slope(qs[2]), lin = coef(f1)[2]) }
for (g in c("pretrained","finetuned")) { s0 <- A$grp == g; d <- list(Z = Z[s0, ]); cl <- A$base[s0]; uu <- unique(cl)
  cat_("%s (n = %d):\n", g, sum(s0))
  for (i in seq_along(vl)) { e0 <- lin(d, i)
    bb <- t(replicate(500, { s <- sample(uu, length(uu), replace = TRUE); ii <- unlist(lapply(s, function(z) which(cl == z))); lin(list(Z = d$Z[ii, ]), i) }))
    q <- apply(bb, 2, quantile, c(.025, .975))
    cat_("  %-8s linear slope %.2f; R2 gain from curvature %.3f [%.3f, %.3f]; local slope at 10th / 90th percentile %.2f / %.2f (ratio %.2f [%.2f, %.2f])\n",
      v[i], e0["lin.r"], e0["dR2"], q[1, "dR2"], q[2, "dR2"], e0["s10.r"], e0["s90.r"], e0["s90.r"] / e0["s10.r"],
      quantile(bb[, "s90.r"] / bb[, "s10.r"], .025), quantile(bb[, "s90.r"] / bb[, "s10.r"], .975)) } }
# Common support: training-role invariance on the models whose rest score on all six benchmarks lies inside the range shared by both groups
r_all <- as.numeric(Z %*% k0$lam) / sum(k0$lam); lo <- max(tapply(r_all, A$grp, quantile, .05)); hi <- min(tapply(r_all, A$grp, quantile, .95))
Acs <- A[r_all >= lo & r_all <= hi, ]; Acs$grp <- factor(Acs$grp, levels = c("pretrained","finetuned"))
dl <- tapply(r_all[r_all >= lo & r_all <= hi], Acs$grp, mean)
cat_("\nCommon support (overall score between %.2f and %.2f): n = %s; group difference in mean overall score %.2f SD (full sample %.2f)\n  training role: %s\n",
  lo, hi, paste(table(Acs$grp), collapse = "/"), diff(dl), diff(tapply(r_all, factor(A$grp, levels = c("pretrained","finetuned")), mean)), inv(Acs, "base"))

cat("\n=== (L2) Curvature combined with different latent levels ===\n")
# Matched samples: groups given the same distribution of the overall score, so a shared curved measurement function cannot
# produce different linear loadings. Greedy nearest-neighbour matching without replacement, caliper 0.05 SD of the overall score.
set.seed(20261117); cal <- .05 * sd(r_all)
nn <- function(ref, pool) { used <- rep(FALSE, length(pool)); out <- rep(NA_integer_, length(ref))   # 1:1 match for each reference unit, NA if none within the caliper
  for (i in sample(seq_along(ref))) { d <- abs(pool - ref[i]); d[used] <- Inf; m <- which.min(d); if (d[m] <= cal) { used[m] <- TRUE; out[i] <- m } }
  out }
iP <- which(A$grp == "pretrained"); iF <- which(A$grp == "finetuned"); mt <- nn(r_all[iP], r_all[iF]); ok <- !is.na(mt)
Am <- A[c(iP[ok], iF[mt[ok]]), ]; Am$grp <- factor(Am$grp, levels = c("pretrained","finetuned")); rm_ <- r_all[c(iP[ok], iF[mt[ok]])]
cat_("Training role, matched pairs (1:1, unmatched pretrained models dropped): n = %s; overall score mean %s, SD %s\n  %s\n", paste(table(Am$grp), collapse = "/"),
  paste(sprintf("%.3f", tapply(rm_, Am$grp, mean)), collapse = " vs "), paste(sprintf("%.3f", tapply(rm_, Am$grp, sd)), collapse = " vs "), inv(Am, "base"))
write.csv(x[match(Am$model, x$model), ], "derived/hf_models_v2_matched_role.csv", row.names = FALSE)
Pc <- x[isP(x$type) & !is.na(x$upload_year), ]; Pc$coh <- cut(Pc$upload_year, c(2021, 2023, 2024, 2025), labels = c("c2022_23","c2024","c2025"))
rP <- r_all[match(Pc$model, A$model)]; i0 <- which(Pc$coh == "c2022_23"); i1 <- which(Pc$coh == "c2024"); i2 <- which(Pc$coh == "c2025")
a1 <- nn(rP[i0], rP[i1]); a2 <- nn(rP[i0], rP[i2]); ok <- !is.na(a1) & !is.na(a2)
keep <- c(i0[ok], i1[a1[ok]], i2[a2[ok]]); Pm <- Pc[keep, ]; Pm$grp <- Pm$coh; Pm$org <- sub("/.*", "", Pm$model)
cat_("Cohorts, matched triplets (one 2024 and one 2025 model per 2022-23 model, unmatched dropped): n = %s; overall score mean %s, SD %s\n  %s\n", paste(table(Pm$grp), collapse = "/"),
  paste(sprintf("%.3f", tapply(rP[keep], Pm$grp, mean)), collapse = " / "), paste(sprintf("%.3f", tapply(rP[keep], Pm$grp, sd)), collapse = " / "), inv(Pm, "org"))
write.csv(x[match(Pm$model, x$model), ], "derived/hf_models_v2_matched_cohort.csv", row.names = FALSE)
s2 <- Pc$coh %in% c("c2022_23","c2025"); P2 <- Pc[s2, ]; P2$grp <- droplevels(P2$coh); P2$org <- sub("/.*", "", P2$model)
cat_("Cohorts 2022-23 and 2025 only (similar overall level, all models): n = %s; overall score mean %s, SD %s\n  %s\n", paste(table(P2$grp), collapse = "/"),
  paste(sprintf("%.3f", tapply(rP[s2], P2$grp, mean)), collapse = " / "), paste(sprintf("%.3f", tapply(rP[s2], P2$grp, sd)), collapse = " / "), inv(P2, "org"))
write.csv(x[match(P2$model, x$model), ], "derived/hf_models_v2_cohort_2223_25.csv", row.names = FALSE)

# Local-slope IFEval contrast: loadings replaced by each benchmark's local slope (cubic fit in pretrained models) at the base model's position
Zp <- Z[A$grp == "pretrained", ]; fits <- lapply(seq_along(vl), function(i) { r <- as.numeric(Zp[, -i] %*% k0$lam[-i]) / sum(k0$lam[-i]); lm(Zp[, i] ~ r + I(r^2) + I(r^3)) })
zB <- sweep(sweep(as.matrix(k0$B[, vl]), 2, mu), 2, sdv, "/"); gz <- k0$G / rep(sdv, each = nrow(k0$G))
Ls <- sapply(seq_along(vl), function(i) { r <- as.numeric(zB[, -i] %*% k0$lam[-i]) / sum(k0$lam[-i]); cf <- coef(fits[[i]]); cf[2] + 2 * cf[3] * r + 3 * cf[4] * r^2 })
cpair <- gz[, 1] / Ls[, 1] - rowSums(Ls[, -1] * gz[, -1]) / rowSums(Ls[, -1]^2)
cat_("Local slopes at the base models (median across pairs): %s; share of pairs with any non-positive local slope %.3f\n", paste(sprintf("%s %.2f", v, apply(Ls, 2, median)), collapse = "; "), mean(apply(Ls, 1, min) <= 0))
okp <- Ls[, 1] > .1
for (lev in c("lineage","organisation")) { cl <- if (lev == "lineage") k0$C$base[okp] else sub("/.*", "", k0$C$base[okp])
  cm <- tapply(tapply(cpair[okp], k0$C$base[okp], mean), if (lev == "lineage") names(tapply(cpair[okp], k0$C$base[okp], mean)) else sub("/.*", "", names(tapply(cpair[okp], k0$C$base[okp], mean))), mean)
  tt <- t.test(cm); cat_("Local-slope IFEval contrast, %s units (n = %d, pairs with IFEval local slope > 0.1: %d): %+.2f, t(%d) = %.2f, p = %.2g\n", lev, length(cm), sum(okp), mean(cm), length(cm) - 1, tt$statistic, tt$p.value) }
sink()
