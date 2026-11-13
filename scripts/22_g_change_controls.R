suppressMessages({library(lavaan); library(MASS)})
# Strengthening analyses: g-change estimation with equivalence test, subgroup/depth replication,
# precision-duplicate negative control, pure-g positive control.
set.seed(20261113)
x <- read.csv("derived/hf_models_v2.csv"); v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
Vlam <- vcov(fp)[paste0("g=~", vl), paste0("g=~", vl)]
Theta <- lavInspect(fp, "est")$theta[vl, vl]; nu <- lavInspect(fp, "est")$nu[vl, 1]
# Bartlett factor-score weights from the pretrained model: w = (L' T^-1 L)^-1 L' T^-1
Ti <- solve(Theta); wB <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti); names(wB) <- vl
sdP <- sapply(P[, vl], sd)
idx <- match(x$base, x$model)
ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl]); colnames(G) <- vl

est <- function(G, cl, lm = lam) {
  L <- aggregate(as.data.frame(G), by = list(cl = cl), FUN = mean); M <- as.matrix(L[, vl]); n <- nrow(M)
  mbar <- colMeans(M); V <- cov(M) / n; Vi <- solve(V)
  mu <- as.numeric(t(lm) %*% Vi %*% mbar / (t(lm) %*% Vi %*% lm)); r <- mbar - lm * mu
  c(mu = mu, se_model = 1 / sqrt(as.numeric(t(lm) %*% Vi %*% lm)), Q = as.numeric(t(r) %*% Vi %*% r),
    bart = sum(wB * mbar), bart_se = sqrt(as.numeric(t(wB) %*% V %*% wB)), n = n) }
boot <- function(G, cl, R = 2000) {
  u <- unique(cl)
  t(replicate(R, { s <- sample(u, length(u), replace = TRUE)
    ii <- unlist(lapply(seq_along(s), function(j) which(cl == s[j]))); cc <- unlist(lapply(seq_along(s), function(j) rep(j, sum(cl == s[j]))))
    lm <- mvrnorm(1, lam, Vlam); e <- tryCatch(est(G[ii, , drop = FALSE], cc, lm)[c("mu","bart")], error = function(z) c(mu = NA, bart = NA)); e })) }
rep_block <- function(lab, G, cl, R = 2000) {
  nl <- length(unique(cl)); if (nl < 12) { cat(sprintf("\n[%s] pairs %d, lineages %d: too few lineages for a 6x6 lineage covariance; gains only\n", lab, nrow(G), nl)); g <- colMeans(G) / sdP; cat("  benchmark gains (pretrained SD, pair means):", paste(sprintf("%s %+.2f", v, g), collapse = "; "), "\n"); return(invisible(NULL)) }
  e <- est(G, cl); b <- boot(G, cl, R)
  cat(sprintf("\n  (usable bootstrap replicates: %d of %d)", sum(complete.cases(b)), nrow(b)))
  q90 <- apply(b, 2, quantile, c(.05, .95), na.rm = TRUE); q95 <- apply(b, 2, quantile, c(.025, .975), na.rm = TRUE)
  cat(sprintf("\n[%s] pairs %d, lineages %d\n", lab, nrow(G), e["n"]))
  cat(sprintf("  Q(5) = %.1f, p = %.2g\n", e["Q"], pchisq(e["Q"], 5, lower.tail = FALSE)))
  cat(sprintf("  g change, GLS projection: %.3f (model SE %.3f); bootstrap 95%% CI [%.3f, %.3f]; 90%% CI [%.3f, %.3f]\n", e["mu"], e["se_model"], q95[1,"mu"], q95[2,"mu"], q90[1,"mu"], q90[2,"mu"]))
  cat(sprintf("  g change, Bartlett scores:  %.3f (SE %.3f); bootstrap 95%% CI [%.3f, %.3f]; 90%% CI [%.3f, %.3f]\n", e["bart"], e["bart_se"], q95[1,"bart"], q95[2,"bart"], q90[1,"bart"], q90[2,"bart"]))
  cat(sprintf("  TOST: smallest symmetric bound with 90%% CI inside: GLS %.3f, Bartlett %.3f (pretrained g SD units)\n", max(abs(q90[, "mu"])), max(abs(q90[, "bart"]))))
  g <- colMeans(as.matrix(aggregate(as.data.frame(G), by = list(cl), mean)[, vl])) / sdP
  cat("  benchmark gains (pretrained SD):", paste(sprintf("%s %+.2f", v, g), collapse = "; "), "\n")
  invisible(c(e, q90[, "mu"], q90[, "bart"])) }

sink("out/p_strong_v1.txt")
cat("Bartlett weights (pretrained M4, var(g)=1):", paste(sprintf("%s %.3f", v, wB), collapse = "; "), "\n")
cat("Check: sum(wB * lam) =", round(sum(wB * lam), 4), "(must be 1)\n")
cat("\n=== 1. Strict pair set: equivalence test ===\n")
rep_block("strict, all", G, C$base)
cat("\n=== 2. Subgroups and depth ===\n")
for (tp in c("chat models (RLHF, DPO, IFT, ...)", "fine-tuned on domain-specific datasets")) { s <- C$type == tp; rep_block(tp, G[s, ], C$base[s]) }
for (dp in list(1, 2, 3:8)) { s <- C$gen %in% dp; rep_block(paste("depth", paste(range(dp), collapse = "-")), G[s, ], C$base[s]) }
cat("\nDose-response across depth (child-level regression on depth, lineage-clustered by lineage means of depth groups):\n")
dd <- data.frame(G, gen = pmin(C$gen, 3), base = C$base); fsd <- as.numeric(G %*% wB)
for (j in c(vl)) { f <- lm(dd[[j]] / sdP[j] ~ gen, data = dd); cat(sprintf("  %-11s slope per depth step %+.3f\n", j, coef(f)[2])) }
fit <- lm(fsd ~ gen, data = dd); cat(sprintf("  Bartlett g   slope per depth step %+.3f\n", coef(fit)[2]))
# cluster bootstrap of slopes
u <- unique(dd$base); bs <- t(replicate(2000, { s <- sample(u, length(u), replace = TRUE); ii <- unlist(lapply(s, function(z) which(dd$base == z)))
  c(sapply(vl, function(j) coef(lm(dd[[j]][ii] / sdP[j] ~ dd$gen[ii]))[2]), g = coef(lm(fsd[ii] ~ dd$gen[ii]))[2]) }))
print(round(apply(bs, 2, quantile, c(.025, .975), na.rm = TRUE), 3))
cat("\n=== 3. Negative control: bfloat16 vs float16 evaluations of the same model ===\n")
pdif <- read.csv("out/precision_diffs_v1.csv"); Gp <- as.matrix(pdif[, paste0(v, "_d")]); colnames(Gp) <- vl
rep_block("precision duplicates (fp16 - bf16), each model its own cluster", Gp, pdif$fullname)
cat("\n=== 4. Positive control: pure g gain injected into the real pair structure ===\n")
L <- aggregate(as.data.frame(G), by = list(cl = C$base), FUN = mean); M <- as.matrix(L[, vl]); mbar <- colMeans(M)
for (delta in c(0, 0.2, 0.5)) {
  Ms <- sweep(M, 2, mbar) ; Ms <- sweep(Ms, 2, lam * delta, "+")
  res <- t(replicate(1000, { s <- sample(nrow(Ms), nrow(Ms), replace = TRUE); lm <- mvrnorm(1, lam, Vlam); e <- est(Ms[s, ], seq_along(s), lm); e[c("mu","Q","bart")] }))
  cat(sprintf("  injected g gain %.1f SD: GLS mu mean %.3f (2.5-97.5%%: %.3f, %.3f); Bartlett mean %.3f; Q rejection rate at .05 = %.3f\n",
      delta, mean(res[, "mu"]), quantile(res[, "mu"], .025), quantile(res[, "mu"], .975), mean(res[, "bart"]), mean(pchisq(res[, "Q"], 5, lower.tail = FALSE) < .05)))
}
sink()
