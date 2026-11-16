suppressMessages({library(lavaan); library(mice)})
# Frontier replication, missing data and cohort split.
# (M1) what predicts a missing score; (M2) complete cases only; (M3) missing-not-at-random sensitivity: multiple imputation within each cohort with the
# imputed values shifted by delta standard deviations, each imputed data set tested by permutation; (S) cohort split at the median release date.
# args: [block: all|M|D-1|D-0.5|D0|D0.5|S] [outfile]
args <- commandArgs(TRUE); part <- if (length(args) >= 1) args[1] else "all"; outf <- if (length(args) >= 2) args[2] else "out/u_epoch_missing_v1.txt"
vl <- paste0(c("gpqa","aime","arc1","arc2","dtb"), "_lg"); m <- paste("g =~", paste(vl, collapse = " + "))
x <- read.csv("derived/epoch_wide_v1.csv"); x$grp <- factor(ifelse(x$year <= 2025, "le2025", "y2026")); x <- x[rowSums(!is.na(x[, vl])) >= 3, ]
eqs <- list("", "loadings", c("loadings","intercepts"))
stat <- function(d, miss) { f <- lapply(eqs, function(e) cfa(m, data = d, group = "grp", estimator = "MLR", missing = miss, group.equal = e))
  lr <- lavTestLRT(f[[1]], f[[2]], f[[3]], method = "satorra.bentler.2001"); c(metric = lr[2, "Chisq diff"], scalar = lr[3, "Chisq diff"]) }
perm <- function(d, miss, R, seed) { set.seed(seed); o <- stat(d, miss)
  p <- t(replicate(R, { dp <- d; dp$grp <- sample(dp$grp); tryCatch(stat(dp, miss), error = function(e) c(NA, NA)) })); p <- p[complete.cases(p), , drop = FALSE]
  c(o, p_metric = (1 + sum(p[, 1] >= o[1])) / (nrow(p) + 1), p_scalar = (1 + sum(p[, 2] >= o[2])) / (nrow(p) + 1), R = nrow(p)) }
sink(outf)
if (part %in% c("all","M")) {
  cat("=== (M1) Predictors of a missing score (logistic regression on cohort and the observed GPQA and AIME logits) ===\n")
  for (b in c("arc1_lg","arc2_lg","dtb_lg")) { d <- x[!is.na(x$gpqa_lg) & !is.na(x$aime_lg), ]; d$miss <- is.na(d[[b]])
    f <- glm(miss ~ grp + gpqa_lg + aime_lg, data = d, family = binomial); cf <- summary(f)$coefficients
    cat(sprintf("  %s missing %d of %d: %s\n", b, sum(d$miss), nrow(d), paste(sprintf("%s OR %.2f (p = %.2g)", rownames(cf)[-1], exp(cf[-1, 1]), cf[-1, 4]), collapse = "; "))) }
  cat("\n=== (M2) Complete cases (all five scores observed) ===\n")
  cc <- x[complete.cases(x[, vl]), ]; cat("  n =", nrow(cc), "(", paste(table(cc$grp), collapse = "/"), ")\n")
  r <- tryCatch(perm(cc, "listwise", 500, 31), error = function(e) { cat("  failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) cat(sprintf("  metric dchi2 %.1f, permutation p = %.3f; scalar dchi2 %.1f, permutation p = %.3f (%d permutations)\n", r[1], r[3], r[2], r[4], r[5])) }
for (delta in c(-1, -0.5, 0, 0.5)) { if (!(part %in% c("all", paste0("D", delta)))) next
  set.seed(3100 + round(10 * delta)); sdv <- sapply(x[, vl], sd, na.rm = TRUE)
  # Imputation within each cohort, so that the imputation model allows the cohorts different means and covariances (as the configural model does)
  imps <- lapply(levels(x$grp), function(g) mice(x[x$grp == g, vl], m = 10, method = "norm", printFlag = FALSE))
  res <- t(sapply(1:10, function(i) { d <- do.call(rbind, lapply(seq_along(imps), function(j) data.frame(complete(imps[[j]], i), grp = levels(x$grp)[j])))
    d$grp <- factor(d$grp, levels = levels(x$grp)); xo <- x[order(x$grp), ]
    for (b in vl) { mm <- is.na(xo[[b]]); d[[b]][mm] <- d[[b]][mm] + delta * sdv[b] }
    perm(d, "listwise", 200, 3200 + i) }))
  cat(sprintf("\n=== (M3) Imputed values shifted by %+.1f SD (10 imputations, 200 permutations each) ===\n  metric dchi2 median %.1f (range %.1f to %.1f); permutation p median %.3f, max %.3f\n  scalar dchi2 median %.1f; permutation p median %.3f, min %.3f\n",
    delta, median(res[, 1]), min(res[, 1]), max(res[, 1]), median(res[, 3]), max(res[, 3]), median(res[, 2]), median(res[, 4]), min(res[, 4]))) }
if (part %in% c("all","S")) {
  med <- median(as.Date(x$release)); d <- x; d$grp <- factor(ifelse(as.Date(d$release) <= med, "early", "late"))
  cat(sprintf("\n=== (S) Cohort split at the median release date (%s): n = %s ===\n", med, paste(table(d$grp), collapse = "/")))
  r <- perm(d, "ml", 500, 33); cat(sprintf("  metric dchi2 %.1f, permutation p = %.3f; scalar dchi2 %.1f, permutation p = %.3f\n", r[1], r[3], r[2], r[4])) }
sink()
