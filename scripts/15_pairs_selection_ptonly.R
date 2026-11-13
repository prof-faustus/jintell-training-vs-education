suppressMessages({library(lavaan)})
x <- read.csv("derived/hf_models_v2.csv"); v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained"); isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
sink("out/j_extra_v1.txt")
idx <- match(x$base, x$model); ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
L <- aggregate(G, by = list(base = C$base), FUN = mean)
cat("Strict pairs: lineage-level one-sample t-tests (df = n_lineages - 1)\n")
for (j in vl) { t <- t.test(L[[j]]); cat(sprintf("%-11s gain %.3f  t(%d) = %.2f  p = %.2g\n", j, mean(L[[j]]), t$parameter, t$statistic, t$p.value)) }
P <- x[isP(x$type), ]; sdP <- sapply(P[, vl], sd)
ub <- unique(B$model); Bu <- x[match(ub, x$model), ]; others <- P[!(P$model %in% ub), ]
cat("\nBase selection: bases used by strict fine-tunes (n =", length(ub), ") vs other pretrained models (n =", nrow(others), "), SD units of the pretrained pool; Welch t\n")
for (j in vl) { t <- t.test(Bu[[j]], others[[j]]); cat(sprintf("%-11s diff %.3f SD  t = %.2f  p = %.2g\n", j, (mean(Bu[[j]]) - mean(others[[j]])) / sdP[j], t$statistic, t$p.value)) }
fp <- cfa(m4, data = P, estimator = "MLR", std.lv = TRUE); fsc <- lavPredict(fp)[, 1]; names(fsc) <- P$model
t <- t.test(fsc[ub], fsc[!(names(fsc) %in% ub)]); cat(sprintf("factor score: bases %.3f vs others %.3f, t = %.2f, p = %.2g\n", mean(fsc[ub]), mean(fsc[!(names(fsc) %in% ub)]), t$statistic, t$p.value))
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
fm <- c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr")
A <- x[x$type == "pretrained" | isF(x$type), ]; A$grp <- factor(ifelse(A$type == "pretrained", "pretrained", "finetuned"), levels = c("pretrained","finetuned"))
fits <- lapply(eqs, function(e) cfa(m4, data = A, group = "grp", estimator = "MLR", cluster = "base", group.equal = e))
tab <- t(sapply(fits, function(f) round(fitMeasures(f, fm), 3))); cat("\nRole sequence excluding continuously pretrained | n:", paste(table(A$grp), collapse = "/"), "\n"); print(tab); cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), "\n")
print(lavTestLRT(fits$configural, fits$metric, fits$scalar, fits$strict, method = "satorra.bentler.2001"))
Pc <- x[isP(x$type) & !is.na(x$upload_year), ]; cat("\nCohort sample orgs:", length(unique(sub("/.*", "", Pc$model))), " role clusters:", length(unique(x$base[isP(x$type) | isF(x$type)])), "\n")
sink()

# Base-selection gap in factor-score standard deviations
sink("out/j_extra_v1.txt", append = TRUE)
cat(sprintf("\nFactor-score SD %.3f; gap in SD units %.3f\n", sd(fsc), (mean(fsc[ub]) - mean(fsc[!(names(fsc) %in% ub)])) / sd(fsc)))
sink()
