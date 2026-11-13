suppressMessages({library(lavaan); library(psych)})
x <- read.csv("derived/hf_models_v2.csv"); ind <- c("ifeval_lg","bbh_lg","math_lg","gpqa_lg","musr_lg","mmlupro_lg")
set.seed(20261113); lin <- unique(x$base); h <- sample(c(1,2), length(lin), replace = TRUE); names(h) <- lin
x$half <- h[x$base]; D <- x[x$half == 1, ]; H <- x[x$half == 2, ]
fm <- c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr","bic")
sink("out/g2_twofactor_v1.txt")
for (nm in c("discovery","holdout")) { dd <- if (nm == "discovery") D else H
  e <- fa(dd[, ind], nfactors = 2, fm = "ml", rotate = "oblimin"); cat("\nEFA 2 factors (ML, oblimin) |", nm, "\n"); print(round(unclass(e$loadings), 3)); cat("factor r:", round(e$Phi[1,2], 3), "\n") }
# CFA translation of the discovery EFA: indicators to the factor with the largest discovery loading; cross-loadings > .30 kept
e <- fa(D[, ind], nfactors = 2, fm = "ml", rotate = "oblimin"); L <- unclass(e$loadings)
f1 <- ind[abs(L[,1]) > .30]; f2 <- ind[abs(L[,2]) > .30]
m2 <- paste0("F1 =~ ", paste(f1, collapse = "+"), "\nF2 =~ ", paste(f2, collapse = "+"))
cat("\nCFA 2-factor syntax from discovery EFA:\n", m2, "\n"); writeLines(m2, "out/m2f_syntax_v1.txt")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
for (mod in list(c("M4", m4), c("2F", m2))) { f <- try(cfa(mod[2], data = H, estimator = "MLR", cluster = "base", std.lv = TRUE), silent = TRUE)
  cat("\nHoldout", mod[1], "\n"); if (inherits(f, "try-error")) { cat("failed\n"); next }; print(round(fitMeasures(f, fm), 3))
  if (mod[1] == "2F") print(standardizedSolution(f)[standardizedSolution(f)$op %in% c("=~","~~") & standardizedSolution(f)$lhs %in% c("F1","F2"), c("lhs","op","rhs","est.std","ci.lower","ci.upper")], digits = 3) }
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
x$role <- ifelse(x$type %in% c("pretrained","continuously pretrained"), "pretrained", ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"), "finetuned", NA))
A <- x[!is.na(x$role), ]; A$grp <- factor(A$role, levels = c("pretrained","finetuned")); A$cl <- A$base
P <- x[x$type %in% c("pretrained","continuously pretrained") & !is.na(x$upload_year), ]
P$grp <- factor(ifelse(P$upload_year <= 2023, "c2022_23", ifelse(P$upload_year == 2024, "c2024", "c2025"))); P$cl <- sub("/.*", "", P$model)
for (nm in c("role","cohort")) { dd <- if (nm == "role") A else P
  tab <- t(sapply(names(eqs), function(k) { f <- try(cfa(m2, data = dd, group = "grp", estimator = "MLR", cluster = "cl", group.equal = eqs[[k]]), silent = TRUE)
    if (inherits(f, "try-error") || !lavInspect(f, "converged")) return(rep(NA, 5)); round(fitMeasures(f, fm[-6]), 3) }))
  cat("\n=== 2-factor invariance |", nm, "===\n"); print(tab); cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), "\n") }
sink()

# Admissibility of the two-factor model in the cohort contrast
sink("out/g2_twofactor_cohort_admissibility_v1.txt")
for (e in list("", "loadings", c("loadings","intercepts"))) {
  f <- cfa(m2, data = P, group = "grp", estimator = "MLR", cluster = "cl", group.equal = e)
  pe <- parameterEstimates(f); neg <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0, c("group","lhs","est")]
  cat(paste(e, collapse = "+"), "| converged", lavInspect(f, "converged"), "| negative variances:", nrow(neg), "\n"); if (nrow(neg)) print(neg)
  cat("factor correlations by group:", round(sapply(lavInspect(f, "cov.lv"), function(z) z[1, 2] / sqrt(z[1, 1] * z[2, 2])), 3), "\n")
  ev <- min(eigen(vcov(f), symmetric = TRUE, only.values = TRUE)$values); cat("smallest eigenvalue of vcov:", ev, "\n") }
sink()
