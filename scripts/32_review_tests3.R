suppressMessages(library(lavaan))
# Third set of review tests: (S) selection of which descendants reach the leaderboard; (H) multiplicity across the confirmatory tests.
set.seed(20261118)
x <- read.csv("derived/hf_models_v2.csv")
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
Th <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti); sdP <- sapply(P[, vl], sd)
idx <- match(x$base, x$model)
ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
lmean <- function(G, cl) as.matrix(aggregate(as.data.frame(G), by = list(cl = cl), FUN = mean)[, -1])
sink("out/v_review_tests3_v1.txt")
cat("=== (S) Selection of descendants: releases that do not depend on the leaderboard ===\n")
first <- sub("/.*", "", C$model) == sub("/.*", "", C$base); off <- C$official %in% c("True", TRUE)
for (s in list(list("first-party releases (descendant from the base model's own organisation)", first),
               list("official-provider descendants", off),
               list("community descendants (not official, other organisation)", !off & !first))) {
  sel <- s[[2]]; M <- lmean(G[sel, , drop = FALSE], C$base[sel]); n <- nrow(M)
  if (n < 5) { cat(sprintf("[%s] lineages %d: too few\n", s[[1]], n)); next }
  bart <- as.numeric(M %*% w); bb <- replicate(2000, mean(bart[sample(n, n, replace = TRUE)]))
  cc <- M[, 1] / lam[1] - as.numeric(M[, -1] %*% lam[-1]) / sum(lam[-1]^2); tt <- t.test(cc)
  gi <- sapply(1:6, function(j) { t0 <- t.test(M[, j] / sdP[j]); c(t0$estimate, t0$conf.int) })
  cat(sprintf("[%s] pairs %d, lineages %d\n  Bartlett factor change %+.3f, lineage-bootstrap 95%% CI [%+.3f, %+.3f]\n  IFEval contrast %+.2f, t(%d) = %.2f, p = %.2g\n  gains (pretrained SD, 95%% CI): %s\n",
    s[[1]], sum(sel), n, mean(bart), quantile(bb, .025), quantile(bb, .975), mean(cc), n - 1, tt$statistic, tt$p.value,
    paste(sprintf("%s %+.2f [%+.2f, %+.2f]", v, gi[1, ], gi[2, ], gi[3, ]), collapse = "; "))) }

cat("\n=== (H) Holm adjustment across the five confirmatory tests ===\n")
# Empirical p = (1 + replicates at least as extreme) / (1 + replicates); for the cohorts the most conservative null (ext) is used.
ep <- function(txt, reps, col) { o <- as.numeric(strsplit(trimws(sub("^observed", "", grep("^observed", readLines(txt), value = TRUE)[1])), " +")[[1]])
  r <- read.csv(reps); k <- match(col, names(r)); (1 + sum(r[[col]] <= o[k], na.rm = TRUE)) / (1 + sum(!is.na(r[[col]]))) }
pr <- c(role_metric = ep("out/e_mcnull_role_exact_v1.txt", "out/e_mcnull_role_exact_v1_reps.csv", "dCFI_m"),
        role_scalar = ep("out/e_mcnull_role_exact_v1.txt", "out/e_mcnull_role_exact_v1_reps.csv", "dCFI_s"),
        cohort_metric = ep("out/e_mcnull_cohort_ext_v1.txt", "out/e_mcnull_cohort_ext_v1_reps.csv", "dCFI_m"),
        cohort_scalar = ep("out/e_mcnull_cohort_ext_v1.txt", "out/e_mcnull_cohort_ext_v1_reps.csv", "dCFI_s"),
        ifeval_contrast_org = as.numeric(sub(".*p = ([0-9.e-]+);.*", "\\1", grep("organisation units", readLines("out/r_review_tests_v1.txt"), value = TRUE)[1])))
print(round(rbind(p = pr, holm = p.adjust(pr, "holm")), 4))
sink()
