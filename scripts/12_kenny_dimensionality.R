suppressMessages({library(lavaan); library(psych)})
args <- commandArgs(TRUE); x <- read.csv(args[1]); outf <- args[2]
ind <- c("ifeval_lg","bbh_lg","math_lg","gpqa_lg","musr_lg","mmlupro_lg")
fm <- c("chisq.scaled","df","pvalue.scaled","cfi.robust","rmsea.robust","srmr")
set.seed(20261113)  # identical split to b_select_logit.R
lin <- unique(x$base); h <- sample(c(1,2), length(lin), replace = TRUE); names(h) <- lin
x$half <- h[x$base]; D <- x[x$half == 1, ]; H <- x[x$half == 2, ]
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
sink(outf)
cat("Kenny, Kaniskan & McCoach (2015) Eq. 2: P(chi2_df > df[(N-1)c^2 + 1]) under a correct model\n")
for (cc in c(.05,.08,.084,.10)) cat(sprintf("df=6 N=2460 c=%.3f: threshold=%.1f  P=%.3g\n", cc, 6*((2460-1)*cc^2+1), pchisq(6*((2460-1)*cc^2+1), 6, lower.tail = FALSE)))
# extend M4 on discovery half by largest MI, down to df = 3
m <- m4; added <- c()
for (i in 1:3) { f <- cfa(m, data = D, estimator = "MLR", cluster = "base", std.lv = TRUE)
  mi <- modindices(f); mi <- mi[mi$op == "~~" & mi$lhs != mi$rhs, ]; mi <- mi[order(-mi$mi), ]
  t <- paste(mi$lhs[1], "~~", mi$rhs[1]); cat("\nDiscovery step", i, "add", t, "MI =", round(mi$mi[1], 1), "\n"); added <- c(added, t); m <- paste(m, "\n", t) }
cat("\n=== Holdout: M4 and each extension ===\n")
mm <- m4
for (i in 0:3) { if (i > 0) mm <- paste(mm, "\n", added[i])
  f <- cfa(mm, data = H, estimator = "MLR", cluster = "base", std.lv = TRUE)
  cat("\nModel M4 +", i, "\n"); print(round(fitMeasures(f, fm), 3))
  pe <- parameterEstimates(f, standardized = TRUE); pe <- pe[pe$op == "~~" & pe$lhs != pe$rhs, c("lhs","rhs","est","ci.lower","ci.upper","pvalue","std.all")]
  print(pe, digits = 3)
  ld <- standardizedSolution(f); print(ld[ld$op == "=~", c("rhs","est.std","ci.lower","ci.upper")], digits = 3) }
writeLines(mm, "out/m5_syntax_logit_v1.txt")
cat("\n=== Parallel analysis (Horn), full sample and holdout ===\n")
for (nm in c("full","holdout")) { dd <- if (nm == "full") x else H
  options(mc.cores = 1); set.seed(20261113); pa <- fa.parallel(dd[, ind], fm = "ml", fa = "fa", n.iter = 200, plot = FALSE)
  cat("\n", nm, ": n =", nrow(dd), " suggested factors =", pa$nfact, "\n observed eigen (FA):", round(pa$fa.values, 3), "\n simulated mean:", round(pa$fa.sim, 3), "\n") }
# invariance under the extended model (M4 + 1 and M4 + 3) for role and cohort
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
x$role <- ifelse(x$type %in% c("pretrained","continuously pretrained"), "pretrained",
          ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"), "finetuned", NA))
A <- x[!is.na(x$role), ]; A$grp <- factor(A$role, levels = c("pretrained","finetuned")); A$cl <- A$base
P <- x[x$type %in% c("pretrained","continuously pretrained") & !is.na(x$upload_year), ]
P$grp <- factor(ifelse(P$upload_year <= 2023, "c2022_23", ifelse(P$upload_year == 2024, "c2024", "c2025"))); P$cl <- sub("/.*", "", P$model)
for (k in c(1, 3)) { mk <- paste(c(m4, added[1:k]), collapse = "\n")
  for (nm in c("role","cohort")) { dd <- if (nm == "role") A else P
    tab <- t(sapply(names(eqs), function(e) { f <- cfa(mk, data = dd, group = "grp", estimator = "MLR", cluster = "cl", group.equal = eqs[[e]]); round(fitMeasures(f, fm[-3]), 3) }))
    cat("\n=== Invariance under M4 +", k, "|", nm, "===\n"); print(tab); cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), " Delta RMSEA:", round(diff(tab[, "rmsea.robust"]), 3), "\n") } }
sink()
