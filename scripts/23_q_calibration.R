suppressMessages({library(lavaan); library(MASS)}); set.seed(20261114)
x <- read.csv("derived/hf_models_v2.csv"); v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained"); isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model); fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]; Vlam <- vcov(fp)[paste0("g=~", vl), paste0("g=~", vl)]
idx <- match(x$base, x$model); ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; G <- as.matrix(C[, vl]) - as.matrix(x[idx[ok], vl])
M <- as.matrix(aggregate(as.data.frame(G), by = list(C$base), mean)[, vl]); mbar <- colMeans(M)
Q <- function(M, lm) { V <- cov(M) / nrow(M); Vi <- solve(V); mb <- colMeans(M); mu <- as.numeric(t(lm) %*% Vi %*% mb / (t(lm) %*% Vi %*% lm)); r <- mb - lm * mu; as.numeric(t(r) %*% Vi %*% r) }
q0 <- Q(M, lam)
for (delta in c(0, 0.2, 0.5)) { Ms <- sweep(sweep(M, 2, mbar), 2, lam * delta, "+")
  qs <- replicate(5000, { s <- sample(nrow(Ms), nrow(Ms), replace = TRUE); tryCatch(Q(Ms[s, ], mvrnorm(1, lam, Vlam)), error = function(e) NA) })
  cat(sprintf("pure-g null delta %.1f: Q quantiles 95%% %.1f, 99%% %.1f, max %.1f; observed Q %.1f; empirical p = %.4f (n = %d)\n",
      delta, quantile(qs, .95, na.rm = TRUE), quantile(qs, .99, na.rm = TRUE), max(qs, na.rm = TRUE), q0, mean(qs >= q0, na.rm = TRUE), sum(!is.na(qs)))) }
