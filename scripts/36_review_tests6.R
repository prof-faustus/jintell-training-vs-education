suppressMessages(library(lavaan))
# Sixth set of review tests:
# (D) decomposition of the Bartlett factor change into benchmark contributions, weights per logit and per pretrained SD;
# (S) metric dependence of the linking estimators: unweighted least squares in the logit metric and in the
#     pretrained-SD metric (indicators standardised), with calibrated upper bounds; Bartlett, GLS and median are scale-invariant;
# (L) leave-one-organisation-out stability of every estimator and of the IFEval contrast;
# (F) scoring-format contrast: IFEval-implied against MATH-implied factor change (both scored from generated text);
# (K) Kenny, Kaniskan & McCoach (2015) Eq. 2 with the uncorrected ML RMSEA;
# (T) median submission date in the training-role sample and in all models.
# args: [outfile] [seed]
args <- commandArgs(TRUE); outf <- if (length(args) >= 1) args[1] else "out/x_review_tests6_v1.txt"
seed <- if (length(args) >= 2) as.integer(args[2]) else 20261121L
x <- read.csv("derived/hf_models_v2.csv")
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
sink(outf)
set.seed(seed); cat("Seed:", seed, "\n")
P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
Th <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti)
sdP <- sapply(P[, vl], sd)
idx <- match(x$base, x$model)
ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
Ml <- as.matrix(aggregate(as.data.frame(G), by = list(cl = C$base), FUN = mean)[, -1]); lin <- sort(unique(C$base))
og <- sub("/.*", "", lin)
ag <- aggregate(as.data.frame(Ml), by = list(o = og), FUN = mean); Mo <- as.matrix(ag[, -1]); rownames(Mo) <- ag$o
lamS <- lam / sdP   # loadings in the pretrained-SD metric of each benchmark
ests <- function(M) { mb <- colMeans(M); V <- cov(M) / nrow(M); Vi <- tryCatch(solve(V), error = function(e) NULL)
  gls <- if (is.null(Vi)) NA else as.numeric(t(lam) %*% Vi %*% mb / (t(lam) %*% Vi %*% lam))
  mS <- mb / sdP
  c(bartlett = sum(w * mb), gls = gls, uls_logit = sum(lam * mb) / sum(lam^2), uls_sd = sum(lamS * mS) / sum(lamS^2), median = median(mb / lam)) }
ifc <- function(M, metric = c("logit","sd")) { metric <- match.arg(metric)
  l <- if (metric == "logit") lam else lamS; MM <- if (metric == "logit") M else sweep(M, 2, sdP, "/")
  MM[, 1] / l[1] - as.numeric(MM[, -1] %*% l[-1]) / sum(l[-1]^2) }

cat("=== (D) Decomposition of the Bartlett factor change (strict pair set, 48 lineages, lineage-weighted) ===\n")
mb <- colMeans(Ml); mo <- colMeans(Mo)
D <- data.frame(benchmark = v, gain_logit = mb, gain_sd = mb / sdP, w_per_logit = w, w_per_sd = w * sdP,
  contribution = w * mb, contribution_org = w * mo, loading_std = lam / sqrt(lam^2 + diag(Th)))
print(format(D, digits = 3), row.names = FALSE)
cat(sprintf("Sum of contributions: lineage-weighted %+.4f; organisation-weighted %+.4f\n", sum(D$contribution), sum(D$contribution_org)))
cat(sprintf("Share of the lineage-weighted estimate carried by BBH and MMLU-Pro: %.3f\n", sum(D$contribution[c(2, 6)]) / sum(D$contribution)))
cat(sprintf("Weights per pretrained SD sum to %.3f (the factor-score metric); BBH + MMLU-Pro share of that sum %.3f\n", sum(D$w_per_sd), sum(D$w_per_sd[c(2, 6)]) / sum(D$w_per_sd)))

cat("\n=== (S) Metric dependence of the linking estimators ===\n")
cat("Rescaling every indicator by a constant leaves Bartlett, GLS and median estimates unchanged; ULS changes. Check on the data:\n")
sc <- runif(6, .5, 2); lam0 <- lam; sd0 <- sdP
lam <- lam0 * sc; Tsc <- diag(sc) %*% Th %*% diag(sc); Ti <- solve(Tsc); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti); sdP <- sd0 * sc
e_sc <- ests(sweep(Ml, 2, sc, "*")); lam <- lam0; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti); sdP <- sd0
e0l <- ests(Ml); e0o <- ests(Mo)
cat("  lineage-weighted, original scale:   ", paste(sprintf("%s %+.4f", names(e0l), e0l), collapse = "; "), "\n")
cat("  lineage-weighted, indicators x random constants:", paste(sprintf("%s %+.4f", names(e_sc), e_sc), collapse = "; "), "\n")
cat("  (uls_sd is defined on standardised indicators and is therefore unchanged by the rescaling; uls_logit is not)\n")
cat("  organisation-weighted:", paste(sprintf("%s %+.4f", names(e0o), e0o), collapse = "; "), "\n")
cat(sprintf("  ULS weights per benchmark gain in SD (lambda_SD^2 share): %s\n", paste(sprintf("%s %.3f", v, lamS^2 / sum(lamS^2)), collapse = "; ")))
cat(sprintf("  ULS weights per benchmark gain in logits (lambda^2 share): %s\n", paste(sprintf("%s %.3f", v, lam^2 / sum(lam^2)), collapse = "; ")))
qs <- c(.975, .99, .995, .999)
rs <- function(lab) { u <- unique(lab); s <- sample(u, length(u), replace = TRUE)
  list(idx = unlist(lapply(s, function(z) which(lab == z))), lab = unlist(lapply(seq_along(s), function(k) rep(k, sum(lab == s[k]))))) }
calib <- function(lab, M, resamp) {   # resamp: "unit" (rows of M) or "cluster" (organisations of lineages)
  n <- nrow(M); Mc <- sweep(M, 2, colMeans(M))
  draw <- function(MM, L) if (resamp == "unit") MM[sample(nrow(MM), nrow(MM), replace = TRUE), , drop = FALSE] else { r <- rs(L); MM[r$idx, , drop = FALSE] }
  cv <- t(replicate(500, { if (resamp == "unit") { Ms <- Mc[sample(n, n, replace = TRUE), , drop = FALSE]; L2 <- NULL } else { r1 <- rs(og); Ms <- Mc[r1$idx, , drop = FALSE]; L2 <- r1$lab }
    bm <- t(replicate(1000, ests(draw(Ms, L2))))
    as.vector(sapply(1:5, function(k) sapply(qs, function(q) quantile(bm[, k], q, na.rm = TRUE) >= 0))) }))
  cov <- matrix(colMeans(cv, na.rm = TRUE), length(qs), 5, dimnames = list(qs, names(e0l)))
  bm <- t(replicate(5000, ests(draw(M, og)))); e0 <- ests(M)
  cat(sprintf("\n%s: coverage of the upper bound by percentile q\n", lab)); print(round(cov, 3))
  for (k in 1:5) { qk <- qs[which(cov[, k] >= .975)[1]]; if (is.na(qk)) qk <- max(qs)
    cat(sprintf("  %-10s estimate %+.3f; 2.5th %+.3f; 97.5th %+.3f; calibrated upper bound (q = %.3f) %+.3f\n", colnames(cov)[k], e0[k],
      quantile(bm[, k], .025, na.rm = TRUE), quantile(bm[, k], .975, na.rm = TRUE), qk, quantile(bm[, k], qk, na.rm = TRUE))) } }
calib("Lineage-weighted, lineages resampled (48)", Ml, "unit")
calib("Lineage-weighted, organisations resampled (16 clusters)", Ml, "cluster")
calib("Organisation-weighted, organisations resampled (16)", Mo, "unit")
lobo <- t(sapply(1:6, function(j) { l <- lamS[-j]; mS <- colMeans(Ml)[-j] / sdP[-j]; lo <- lam[-j]; mo2 <- colMeans(Ml)[-j]
  c(uls_sd = sum(l * mS) / sum(l^2), uls_logit = sum(lo * mo2) / sum(lo^2)) })); rownames(lobo) <- v
cat("\nLeave-one-benchmark-out (lineage-weighted), ULS in both metrics:\n"); print(round(lobo, 3))
cat(sprintf("Kenny check with the scaled statistic is reported in (K).\n"))
for (mt in c("logit","sd")) { cl <- ifc(Ml, mt); co <- ifc(Mo, mt); tl <- t.test(cl); to <- t.test(co)
  cat(sprintf("\nIFEval contrast (%s metric): lineages %+.2f, t(%d) = %.2f, p = %.2g; organisations %+.2f, t(%d) = %.2f, p = %.2g; positive in %d of %d organisations\n",
    mt, mean(cl), length(cl) - 1, tl$statistic, tl$p.value, mean(co), length(co) - 1, to$statistic, to$p.value, sum(co > 0), length(co))) }

cat("\n=== (L) Leave-one-organisation-out (lineage-weighted) ===\n")
cnt <- data.frame(org = og, pairs = as.numeric(table(C$base)[lin]))
tab <- aggregate(pairs ~ org, data = cnt, FUN = sum); tab$lineages <- as.numeric(table(og)[tab$org]); tab <- tab[order(-tab$pairs), ]
print(tab, row.names = FALSE)
L <- t(sapply(tab$org, function(o) { keep <- og != o; M <- Ml[keep, , drop = FALSE]; e <- ests(M); cl <- ifc(M); tt <- t.test(cl)
  c(lineages = nrow(M), e, ifeval_sd = mean(M[, 1]) / sdP[1], ifc = mean(cl), ifc_t = unname(tt$statistic), ifc_p = tt$p.value) }))
print(round(L, 3))
cat("Ranges over the 16 deletions:\n"); print(round(apply(L[, 2:10], 2, range), 3))
cat("\nOrganisation-weighted, leave one organisation out (15 units):\n")
Lo <- t(sapply(rownames(Mo), function(o) { M <- Mo[rownames(Mo) != o, , drop = FALSE]; e <- ests(M); co <- ifc(M); tt <- t.test(co)
  c(e, ifc = mean(co), ifc_t = unname(tt$statistic), ifc_p = tt$p.value) }))
print(round(apply(Lo, 2, range), 3))

cat("\n=== (F) Scoring-format contrast: IFEval-implied minus MATH-implied factor change ===\n")
f1 <- Ml[, 1] / lam[1] - Ml[, 3] / lam[3]; f2 <- Mo[, 1] / lam[1] - Mo[, 3] / lam[3]
t1 <- t.test(f1); t2 <- t.test(f2)
cat(sprintf("lineages: %+.2f g-SD, t(%d) = %.2f, p = %.2g; organisations: %+.2f, t(%d) = %.2f, p = %.2g\n", mean(f1), length(f1) - 1, t1$statistic, t1$p.value, mean(f2), length(f2) - 1, t2$statistic, t2$p.value))
g1 <- Ml[, 1] / sdP[1] - Ml[, 3] / sdP[3]; tg <- t.test(g1)
cat(sprintf("gain difference IFEval - MATH in pretrained SD: %+.2f, t(%d) = %.2f, p = %.2g\n", mean(g1), length(g1) - 1, tg$statistic, tg$p.value))

cat("\n=== (K) Kenny et al. (2015) Eq. 2 on the holdout fit of M4 ===\n")
set.seed(20261113); lins <- unique(x$base); h <- sample(c(1, 2), length(lins), replace = TRUE); names(h) <- lins; H <- x[h[x$base] == 2, ]
fh <- cfa(m4, data = H, estimator = "MLR", cluster = "base", std.lv = TRUE)
fm <- fitMeasures(fh, c("chisq","df","rmsea","chisq.scaled","rmsea.scaled","rmsea.robust","ntotal"))
print(round(fm, 4))
for (cc in c(fm["rmsea"], fm["rmsea.robust"])) { thr <- fm["df"] * ((fm["ntotal"] - 1) * cc^2 + 1)
  cat(sprintf("c = %.4f: threshold chi2 = %.1f, P = %.3g\n", cc, thr, pchisq(thr, fm["df"], lower.tail = FALSE))) }
cat(sprintf("Scaled statistic (RMSEA.scaled %.4f): P(chi2_%d > %.2f) = %.3g\n", fm["rmsea.scaled"], fm["df"], fm["chisq.scaled"], pchisq(fm["chisq.scaled"], fm["df"], lower.tail = FALSE)))

cat("\n=== (T) Median submission dates ===\n")
A <- x[isP(x$type) | isF(x$type), ]; s1 <- as.Date(substr(A$submit, 1, 10)); s2 <- as.Date(substr(x$submit, 1, 10))
cat(sprintf("Training-role models (n = %d with date): median %s; all models (n = %d with date): median %s\n", sum(!is.na(s1)), median(s1, na.rm = TRUE), sum(!is.na(s2)), median(s2, na.rm = TRUE)))
sink()
