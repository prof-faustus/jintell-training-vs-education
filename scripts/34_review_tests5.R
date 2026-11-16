suppressMessages(library(lavaan))
# Fifth set of review tests.
# (A) Sample construction: the within-lineage estimators and the IFEval contrast recomputed on the files written by 34_prep_alt.py
#     (main = analysis file; flagged = community-flagged model kept; f16 = float16 preferred over bfloat16). The main file must
#     reproduce the published values. The invariance sequences for the alternative files are run by 10_invariance_null.R.
# (S) Stability of the measurement-model selection: the selection procedure of 06_select_logit.R repeated on 100 random lineage splits.
# args: [block: all|A|S] [outfile]
args <- commandArgs(TRUE); part <- if (length(args) >= 1) args[1] else "all"; outf <- if (length(args) >= 2) args[2] else "out/x_review_tests5_v1.txt"
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
sink(outf)
if (part %in% c("all","A")) {
  cat("=== (A) Within-lineage results under alternative sample constructions ===\n")
  for (variant in c("main","flagged","f16")) {
    x <- read.csv(sprintf("derived/hf_models_v2_%s.csv", variant))
    P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
    fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
    pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
    Th <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti)
    idx <- match(x$base, x$model)
    ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
    C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
    Ml <- as.matrix(aggregate(as.data.frame(G), by = list(cl = C$base), FUN = mean)[, -1]); lin <- sort(unique(C$base))
    Mo <- as.matrix(aggregate(as.data.frame(Ml), by = list(o = sub("/.*", "", lin)), FUN = mean)[, -1])
    mb <- colMeans(Ml); V <- cov(Ml) / nrow(Ml); Vi <- solve(V)
    e <- c(bartlett = sum(w * mb), gls = as.numeric(t(lam) %*% Vi %*% mb / (t(lam) %*% Vi %*% lam)), uls = sum(lam * mb) / sum(lam^2), median = median(mb / lam))
    sdP <- sapply(P[, vl], sd)
    con <- function(Mx, l) Mx[, 1] / l[1] - as.numeric(Mx[, -1] %*% l[-1]) / sum(l[-1]^2)
    tl <- t.test(con(Ml, lam)); to <- t.test(con(Mo, lam))
    cat(sprintf("\n%s: %d models; %d pretrained; %d pairs in %d lineages, %d organisations\n", variant, nrow(x), nrow(P), nrow(C), nrow(Ml), nrow(Mo)))
    cat("  gains in pretrained SD:", paste(sprintf("%s %+.2f", v, mb / sdP), collapse = "; "), "\n")
    cat("  factor change:", paste(sprintf("%s %+.3f", names(e), e), collapse = "; "), "\n")
    cat(sprintf("  IFEval contrast: lineages t(%d) = %.2f, p = %.2g; organisations t(%d) = %.2f, p = %.2g\n",
      nrow(Ml) - 1, tl$statistic, tl$p.value, nrow(Mo) - 1, to$statistic, to$p.value)) } }
if (part %in% c("all","S")) {
  set.seed(20261121); x <- read.csv("derived/hf_models_v2.csv")
  fm <- c("cfi.robust","rmsea.robust","srmr"); base1 <- paste("g =~", paste(vl, collapse = "+"))
  key <- function(a, b) paste(sort(c(a, b)), collapse = "~~")
  m4set <- unname(sort(sapply(regmatches(m4, gregexpr("[a-z]+_lg ~~ [a-z]+_lg", m4))[[1]], function(s) { p <- strsplit(s, " ~~ ")[[1]]; key(p[1], p[2]) })))
  res <- t(replicate(100, { lin <- unique(x$base); h <- sample(c(1, 2), length(lin), replace = TRUE); names(h) <- lin
    D <- x[h[x$base] == 1, ]; H <- x[h[x$base] == 2, ]; m <- base1; added <- c()
    for (i in 1:3) { f <- cfa(m, data = D, estimator = "MLR", cluster = "base", std.lv = TRUE); mi <- modindices(f); mi <- mi[mi$op == "~~", ]; mi <- mi[order(-mi$mi), ]
      if (mi$mi[1] < 50) break; added <- c(added, key(mi$lhs[1], mi$rhs[1])); m <- paste(m, "\n", mi$lhs[1], "~~", mi$rhs[1]) }
    fs <- fitMeasures(cfa(m, data = H, estimator = "MLR", cluster = "base", std.lv = TRUE), fm)
    f4 <- fitMeasures(cfa(m4, data = H, estimator = "MLR", cluster = "base", std.lv = TRUE), fm)
    c(set = paste(sort(added), collapse = " | "), n_added = length(added), same = identical(sort(added), m4set), cfi_sel = fs[1], rmsea_sel = fs[2], cfi_m4 = f4[1], rmsea_m4 = f4[2], srmr_m4 = f4[3]) }))
  cat("\n=== (S) Measurement-model selection repeated on 100 random lineage splits (procedure of 06_select_logit.R) ===\n")
  cat("Published M4 covariances:", paste(m4set, collapse = " | "), "\n")
  cat("Selections identical to M4:", sum(res[, "same"] == "TRUE"), "of 100\n"); cat("Number of covariances added:\n"); print(table(res[, "n_added"]))
  cat("Selected sets:\n"); print(sort(table(res[, "set"]), decreasing = TRUE))
  cat("Splits whose selection contains each M4 covariance:", paste(sapply(m4set, function(k) sum(grepl(k, res[, "set"], fixed = TRUE))), collapse = " / "), "\n")
  q <- function(z) { z <- as.numeric(z); sprintf("median %.3f, range %.3f to %.3f", median(z), min(z), max(z)) }
  cat("Holdout fit of the selected model: robust CFI", q(res[, grep("^cfi_sel", colnames(res))]), "; robust RMSEA", q(res[, grep("^rmsea_sel", colnames(res))]), "\n")
  cat("Holdout fit of M4 itself:          robust CFI", q(res[, grep("^cfi_m4", colnames(res))]), "; robust RMSEA", q(res[, grep("^rmsea_m4", colnames(res))]), "; SRMR", q(res[, grep("^srmr_m4", colnames(res))]), "\n") }
sink()
