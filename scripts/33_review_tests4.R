suppressMessages(library(lavaan))
# Fourth set of review tests: (C) coverage and calibration of the bootstrap upper bounds for the within-lineage factor change;
# (P) invariance across evaluation periods (submission date), calibrated by permuting the period labels.
# args: [block: all|C|P] [outfile]
args <- commandArgs(TRUE); part <- if (length(args) >= 1) args[1] else "all"; outf <- if (length(args) >= 2) args[2] else "out/w_review_tests4_v1.txt"
x <- read.csv("derived/hf_models_v2.csv")
v <- c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl <- paste0(v, "_lg")
m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained")
isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
sink(outf)
if (part %in% c("all","C")) {
  set.seed(20261119)
  P <- x[isP(x$type), ]; P$org <- sub("/.*", "", P$model)
  fp <- cfa(m4, data = P, estimator = "MLR", cluster = "org", std.lv = TRUE)
  pe <- parameterEstimates(fp); lam <- setNames(pe$est[pe$op == "=~"], pe$rhs[pe$op == "=~"])[vl]
  Th <- lavInspect(fp, "est")$theta[vl, vl]; Ti <- solve(Th); w <- as.numeric(solve(t(lam) %*% Ti %*% lam) %*% t(lam) %*% Ti)
  idx <- match(x$base, x$model)
  ok <- x$gen >= 1 & !is.na(idx) & x$model != x$base & isF(x$type) & x$merged == "False" & isP(x$type[idx]); ok[is.na(ok)] <- FALSE
  C <- x[ok, ]; B <- x[idx[ok], ]; G <- as.matrix(C[, vl]) - as.matrix(B[, vl])
  Ml <- as.matrix(aggregate(as.data.frame(G), by = list(cl = C$base), FUN = mean)[, -1]); lin <- sort(unique(C$base))
  Mo <- as.matrix(aggregate(as.data.frame(Ml), by = list(o = sub("/.*", "", lin)), FUN = mean)[, -1])
  ests <- function(M) { mb <- colMeans(M); V <- cov(M) / nrow(M); Vi <- tryCatch(solve(V), error = function(e) NULL)
    gls <- if (is.null(Vi)) NA else as.numeric(t(lam) %*% Vi %*% mb / (t(lam) %*% Vi %*% lam))
    c(bartlett = sum(w * mb), gls = gls, uls = sum(lam * mb) / sum(lam^2), median = median(mb / lam)) }
  sdP <- sapply(P[, vl], sd)
  cat("Benchmark gains in pretrained SD, lineage-weighted:", paste(sprintf("%s %+.2f", v, colMeans(Ml) / sdP), collapse = "; "), "\n")
  cat("Benchmark gains in pretrained SD, organisation-weighted:", paste(sprintf("%s %+.2f", v, colMeans(Mo) / sdP), collapse = "; "), "\n\n")
  qs <- c(.975, .99, .995, .999)
  cat("=== (C) Calibration of the bootstrap upper bounds for the within-lineage factor change ===\n")
  cat("Double bootstrap: 500 data sets simulated by resampling units from the recentred unit means plus a pure g change of zero\n(loadings x 0), each with 1000 inner resamples. Coverage = share of data sets whose upper bound at percentile q is at or above zero.\n")
  cat("The calibrated upper bound uses the smallest q with coverage of at least .975 for that estimator and level.\n")
  for (lev in c("lineage","organisation")) { M <- if (lev == "lineage") Ml else Mo; n <- nrow(M); Mc <- sweep(M, 2, colMeans(M))
    cv <- t(replicate(500, { Ms <- Mc[sample(n, n, replace = TRUE), ]
      bm <- t(replicate(1000, ests(Ms[sample(n, n, replace = TRUE), , drop = FALSE])))
      as.vector(sapply(1:4, function(k) sapply(qs, function(q) quantile(bm[, k], q, na.rm = TRUE) >= 0))) }))
    cov <- matrix(colMeans(cv, na.rm = TRUE), length(qs), 4, dimnames = list(qs, c("bartlett","gls","uls","median")))
    bm <- t(replicate(5000, ests(M[sample(n, n, replace = TRUE), , drop = FALSE]))); e0 <- ests(M)
    cat(sprintf("\n%s units (n = %d): coverage of the upper bound by percentile q\n", lev, n)); print(round(cov, 3))
    for (k in 1:4) { qk <- qs[which(cov[, k] >= .975)[1]]; if (is.na(qk)) qk <- max(qs)
      cat(sprintf("  %-9s estimate %+.3f; 97.5th percentile %+.3f; calibrated upper bound (q = %.3f) %+.3f\n", colnames(cov)[k], e0[k],
        quantile(bm[, k], .975, na.rm = TRUE), qk, quantile(bm[, k], qk, na.rm = TRUE))) } }
  # Lineage-weighted estimand (as in the article) with organisations as the resampled units
  og <- sub("/.*", "", lin); Mc <- sweep(Ml, 2, colMeans(Ml))
  rs <- function(lab) { u <- unique(lab); s <- sample(u, length(u), replace = TRUE)   # cluster resample; returns rows and fresh cluster labels
    list(idx = unlist(lapply(s, function(z) which(lab == z))), lab = unlist(lapply(seq_along(s), function(k) rep(k, sum(lab == s[k]))))) }
  cv <- t(replicate(500, { r1 <- rs(og); Ms <- Mc[r1$idx, , drop = FALSE]
    bm <- t(replicate(1000, { r2 <- rs(r1$lab); ests(Ms[r2$idx, , drop = FALSE]) }))
    as.vector(sapply(1:4, function(k) sapply(qs, function(q) quantile(bm[, k], q, na.rm = TRUE) >= 0))) }))
  cov <- matrix(colMeans(cv, na.rm = TRUE), length(qs), 4, dimnames = list(qs, c("bartlett","gls","uls","median")))
  bm <- t(replicate(5000, { r <- rs(og); ests(Ml[r$idx, , drop = FALSE]) })); e0 <- ests(Ml)
  cat("\nLineage-weighted estimand, organisations resampled (16 clusters of 48 lineages): coverage of the upper bound by percentile q\n"); print(round(cov, 3))
  for (k in 1:4) { qk <- qs[which(cov[, k] >= .975)[1]]; if (is.na(qk)) qk <- max(qs)
    cat(sprintf("  %-9s estimate %+.3f; 97.5th percentile %+.3f; calibrated upper bound (q = %.3f) %+.3f\n", colnames(cov)[k], e0[k],
      quantile(bm[, k], .975, na.rm = TRUE), qk, quantile(bm[, k], qk, na.rm = TRUE))) } }
if (part %in% c("all","P")) {
  set.seed(20261120)
  A <- x[isP(x$type) | isF(x$type), ]; A$sub <- as.Date(substr(A$submit, 1, 10)); cat("\nModels without a submission date (excluded from (P)):", sum(is.na(A$sub)), "\n")
  A <- A[!is.na(A$sub), ]; med <- median(A$sub)
  A$per <- factor(ifelse(A$sub <= med, "early", "late")); A$role <- ifelse(isP(A$type), "pretrained", "finetuned")
  A$cl <- ifelse(A$role == "pretrained", sub("/.*", "", A$model), A$base)   # pretrained: organisation (each is its own base); fine-tuned: base model
  st <- function(d, cl) { f <- lapply(list("", "loadings", c("loadings","intercepts")), function(e) cfa(m4, data = d, group = "grp", estimator = "MLR", cluster = cl, group.equal = e))
    cf <- sapply(f, fitMeasures, "cfi.robust"); c(m = cf[2] - cf[1], s = cf[3] - cf[2]) }
  cat(sprintf("\n=== (P) Evaluation period: models submitted on or before %s (early) versus later (late) ===\n", med))
  for (r in c("pretrained","finetuned")) { d <- A[A$role == r, ]; d$grp <- d$per; cl <- "cl"
    o <- st(d, cl); pr <- t(replicate(200, { dp <- d; dp$grp <- sample(dp$grp); tryCatch(st(dp, cl), error = function(e) c(NA, NA)) })); pr <- pr[complete.cases(pr), , drop = FALSE]
    cat(sprintf("  %s, early vs late (n = %s): dCFI metric %+.4f (permutation p = %.3f), scalar %+.4f (p = %.3f); %d permutations\n", r, paste(table(d$grp), collapse = "/"),
      o[1], (1 + sum(pr[, 1] <= o[1])) / (1 + nrow(pr)), o[2], (1 + sum(pr[, 2] <= o[2])) / (1 + nrow(pr)), nrow(pr))) }
  # Drift test separating evaluation date from model vintage: MIMIC with upload date and submission date (years) as covariates of the factor
  # and direct effects of both on every benchmark except a referent; a Wald test of the five submission-date direct effects.
  base_syn <- readLines("out/m4_syntax_logit_v1.txt")
  for (r in c("pretrained","finetuned")) { up <- as.Date(substr(A$upload, 1, 10)); d <- A[A$role == r & !is.na(up), ]
    d$upy <- as.numeric(as.Date(substr(d$upload, 1, 10)) - as.Date("2024-01-01")) / 365.25; d$suby <- as.numeric(d$sub - as.Date("2024-01-01")) / 365.25
    for (ref in c("bbh_lg","mmlupro_lg")) { oth <- setdiff(vl, ref)
      syn <- c(base_syn, "g ~ upy + suby", paste0(oth, " ~ du_", sub("_lg", "", oth), "*upy + ds_", sub("_lg", "", oth), "*suby"))
      f <- tryCatch(sem(paste(syn, collapse = "\n"), data = d, estimator = "MLR", cluster = "cl"), error = function(e) NULL)
      if (is.null(f)) { cat(sprintf("  %s, referent %s: model failed\n", r, ref)); next }
      pe <- parameterEstimates(f); ds <- pe[pe$label %in% paste0("ds_", sub("_lg", "", oth)), ]; sdv <- sapply(d[, oth], sd)
      wt <- lavTestWald(f, paste(paste0("ds_", sub("_lg", "", oth), " == 0"), collapse = "\n"))
      cat(sprintf("  %s (n = %d, upload-submission correlation %.2f), referent %s: submission-date direct effects (SD per year) %s; Wald chi2(5) = %.1f, p = %.2g\n",
        r, nrow(d), cor(d$upy, d$suby), sub("_lg", "", ref), paste(sprintf("%s %+.2f", sub("_lg", "", ds$lhs), ds$est / sdv[ds$lhs]), collapse = "; "), wt$stat, wt$p.value)) } }
  # Subsets for exact nulls (10_invariance_null.R): training role within each evaluation period; cohorts within each period
  for (p in levels(A$per)) write.csv(x[match(A$model[A$per == p], x$model), ], sprintf("derived/hf_models_v2_period_%s.csv", p), row.names = FALSE)
  Pe <- A[A$role == "pretrained" & A$per == "early" & A$upload_year %in% c(2022, 2023, 2024), ]
  Pl <- A[A$role == "pretrained" & A$per == "late" & A$upload_year %in% c(2024, 2025), ]
  cat(sprintf("  Cohorts within the early period (2022-23 vs 2024): n = %s; within the late period (2024 vs 2025): n = %s\n",
    paste(table(Pe$upload_year >= 2024), collapse = "/"), paste(table(Pl$upload_year), collapse = "/")))
  write.csv(x[match(Pe$model, x$model), ], "derived/hf_models_v2_cohort_early.csv", row.names = FALSE)
  write.csv(x[match(Pl$model, x$model), ], "derived/hf_models_v2_cohort_late.csv", row.names = FALSE)
  for (p in levels(A$per)) { d <- A[A$per == p, ]; d$grp <- factor(d$role, levels = c("pretrained","finetuned")); o <- st(d, "base")
    cat(sprintf("  training role within the %s period (n = %s): dCFI metric %+.4f, scalar %+.4f\n", p, paste(table(d$grp), collapse = "/"), o[1], o[2])) } }
sink()
