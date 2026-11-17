suppressMessages({library(lavaan)})
# Step-specific Bollen-Stine cluster-bootstrap nulls with joint cluster resampling.
# For each invariance step the data are transformed, group by group, to the moments implied by the model that
# defines that step's null hypothesis (metric step: equal loadings; scalar step: equal loadings and intercepts;
# strict step: equal loadings, intercepts and residual variances), so every nuisance parameter left free under that
# null keeps its observed value. Clusters are then resampled jointly across groups: a sampled lineage or organisation
# carries all of its models, in every group, into the replicate, preserving dependence between groups.
# Only the fit change of the step whose null generated the data is recorded.
# args: data.csv  syntax.txt  which(role|cohort)  R  outfile  seed  [resample joint|within]
args <- commandArgs(TRUE)
x <- read.csv(args[1]); m <- paste(readLines(args[2]), collapse = "\n")
which <- args[3]; R <- as.integer(args[4]); outf <- args[5]; seed <- as.integer(args[6])
rs <- if (length(args) >= 7) args[7] else "joint"
ind <- c("ifeval_lg","bbh_lg","math_lg","gpqa_lg","musr_lg","mmlupro_lg")
if (which == "role") {
  x$grp <- ifelse(x$type %in% c("pretrained","continuously pretrained"), "pretrained",
           ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"), "finetuned", NA))
  d <- x[!is.na(x$grp), ]; d$cl <- d$base
} else {
  d <- x[x$type %in% c("pretrained","continuously pretrained") & !is.na(x$upload_year), ]
  d$grp <- ifelse(d$upload_year <= 2023, "c2022_23", ifelse(d$upload_year == 2024, "c2024", "c2025"))
  d$cl <- sub("/.*", "", d$model)
}
d$grp <- factor(d$grp); d <- d[complete.cases(d[, ind]), ]; lv <- levels(d$grp)
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
fm <- c("cfi.robust","rmsea.robust","srmr")
fit1 <- function(dd, k) { f <- try(cfa(m, data = dd, group = "grp", estimator = "MLR", cluster = "cl", group.equal = eqs[[k]]), silent = TRUE)
  if (inherits(f, "try-error") || !lavInspect(f, "converged")) return(rep(NA, length(fm))); as.numeric(fitMeasures(f, fm)) }
msqrt <- function(S, inv = FALSE) { e <- eigen(S, symmetric = TRUE); v <- if (inv) 1/sqrt(e$values) else sqrt(e$values); e$vectors %*% diag(v) %*% t(e$vectors) }
transform_to <- function(k) {
  fs <- cfa(m, data = d, group = "grp", estimator = "MLR", cluster = "cl", group.equal = eqs[[k]]); imp <- lavInspect(fs, "implied"); dt <- d
  for (gi in seq_along(lv)) { r <- dt$grp == lv[gi]; X <- as.matrix(dt[r, ind]); S <- cov(X) * (nrow(X) - 1) / nrow(X)
    Z <- sweep(X, 2, colMeans(X)) %*% msqrt(S, TRUE) %*% msqrt(imp[[gi]]$cov[ind, ind]); dt[r, ind] <- sweep(Z, 2, imp[[gi]]$mean[ind], "+") }
  dt }
resample <- function(dt) {
  if (rs == "joint") { cls <- unique(dt$cl); pick <- sample(cls, length(cls), replace = TRUE)
    bb <- do.call(rbind, lapply(seq_along(pick), function(j) { z <- dt[dt$cl == pick[j], ]; z$cl <- paste0("c", j); z }))
  } else bb <- do.call(rbind, lapply(lv, function(g) { dg <- dt[dt$grp == g, ]; cls <- unique(dg$cl); pick <- sample(cls, length(cls), replace = TRUE)
    do.call(rbind, lapply(seq_along(pick), function(j) { z <- dg[dg$cl == pick[j], ]; z$cl <- paste0(g, "_", j); z })) }))
  bb$grp <- factor(bb$grp, levels = lv); bb }
steps <- c(metric = "configural", scalar = "metric", strict = "scalar")   # step -> its baseline model
obs <- sapply(names(eqs), function(k) fit1(d, k)); rownames(obs) <- fm
od <- sapply(names(steps), function(s) obs[, s] - obs[, steps[[s]]])
t0 <- Sys.time(); out <- list(); chk <- list()
for (s in names(steps)) {
  set.seed(seed + match(s, names(steps)))
  dt <- transform_to(s); chk[[s]] <- c(fit1(dt, steps[[s]]), fit1(dt, s))
  res <- t(sapply(1:R, function(b) { bb <- resample(dt); fit1(bb, s) - fit1(bb, steps[[s]]) })); colnames(res) <- paste0("d", c("CFI","RMSEA","SRMR"))
  out[[s]] <- res; write.csv(res, sub("\\.txt$", paste0("_", s, "_reps.csv"), outf), row.names = FALSE) }
sink(outf)
cat("Step-specific nulls | design:", which, "| resampling:", rs, "| R =", R, "| seed:", seed, "+ step index\n")
cat("Groups:", paste(lv, table(d$grp), collapse = "; "), "| clusters:", length(unique(d$cl)), "| clusters spanning groups:",
    sum(rowSums(table(d$cl, d$grp) > 0) > 1), "| elapsed:", format(Sys.time() - t0), "\n\nObserved fit (robust CFI, robust RMSEA, SRMR)\n"); print(round(obs, 4))
for (s in names(steps)) { res <- out[[s]]; n <- sum(complete.cases(res)); o <- od[, s]
  k <- c(sum(res[, 1] <= o[1], na.rm = TRUE), sum(res[, 2] >= o[2], na.rm = TRUE), sum(res[, 3] >= o[3], na.rm = TRUE))
  cat(sprintf("\n%s step: data transformed to the %s model (check: baseline and constrained fit CFI %.4f / %.4f); converged %d of %d\n", s, s, chk[[s]][1], chk[[s]][4], n, R))
  print(round(rbind(apply(res, 2, quantile, c(.01, .05, .5, .95, .99), na.rm = TRUE), observed = o), 4))
  cat("k (replicates at least as unfavourable):", paste(colnames(res), k, collapse = "; "), "\n")
  cat("p = k/R:", paste(colnames(res), sprintf("%.4f", k / n), collapse = "; "), " | conservative (1+k)/(1+R):", paste(colnames(res), sprintf("%.4f", (1 + k) / (1 + n)), collapse = "; "), "\n") }
sink()
