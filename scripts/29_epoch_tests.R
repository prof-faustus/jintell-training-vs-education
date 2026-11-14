suppressMessages(library(lavaan))
# Frontier replication: tests that do not rely on a correctly specified configural model.
# (1) probability of the configural chi-square under a correct model; (2) permutation test of the cohort labels
# for the metric and scalar steps; (3) the same analysis with model versions collapsed by their median instead of their best score.
# args: [R for the main block, default 1000] [R for the three robustness blocks, default 500] [block: all|1|2|3|4] [outfile]
# Each block has its own seed, so a block gives the same result whether it is run alone or with the others.
args <- commandArgs(TRUE); Rmain <- if (length(args) >= 1) as.integer(args[1]) else 1000; Rrob <- if (length(args) >= 2) as.integer(args[2]) else 500
part <- if (length(args) >= 3) args[3] else "all"; outf <- if (length(args) >= 4) args[4] else "out/s_epoch_tests_v1.txt"
vl <- paste0(c("gpqa","aime","arc1","arc2","dtb"), "_lg"); m <- paste("g =~", paste(vl, collapse = " + "))
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"))
seqfit <- function(d) { f <- lapply(eqs, function(e) cfa(m, data = d, group = "grp", estimator = "MLR", missing = "ml", group.equal = e))
  lr <- lavTestLRT(f$configural, f$metric, f$scalar, method = "satorra.bentler.2001")
  c(conf_chisq = fitMeasures(f$configural, "chisq.scaled"), conf_p = fitMeasures(f$configural, "pvalue.scaled"),
    d_metric = lr[2, "Chisq diff"], d_scalar = lr[3, "Chisq diff"],
    dcfi_m = fitMeasures(f$metric, "cfi.robust") - fitMeasures(f$configural, "cfi.robust"),
    dcfi_s = fitMeasures(f$scalar, "cfi.robust") - fitMeasures(f$metric, "cfi.robust")) }
run <- function(file, lab, k0) {
  x <- read.csv(file); x$grp <- factor(ifelse(x$year <= 2025, "le2025", "y2026"))
  for (nm in list(c("all scores", "none"), c("floor-clipped scores as missing", "clip"))) {
    k0 <- k0 + 1; if (part != "all" && part != as.character(k0)) next
    set.seed(20261113 + k0); R <- if (k0 == 1) Rmain else Rrob
    d <- x; if (nm[2] == "clip") for (b in sub("_lg", "", vl)) d[[paste0(b, "_lg")]][d[[paste0(b, "_clip")]] %in% 1] <- NA
    d <- d[rowSums(!is.na(d[, vl])) >= 3, ]
    o <- seqfit(d)
    perm <- t(replicate(R, { dp <- d; dp$grp <- sample(dp$grp); tryCatch(seqfit(dp), error = function(e) rep(NA, 6)) }))
    ok <- complete.cases(perm); perm <- perm[ok, , drop = FALSE]; k <- nrow(perm)
    cat(sprintf("\n[%s | %s] n = %d (%s); usable permutations %d of %d\n", lab, nm[1], nrow(d), paste(table(d$grp), collapse = "/"), k, R))
    cat(sprintf("  configural scaled chi2(10) = %.2f, p = %.3f (probability of a chi-square this large under a correct model)\n", o[1], o[2]))
    cat(sprintf("  metric: scaled dchi2(4) = %.1f, asymptotic p = %.2g, permutation p = %.4f; dCFI = %+.3f, permutation p = %.4f\n",
      o[3], pchisq(o[3], 4, lower.tail = FALSE), (1 + sum(perm[, 3] >= o[3])) / (k + 1), o[5], (1 + sum(perm[, 5] <= o[5])) / (k + 1)))
    cat(sprintf("  scalar: scaled dchi2(4) = %.1f, asymptotic p = %.2g, permutation p = %.4f; dCFI = %+.3f, permutation p = %.4f\n",
      o[4], pchisq(o[4], 4, lower.tail = FALSE), (1 + sum(perm[, 4] >= o[4])) / (k + 1), o[6], (1 + sum(perm[, 6] <= o[6])) / (k + 1)))
    cat(sprintf("  permutation 95th/99th percentiles of metric dchi2: %.1f / %.1f\n", quantile(perm[, 3], .95), quantile(perm[, 3], .99))) } }
sink(outf)
run("derived/epoch_wide_v1.csv", "best score per model group", 0)
run("derived/epoch_wide_median_v1.csv", "median score per model group", 2)
sink()
