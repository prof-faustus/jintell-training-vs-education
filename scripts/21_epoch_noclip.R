suppressMessages({library(lavaan)})
x <- read.csv("derived/epoch_wide_v1.csv"); b <- c("gpqa","aime","arc1","arc2","dtb"); vl <- paste0(b, "_lg")
for (k in b) x[[paste0(k, "_lg")]][x[[paste0(k, "_clip")]] %in% 1] <- NA
x <- x[rowSums(!is.na(x[, vl])) >= 3, ]
x$grp <- factor(ifelse(x$year <= 2025, "le2025", "y2026")); fm <- c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr")
m <- paste("g =~", paste(vl, collapse = " + "))
sink("out/i_epoch_noclip_v1.txt"); cat("Clipped values set to missing | n:", nrow(x), "\n"); print(table(x$grp))
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"), strict = c("loadings","intercepts","residuals"))
fits <- lapply(eqs, function(e) cfa(m, data = x, group = "grp", estimator = "MLR", missing = "ml", group.equal = e))
tab <- t(sapply(fits, function(f) round(fitMeasures(f, fm), 3))); print(tab); cat("Delta CFI:", round(diff(tab[, "cfi.robust"]), 3), "\n")
print(lavTestLRT(fits$configural, fits$metric, fits$scalar, fits$strict, method = "satorra.bentler.2001"))
sink()
