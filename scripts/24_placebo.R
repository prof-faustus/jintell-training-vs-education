suppressMessages(library(lavaan)); set.seed(20261115)
x <- read.csv("derived/hf_models_v2.csv"); m4 <- paste(readLines("out/m4_syntax_logit_v1.txt"), collapse = "\n")
isP <- function(t) t %in% c("pretrained","continuously pretrained"); isF <- function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
eqs <- list(configural = "", metric = "loadings", scalar = c("loadings","intercepts"))
dc <- function(dd) { c <- sapply(eqs, function(e) { f <- try(cfa(m4, data = dd, group = "grp", estimator = "MLR", cluster = "cl", group.equal = e), silent = TRUE)
  if (inherits(f, "try-error") || !lavInspect(f, "converged")) NA else fitMeasures(f, "cfi.robust") }); diff(c) }
run <- function(d, R) t(replicate(R, { u <- unique(d$cl); h <- sample(u, floor(length(u) / 2)); d$grp <- factor(ifelse(d$cl %in% h, "a", "b")); dc(d) }))
F <- x[isF(x$type), ]; F$cl <- F$base
P <- x[isP(x$type), ]; P$cl <- sub("/.*", "", P$model)
rf <- run(F, 100); rp <- run(P, 100)
sink("out/q_placebo_v1.txt")
for (nm in c("fine-tuned split by lineage", "pretrained split by organisation")) { r <- if (grepl("fine", nm)) rf else rp
  cat("\n", nm, "| converged:", sum(complete.cases(r)), "\n"); print(round(apply(r, 2, quantile, c(0, .01, .05, .5, .95), na.rm = TRUE), 4))
  cat("share with dCFI <= -.012 (metric):", mean(r[, 1] <= -.0122, na.rm = TRUE), " share with dCFI <= -.015 (scalar):", mean(r[, 2] <= -.0150, na.rm = TRUE), "\n") }
sink(); write.csv(rbind(data.frame(set = "ft", rf), data.frame(set = "pt", rp)), "out/q_placebo_v1_reps.csv", row.names = FALSE)
