suppressMessages({library(lavaan)})
# dMACS (Nye & Drasgow 2011) per indicator, every referent; cluster-bootstrap CI for the primary referent.
args <- commandArgs(TRUE); x <- read.csv(args[1]); outf <- args[2]; B <- as.integer(args[3]); set.seed(20261113)
ind <- c("ifeval_lg","bbh_lg","math_lg","gpqa_lg","musr_lg","mmlupro_lg")
rc <- "ifeval_lg ~~ math_lg\ngpqa_lg ~~ musr_lg\nbbh_lg ~~ musr_lg"
mk <- function(r) paste0("g =~ 1*", r, " + ", paste(setdiff(ind, r), collapse = " + "), "\n", r, " ~ 0*1\ng ~ NA*1\n", rc)
prep <- function(which) {
  if (which == "role") {
    x$grp <- ifelse(x$type %in% c("pretrained","continuously pretrained"), "pretrained",
             ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"), "finetuned", NA))
    d <- x[!is.na(x$grp), ]; d$cl <- d$base; d$grp <- factor(d$grp, levels = c("pretrained","finetuned"))
  } else {
    d <- x[x$type %in% c("pretrained","continuously pretrained") & !is.na(x$upload_year), ]
    d$grp <- factor(ifelse(d$upload_year <= 2023, "c2022_23", ifelse(d$upload_year == 2024, "c2024", "c2025")), levels = c("c2022_23","c2024","c2025"))
    d$cl <- sub("/.*", "", d$model) }
  d[complete.cases(d[, ind]), ] }
dmacs <- function(d, r) {
  f <- cfa(mk(r), data = d, group = "grp", estimator = "ML")
  if (!lavInspect(f, "converged")) return(NULL)
  est <- lavInspect(f, "est"); lv <- levels(d$grp)
  ng <- table(d$grp); sdp <- sapply(ind, function(v) sum(sapply(lv, function(g) (ng[g] - 1) * sd(d[d$grp == g, v]))) / sum(ng - 1))
  R <- est[[1]]; out <- list()
  for (gi in 2:length(lv)) { Fp <- est[[gi]]
    dl <- R$lambda[ind, "g"] - Fp$lambda[ind, "g"]; dt <- R$nu[ind, 1] - Fp$nu[ind, 1]
    kF <- Fp$alpha[1, 1]; pF <- Fp$psi[1, 1]
    dm <- sqrt((dt + dl * kF)^2 + dl^2 * pF) / sdp
    sg <- (dt + dl * kF) / sdp   # signed: positive = focal scores LOWER than reference at equal g
    out[[lv[gi]]] <- rbind(dMACS = dm, signed = sg) }
  out }
sink(outf)
for (which in c("role","cohort")) {
  d <- prep(which); cat("\n########", which, "| n:", paste(names(table(d$grp)), table(d$grp), collapse = "; "), "\n")
  cat("Pooled-SD rule: Nye & Drasgow Eq. 4 over all groups. Signed index > 0 means the focal group scores lower than the reference at equal g.\n")
  allr <- list()
  for (r in ind) { o <- dmacs(d, r); if (is.null(o)) { cat("referent", r, "did not converge\n"); next }
    allr[[r]] <- o; for (g in names(o)) { cat("\nReferent", r, "| focal", g, "\n"); print(round(o[[g]], 3)) } }
  # range across referents (excluding the referent itself, whose dMACS is 0 by construction)
  for (g in names(allr[[1]])) { M <- sapply(names(allr), function(r) { v <- allr[[r]][[g]]["dMACS", ]; v[r] <- NA; v })
    cat("\nRange of dMACS across referents | focal", g, "\n"); print(round(t(apply(M, 1, range, na.rm = TRUE)), 3)) }
  pr <- if (which == "role") "bbh_lg" else "mmlupro_lg"
  cat("\nPrimary referent:", pr, "| cluster bootstrap B =", B, "\n")
  bs <- list(); cls <- unique(d$cl)
  for (b in 1:B) { pick <- sample(cls, length(cls), replace = TRUE)
    bb <- do.call(rbind, lapply(seq_along(pick), function(j) d[d$cl == pick[j], ]))
    if (any(table(bb$grp) < 20)) next
    o <- try(dmacs(bb, pr), silent = TRUE); if (!inherits(o, "try-error") && !is.null(o)) bs[[length(bs) + 1]] <- o }
  cat("usable bootstrap replicates:", length(bs), "\n")
  for (g in names(allr[[pr]])) { A <- sapply(bs, function(o) o[[g]]["dMACS", ]); S <- sapply(bs, function(o) o[[g]]["signed", ])
    tab <- cbind(dMACS = allr[[pr]][[g]]["dMACS", ], lo = apply(A, 1, quantile, .025), hi = apply(A, 1, quantile, .975),
                 signed = allr[[pr]][[g]]["signed", ], s_lo = apply(S, 1, quantile, .025), s_hi = apply(S, 1, quantile, .975))
    cat("\nFocal", g, "(referent row is 0 by construction)\n"); print(round(tab, 3)) } }
sink()
