suppressMessages({library(lavaan)})
args<-commandArgs(TRUE); x<-read.csv(args[1]); od<-args[2]; m<-paste(readLines(args[3]),collapse="\n")
fm<-c("chisq.scaled","df","cfi.robust","tli.robust","rmsea.robust","srmr")
inv<-function(dat,grp,label,cl="base"){
  cat("\n\n########",label,"########\n"); print(table(dat[[grp]]))
  fits<-list(); eqs<-list(configural=NULL,metric="loadings",scalar=c("loadings","intercepts"),strict=c("loadings","intercepts","residuals"))
  for(k in names(eqs)) fits[[k]]<-cfa(m,data=dat,group=grp,estimator="MLR",cluster=cl,std.lv=FALSE,group.equal=if(is.null(eqs[[k]])) "" else eqs[[k]])
  tab<-t(sapply(fits,function(f) round(fitMeasures(f,fm),3))); print(tab)
  cat("\nDelta CFI (robust):",round(diff(tab[,"cfi.robust"]),4),"\nDelta RMSEA:",round(diff(tab[,"rmsea.robust"]),4),"\n")
  print(lavTestLRT(fits$configural,fits$metric,fits$scalar,fits$strict,method="satorra.bentler.2001"))
  sc<-lavTestScore(fits$scalar); u<-sc$uni; pt<-parTable(fits$scalar); u$lab<-paste(pt$lhs[match(u$lhs,pt$plabel)],pt$op[match(u$lhs,pt$plabel)],pt$rhs[match(u$lhs,pt$plabel)])
  cat("\nScore tests on scalar constraints (top)\n"); print(head(u[order(-u$X2),c("lab","X2","p.value")],8))
  # latent mean difference under scalar model
  pe<-parameterEstimates(fits$scalar); cat("\nLatent means/variances (scalar)\n"); print(pe[pe$op %in% c("~1","~~") & pe$lhs=="g" & pe$rhs %in% c("","g"),c("group","op","est","se","pvalue")])
  invisible(fits)}
sink(file.path(od,"c_invariance_cohort_v1.txt"))
x$role<-ifelse(x$type %in% c("pretrained","continuously pretrained"),"pretrained",ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"),"finetuned",NA))

P<-x[x$type %in% c("pretrained","continuously pretrained") & !is.na(x$upload_year),]
P$cohort<-ifelse(P$upload_year<=2023,"c2022_23",ifelse(P$upload_year==2024,"c2024","c2025"))
P$org<-sub("/.*","",P$model); cat("orgs",length(unique(P$org)),"\n"); inv(P,"cohort","(ii) pretrained cohorts by upload year",cl="org")
sink()
