suppressMessages({library(lavaan)})
args<-commandArgs(TRUE); x<-read.csv(args[1]); pr<-read.csv(args[2]); od<-args[3]; m<-paste(readLines(args[4]),collapse="\n")
v<-c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl<-paste0(v,"_lg")
rel<-setNames(pr$retest_r,pr$bench)[v]
sink(file.path(od,"d_mcv_partial_v1.txt"))
cat("Retest reliabilities (logit scale)\n"); print(round(rel,3))
# loadings from pretrained sample
P<-x[x$type %in% c("pretrained","continuously pretrained"),]; P$org<-sub("/.*","",P$model)
fp<-cfa(m,data=P,estimator="MLR",cluster="org",std.lv=TRUE); ld<-standardizedSolution(fp); ld<-setNames(ld$est.std[ld$op=="=~"],ld$rhs[ld$op=="=~"])[vl]
cat("\nPretrained g loadings\n"); print(round(ld,3)); sdP<-sapply(P[,vl],sd)
# pairs: base and child logits
idx<-match(x$base,x$model); pr_ok<-x$gen>=1 & !is.na(idx) & x$model!=x$base
C<-x[pr_ok,]; Bm<-x[idx[pr_ok],]
G<-as.matrix(C[,vl])-as.matrix(Bm[,vl]); colnames(G)<-v
cat("\nPairs:",nrow(G),"bases:",length(unique(C$base)),"\n")
dgain<-colMeans(G)/sdP; cat("Standardised mean gains (pooled pretrained SD units)\n"); print(round(dgain,3))
cat("Paired t-tests (cluster-robust via lineage means)\n")
lm_<-aggregate(G,by=list(base=C$base),FUN=mean); for(j in v){tt<-t.test(lm_[[j]]); cat(j,"lineage-mean gain",round(mean(lm_[[j]]),3),"t",round(tt$statistic,2),"p",signif(tt$p.value,3),"\n")}
mcv<-function(g,l,r) cor(l/sqrt(r), g/sqrt(r))
r_unc<-cor(ld,dgain); r_cor<-mcv(dgain,ld,rel)
cat("\nMCV r (uncorrected):",round(r_unc,3)," r (reliability-corrected):",round(r_cor,3),"\n")
# cluster bootstrap over lineages
set.seed(20261113); bs<-unique(C$base); B<-2000; rb<-numeric(B)
for(b in 1:B){ s<-sample(bs,replace=TRUE); ii<-unlist(lapply(s,function(z) which(C$base==z))); rb[b]<-mcv(colMeans(G[ii,,drop=FALSE])/sdP,ld,rel)}
cat("Cluster-bootstrap 95% CI (corrected):",round(quantile(rb,c(.025,.975),na.rm=TRUE),3),"\n")
# permutation of loadings (6! = 720 exact)
perms<-function(n){ if(n==1) return(matrix(1)); p<-perms(n-1); do.call(rbind,lapply(1:n,function(i) cbind(i,ifelse(p>=i,p+1,p))))}
Pm<-perms(6); null<-apply(Pm,1,function(o) mcv(dgain,ld[o],rel[o])); cat("Exact permutation p (two-sided):",mean(abs(null)>=abs(r_cor)-1e-12),"\n")
# partial invariance with proper labels
partial<-function(dat,grp,cl,label,free){
  f<-cfa(m,data=dat,group=grp,estimator="MLR",cluster=cl,group.equal=c("loadings","intercepts"),group.partial=free)
  cat("\n####",label,"| groups:",paste(lavInspect(f,"group.label"),collapse=", "),"| freed:",paste(free,collapse="; "),"\n")
  print(round(fitMeasures(f,c("chisq.scaled","df","cfi.robust","rmsea.robust","srmr")),3))
  pe<-parameterEstimates(f); print(pe[pe$lhs=="g" & pe$op=="~1",c("group","est","se","pvalue")])
  pe2<-pe[pe$lhs=="g"&pe$op=="~~",c("group","est")]; print(pe2); f}
x$role<-ifelse(x$type %in% c("pretrained","continuously pretrained"),"pretrained",ifelse(x$type %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)"),"finetuned",NA))
A<-x[!is.na(x$role),]; A$role<-factor(A$role,levels=c("pretrained","finetuned")); A<-A[order(A$role),]
fa1<-partial(A,"role","base","(i) role, full scalar",character(0))
fa2<-partial(A,"role","base","(i) role, partial",c("g=~musr_lg","g=~mmlupro_lg","ifeval_lg~1","musr_lg~1"))
P$cohort<-factor(ifelse(P$upload_year<=2023,"c2022_23",ifelse(P$upload_year==2024,"c2024","c2025"))); P<-P[!is.na(P$cohort),]; P<-P[order(P$cohort),]
fb1<-partial(P,"cohort","org","(ii) cohort, full scalar",character(0))
fb2<-partial(P,"cohort","org","(ii) cohort, partial",c("gpqa_lg~1","math_lg~1","g=~gpqa_lg","bbh_lg~1"))
sink()
