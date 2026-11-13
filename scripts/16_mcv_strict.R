suppressMessages(library(lavaan))
x<-read.csv("derived/hf_models_v2.csv"); pr<-read.csv("out/retest_logit_v1.csv"); m<-paste(readLines("out/m4_syntax_logit_v1.txt"),collapse="\n")
v<-c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl<-paste0(v,"_lg"); rel<-setNames(pr$retest_r,pr$bench)[v]
isP<-function(t) t %in% c("pretrained","continuously pretrained"); isF<-function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
P<-x[isP(x$type),]; P$org<-sub("/.*","",P$model)
fp<-cfa(m,data=P,estimator="MLR",cluster="org",std.lv=TRUE); ld<-standardizedSolution(fp); ld<-setNames(ld$est.std[ld$op=="=~"],ld$rhs[ld$op=="=~"])[vl]
idx<-match(x$base,x$model); ok<-x$gen>=1 & !is.na(idx) & x$model!=x$base & isF(x$type) & x$merged=="False" & isP(x$type[idx]); ok[is.na(ok)]<-FALSE
C<-x[ok,]; G<-as.matrix(C[,vl])-as.matrix(x[idx[ok],vl]); sdP<-sapply(P[,vl],sd)
mcv<-function(g) cor(ld/sqrt(rel), g/sqrt(rel))
dg<-colMeans(G)/sdP; r0<-mcv(dg); set.seed(20261113); bs<-unique(C$base)
rb<-replicate(2000,{s<-sample(bs,replace=TRUE); ii<-unlist(lapply(s,function(z) which(C$base==z))); mcv(colMeans(G[ii,,drop=FALSE])/sdP)})
perms<-function(n){ if(n==1) return(matrix(1)); p<-perms(n-1); do.call(rbind,lapply(1:n,function(i) cbind(i,ifelse(p>=i,p+1,p))))}
Pm<-perms(6); null<-apply(Pm,1,function(o) cor(ld[o]/sqrt(rel[o]), dg/sqrt(rel)))
cat("strict pairs",nrow(G),"MCV corrected r",round(r0,3),"uncorrected",round(cor(ld,dg),3),"CI",round(quantile(rb,c(.025,.975)),3),"perm p",mean(abs(null)>=abs(r0)-1e-12),"\n")
cat("logit correlations (all 4496)\n"); print(round(cor(x[,vl]),2))
