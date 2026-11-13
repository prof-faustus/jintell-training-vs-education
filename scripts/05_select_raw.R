suppressMessages({library(lavaan)})
args<-commandArgs(TRUE); x<-read.csv(args[1]); od<-args[2]
set.seed(20261113)
# split by lineage so related models do not straddle halves
lin<-unique(x$base); h<-sample(c(1,2),length(lin),replace=TRUE); names(h)<-lin
x$half<-h[x$base]; D<-x[x$half==1,]; H<-x[x$half==2,]
sink(file.path(od,"b_select_v1.txt"))
cat("discovery n=",nrow(D)," holdout n=",nrow(H),"\n")
fm<-c("chisq.scaled","df","pvalue.scaled","cfi.robust","tli.robust","rmsea.robust","srmr","bic")
M<-list(
 M1='g =~ ifeval+bbh+math+gpqa+musr+mmlupro',
 M2='g =~ bbh+math+gpqa+musr+mmlupro',
 M3='g =~ bbh+gpqa+musr+mmlupro+math
     p =~ ifeval+math
     ifeval ~~ 0*ifeval',
 M3b='k =~ bbh+gpqa+musr+mmlupro
      q =~ math+ifeval+mmlupro')
for(n in names(M)){ f<-try(cfa(M[[n]],data=D,estimator="MLR",cluster="base",std.lv=TRUE),silent=TRUE)
 cat("\n==",n,"\n"); if(inherits(f,"try-error")) {cat("FAILED\n");next}; print(round(fitMeasures(f,fm),3)); cat("converged",lavInspect(f,"converged"),"\n")}
# M4: residual covariances added sequentially on discovery half only, stop when MI<50 or 3 added
m<-M$M1; added<-c()
for(i in 1:3){ f<-cfa(m,data=D,estimator="MLR",cluster="base",std.lv=TRUE); mi<-modindices(f); mi<-mi[mi$op=="~~",]; mi<-mi[order(-mi$mi),]
 if(mi$mi[1]<50) break; t<-paste(mi$lhs[1],"~~",mi$rhs[1]); added<-c(added,t); m<-paste(m,"\n",t)}
cat("\n== M4 added:",added,"\n"); f4<-cfa(m,data=D,estimator="MLR",cluster="base",std.lv=TRUE); print(round(fitMeasures(f4,fm),3))
cat("\n=== HOLDOUT CONFIRMATION ===\n")
for(n in c("M2")){ f<-cfa(M[[n]],data=H,estimator="MLR",cluster="base",std.lv=TRUE); cat("\n",n,"\n"); print(round(fitMeasures(f,fm),3)); print(standardizedSolution(f)[1:5,c("rhs","est.std","se")])}
f<-cfa(m,data=H,estimator="MLR",cluster="base",std.lv=TRUE); cat("\nM4 holdout\n"); print(round(fitMeasures(f,fm),3)); print(standardizedSolution(f)[1:6,c("rhs","est.std","se")])
writeLines(m,file.path(od,"m4_syntax_v1.txt"))
sink()
