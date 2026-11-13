suppressMessages(library(lavaan))
x<-read.csv("derived/hf_models_v2.csv"); v<-c("ifeval","bbh","math","gpqa","musr","mmlupro"); vl<-paste0(v,"_lg")
m4<-paste(readLines("out/m4_syntax_logit_v1.txt"),collapse="\n"); m2f<-paste(readLines("out/m2f_syntax_v1.txt"),collapse="\n")
add<-c("ifeval_lg ~~ gpqa_lg","musr_lg ~~ mmlupro_lg","bbh_lg ~~ math_lg")
isP<-function(t) t %in% c("pretrained","continuously pretrained"); isF<-function(t) t %in% c("fine-tuned on domain-specific datasets","chat models (RLHF, DPO, IFT, ...)")
eqs<-list(configural="",metric="loadings",scalar=c("loadings","intercepts"))
mkrole<-function(d,ptonly=FALSE){ pp<- if(ptonly) d$type=="pretrained" else isP(d$type); d<-d[pp|isF(d$type),]; pp<- if(ptonly) d$type=="pretrained" else isP(d$type); d$grp<-factor(ifelse(pp,"pretrained","finetuned"),levels=c("pretrained","finetuned")); d$cl<-d$base; d}
mkcoh<-function(d){ d<-d[isP(d$type)&!is.na(d$upload_year),]; d$grp<-factor(ifelse(d$upload_year<=2023,"c2022_23",ifelse(d$upload_year==2024,"c2024","c2025"))); d$cl<-sub("/.*","",d$model); d}
dcfi<-function(dd,m){ c<-sapply(eqs,function(e){ f<-cfa(m,data=dd,group="grp",estimator="MLR",cluster="cl",group.equal=e); fitMeasures(f,"cfi.robust")}); round(diff(c),6)}
dup<-duplicated(round(x[,v],4))
sets<-list(baseline=x,no_clip=x[x$any_clip==0,],no_merged_flag=x[x$merged=="False",],no_dup=x[!dup,],all_three=x[x$any_clip==0 & x$merged=="False" & !dup,])
out<-list()
for(n in names(sets)){ out[[paste("role",n)]]<-dcfi(mkrole(sets[[n]]),m4); out[[paste("cohort",n)]]<-dcfi(mkcoh(sets[[n]]),m4)}
out[["role pretrained_only"]]<-dcfi(mkrole(x,TRUE),m4)
b<-mkcoh(x); b<-b[b$params>=0.5 & b$params<=15,]; out[["cohort size_band"]]<-dcfi(b,m4)
for(k in c(1,3)){ mk<-paste(c(m4,add[1:k]),collapse="\n"); out[[paste("role M4+",k)]]<-dcfi(mkrole(x),mk); out[[paste("cohort M4+",k)]]<-dcfi(mkcoh(x),mk)}
out[["role 2F"]]<-dcfi(mkrole(x),m2f); out[["cohort 2F"]]<-dcfi(mkcoh(x),m2f)
print(do.call(rbind,out))
r<-read.csv("out/e_mcnull_cohort_exact_v1_reps.csv"); cat("cohort exact null: share dCFI_m <= -.01:",mean(r$dCFI_m<=-.01,na.rm=TRUE)," dCFI_s <= -.01:",mean(r$dCFI_s<=-.01,na.rm=TRUE),"\n")
r<-read.csv("out/e_mcnull_role_exact_v1_reps.csv"); cat("role exact null: min dCFI_m",min(r$dCFI_m),"min dCFI_s",min(r$dCFI_s),"\n")
