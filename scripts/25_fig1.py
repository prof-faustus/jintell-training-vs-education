import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({'font.family':'serif','font.serif':['TeX Gyre Pagella','Palatino','DejaVu Serif'],'font.size':9})
x=pd.read_csv('derived/hf_models_v2.csv'); v=['ifeval','bbh','math','gpqa','musr','mmlupro']; vl=[c+'_lg' for c in v]
lab=['IFEval','BBH','MATH Lvl 5','GPQA','MuSR','MMLU-Pro']
isP=x.type.isin(['pretrained','continuously pretrained']); isF=x.type.isin(['fine-tuned on domain-specific datasets','chat models (RLHF, DPO, IFT, ...)'])
pos={m:i for i,m in enumerate(x.model)}; bi=x.base.map(pos)
ok=(x.gen>=1)&bi.notna()&(x.model!=x.base)&isF&(x.merged.astype(str)=='False')
ok&=bi.fillna(0).astype(int).map(lambda i: isP.iloc[i])
C=x[ok]; B=x.iloc[bi[ok].astype(int)]
G=pd.DataFrame(C[vl].values-B[vl].values,columns=vl); G['base']=C.base.values; L=G.groupby('base')[vl].mean()
sdP=x.loc[isP,vl].std(); n=len(L); tc=stats.t.ppf(.975,n-1)
m=L.mean()/sdP; se=L.std()/np.sqrt(n)/sdP
lam=[0.622,0.979,0.807,0.838,0.638,0.979]   # M4 holdout standardised loadings (out/b_select_logit_v1.txt)
fig,ax=plt.subplots(1,2,figsize=(6.6,2.6))
ax[0].bar(lab,lam,color='0.55',edgecolor='black',linewidth=.6); ax[0].set_ylim(0,1.05); ax[0].set_ylabel('Standardised loading on g'); ax[0].set_title('(a) g loadings (holdout half)')
ax[1].bar(lab,m.values,color='0.8',edgecolor='black',linewidth=.6,yerr=tc*se.values,capsize=3,error_kw={'linewidth':.8})
ax[1].axhline(0,color='black',linewidth=.6); ax[1].set_ylabel('Gain (pretrained SD units)'); ax[1].set_title('(b) Within-lineage fine-tuning gain')
for a in ax: a.tick_params(axis='x',rotation=40); [s.set_linewidth(.6) for s in a.spines.values()]
plt.tight_layout(); plt.savefig('fig/fig1_v1.pdf'); plt.savefig('fig/fig1_v1.png',dpi=200)
print('pairs',len(C),'lineages',n); print(pd.DataFrame({'gain_sd':m,'lo':m-tc*se,'hi':m+tc*se}).round(3))
