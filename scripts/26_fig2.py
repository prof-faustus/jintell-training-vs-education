import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
plt.rcParams.update({'font.family':'serif','font.serif':['TeX Gyre Pagella','Palatino','DejaVu Serif'],'font.size':8})
R=pd.read_csv('out/e_mcnull_role_exact_v1_reps.csv'); C=pd.read_csv('out/e_mcnull_cohort_exact_v1_reps.csv')
obs={('role','m'):-0.0122,('role','s'):-0.0150,('coh','m'):-0.0267,('coh','s'):-0.0392}
fig,ax=plt.subplots(2,2,figsize=(6.6,4.0))
for i,(nm,D) in enumerate([('role',R),('coh',C)]):
    for j,(st,col) in enumerate([('m','dCFI_m'),('s','dCFI_s')]):
        a=ax[i,j]; x=D[col].dropna(); o=obs[(nm,st)]
        lo=min(x.min(),o,-0.01)*1.15; hi=max(x.max(),0.002)*1.15
        a.hist(x,bins=np.linspace(lo,hi,45),color='0.7',edgecolor='0.3',linewidth=.4)
        a.axvline(o,color='black',linewidth=1.4); a.axvline(-0.01,color='black',linewidth=.8,linestyle='--')
        a.set_title(('(a) Training role' if nm=='role' else '(c) Cohort' if st=='m' else '(d) Cohort') if False else
                    {('role','m'):'(a) Training role, metric',('role','s'):'(b) Training role, scalar',('coh','m'):'(c) Cohort, metric',('coh','s'):'(d) Cohort, scalar'}[(nm,st)])
        a.set_xlabel(r'$\Delta$CFI'); a.set_ylabel('Replicates'); [s.set_linewidth(.6) for s in a.spines.values()]
        p=(x<=o).mean(); a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5)); a.set_title(a.get_title()+(f' (null p = {p:.3f})' if p>0 else f' (null p < {1/len(x):.3f})'),fontsize=8)
plt.tight_layout(); plt.savefig('fig/fig2_v1.pdf'); plt.savefig('fig/fig2_v1.png',dpi=200)
print(len(R.dropna()),len(C.dropna()))
