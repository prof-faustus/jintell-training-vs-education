import sys,pandas as pd,numpy as np
d=pd.read_parquet(sys.argv[1]); eps=0.005
cols=[('IFEval','ifeval'),('BBH','bbh'),('MATH Lvl 5','math'),('GPQA','gpqa'),('MUSR','musr'),('MMLU-PRO','mmlupro')]
dup=d.groupby('fullname').filter(lambda g: len(g)==2).sort_values(['fullname','Precision'])
A=dup[dup.Precision=='bfloat16'].set_index('fullname'); B=dup[dup.Precision=='float16'].set_index('fullname'); B=B.loc[A.index]
out=[]
for n,c in cols:
    a=np.log(np.clip(A[n]/100,eps,1-eps)/(1-np.clip(A[n]/100,eps,1-eps))); b=np.log(np.clip(B[n]/100,eps,1-eps)/(1-np.clip(B[n]/100,eps,1-eps)))
    out.append((c,np.corrcoef(a,b)[0,1],len(a)))
pd.DataFrame(out,columns=['bench','retest_r','n']).to_csv(sys.argv[2],index=False); print(out)

# Difference vectors (float16 minus bfloat16) of the same models, used as the negative control in 22_g_change_controls.R
lg = lambda s: np.log(np.clip(s/100, eps, 1-eps) / (1 - np.clip(s/100, eps, 1-eps)))
D = pd.DataFrame({c + '_d': lg(B[n]).values - lg(A[n]).values for n, c in cols}, index=A.index)
D.index.name = 'fullname'; D['base'] = A['Base Model'].values
D.to_csv(sys.argv[3])
