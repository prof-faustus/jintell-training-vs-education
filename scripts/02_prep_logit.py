import sys,pandas as pd,numpy as np
d=pd.read_parquet(sys.argv[1]); x=pd.read_csv(sys.argv[2]); out=sys.argv[3]
pairs=[('IFEval','ifeval'),('BBH','bbh'),('MATH Lvl 5','math'),('GPQA','gpqa'),('MUSR','musr'),('MMLU-PRO','mmlupro')]
d['prec_rank']=d.Precision.map({'bfloat16':0,'float16':1,'4bit':2})
d=d.sort_values('prec_rank').drop_duplicates('fullname',keep='first').set_index('fullname')
eps=0.005
for n,c in pairs:
    p=d.loc[x.model,n].values/100.0   # leaderboard chance-corrected proportion
    x[c+'_clip']=((p<=eps)|(p>=1-eps)).astype(int)
    p=np.clip(p,eps,1-eps); x[c+'_lg']=np.log(p/(1-p))
x['any_clip']=x[[c+'_clip' for _,c in pairs]].max(axis=1)
x.to_csv(f'{out}/hf_models_v2.csv',index=False)
print('any clipped share %.4f'%x.any_clip.mean()); print(x[[c+'_clip' for _,c in pairs]].mean().round(4).to_dict())
