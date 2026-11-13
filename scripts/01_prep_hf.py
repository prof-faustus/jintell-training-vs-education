import sys,pandas as pd,numpy as np
d=pd.read_parquet(sys.argv[1]); out=sys.argv[2]
raw={'IFEval Raw':'ifeval','BBH Raw':'bbh','MATH Lvl 5 Raw':'math','GPQA Raw':'gpqa','MUSR Raw':'musr','MMLU-PRO Raw':'mmlupro'}
d['prec_rank']=d.Precision.map({'bfloat16':0,'float16':1,'4bit':2})
d=d.sort_values('prec_rank').drop_duplicates('fullname',keep='first')
d=d[~d.Flagged]
x=d.rename(columns=raw)[['fullname','Base Model','Generation','Type','Merged','Official Providers','#Params (B)','MoE','Upload To Hub Date','Submission Date','Architecture']+list(raw.values())]
x.columns=['model','base','gen','type','merged','official','params','moe','upload','submit','arch']+list(raw.values())
x['type']=x.type.str.replace(r'^\W+\s*','',regex=True)
x['upload_year']=pd.to_datetime(x.upload,errors='coerce').dt.year
x.to_csv(f'{out}/hf_models_v1.csv',index=False)
# pairs
names=x.set_index('model')
ft=x[(x.gen>=1)&x.base.isin(names.index)&(x.model!=x.base)].copy()
b=names.loc[ft.base, list(raw.values())].reset_index(drop=True)
ft=ft.reset_index(drop=True)
for c in raw.values(): ft[c+'_base']=b[c].values
ft.to_csv(f'{out}/hf_pairs_v1.csv',index=False)
print('models',len(x),'pairs',len(ft),'bases',ft.base.nunique())
print(x.type.value_counts().to_dict())
pt=x[x.type.isin(['pretrained','continuously pretrained'])]
print('pretrained by year',pt.upload_year.value_counts().sort_index().to_dict())
print('all by year',x.upload_year.value_counts().sort_index().to_dict())
print('pairs per base', ft.groupby('base').size().describe().to_dict())
