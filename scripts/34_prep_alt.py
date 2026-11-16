import sys, io, pandas as pd, numpy as np
# Alternative constructions of the analysis file, to test the two sample-construction choices of 01_prep_hf.py:
#   main    : one record per model preferring bfloat16, community-flagged models removed (reproduces derived/hf_models_v2.csv)
#   flagged : as main, but community-flagged models kept
#   f16     : as main, but float16 preferred over bfloat16 for models evaluated at both precisions
# args: parquet outdir variant
d = pd.read_parquet(sys.argv[1]); out = sys.argv[2]; variant = sys.argv[3]
raw = {'IFEval Raw':'ifeval','BBH Raw':'bbh','MATH Lvl 5 Raw':'math','GPQA Raw':'gpqa','MUSR Raw':'musr','MMLU-PRO Raw':'mmlupro'}
rank = {'bfloat16':0,'float16':1,'4bit':2} if variant != 'f16' else {'float16':0,'bfloat16':1,'4bit':2}
d['prec_rank'] = d.Precision.map(rank)
d = d.sort_values('prec_rank').drop_duplicates('fullname', keep='first')
if variant != 'flagged': d = d[~d.Flagged]
x = d.rename(columns=raw)[['fullname','Base Model','Generation','Type','Merged','Official Providers','#Params (B)','MoE','Upload To Hub Date','Submission Date','Architecture']+list(raw.values())]
x.columns = ['model','base','gen','type','merged','official','params','moe','upload','submit','arch']+list(raw.values())
x['type'] = x.type.str.replace(r'^\W+\s*','',regex=True)
x['upload_year'] = pd.to_datetime(x.upload, errors='coerce').dt.year
x = pd.read_csv(io.StringIO(x.to_csv(index=False)))   # same write-and-reread as 01_prep_hf.py -> 02_prep_logit.py
pairs = [('IFEval','ifeval'),('BBH','bbh'),('MATH Lvl 5','math'),('GPQA','gpqa'),('MUSR','musr'),('MMLU-PRO','mmlupro')]
dd = d.set_index('fullname'); eps = 0.005
for n, c in pairs:
    p = dd.loc[x.model, n].values/100.0
    x[c+'_clip'] = ((p <= eps) | (p >= 1-eps)).astype(int)
    p = np.clip(p, eps, 1-eps); x[c+'_lg'] = np.log(p/(1-p))
x['any_clip'] = x[[c+'_clip' for _, c in pairs]].max(axis=1)
x.to_csv(f'{out}/hf_models_v2_{variant}.csv', index=False)
print(variant, 'models', len(x), 'flagged kept', int(d.Flagged.sum()), 'precision', d.Precision.value_counts().to_dict())
