import sys, pandas as pd
# Normalised (chance-corrected) leaderboard scores, 0-100, for the analysed models; used by 30_review_tests2.R to vary the logit clipping constant.
d = pd.read_parquet(sys.argv[1]); x = pd.read_csv(sys.argv[2]); out = sys.argv[3]
pairs = [('IFEval','ifeval'),('BBH','bbh'),('MATH Lvl 5','math'),('GPQA','gpqa'),('MUSR','musr'),('MMLU-PRO','mmlupro')]
d['prec_rank'] = d.Precision.map({'bfloat16': 0, 'float16': 1, '4bit': 2})
d = d.sort_values('prec_rank').drop_duplicates('fullname', keep='first').set_index('fullname')
y = pd.DataFrame({'model': x.model})
for n, c in pairs: y[c + '_norm'] = d.loc[x.model, n].values
y.to_csv(f'{out}/hf_norm_v1.csv', index=False); print(len(y))
