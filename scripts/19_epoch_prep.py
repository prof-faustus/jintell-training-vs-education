import sys, os, pandas as pd, numpy as np
d, out = sys.argv[1], sys.argv[2]
meta = pd.read_csv(os.path.join(d, 'benchmark_metadata.csv')).set_index('benchmark')
mm = pd.read_csv(os.path.join(d, 'model_metadata.csv'))
L = pd.read_csv(out + '/long_v1.csv')
B = {'GPQA diamond': 'gpqa', 'OTIS Mock AIME 2024-2025': 'aime', 'ARC-AGI': 'arc1', 'ARC-AGI-2': 'arc2', 'DTBench': 'dtb'}
L = L[L.benchmark.isin(B)].copy()
L['p'] = [(r - meta.loc[b, 'random_baseline']) / (meta.loc[b, 'score_ceiling'] - meta.loc[b, 'random_baseline']) for r, b in zip(L.raw, L.benchmark)]
g = mm.dropna(subset=['model_version','model_group']).drop_duplicates('model_version').set_index('model_version')['model_group']
L['group'] = L.model.map(g).fillna(L.model)
W = L.groupby(['group', 'benchmark']).p.max().unstack().rename(columns=B)
info = L.groupby('group').agg(release=('release', 'min'), org=('org', 'first'))
W = W.join(info); W['year'] = pd.to_datetime(W.release, errors='coerce').dt.year
eps = .005; cols = list(B.values())
for c in cols:
    W[c + '_clip'] = ((W[c] <= eps) | (W[c] >= 1 - eps)).astype(float).where(W[c].notna())
    W[c + '_lg'] = np.log(np.clip(W[c], eps, 1 - eps) / (1 - np.clip(W[c], eps, 1 - eps)))
W['nobs'] = W[cols].notna().sum(axis=1)
W = W[(W.nobs >= 3) & W.year.notna()].copy()
W['org'] = W.org.fillna('unknown').str.split(',').str[0]
W.to_csv(out + '/epoch_wide_v1.csv')
print('model groups', len(W), 'orgs', W.org.nunique()); print(W.year.value_counts().sort_index().to_dict())
print('coverage', W[cols].notna().sum().to_dict()); print('clipped share', {c: round(W[c + '_clip'].mean(), 3) for c in cols})
