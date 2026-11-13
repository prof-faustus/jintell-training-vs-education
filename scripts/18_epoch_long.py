import sys, os, pandas as pd, numpy as np
d, out = sys.argv[1], sys.argv[2]
meta = pd.read_csv(os.path.join(d,'benchmark_metadata.csv'))
meta = meta[meta.source_file.notna()]
rows=[]
for _,m in meta.iterrows():
    f=os.path.join(d,m.source_file)
    if not os.path.exists(f): print('MISSING',m.source_file); continue
    df=pd.read_csv(f)
    if m.score_column not in df.columns: print('NOCOL',m.benchmark,list(df.columns)[:8]); continue
    s=pd.to_numeric(df[m.score_column],errors='coerce')*m.scale
    t=pd.DataFrame({'model':df['Model version'],'benchmark':m.benchmark,'raw':s,
                    'release':df.get('Release date'),'org':df.get('Organization'),
                    'compute':pd.to_numeric(df.get('Training compute (FLOP)'),errors='coerce'),
                    'stderr':pd.to_numeric(df['stderr'],errors='coerce') if 'stderr' in df else np.nan})
    rows.append(t.dropna(subset=['raw']))
long=pd.concat(rows)
long.to_csv(os.path.join(out,'long_v1.csv'),index=False)
print('rows',len(long),'models',long.model.nunique())
cov=long.groupby('benchmark').model.nunique().sort_values(ascending=False)
print(cov.to_string())
