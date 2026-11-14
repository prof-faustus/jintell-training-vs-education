import sys, pandas as pd
# Extra leaderboard fields used by the review tests (28_review_tests.R): evaluation precision,
# chat-template flag and weight type, attached to each model under the same de-duplication rule as 01_prep_hf.py.
d = pd.read_parquet(sys.argv[1]); out = sys.argv[2]
d['prec_rank'] = d.Precision.map({'bfloat16': 0, 'float16': 1, '4bit': 2})
d = d.sort_values('prec_rank').drop_duplicates('fullname', keep='first')
d = d[~d.Flagged]
x = d[['fullname', 'Precision', 'Chat Template', 'Weight type', 'Hub License']].copy()
x.columns = ['model', 'precision', 'chat_template', 'weight_type', 'license']
x.to_csv(f'{out}/hf_extra_v1.csv', index=False)
print('models', len(x)); print(x.precision.value_counts().to_dict()); print(x.weight_type.value_counts().to_dict())
