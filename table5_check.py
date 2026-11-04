"""Quick regeneration of Table 5 (seed 42): row-level LOOCV and leave-one-drug-out for RF, GB and Dummy."""
from pipe import *
d = load(); X = d[FEATS].values.astype(float); y = d.y.values
rows = []
for cv, groups in (('LOOCV', np.arange(len(d)).astype(str)), ('LODO', d.drug.values)):
    for name in ('RF', 'GB', 'Dummy'):
        m = metrics(y, oof(X, y, groups, name, 42, 300))
        rows.append(dict(cv=cv, model=name, **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()}))
out = pd.DataFrame(rows); print(out.to_string(index=False)); out.to_csv('table5_check.csv', index=False)
