import json, numpy as np
from scipy.stats import mannwhitneyu
from pipe import *
d = load(); d['Endpoint_category'] = d['Endpoint_Category']; d['Formulation_type'] = d['Formulation_Type']
# Ra recomputation check
def ra(df): return np.sqrt(4*(df.Drug_Hansen_dD-df.Stabilizer_Hansen_dD)**2+(df.Drug_Hansen_dP-df.Stabilizer_Hansen_dP)**2+(df.Drug_Hansen_dH-df.Stabilizer_Hansen_dH)**2)
print('max |Ra recomputed - Ra file| =', float((ra(d)-d.Ra_Hansen_Distance).abs().max()))
out = {}
y = d.y.values; g = d.drug.values; rid = d.ID.values
def run_models(df, feats, seeds=(42,), label=''):
    X = df[feats].values.astype(float); yy = df.y.values; gg = df.drug.values; rr = np.arange(len(df))
    res = {}
    for cvn, grp in (('LOOCV', rr), ('LODO', gg)):
        for nm in ('RF', 'GB', 'Dummy'):
            if nm == 'Dummy':
                m = metrics(yy, np.zeros(len(yy))); m['AUC'] = 0.5; res[(cvn, nm)] = m; continue
            ms = []
            for sd in (seeds if nm == 'RF' else (42,)):
                p = oof(X, yy, grp, nm, sd); ms.append(metrics(yy, p)); 
            m = ms[0].copy()
            if len(ms) > 1:
                for k in ('acc','F1','AUC','sens','spec'):
                    v = [q[k] for q in ms]; m[k+'_min'] = min(v); m[k+'_max'] = max(v); m[k+'_mean'] = float(np.mean(v))
            res[(cvn, nm)] = m
    return res
# ---- Table 5 (final pipeline)
t5 = run_models(d, FEATS, seeds=range(42, 52))
rows = []
for (cvn, nm), m in t5.items():
    rows.append(dict(cv=cvn, model=nm, **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()}))
pd.DataFrame(rows).to_csv('table5_final.csv', index=False); print(pd.DataFrame(rows).to_string())
# ---- cluster bootstrap CI for LODO AUC (seed 42 predictions)
X = d[FEATS].values.astype(float)
for nm in ('RF', 'GB'):
    p = oof(X, y, g, nm, 42); rng = np.random.default_rng(7); ds = pd.unique(g); idx = {k: np.where(g == k)[0] for k in ds}; a = []
    for _ in range(4000):
        gs = rng.choice(ds, len(ds)); ii = np.concatenate([idx[k] for k in gs])
        if y[ii].sum() == 0 or (y[ii] == 0).sum() == 0: continue
        a.append(roc_auc_score(y[ii], p[ii]))
    print(nm, 'LODO AUC drug-cluster bootstrap 95% CI', np.percentile(a, [2.5, 97.5]).round(3), 'n valid', len(a))
# ---- feature-set hierarchy (LODO and LOOCV AUC)
H = {'A Ra only': ['Ra_Hansen_Distance'], 'B HSP components (6) + Ra': HSPC + ['Ra_Hansen_Distance'], 'C conventional descriptors': CONV, 'D combined (all 15)': FEATS}
hr = []
for k, f in H.items():
    r = run_models(d, f)
    for nm in ('RF', 'GB'):
        hr.append(dict(set=k, model=nm, LOOCV_AUC=round(r[('LOOCV', nm)]['AUC'], 3), LODO_AUC=round(r[('LODO', nm)]['AUC'], 3), LODO_found=f"{r[('LODO', nm)]['TP']}/11", LODO_sens=round(r[('LODO', nm)]['sens'], 3), LODO_spec=round(r[('LODO', nm)]['spec'], 3)))
pd.DataFrame(hr).to_csv('table_feature_hierarchy.csv', index=False); print(pd.DataFrame(hr).to_string())
# ---- direct Ra statistics
rng = np.random.default_rng(1)
def rb(u, s): return 2*mannwhitneyu(u, s, alternative='two-sided').statistic/(len(u)*len(s)) - 1
def direct(df, col='Ra_Hansen_Distance', label=''):
    u = df[df.y == 1][col].values; s = df[df.y == 0][col].values
    t = mannwhitneyu(u, s, alternative='two-sided'); gr = mannwhitneyu(u, s, alternative='greater')
    ds = df.drug.unique(); idx = {k: df.index[df.drug == k].values for k in ds}; bs = []
    for _ in range(4000):
        gs = rng.choice(ds, len(ds)); b = df.loc[np.concatenate([idx[k] for k in gs])]
        if b.y.sum() == 0 or (b.y == 0).sum() == 0: continue
        bs.append(rb(b[b.y == 1][col].values, b[b.y == 0][col].values))
    dm = df.groupby(['drug', 'y'])[col].mean().reset_index(); du = dm[dm.y == 1][col].values; dsv = dm[dm.y == 0][col].values
    dl = mannwhitneyu(du, dsv, alternative='two-sided')
    unst_means_vs_rows = mannwhitneyu(du, s, alternative='two-sided').pvalue
    return dict(label=label, nU=len(u), nS=len(s), medU=round(float(np.median(u)), 2), medS=round(float(np.median(s)), 2), U=float(t.statistic), p_two=round(t.pvalue, 3), p_one_U_gt_S=round(gr.pvalue, 3), r_rb=round(rb(u, s), 2), CI_lo=round(float(np.percentile(bs, 2.5)), 2), CI_hi=round(float(np.percentile(bs, 97.5)), 2), nUdrugs=len(du), nSdrugs=len(dsv), p_drug_level=round(dl.pvalue, 3), r_drug_level=round(rb(du, dsv), 2), p_4means_vs_43rows=round(unst_means_vs_rows, 3), Unstable_drug_means=[round(float(x), 2) for x in du])
res = [direct(d, label='Baseline (54 rows)')]
print(res[0])
# ---- sensitivity analyses
DHP = ['Nifedipine','Felodipine','Nitrendipine','Nicardipine','Lercanidipine']
def with_ra(df, newra):
    e = df.copy(); e['Ra_Hansen_Distance'] = newra; return e
def sens(label, df, shift_note=''):
    r = direct(df, label=label)
    m = run_models(df.reset_index(drop=True), FEATS)
    r.update(RF_LODO_AUC=round(m[('LODO','RF')]['AUC'],3), RF_LODO_found=f"{m[('LODO','RF')]['TP']}/{int(df.y.sum())}", GB_LODO_AUC=round(m[('LODO','GB')]['AUC'],3), GB_LODO_found=f"{m[('LODO','GB')]['TP']}/{int(df.y.sum())}", RF_LOOCV_found=f"{m[('LOOCV','RF')]['TP']}/{int(df.y.sum())}")
    return r
S = []
S.append(sens('Baseline', d.reset_index(drop=True)))
for sh in (-1.0, -1.5):
    e = d.copy(); mk = e.drug.isin(DHP); e.loc[mk, 'Ra_Hansen_Distance'] += sh; S.append(sens(f'DHP amide/lactam alternative, uniform shift {sh}', e.reset_index(drop=True)))
e = d.copy(); e.loc[e.ID == '11', 'Ra_Hansen_Distance'] = 5.18; e.loc[e.ID == '19', 'Ra_Hansen_Distance'] = 6.98
mk = e.drug.isin(['Nitrendipine','Nicardipine','Lercanidipine']); e.loc[mk, 'Ra_Hansen_Distance'] += -1.25
S.append(sens('DHP alternative: reported nifedipine/felodipine values; other three -1.25', e.reset_index(drop=True)))
S.append(sens('Excluding glimepiride (processing-stage endpoint)', d[d.Drug_Name != 'Glimepiride'].reset_index(drop=True)))
S.append(sens('Excluding Felodipine/PVA (ID 19, relative-comparison record)', d[d.ID != '19'].reset_index(drop=True)))
S.append(sens('Excluding glimepiride and ID 19', d[(d.Drug_Name != 'Glimepiride') & (d.ID != '19')].reset_index(drop=True)))
st = d[(d.y == 1) | d.Endpoint_category.fillna('').str.startswith('Stable on stated storage')].reset_index(drop=True)
S.append(sens('Stable class restricted to entries with stated storage conditions', st))
S.append(sens('Excluding glimepiride and the four ibuprofen screening rows (indomethacin and UDCA failures only)', d[~d.drug.isin(['Glimepiride']) & ~((d.drug=='Ibuprofen') & (d.y==1))].reset_index(drop=True)))
S.append(sens('Restricted to entries with a recorded formulation route', d[~d.Formulation_type.fillna('').str.contains('Not recorded')].reset_index(drop=True)))
# source substitution: copovidone (celecoxib, ID 8) and Soluplus (ID 67) polar component, manufacturer-type value instead of measured
e = d.copy()
e.loc[e.ID == '8', 'Stabilizer_Hansen_dP'] = 0.5/2 + 9.11 - 10.2/2   # 1:1 mix: swap the copovidone polar component 10.2 -> 0.5
e.loc[e.ID == '67', 'Stabilizer_Hansen_dP'] = 0.3
e['Ra_Hansen_Distance'] = np.where(e.ID.isin(['8','67']), ra(e), e.Ra_Hansen_Distance)
print('substituted Ra:', e[e.ID.isin(['8','67'])][['ID','Ra_Hansen_Distance']].round(2).to_dict('records'), 'orig', d[d.ID.isin(['8','67'])][['ID','Ra_Hansen_Distance']].to_dict('records'))
S.append(sens('Alternative polar component for copovidone (ID 8) and Soluplus (ID 67)', e.reset_index(drop=True)))
sd = pd.DataFrame(S); sd.to_csv('table_sensitivity.csv', index=False); print(sd.drop(columns=['Unstable_drug_means']).to_string())
json.dump(dict(S=S), open('sens.json', 'w'), default=str)
