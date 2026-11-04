"""Shared analysis pipeline (final). Python 3, scikit-learn, pandas, numpy, scipy."""
import re, warnings, numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.dummy import DummyClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score
warnings.filterwarnings("ignore")
FEATS = ['MW_g_mol','LogP','pKa','HBD','HBA','TPSA_A2','MP_C','Drug_Hansen_dD','Drug_Hansen_dP','Drug_Hansen_dH',
         'Stabilizer_Hansen_dD','Stabilizer_Hansen_dP','Stabilizer_Hansen_dH','Ra_Hansen_Distance','charged']
CONV = ['MW_g_mol','LogP','pKa','HBD','HBA','TPSA_A2','MP_C','charged']
HSPC = ['Drug_Hansen_dD','Drug_Hansen_dP','Drug_Hansen_dH','Stabilizer_Hansen_dD','Stabilizer_Hansen_dP','Stabilizer_Hansen_dH']
def load(path='Supplementary_Dataset_S1_v2.csv'):
    d = pd.read_csv(path, dtype={'ID': str})
    d['y'] = (d.Stability_Outcome == 'Unstable').astype(int)
    d['drug'] = d.Drug_Name.str.replace(r' \(.*\)', '', regex=True)
    d['charged'] = d.Stabilizer_Charge.str.lower().str.contains('anionic|cationic').astype(int)
    for c in ['pKa', 'MP_C']:   # free-text fields: first numeric value is used; non-numeric -> missing (imputed in-fold)
        d[c] = pd.to_numeric(d[c].astype(str).str.extract(r'^\s*(-?\d+\.?\d*)')[0], errors='coerce')
    return d
def model(name, seed=42, trees=300):
    if name == 'RF': return RandomForestClassifier(n_estimators=trees, max_depth=5, class_weight='balanced', random_state=seed, n_jobs=1)
    if name == 'GB': return GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, random_state=seed)
    return DummyClassifier(strategy='most_frequent')
def oof(X, y, groups, name, seed=42, trees=300):
    """Out-of-fold probabilities. groups = row ids (LOOCV) or drug names (leave-one-drug-out).
    Median imputation is fitted inside each training fold; no other preprocessing or feature selection is used."""
    p = np.zeros(len(y))
    for g in pd.unique(groups):
        te = groups == g; tr = ~te
        if name == 'Dummy' or y[tr].sum() == 0:
            p[te] = 0.0; continue
        m = make_pipeline(SimpleImputer(strategy='median'), model(name, seed, trees)).fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p
def metrics(y, p, thr=0.5):
    pred = (p >= thr).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
    auc = roc_auc_score(y, p) if len(np.unique(p)) > 1 else 0.5
    return dict(acc=(tp+tn)/len(y), F1=f1_score(y, pred, zero_division=0), AUC=auc, sens=tp/(tp+fn), spec=tn/(tn+fp), TP=tp, FN=fn, FP=fp, TN=tn)
