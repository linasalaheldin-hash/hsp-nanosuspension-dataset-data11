import sys, numpy as np
from pipe import *
mode, seed, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = load(); X = d[FEATS].values.astype(float); g = d.drug.values; y = d.y.values
drugs = pd.unique(g); prof = [4, 4, 2, 1]   # observed Unstable-row counts of the four failing drugs
rng = np.random.default_rng(seed)
def auc(yy): return roc_auc_score(yy, oof(X, yy, g, 'RF', 42, 300))
if mode == 'obs': print(auc(y)); sys.exit()
for i in range(N):
    if mode == 'group2':   # drug-preserving: exactly 11 Unstable rows on 4 drugs, with the observed profile 4, 4, 2, 1
        yy = np.zeros(len(y), int); used = set()
        for c in prof:   # largest first; a drug is eligible only if it has at least c formulations
            elig = [k for k in drugs if k not in used and (g == k).sum() >= c]
            k = rng.choice(elig); used.add(k)
            yy[rng.choice(np.where(g == k)[0], c, replace=False)] = 1
    else:
        yy = rng.permutation(y)
    print(auc(yy), int(yy.sum()), flush=True)
