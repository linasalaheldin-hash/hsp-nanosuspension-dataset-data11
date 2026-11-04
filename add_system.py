import pandas as pd,re,sys
M={'tpgs':'TPGS','d-alpha-tocopheryl polyethylene glycol':'TPGS','pvp k30':'PVP','pvp (10%)':'PVP','pvp k25':'PVP','polyvinylpyrrolidone':'PVP','pluronic f127':'Poloxamer 407','poloxamer 407':'Poloxamer 407','poloxamer p407':'Poloxamer 407','poloxamer 188':'Poloxamer 188','polysorbate 80':'Polysorbate 80','tween 80':'Polysorbate 80','sodium dodecyl sulfate':'SDS','sds':'SDS','hpmc':'HPMC','polyvinyl alcohol':'PVA'}
FIX=[('hpc + pluronic',['HPC','Poloxamer 188','SDS']),('hpc ssl',['HPC','Poloxamer 407']),('kollicoat',['Kollicoat MAE 100P','Polysorbate 80']),('pvp va64',['PVP VA64','SDS']),('tween 80 + span',['Polysorbate 80','Span 80']),('hydroxypropyl beta',['HP-beta-CD']),('ethyl',['Ethyl cellulose']),('peg',['PEG 400']),('soluplus',['Soluplus']),('tetronic',['Poloxamine 1107']),('poloxamer 188 + chitosan',['Poloxamer 188','Chitosan'])]
def system(s):
    s0=s.lower()
    for k,v in FIX:
        if s0.startswith(k): return ' + '.join(sorted(v))
    out=[]
    for p in re.split(r'\s\+\s',s0):
        p2=re.sub(r'\s*\(.*$','',p) if not p.startswith('pluronic') else p
        for k,v in M.items():
            if p2.startswith(k) or p.startswith(k): out.append(v); break
        else: raise SystemExit('unmapped '+p)
    return ' + '.join(sorted(set(out)))
f=sys.argv[1]; d=pd.read_csv(f,dtype=str,keep_default_na=False)
if 'Stabilizer_System' in d: d=d.drop(columns='Stabilizer_System')
d.insert(d.columns.get_loc('Stabilizer')+1,'Stabilizer_System',d.Stabilizer.map(system))
d.to_csv(f,index=False); print(f,d.shape,d.Stabilizer_System.nunique(),d[d.ID.astype(int) <= 62].Stabilizer_System.nunique())
