# Q26 injection-recovery in the uv plane: sky model (N+S points + 10" Gaussian) sampled at the REAL uv points/weights of
# Xb55 (2-1), Xb4f (3-2), Xb4d (3-2), noise = real sigma (x noise_scale), then the SAME binning (q26_uvana.binstat/comb) and fit.
# Truth: R_ext(2-1)=0.54; R_ext(3-2)= 0.54 (no fall) or 0.46 (15% fall). Core ratios from the data fits (Q22-like): N .54/.40, S .93/.73.
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts'); import q26_uvana as U
ARC=U.ARC; rng=np.random.default_rng(26)
FL={'2-1':dict(HCO=[18.6,15.7,80.4],RN=.54,RS=.93),'3-2':dict(HCO=[25.5,26.8,85.8],RN=.40,RS=.73)}
def model(d,fw,fl,Rext):
    X=U.design(d,fw); hco=np.array(fl['HCO']); hcn=hco*np.array([fl['RN'],fl['RS'],Rext])
    return X@hcn, X@hco
def sim(sets,Rext,fwt,scale,fwfit_list=(6,10,15)):
    out={}
    for o,J in sets:
        Ld={sp:U.load(o,sp) for sp in('HCN','HCOp')}; sims={'HCN':[],'HCOp':[]}
        for i in range(len(Ld['HCN'])):
            for sp,(k) in (('HCN',0),('HCOp',1)):
                d=Ld[sp][i]; m=model(d,fwt,FL[J],Rext[J])[k]
                e=dict(d); e['re']=m.real+rng.normal(0,d['sg'])*scale; e['im']=m.imag+rng.normal(0,d['sg'])*scale
                sims[sp].append(e)
        out[o]={}
        for lo,hi in zip(U.BINS[:-1],U.BINS[1:]):
            if o=='Xb4d' and lo<10: continue
            h,_=U.comb(sims['HCN'],lo,hi); q,_=U.comb(sims['HCOp'],lo,hi)
            if h and q: out[o][f'{lo}-{hi}']=h['re']/q['re']
        for fw in fwfit_list:
            ul=8 if o!='Xb4d' else 10
            ph=U.fit(sims['HCN'],fw,ul)[0]; pq=U.fit(sims['HCOp'],fw,ul)[0]
            out[o][f'fit{fw}']=ph[2]/pq[2]
    return out
if __name__=='__main__':
    sets=[('Xb55','2-1'),('Xb4f','3-2')]; NR=int(sys.argv[1]) if len(sys.argv)>1 else 100
    res={}
    for fwt in (10,):
      for name,R32 in (('noFall',0.54),('fall15',0.459)):
        for scale in (1.0,):
            acc=[sim(sets,{'2-1':0.54,'3-2':R32},fwt,scale) for _ in range(NR)]
            keys=acc[0]['Xb55'].keys()
            r={}
            for k in acc[0]['Xb55']:
                if k in acc[0]['Xb4f']:
                    a=np.array([x['Xb55'][k] for x in acc]); b=np.array([x['Xb4f'][k] for x in acc]); q=b/a
                    r[k]=(float(a.mean()),float(b.mean()),float(q.mean()),float(q.std()))
                    print(f'truth={name} fwTrue={fwt} {k:8s} R21={a.mean():.3f} R32={b.mean():.3f} R32/R21={q.mean():.3f}+-{q.std():.3f}')
            res[f'{name}|{fwt}']=r
    json.dump(res,open('/workspace/out/q26/q26_inject.json','w'),indent=1)
