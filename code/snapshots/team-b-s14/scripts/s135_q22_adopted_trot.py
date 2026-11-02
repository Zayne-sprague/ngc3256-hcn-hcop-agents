# Q22: adopted blend integrals with the empirical edge-ramp baseline (s134) and propagated Trot (same thin-LTE model/254 inputs as s120c).
import json, numpy as np
rng=np.random.default_rng(2)
E=json.load(open('/workspace/out/q22/s134.json')); B=json.load(open('/workspace/out/q9/s120b.json'))
S21=json.load(open('/workspace/out/s21_ch3oh254.json')); sm=S21['smoothed']; mad=S21['smoothed_mad']
L254=(254.015377,20.09,-4.72074,5)
BL=[(165.050175,23.37,-4.62911,3),(165.061130,28.01,-4.63041,5),(165.099240,34.98,-4.63216,7)]
def W(l,T): f,Ee,la,g=l; return g*10**la*np.exp(-Ee/T)/f**2
Tg=np.linspace(3,60,1141); Rg=W(L254,Tg)/sum(W(b,Tg) for b in BL); nuc=(163.53/251.68)**2
out={}
for reg in ('N','S','Bridge','W'):
    o=E[f'{reg}_ap']; vals=np.array([[o[w][m] for m in 'RLQ'] for w in ('150','200','250','300')])  # W x method
    a200=np.array([o['200'][m] for m in 'RLQ']); stat=o['200']['stat']
    adopted=a200.mean(); sysm=0.5*(vals.max()-vals.min())
    old=B[f'{reg}_ap']; h4=[r[2] for r in old['rows'] if r[0]==400]; oldlo,oldhi,oldmed=min(h4),max(h4),old['med']   # old 5-baseline range at the +-400 km/s window
    d=dict(I_adopted=float(adopted),stat=float(stat),sys_half_range=float(sysm),range=[float(vals.min()),float(vals.max())],
           zero_baseline_bound=float(o['200']['Z']),old_range_s120b=[oldlo,oldhi],old_median=oldmed,
           old_range_ratio=float(oldhi/oldlo),new_range_ratio=float(vals.max()/vals.min()),scale_vs_old_median=float(adopted/oldmed))
    s=sm.get(reg)
    if s is not None:
        a=s['I254_55']; ea=np.hypot(mad,0.05*a); n=20000
        draws=rng.choice(vals.ravel(),n)+rng.normal(0,stat,n); as_=a+rng.normal(0,ea,n); r=as_/draws*nuc; ok=(draws>0)&(r>0)
        T=np.interp(r[ok],Rg[::-1],Tg[::-1]); d.update(I254=float(a),eI254=float(ea),R_K=float(np.median(r[ok])),R_K_16_84=[float(x) for x in np.percentile(r,[16,84])],Trot_16_50_84=[float(x) for x in np.percentile(T,[16,50,84])],frac_at_60K_cap=float((T>=59.9).mean()))
    out[reg]=d; print(reg,json.dumps({k:(np.round(v,3).tolist() if not isinstance(v,float) else round(v,3)) for k,v in d.items()}))
json.dump(out,open('/workspace/out/q22/s135.json','w'),indent=1)
