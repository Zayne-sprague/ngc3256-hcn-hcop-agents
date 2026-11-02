# Q2 round 2, fix 2(i): lobe (and red-feature) wing flux in independent imagings of the same data: standard Xb4d (3.9 MHz), hand-made Xb4d tapered (spw27 HCN, spw25 HCO+),
# standard Xb4b smoothed, Xb4b manual pbcor / manual flat-noise HCO+ (26.6 MHz). All convolved to a 1.05" circular beam (flux rescaled by beam area), r=1.0"/1.5" apertures, window edges
# made identical by fractional-channel weights. Baselines: (A) none, (B) mean of off-source random apertures (field) per channel, (C) per-pixel median(-900..-450) (std cubes only).
# Output out/q2r2/s77a.txt/json. usage: python3 -I s77a_flux_systematics.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(11)
def wflux(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
WIN={'blue':(-405.,-150.),'redC':(180.,320.)}
PRODS=[('std_Xb4d','HCN'),('std_Xb4d','HCO+'),('tap_Xb4d','HCN'),('tap_Xb4d','HCO+'),('std_Xb4b','HCN'),('std_Xb4b','HCO+'),('man_Xb4b','HCO+'),('manflat_Xb4b','HCO+')]
res={}
for prod,mol in PRODS:
    R=load_any(prod,mol,1.05,12.0 if 'Xb4b' not in prod else 9.0); v=R['v']; dv=R['dv']
    P(f'## {prod} {mol}: nchan {len(v)}, v {v.min():.0f}..{v.max():.0f}, dv {dv:.1f}, ps {R["ps"]:.3f}", n {R["n"]}')
    fm=field_pixels(R); ys,xs=np.where(fm); Y,X=np.mgrid[:R['n'],:R['n']]
    def apspec(pos,rad,c=R['c']): return aspec(c,R['ba'],apmask(R['n'],R['ps'],R['pxy'](*pos),rad))
    for rad in (1.0,1.5):
        # field apertures
        fs=[]
        while len(fs)<60:
            i=rng.integers(len(ys)); m=np.hypot(X-xs[i],Y-ys[i])*R['ps']<=rad
            if fm[m].all(): fs.append(aspec(R['c'],R['ba'],m))
        fs=np.array(fs); base_field=fs.mean(0)
        for bl in ('A','B','C'):
            if bl=='C' and not prod.startswith('std') : continue
            if bl=='C': cc=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); sub=lambda pos,rr: apspec(pos,rr,cc); fsB=np.array([aspec(cc,R['ba'],np.hypot(X-xs[i],Y-ys[i])*R['ps']<=rad) for i in rng.integers(len(ys),size=60) if fm[np.hypot(X-xs[i],Y-ys[i])*R['ps']<=rad].all()])
            elif bl=='B': sub=lambda pos,rr: apspec(pos,rr)-base_field; fsB=fs-base_field
            else: sub=lambda pos,rr: apspec(pos,rr); fsB=fs
            for wn,(a,b) in WIN.items():
                for pn,pos in (('lobe',LOBE),('red',RED)):
                    if (wn=='blue')!=(pn=='lobe') and not (wn=='redC' and pn=='lobe'): continue
                    sp=sub(pos,rad); F=wflux(v,dv,sp,a,b); sg=np.std([wflux(v,dv,f,a,b) for f in fsB])
                    res[f'{prod}|{mol}|r{rad}|{bl}|{wn}|{pn}']=dict(F=F,sig=float(sg))
                    P(f'{prod:13s} {mol:5s} r={rad} base {bl} {wn:5s}({a:.0f}..{b:.0f}) at {pn:4s}: {F:+.3f} +- {sg:.3f} Jy km/s ({F/sg:+.1f} sigma, {len(fsB)} field apertures)')
json.dump(res,open(OUT+'/s77a.json','w'),indent=1); open(OUT+'/s77a.txt','w').write('\n'.join(txt))
