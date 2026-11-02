# Q40 part 2: nuclear-wing decontamination at the native Xb55 beam (8.64x6.00) for 2-1 and the native Xb4f beam (5.72x3.54) for 1-0 and 3-2; compact components are beam independent. Sets: 1-0 (Harada, Michiyama), 2-1 (Xb55), 3-2 (Xb4f, Xb4d).
# Model per map (pixels r<RMAX from field centre): non-negative beam-shaped components at N, S, [Bridge], E1, E2, W (+ optional broad halo) + free-sign constant.
# Noise realisations: the line-free (null) windows smoothed identically are ADDED to the data and refit (correlated noise). Variants: bridge/halo/RMAX. R_K as q40_clump_nucleus.py.
# usage: python3 -I q40_decontam.py -> out/q40/q40_decontam.json/.log
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import q25_ratio10 as Q
from scipy.optimize import nnls
import json, itertools
BEAMS={'Xb55':(8.6426,5.9954,-89.72)}; B32=(5.7222,3.5423,82.02)
TB=None   # set per map set (native Xb55 beam for 2-1; native Xb4f beam for 1-0 and 3-2)
NU={1:(88.631602,89.188526),2:(177.261111,178.375056),3:(265.886434,267.557626)}
KF={J:(b/a)**2 for J,(a,b) in NU.items()}
POS={'N':(156.96346,-43.90390),'S':(156.96340,-43.90532),'Bridge':(156.96357,-43.90457),'E1':(156.96606,-43.90353),'E2':(156.96548,-43.90540),'W':(156.96114,-43.90590)}
log=[]
def P(s=''): print(s,flush=True); log.append(s)
def smooth_tb(d,arr):
    if abs(d['beam'][0]/TB[0]-1)<0.02 and abs(d['beam'][1]/TB[1]-1)<0.02: return np.nan_to_num(arr)   # native beam = target (Xb55)
    k=kernel(d['hd'],TB,d['beam']); a=np.nan_to_num(arr); return fftconvolve(a,k,mode='same')*OM(TB)/OM(d['beam'])
def getmaps(name,ln):
    if name in ('Harada','Michiyama'):
        d=Q.load10(Q.P[name],list(Q.REST.values())[ln]); d['x0']=0;d['y0']=0; nulls=d['nulls']
    else:
        key={'Xb55':['HCN 2-1','HCO+ 2-1'],'Xb4f':['HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'],'Xb4d':['HCN 3-2@Xb4d','HCO+ 3-2@Xb4d']}[name][ln]
        d=load(key); nulls=[n[1] for n in d['nulls']]
    return d,smooth_tb(d,d['m0']),[smooth_tb(d,n) for n in nulls]
def comps(d,names,halo,rmax):
    dx,dy=d['dx'],d['dy']; r=np.hypot(dx,dy); sel=(r<rmax)
    sg=[TB[0]/2.3548,TB[1]/2.3548]   # beam-shaped (unresolved), PA ~ -90 => major axis E-W
    cols=[]
    for n in names:
        ra,dec=POS[n]; ex=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; ey=(dec-DE0)*3600
        cols.append(np.exp(-0.5*(((dx-ex)/sg[0])**2+((dy-ey)/sg[1])**2))[sel])
    if halo: cols.append(np.exp(-0.5*((dx**2+dy**2)/(halo/2.3548)**2))[sel])
    return sel,cols
def fit(y,cols,sel,sgn=True):
    A=np.array(cols).T; 
    A2=np.hstack([A,np.ones((A.shape[0],1)),-np.ones((A.shape[0],1))]) if sgn else A
    x,_=nnls(A2,y[sel]); amp=x[:len(cols)]; 
    return amp,(x[len(cols)]-x[len(cols)+1] if sgn else 0.)
VAR=[(('N','S','E1','E2','W'),0,14),(('N','S','Bridge','E1','E2','W'),0,14),(('N','S','Bridge','E1','E2','W'),0,17),(('N','S','E1','E2','W'),20,14),(('N','S','Bridge','E1','E2','W'),20,17),(('N','S','Bridge','E1','E2','W'),0,11)]
SETS=[(1,'Harada'),(1,'Michiyama'),(2,'Xb55'),(3,'Xb4f'),(3,'Xb4d')]
OUT={}
for J,name in SETS:
    TB=BEAMS.get(name,B32)
    for ln,lab in ((0,'HCN'),(1,'HCO+')):
        d,s,ns=getmaps(name,ln)
        for vi,(names,halo,rmax) in enumerate(VAR):
            sel,cols=comps(d,names,halo,rmax)
            amp,c=fit(s,cols,sel)
            noise=[]
            for n in ns:
                ok=np.isfinite(n)&(n!=0)
                if ok[sel].mean()<0.9: continue
                a2,_=fit(s+n,cols,sel); noise.append(a2-amp)
            e=np.std(noise,axis=0) if len(noise)>=3 else np.full(len(cols),np.nan)
            OUT[(name,lab,vi)]=dict(names=list(names),amp=amp.tolist(),c=c,err=e.tolist(),nn=len(noise))
        P(f'{name} {lab} J={J}: variant0 amps '+' '.join(f'{n}={a:.2f}' for n,a in zip(VAR[0][0],OUT[(name,lab,0)]['amp']))+f' const {OUT[(name,lab,0)]["c"]:.2f}; err '+' '.join(f'{e:.2f}' for e in OUT[(name,lab,0)]['err'])+f' (n_noise {OUT[(name,lab,0)]["nn"]})')
# ratios
def val(name,lab,reg,vi):
    o=OUT[(name,lab,vi)]; 
    if reg=='N+S': return sum(o['amp'][o['names'].index(k)] for k in ('N','S')),np.hypot(*[o['err'][o['names'].index(k)] for k in ('N','S')])
    i=o['names'].index(reg); return o['amp'][i],o['err'][i]
res={}
P('\n## Decontaminated component ratios R_K (HCN/HCO+) per variant (vi=0..5), value +- noise; spread over variants = model systematic')
def R(name,J,reg,vi):
    h=val(name,'HCN',reg,vi); o=val(name,'HCO+',reg,vi)
    if h[0]<=0 or o[0]<=0: return np.nan,np.nan,h,o
    r=h[0]/o[0]*KF[J]; return r,r*np.hypot(h[1]/h[0],o[1]/o[0]),h,o
for reg in ('N','S','N+S','E1','E2','W'):
    for J,name in SETS:
        rs=[R(name,J,reg,vi) for vi in range(len(VAR))]
        # variants where this component exists
        rr=[x for vi,x in enumerate(rs) if (reg in VAR[vi][0] or reg=='N+S')]
        v0=rs[0]; spread=np.nanstd([x[0] for x in rr]); res[(reg,name)]=(v0[0],v0[1],spread,[x[0] for x in rr])
        h,o=v0[2],v0[3]
        P(f' {reg:4s} {name:9s} HCN {h[0]:.2f}+-{h[1]:.2f} HCO+ {o[0]:.2f}+-{o[1]:.2f} Jy km/s  R={v0[0]:.3f}+-{v0[1]:.3f} (var spread {spread:.3f}; range {np.nanmin([x[0] for x in rr]):.3f}-{np.nanmax([x[0] for x in rr]):.3f})')
P('\n## Ladders / double ratios (error = noise + variant spread + 10% cal per ratio of two bands)')
def dr(a,b,cal=True):
    q=a[0]/b[0]; ea=np.hypot(a[1],a[2])/a[0]; eb=np.hypot(b[1],b[2])/b[0]; return q,q*np.hypot(np.hypot(ea,eb),0.10 if cal else 0)
dd={}
for reg in ('N','S','N+S','E1','E2','W'):
    for lab,a,b in (('R32f/R21','Xb4f','Xb55'),('R32d/R21','Xb4d','Xb55'),('R21/R10M','Xb55','Michiyama'),('R21/R10H','Xb55','Harada'),('R10M/R10H','Michiyama','Harada')):
        dd[(reg,lab)]=dr(res[(reg,a)],res[(reg,b)]); 
    P(f' {reg:4s} R10M {res[(reg,"Michiyama")][0]:.3f} R10H {res[(reg,"Harada")][0]:.3f} R21 {res[(reg,"Xb55")][0]:.3f} R32f {res[(reg,"Xb4f")][0]:.3f} R32d {res[(reg,"Xb4d")][0]:.3f} | '+' '.join(f'{l} {dd[(reg,l)][0]:.2f}+-{dd[(reg,l)][1]:.2f}' for l in ('R32f/R21','R32d/R21','R21/R10M','R21/R10H')))
P('\n## clump/nucleus ratio of R per J and of double ratios (nucleus = N+S)')
for c in ('E1','E2','W'):
    P(f' {c}/N+S R per J: '+'; '.join(f'{n} {res[(c,n)][0]/res[("N+S",n)][0]:.2f}' for n in ('Michiyama','Harada','Xb55','Xb4f','Xb4d')))
    for l in ('R32f/R21','R21/R10M','R21/R10H'):
        a=dd[(c,l)]; b=dd[('N+S',l)]; q=a[0]/b[0]; e=q*np.hypot(a[1]/a[0],b[1]/b[0]); P(f'   {c}/N+S {l}: {q:.2f}+-{e:.2f} ({(q-1)/e:+.1f} sigma)')
json.dump({'res':{f'{k[0]}|{k[1]}':v for k,v in res.items()},'dd':{f'{k[0]}|{k[1]}':v for k,v in dd.items()},'fit':{f'{k[0]}|{k[1]}|{k[2]}':v for k,v in OUT.items()}},open('/workspace/out/q40/q40_decontam.json','w'),indent=1,default=float)
open('/workspace/out/q40/q40_decontam.log','w').write('\n'.join(log))
