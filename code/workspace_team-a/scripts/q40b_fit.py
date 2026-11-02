# Q40b: common-beam / position-shift / free-size / injection tests of the clump-vs-nucleus HCN/HCO+ decomposition (reviewer fixes 3, 1, 2).
# Same NNLS decomposition as q40_decontam.py (beam-shaped or beam (+) Gaussian components N,S,E1,E2,W + free-sign constant, r<14"), but
#  (a) ALL J at ONE common target beam (8.64x6.00 PA-89.7, or circular 10"), (b) clump positions shifted by +-1", (c) clump (and nuclear) intrinsic Gaussian FWHM s free,
#  (d) injection/recovery: model(true size s) + real residual map -> refit assuming unresolved.
# usage: python3 -I q40b_fit.py -> out/q40b/q40b_fit.json, q40b_fit.log
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import q25_ratio10 as Q
from scipy.optimize import nnls
import json
NU={1:(88.631602,89.188526),2:(177.261111,178.375056),3:(265.886434,267.557626)}
KF={J:(b/a)**2 for J,(a,b) in NU.items()}
POS={'N':(156.96346,-43.90390),'S':(156.96340,-43.90532),'Bridge':(156.96357,-43.90457),'E1':(156.96606,-43.90353),'E2':(156.96548,-43.90540),'W':(156.96114,-43.90590)}
EXY={k:((ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(dec-DE0)*3600) for k,(ra,dec) in POS.items()}
SETS=[(1,'Harada'),(1,'Michiyama'),(2,'Xb55'),(3,'Xb4f'),(3,'Xb4d')]
TBS={'b86':(8.6426,5.9954,-89.72),'c10':(10.,10.,0.)}
log=[]
def P(s=''): print(s,flush=True); log.append(s)
def _gm(name,ln):
    if name in ('Harada','Michiyama'):
        d=Q.load10(Q.P[name],list(Q.REST.values())[ln]); d['x0']=0;d['y0']=0; nulls=d['nulls']
    else:
        key={'Xb55':['HCN 2-1','HCO+ 2-1'],'Xb4f':['HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'],'Xb4d':['HCN 3-2@Xb4d','HCO+ 3-2@Xb4d']}[name][ln]
        d=load(key); nulls=[n[1] for n in d['nulls']]
    return d,d['m0'],nulls
getmaps=_gm   # (rebindable by q40g)
def smooth_to(d,arr,tb):
    a=np.nan_to_num(arr)
    if abs(d['beam'][0]/tb[0]-1)<0.02 and abs(d['beam'][1]/tb[1]-1)<0.02: return a   # native beam = target
    k=kernel(d['hd'],tb,d['beam'])
    if k is None: raise ValueError(f'target beam smaller than native {d["beam"]}')
    return fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
def cov2(maj,mnr,pa,s=0.):
    th=np.deg2rad(pa); u=np.array([np.sin(th),np.cos(th)]); v=np.array([np.cos(th),-np.sin(th)])
    sm,sn=(maj**2+s**2)/(8*np.log(2)),(mnr**2+s**2)/(8*np.log(2))
    return sm*np.outer(u,u)+sn*np.outer(v,v)
CACHE={}
def maps(name,ln,tbk):
    key=(name,ln,tbk)
    if key not in CACHE:
        d,m,ns=getmaps(name,ln); tb=TBS[tbk]
        s=smooth_to(d,m,tb); nn=[smooth_to(d,n,tb) for n in ns]
        CACHE[key]=(d,s,nn)
    return CACHE[key]
GRID={}
def comp_cols(d,tb,names,sizes,shift,rmax=14.):
    dx,dy=d['dx'],d['dy']; sel=np.hypot(dx,dy)<rmax; cols=[]
    for n in names:
        ex,ey=EXY[n]; sx,sy=shift.get(n,(0,0)); s=sizes.get(n,0.)
        Si=np.linalg.inv(cov2(*tb,s=s)); X=np.stack([dx[sel]-ex-sx,dy[sel]-ey-sy],-1)
        cols.append(np.exp(-0.5*np.einsum('...i,ij,...j->...',X,Si,X)))
    return sel,cols
def fit(y,cols,sel):
    A=np.array(cols).T; A2=np.hstack([A,np.ones((A.shape[0],1)),-np.ones((A.shape[0],1))])
    x,rn=nnls(A2,y[sel]); n=len(cols); res=y[sel]-A2@x
    return x[:n],x[n]-x[n+1],np.sqrt(np.mean(res**2))
NAMES=('N','S','E1','E2','W')
def run(name,ln,tbk,shift,sizes,names=NAMES,noise=True):
    d,s,nn=maps(name,ln,tbk); tb=TBS[tbk]
    sel,cols=comp_cols(d,tb,names,sizes,shift)
    amp,c,rms=fit(s,cols,sel); e=np.full(len(cols),np.nan)
    if noise:
        dd=[]
        for n in nn:
            if n is None or not np.isfinite(n).all() or (n[sel]==0).mean()>0.1: continue
            a2,_,_=fit(s+n,cols,sel); dd.append(a2-amp)
        if len(dd)>=3: e=np.std(dd,axis=0)
    return dict(names=list(names),amp=amp.tolist(),c=float(c),rms=float(rms),err=e.tolist())
def allsets(tbk,shift={},sizes={},names=NAMES,noise=True,sets=SETS):
    out={}
    for J,name in sets:
        out[name]={lab:run(name,ln,tbk,shift,sizes,names,noise) for ln,lab in ((0,'HCN'),(1,'HCO+'))}
    return out
def R(o,J,reg):
    def v(lab):
        a=o[lab]; 
        if reg=='N+S': return sum(a['amp'][a['names'].index(k)] for k in 'NS'),np.hypot(*[a['err'][a['names'].index(k)] for k in 'NS'])
        i=a['names'].index(reg); return a['amp'][i],a['err'][i]
    h,q=v('HCN'),v('HCO+')
    if h[0]<=0 or q[0]<=0: return np.nan,np.nan
    r=h[0]/q[0]*KF[J]; return r,r*np.hypot(h[1]/h[0],q[1]/q[0])
def summ(tag,o,regs=('N','S','N+S','E1','E2','W')):
    P(f'## {tag}')
    for reg in regs:
        P(f' {reg:4s} '+'  '.join(f'{nm}: {R(o[nm],J,reg)[0]:.3f}+-{R(o[nm],J,reg)[1]:.3f}' for J,nm in SETS if nm in o))
    for c in ('E1','W'):
        P(f' {c}/(N+S): '+'  '.join(f'{nm}: {R(o[nm],J,c)[0]/R(o[nm],J,"N+S")[0]:.2f}' for J,nm in SETS if nm in o))
OUT={}
if __name__=='__main__':
    # (a) common beams
    for tbk in TBS:
        OUT[f'common|{tbk}']=allsets(tbk); summ(f'common beam {tbk}',OUT[f'common|{tbk}'])
    # (b) position shifts, one component at a time and all clumps together (b86 beam)
    for comp in ('E1','E2','W','clumps'):
        for sx,sy in ((1,0),(-1,0),(0,1),(0,-1)):
            cs=('E1','E2','W') if comp=='clumps' else (comp,)
            sh={c:(sx,sy) for c in cs}
            k=f'shift|{comp}|{sx},{sy}'; OUT[k]=allsets('b86',sh,noise=False); summ(k,OUT[k],('N+S','E1','W'))
    # free best position of E1 and W from the 3 bands (rms minimisation, 2-1 HCO+ and 3-2 HCO+ and 1-0 as sum of rms^2 normalised)
    for comp in ('E1','W'):
        best=None; tab=[]
        for sx in np.arange(-2,2.01,0.5):
            for sy in np.arange(-2,2.01,0.5):
                tot=0
                for J,nm in (SETS[1],SETS[2],SETS[3]):
                    for ln in (0,1):
                        r=run(nm,ln,'b86',{comp:(sx,sy)},{},noise=False)['rms']; tot+=r**2/ np.mean([run(nm,ln,'b86',{},{},noise=False)['rms']**2])
                tab.append((sx,sy,tot))
        tab=np.array(tab); i=tab[:,2].argmin(); P(f'best shift {comp}: ({tab[i,0]:+.1f},{tab[i,1]:+.1f}) arcsec  sum rel rms^2 {tab[i,2]:.2f} (zero shift: {tab[(tab[:,0]==0)&(tab[:,1]==0),2][0]:.2f}, 3 sets x 2 lines)')
        OUT[f'bestshift|{comp}']=tab.tolist()
    # (c) free size
    for sc in (0,2,4,6,8):
        for sn in (0,3):
            sz={'E1':sc,'E2':sc,'W':sc,'N':sn,'S':sn}
            k=f'size|clump{sc}|nuc{sn}'; OUT[k]=allsets('b86',{},sz); summ(k,OUT[k],('N','S','N+S','E1','W'))
            P('   rms(HCN,HCO+) per set: '+'  '.join(f'{nm}:{OUT[k][nm]["HCN"]["rms"]:.3f},{OUT[k][nm]["HCO+"]["rms"]:.3f}' for J,nm in SETS))
    # (d) injection & recovery: truth = s=0 data-fit amplitudes (so true ratios known) placed as clumps of size s_true (nuclei unresolved) + the REAL residual map of the s=0 fit (real wings, noise); refit assuming unresolved.
    P('## injection-recovery (b86): E1/(N+S), W/(N+S) true (input) -> recovered with the unresolved model')
    INJ={}
    base0=allsets('b86',{},{},noise=False)
    for st in (0,3,5,8):
        sz={'E1':st,'E2':st,'W':st}; rec={}
        for J,nm in SETS:
            for ln,lab in ((0,'HCN'),(1,'HCO+')):
                d,s,nn=maps(nm,ln,'b86'); tb=TBS['b86']
                sel0,cols0=comp_cols(d,tb,NAMES,{},{}); A0=np.array(cols0).T; b=base0[nm][lab]
                resid=s[sel0]-(A0@np.array(b['amp'])+b['c'])
                sel,cols=comp_cols(d,tb,NAMES,sz,{}); A=np.array(cols).T
                syn=np.zeros_like(s); syn[sel]=A@np.array(b['amp'])+b['c']+resid
                amp,c,rms=fit(syn,cols0,sel0)
                rec.setdefault(nm,{})[lab]=dict(names=list(NAMES),amp=amp.tolist(),err=[0]*5)
        INJ[st]=rec
        txt=[]
        for c in ('E1','W'):
            txt.append(f'{c}/(N+S) recov/true: '+' '.join(f'{nm}:{R(rec[nm],J,c)[0]/R(rec[nm],J,"N+S")[0]:.2f}/{R(base0[nm],J,c)[0]/R(base0[nm],J,"N+S")[0]:.2f}' for J,nm in SETS if nm in rec))
        P(f' s_true={st}: '+' | '.join(txt))
    OUT['inj']={str(k):v for k,v in INJ.items()}
    json.dump(OUT,open('/workspace/out/q40b/q40b_fit.json','w'),default=float)
    open('/workspace/out/q40b/q40b_fit.log','w').write('\n'.join(log))
