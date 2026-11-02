# Q40: HCN/HCO+ ladder (1-0,2-1,3-2) at 10" at nuclear complex (N, S, midpoint) and clumps E1, E2, W (team-b positions), with zero-level and PB-width tests.
# Brightness-temperature scale R_K = (S_HCN/S_HCO+)*(nu_HCO+/nu_HCN)^2 (as team-b), maps smoothed to a circular 10" beam, value read at position (bilinear).
# usage: python3 -I q40_clump_nucleus.py  -> out/q40/q40.json, q40.log
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import q25_ratio10 as Q
import json
from scipy.ndimage import map_coordinates
TB=10.
POS={'N':(156.96346,-43.90390),'S':(156.96340,-43.90532),'mid':(156.96357,-43.90457),'E1':(156.96606,-43.90353),'E2':(156.96548,-43.90540),'W':(156.96114,-43.90590)}
NU={1:(88.631602,89.188526),2:(177.261111,178.375056),3:(265.886434,267.557626)}
KF={J:(b/a)**2 for J,(a,b) in NU.items()}
FARR={'Xb55':(30,45),'Xb4f':(26,40),'Xb4d':(17,19),'Harada':(30,45),'Michiyama':(30,45)}
log=[]
def P(s=''): print(s,flush=True); log.append(s)
def sample(hd,x0,y0,arr,ra,dec):
    px,py=[float(a) for a in WCS(hd).celestial.world_to_pixel_values(ra,dec)]
    return map_coordinates(np.nan_to_num(arr),[[py-y0],[px-x0]],order=1)[0]
def info(hd,x0,y0,dx,dy,s,nulls):
    r=np.hypot(dx,dy)
    out={}
    for k,(ra,dec) in POS.items(): out[k]=sample(hd,x0,y0,s,ra,dec)
    # noise: spatial MAD of smoothed null maps within r<12" (quietest quartile), plus null-map mean offset
    sg=[]
    for n in nulls:
        v=n[(r<12)&np.isfinite(n)&(n!=0)]
        if len(v)>20: sg.append(1.4826*np.median(abs(v-np.median(v))))
    sig=np.percentile(sg,25) if sg else np.nan
    far=(r>=FARR_[0])&(r<FARR_[1])&np.isfinite(s)&(s!=0); zl=s[far].mean() if far.sum()>50 else 0.
    return out,sig,zl
res={}
def run_J(J,name,loadfn):
    global FARR_
    FARR_=FARR[name]
    M={}
    for ln in (0,1):
        d,sm_,nul=loadfn(ln)
        hd=d['hd']; s=sm_(d); ns=[sm_(d,n) for n in d['nulls']] if d['nulls'] and not isinstance(d['nulls'][0],tuple) else [sm_(d,n[1]) for n in d['nulls']]
        o,sig,zl=info(hd,d.get('x0',0),d.get('y0',0),d['dx'],d['dy'],s,ns)
        M[ln]=(o,sig,zl)
    return M
def l_xb(key):
    def fn(ln):
        k=key[ln]; d=load(k); return d,(lambda dd,a=None: smooth(dd,TB,a)),None
    return fn
def l_b3(name):
    def fn(ln):
        f0=list(Q.REST.values())[ln]; d=Q.load10(Q.P[name],f0); d['x0']=0; d['y0']=0
        return d,(lambda dd,a=None: Q.sm(dd,a)),None
    return fn
SETS=[(1,'Harada',l_b3('Harada')),(1,'Michiyama',l_b3('Michiyama')),
      (2,'Xb55',l_xb(['HCN 2-1','HCO+ 2-1'])),(3,'Xb4f',l_xb(['HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'])),(3,'Xb4d',l_xb(['HCN 3-2@Xb4d','HCO+ 3-2@Xb4d']))]
VALS={}
for J,name,fn in SETS:
    M=run_J(J,name,fn); VALS[name]=(J,M)
    P(f'## {name} J={J}-{J-1}: sigma(HCN,HCO+) Jy/beam km/s at 10" = {M[0][1]:.3f},{M[1][1]:.3f}; far-field mean (HCN,HCO+) = {M[0][2]:.3f},{M[1][2]:.3f}')
    for p in POS: P(f'   {p:4s} HCN {M[0][0][p]:7.2f}  HCO+ {M[1][0][p]:7.2f}')
def ratio(name,p,f=0.,fixed_off=None,pbscale=1.):
    J,M=VALS[name]; h=M[0][0][p]-f*M[0][2]; o=M[1][0][p]-f*M[1][2]
    return h/o*KF[J] if h>0 and o>0 else np.nan
def err(name,p):
    J,M=VALS[name]; h=M[0][0][p]; o=M[1][0][p]; r=h/o*KF[J]
    return r*np.hypot(M[0][1]/h,M[1][1]/o)
R={}
P('\n## R_K = HCN/HCO+ (K scale) at 10", zero-level offset f x far-field mean (f=-1,0,+1) for each set; err = stat')
for name in VALS:
    for p in POS:
        R[(name,p)]=(ratio(name,p),err(name,p),[ratio(name,p,f) for f in (-1,0,1)])
        P(f'  {name:9s} {p:4s} R={R[(name,p)][0]:.3f}+-{R[(name,p)][1]:.3f}  (f=-1,0,+1: '+' '.join(f'{x:.3f}' for x in R[(name,p)][2])+')')
# PB-width test: ratio bias from PB(nu) mismatch with FWHM x0.7/1/1.4 (pbcor divides by PB(nu); if true PB width differs, HCN/HCO+ ratio changes by PB_H/PB_O)
P('\n## PB-width test (12 m, Gaussian PB, FWHM=1.13 lambda/D): fractional bias of R at x0.7 / x1 / x1.4 PB width')
PB={}
for J,name in ((1,'Harada'),(2,'Xb55'),(3,'Xb4f')):
    nuH,nuO=NU[J]
    for p in POS:
        ra,dec=POS[p]; rr=np.hypot((ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(dec-DE0)*3600)
        out=[]
        for sc in (0.7,1.,1.4):
            fwH=sc*1.13*299792458/(nuH*1e9)/12*206265; fwO=sc*1.13*299792458/(nuO*1e9)/12*206265
            # true PB scales; nominal pbcor assumed the x1 PB at each frequency. bias = PBnom_H/PBtrue_H / (PBnom_O/PBtrue_O)... expressed as ratio of corrections
            nH=1.13*299792458/(nuH*1e9)/12*206265; nO=1.13*299792458/(nuO*1e9)/12*206265
            g=lambda fw,rr: np.exp(-4*np.log(2)*(rr/fw)**2)
            b=(g(fwH,rr)/g(nH,rr))/(g(fwO,rr)/g(nO,rr))-1; out.append(b)
        PB[(name,p)]=out
        P(f'  {name:6s} {p:4s} r={rr:4.1f}": '+' / '.join(f'{100*x:+.1f}%' for x in out))
# ladders
def dr(a,b): # ratio of two R's with errors
    return a[0]/b[0],a[0]/b[0]*np.hypot(a[1]/a[0],b[1]/b[0])
CAL=0.05
P('\n## Ladders and double ratios (stat errors; cal = extra 5% per line per band added in quadrature to double ratios -> 10% per ratio of two bands)')
D={}
pairs=[('R10','Michiyama'),('R10h','Harada'),('R21','Xb55'),('R32f','Xb4f'),('R32d','Xb4d')]
for p in POS:
    P(f' {p}: '+'  '.join(f'{k}={R[(n,p)][0]:.3f}+-{R[(n,p)][1]:.3f}' for k,n in pairs))
    for lab,a,b in (('R32f/R21','Xb4f','Xb55'),('R32d/R21','Xb4d','Xb55'),('R21/R10M','Xb55','Michiyama'),('R21/R10H','Xb55','Harada'),('R10M/R10H','Michiyama','Harada'),('R32f/R32d','Xb4f','Xb4d')):
        v,e=dr(R[(a,p)],R[(b,p)]); ec=v*np.hypot(e/v,0.10); D[(lab,p)]=(v,e,ec)
        P(f'     {lab:10s} {v:.3f} +- {e:.3f} (stat) +- {ec:.3f} (stat+cal)')
P('\n## Clump / nuclear complex double ratios (value >1 = clump higher ratio / shallower decline); nucleus = N and mid')
for ref in ('N','mid'):
    for c in ('E1','E2','W'):
        for lab in ('R32f/R21','R21/R10M','R21/R10H'):
            a=D[(lab,c)]; b=D[(lab,ref)]; q=a[0]/b[0]; e=q*np.hypot(a[2]/a[0],b[2]/b[0]); P(f'  {c}/{ref} {lab:9s} {q:.2f}+-{e:.2f} ({(q-1)/e:+.1f} sigma from 1)')
# direct clump/nucleus R ratios per J
P('\n## R(clump)/R(nucleus) per J (stat only)')
for ref in ('N','mid'):
    for c in ('E1','W'):
        s=[]
        for n in ('Michiyama','Harada','Xb55','Xb4f','Xb4d'):
            q,e=dr(R[(n,c)],R[(n,ref)]); s.append(f'{n} {q:.2f}+-{e:.2f}')
        P(f'  {c}/{ref}: '+'; '.join(s))
json.dump({'R':{f'{k[0]}|{k[1]}':v for k,v in R.items()},'D':{f'{k[0]}|{k[1]}':v for k,v in D.items()},'PB':{f'{k[0]}|{k[1]}':v for k,v in PB.items()},
 'vals':{n:{'J':J,**{ln:{'v':M[ln][0],'sig':M[ln][1],'zl':M[ln][2]} for ln in (0,1)}} for n,(J,M) in VALS.items()}},open('/workspace/out/q40/q40.json','w'),indent=1,default=float)
open('/workspace/out/q40/q40.log','w').write('\n'.join(log))
