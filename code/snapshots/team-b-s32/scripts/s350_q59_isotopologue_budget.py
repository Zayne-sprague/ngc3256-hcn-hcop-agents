"""Q59: 13C/12C line ratios H13CN/HCN and H13CO+/HCO+ (J=1-0) at N and S in ONE common circular beam (2.2" and 2.5"),
2015.1.00993.S (X961: HCN, HCO+ 1-0; X95d: H13CN, H13CO+, SiO 2-1, HN13C, CCH; same night, 1 h apart; 7.8 MHz = 27 km/s channels).
Systematic budget in the style of Q56 (s344/s345): every combination of
  beam {2.2, 2.5"} x aperture {peak pixel, 0.5 theta, theta, S-centred (aperture centred on the HCN moment peak instead of HCO+)}
  x window {+-100, +-150, +-250 km/s about 2780 km/s (box) , template (HCO+ aperture spectrum, |v|<500)} x baseline {linear, constant}
and the weak line amplitude from a JOINT fit of all lines in the spectral window (H13CN + neighbours; H13CO+ + SiO + HN13C + CCH) with local baseline.
Errors: formal (line-free rms) x null-aperture calibration (positions >=7" from the nuclei, same window), floor = null std.
Ratios use the same aperture/window for the 12C line (X961 spw0). Cross-EB calibration systematic 5% added (continuum check: s351).
Outputs: out/q59/s350.json, s350.txt.  usage: python3 -I s350_q59_isotopologue_budget.py"""
import sys, json, os, itertools, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from s182_q29_measure import F0, off
from regions import REGIONS
D='/archive/files/2015.1.00993.S/'
CUB={'X961_spw0':D+'member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',
     'X95d_spw0':D+'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw0_85203MHz.12m.cube.I.pbcor.fits',
     'X95d_spw1':D+'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw1_86692MHz.12m.cube.I.pbcor.fits'}
F0.update(SiO=86.8470,HN13C=87.0908,H13COp=86.7543,H13CN=86.3402,CCHa=87.3169,CCHb=87.4020)
CL={'X95d_spw0':['H13CN','H13COp'],'X95d_spw1':['H13COp','SiO','HN13C','CCHa','CCHb']}
EXTRA=[F0[k] for k in ('SiO','HN13C','H13COp','H13CN','CCHa','CCHb')]+[89.487]
O='/workspace/out/q59/'; os.makedirs(O,exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
def vel(nu,f0): return C*(1-nu/f0)-V0
REG={'N':REGIONS['N'][:2],'S':REGIONS['S'][:2]}
PIX=0.29
def prep(key,theta):
    cut,nu,pix,bm,pxy=load(CUB[key],half=22)
    return to_tb_conv(cut,nu,pix,bm,theta,np.mean(nu)),nu,pix,pxy
def ap(Tb,pix,c,rad):
    ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    if rad<=0.5*pix: return Tb[:,int(round(c[1])),int(round(c[0]))]
    return Tb[:,np.hypot(X-c[0],Y-c[1])*pix<=rad].mean(axis=1)
def lsq(nu,sp,lines,cols,order):
    nuc=(nu-nu.mean())/(nu.max()-nu.min())
    vcl=np.array([vel(nu,F0[l]) for l in lines]); dmin=np.min(np.abs(vcl),axis=0)
    allv=np.array([vel(nu,f) for f in [F0[l] for l in lines]+EXTRA]); near=np.min(np.abs(allv),axis=0)
    use=(dmin<1500)&np.isfinite(sp); free=use&(near>450)
    base=[np.ones_like(nu)]+([nuc] if order==1 else [])
    A=np.array(cols+base).T
    pf=np.polyfit(nuc[free],sp[free],order); sig=float(np.std(sp[free]-np.polyval(pf,nuc[free])))
    cf=np.linalg.lstsq(A[use],sp[use],rcond=None)[0]; cv=np.linalg.inv(A[use].T@A[use])*sig**2
    return cf,np.sqrt(np.diag(cv)),sig
def amp_weak(nu,sp,line,win,order,tmpl,vb):
    """joint fit of all lines in this cube's cluster; amplitude of 'line' as integral K km/s"""
    lines=[l for l in CL[kcube(nu)] if cov(nu,F0[l])]
    v0={l:vel(nu,F0[l]) for l in lines}; dv=abs(v0[lines[0]][1]-v0[lines[0]][0])
    cols=[]
    for l in lines:
        if win=='tmpl': cols.append(np.interp(v0[l],tmpl[0],tmpl[1],left=0,right=0))
        else: cols.append(((np.abs(v0[l])<=win)).astype(float))
    cf,er,sig=lsq(nu,sp,lines,cols,order); i=lines.index(line)
    if win=='tmpl': return cf[i],er[i],sig
    n=cols[i].sum(); return cf[i]*n*dv,er[i]*n*dv,sig
def kcube(nu): return 'X95d_spw0' if nu.min()<84.5 else 'X95d_spw1'
def cov(nu,f0,half=150): v=vel(nu,f0); return v.min()<-half and v.max()>half
def ref_int(nu,sp,f0,win,order,tmpl=None):
    v=vel(nu,f0); o=np.argsort(v); v=v[o]; y=sp[o]; fr=(np.abs(v)>550)&(np.abs(v)<1500)
    y=y-np.polyval(np.polyfit(v[fr],y[fr],order),v); dv=abs(v[1]-v[0])
    if win=='tmpl': return float(y[np.abs(v)<=500].sum()*dv)
    return float(y[np.abs(v)<=win].sum()*dv)
def template(nu,sp):
    v=vel(nu,F0['HCOp']); o=np.argsort(v); v=v[o]; fr=(np.abs(v)>550)&(np.abs(v)<1500)
    t=(sp[o]-np.polyval(np.polyfit(v[fr],sp[o][fr],1),v)); t=np.where(np.abs(v)<=500,np.clip(t,0,None),0); return v,t/np.trapezoid(t,v)
def run(theta):
    DT={k:prep(k,theta) for k in CUB}; T0,nu0,pix,pxy0=DT['X961_spw0']; ny,nx=T0.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    vO=vel(nu0,F0['HCOp']); vN=vel(nu0,F0['HCN']); dv0=abs(vO[1]-vO[0])
    def pos_for(rn,centre):  # sky position of aperture centre
        px,py=pxy0(REG[rn]); line=F0['HCOp'] if centre=='hco' else F0['HCN']
        v=vel(nu0,line); m0=np.tensordot((np.abs(v)<=150).astype(float),T0,axes=(0,0))*dv0
        mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm)
        dx=(float(X.flat[j])-px)*pix; dy=(float(Y.flat[j])-py)*pix
        return (REG[rn][0]-dx/3600/np.cos(np.deg2rad(REG[rn][1])),REG[rn][1]+dy/3600)
    def spectra(pos,rad):
        return {k:ap(DT[k][0],DT[k][2],DT[k][3](pos),rad) for k in CUB}
    # null positions: ring around the Bridge, 12C integral small
    pc=REGIONS['Bridge'][:2]; nul=[]
    for rr in (7.,9.,11.,13.):
        for ang in range(0,360,30):
            p=off(pc,rr*np.sin(np.deg2rad(ang)),rr*np.cos(np.deg2rad(ang))); s=ap(T0,pix,pxy0(p),theta)
            if abs(np.sum(s[np.abs(vO)<350])*dv0)<12: nul.append(p)
    res=[]
    for rn in ('N','S'):
        for centre in ('hco','hcn'):
            pos=pos_for(rn,centre)
            for rname,rad in (('pix',0.0),('half',0.5*theta),('full',theta)):
                sp=spectra(pos,rad); sn=[spectra(p,rad) for p in nul]
                for win,order in itertools.product((100,150,250,'tmpl'),(1,0)):
                    tm=template(nu0,sp['X961_spw0']) if win=='tmpl' else None
                    ref={'HCN':ref_int(nu0,sp['X961_spw0'],F0['HCN'],win,order),'HCOp':ref_int(nu0,sp['X961_spw0'],F0['HCOp'],win,order)}
                    for line,cubes in (('H13CN',['X95d_spw0']),('H13COp',['X95d_spw0','X95d_spw1'])):
                        est=[]
                        for k in cubes:
                            nu=DT[k][1]
                            if not cov(nu,F0[line]): continue
                            a,e,sg=amp_weak(nu,sp[k],line,win,order,tm,None)
                            na=np.array([amp_weak(nu,s_[k],line,win,order,tm,None)[0] for s_ in sn])
                            est.append((a,max(e*max(1,np.sqrt(np.mean((na/np.array([amp_weak(nu,s_[k],line,win,order,tm,None)[1] for s_ in sn]))**2))),na.std()),na.mean()))
                        # H13CO+ combined from both spws (inverse variance) - each already a full-line integral
                        w=np.array([1/x[1]**2 for x in est]); a=float(np.sum([x[0]*wi for x,wi in zip(est,w)])/w.sum()); e=float(1/np.sqrt(w.sum()))
                        r=ref['HCN' if line=='H13CN' else 'HCOp']
                        res.append(dict(theta=theta,reg=rn,centre=centre,ap=rname,win=win,order=order,line=line,amp=a,err=e,ref=r,ratio=a/r,ratio_err=float(np.hypot(e/r,0.05*abs(a/r))),null_mean=float(np.mean([x[2] for x in est]))))
    return res
if __name__=='__main__':
    R=[]
    for th in (2.2,2.5): R+=run(th)
    json.dump(R,open(O+'s350.json','w'),indent=0)
    P('H13CN/HCN and H13CO+/HCO+ 1-0 budget: %d estimates'%len(R))
    def sel(**k): return [r for r in R if all(r[a]==b for a,b in k.items())]
    for line in ('H13CN','H13COp'):
        for rn in ('N','S'):
            x=np.array([r['ratio'] for r in sel(line=line,reg=rn)]); ee=np.array([r['ratio_err'] for r in sel(line=line,reg=rn)])
            P(f'{line:7s} {rn}: ratio mean {x.mean():.4f}  stat err (median) {np.median(ee):.4f}  spread std {x.std():.4f}  range {x.min():.4f}..{x.max():.4f}  n={len(x)}')
    P('\nper-estimator (theta, ap, window, order, centre) tables: standard config = 2.2", full, +-150, linear, HCO+-centre')
    for line in ('H13CN','H13COp'):
        for rn in ('N','S'):
            for r in sel(line=line,reg=rn,theta=2.2,ap='full',win=150,order=1,centre='hco')+sel(line=line,reg=rn,theta=2.5,ap='full',win=150,order=1,centre='hco')+sel(line=line,reg=rn,theta=2.2,ap='pix',win=150,order=1,centre='hco'):
                P(f'  {line:7s} {rn} th{r["theta"]} {r["ap"]:4s} amp {r["amp"]:6.2f}+-{r["err"]:.2f} K km/s  ref {r["ref"]:6.1f}  ratio {r["ratio"]:.4f}+-{r["ratio_err"]:.4f}')
    open(O+'s350.txt','w').write('\n'.join(log))
