"""Q29 s185: weak-line measurement on the SHORT-SPACING-COMPLETE Band 3 cubes (2015.1.00993.S, 12m+ACA-like beams ~2.1x1.7") at ONE circular beam (2.2").
Cubes: X961 spw0 (HCN, HCO+ 1-0 = reference/template), X961 spw1 (HNC, HC3N), X95d spw0 (H13CN, edge of H13CO+), X95d spw1 (H13CO+, SiO 2-1, HN13C, CCH; requested, used if cached).
Method: pbcor Jy/beam -> RJ T_b (native beam) -> unit-area kernel to 2.2" circular -> aperture mean (radius 2.2", as s183/s186) centred on the HCO+ 1-0 +-150 km/s moment peak
within 0.4" of the nominal position (N,S,Bridge; Slobe/Ridge fixed offsets as s182).  Velocity window fixed from HCO+ (N,S,Bridge: HCO+ 1-0 aperture spectrum used as profile template; Slobe -320..-180, Ridge +180..+320 box).
Local baseline: linear fit to line-free channels (450<|v-vline|<1500 km/s, away from every other known line) fitted simultaneously with the templates of all lines in the cluster.
Noise: rms of the aperture-mean spectrum in line-free channels (per channel) -> sigma of the integral = sig*sqrt(sum t^2)... (formal); empirical: NULL apertures (same size, same window/template, 16 offset positions with HCO+ moment < 3% of N)
give the true scatter of the fitted amplitudes; the nominal error is multiplied by max(1,rms(null S/N)) (calibration factor).
Detection classes: >=5 sigma detection, 3-5 marginal, else 3-sigma limit.  Template amplitude = integral over the whole template (|v|<500 km/s) in K km/s.
usage: python3 -I s185_q29_weaklines.py -> out/q29/s185.json, s185.txt"""
import sys, json, os, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from s182_q29_measure import REG, SPECIAL, F0, vel, kern_w, off
D='/archive/files/2015.1.00993.S/'
CUB={'X961_spw0':D+'member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',
     'X961_spw1':D+'member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw1_89515MHz.12m.cube.I.pbcor.fits',
     'X95d_spw0':D+'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw0_85203MHz.12m.cube.I.pbcor.fits',
     'X95d_spw1':D+'member.uid___A001_X2fb_X95d.ari_l.ngc3256_sci.spw1_86692MHz.12m.cube.I.pbcor.fits'}
F0.update(SiO=86.8470,HN13C=87.0908,H13COp=86.7543,H13CN=86.3402,CCHa=87.3169,CCHb=87.4020,HNC=90.6636,HC3N=90.9790)
LINES={'X961_spw0':['HCN','HCOp'],'X961_spw1':['HNC','HC3N'],'X95d_spw0':['H13CN','H13COp'],'X95d_spw1':['H13COp','SiO','HN13C','CCHa','CCHb']}
EXTRA=[F0['SiO'],F0['HN13C'],F0['H13COp'],F0['H13CN'],89.487,F0['CCHa'],F0['CCHb']]   # also keep baseline away from these
TH=2.2; RAD=2.2; VT=500.
O='/workspace/out/q29/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def covered(nu,f0,half=300):
    v=vel(nu,f0); return v.min()<-half and v.max()>half
def prep(key):
    cut,nu,pix,bm,pxy=load(CUB[key],half=22)
    Tb=to_tb_conv(cut,nu,pix,bm,TH,np.mean(nu)); return Tb,nu,pix,bm,pxy
def ap_spec(Tb,pix,cen,rad):
    ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]; m=np.hypot(X-cen[0],Y-cen[1])*pix<=rad; return Tb[:,m].mean(axis=1)
def fit_cluster(nu,sp,lines,tmpl,box=None):
    """joint fit of templates (or boxes) for 'lines' with a linear local baseline; returns {line:(amp,err)}, rms of line-free channels"""
    vs={l:vel(nu,F0[l]) for l in lines}
    allv=[vel(nu,f) for f in list(map(lambda l:F0[l],lines))+EXTRA]
    vref=vs[lines[0]]
    # channels near the cluster: |v - nearest line| 450..1500 km/s and far (>450) from every line
    nearest=np.min(np.abs(np.array(allv)),axis=0)
    lo=min(vel(nu,F0[l]).min() for l in lines)
    # use frequency-independent velocity w.r.t. first line for ordering
    vcl=[vel(nu,F0[l]) for l in lines]; dmin=np.min(np.abs(np.array(vcl)),axis=0)
    free=(nearest>450)&(dmin<1500)&np.isfinite(sp)
    dv=abs(vref[1]-vref[0])
    cols=[]
    for l in lines:
        if box is None: cols.append(np.interp(vs[l],tmpl[0],tmpl[1],left=0,right=0))
        else: cols.append(((vs[l]>=box[0])&(vs[l]<=box[1])).astype(float)/1.0)
    nuc=(nu-nu.mean())/ (nu.max()-nu.min())
    use=(dmin<1500)&((nearest>450)|(np.max(np.abs(np.array(cols)),axis=0)>0)|(dmin<450))
    use=(dmin<1500)
    A=np.array(cols+[np.ones_like(nu),nuc]).T
    # amplitude definition: template normalised to unit integral  (cols from interp of normalised template, units 1/(km/s)); box: unit height -> integral = amp*width*..., handled by caller
    sig=float(np.std(sp[free]-np.polyval(np.polyfit(nuc[free],sp[free],1),nuc[free])))
    Au=A[use]; cf=np.linalg.lstsq(Au,sp[use],rcond=None)[0]; cv=np.linalg.inv(Au.T@Au)*sig**2
    return {l:(float(cf[i]),float(np.sqrt(cv[i,i]))) for i,l in enumerate(lines)},sig,int(free.sum()),A,cf,use
def classify(s): return 'DET' if s>=5 else ('marg' if s>=3 else 'lim')
res={}
CK=[k for k in CUB if os.path.exists(CUB[k])]
P('cubes used: '+', '.join(CK))
DATA={k:prep(k) for k in CK}
Tb0,nu0,pix0,bm0,pxy0=DATA['X961_spw0']
dv0=abs(vel(nu0,F0['HCOp'])[1]-vel(nu0,F0['HCOp'])[0]); vO=vel(nu0,F0['HCOp']); vH=vel(nu0,F0['HCN'])
P('X961 spw0 channel %.1f km/s, beam %.2fx%.2f, convolved to %.1f"'%(dv0,bm0[0],bm0[1],TH))
def positions(rn):
    """aperture centre in each cube's pixel frame (offset determined from the HCO+ +-150 moment peak in spw0)"""
    px,py=pxy0(REG[rn]); ny,nx=Tb0.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    if rn in SPECIAL: return REG[rn]
    k=kern_w(vO,dv0,-150,150); m0=np.tensordot(k,Tb0,axes=(0,0))*dv0
    mm=np.where(np.hypot(X-px,Y-py)*pix0<0.4,m0,-np.inf); j=np.argmax(mm); mp=(float(X.flat[j]),float(Y.flat[j]))
    # back to sky: offset arcsec in pixel axes (x=-RA)
    dxs=(mp[0]-px)*pix0; dys=(mp[1]-py)*pix0
    return (REG[rn][0]-dxs/3600/np.cos(np.deg2rad(REG[rn][1])),REG[rn][1]+dys/3600)
def spectra(pos):
    out={}
    for k,(Tb,nu,pix,bm,pxy) in DATA.items():
        x,y=pxy(pos); out[k]=ap_spec(Tb,pix,(x,y),RAD)
    return out
def template(sp0):
    s=np.interp(vO[::-1],vO[::-1],sp0[::-1])
    o=np.argsort(vO); v=vO[o]; t=sp0[o]; sel=np.abs(v)<=VT
    # linear baseline from line-free ends of HCO+ (|v|>450..700) removed first
    fr=(np.abs(vO)>550)&(np.abs(vO)<1500); c=np.polyfit(vO[fr],sp0[fr],1); t=(sp0-np.polyval(c,vO))[o]
    tt=np.where(sel,t,0.0); tt=np.where(tt<0,0,tt)*sel
    integ=np.trapezoid(tt,v); return (v,tt/integ),integ
def run_region(label,pos,mode,box=None,tmplsrc=None):
    sp=spectra(pos); out={}
    if tmplsrc is None: tm,integ_O=template(sp['X961_spw0'])
    else: tm,integ_O=template(tmplsrc['X961_spw0'])
    def integ(nuv,spec,f0,box):
        v=vel(nuv,f0); o=np.argsort(v); v=v[o]; y=spec[o]; fr=(np.abs(v)>550)&(np.abs(v)<1500); c=np.polyfit(v[fr],y[fr],1); y=y-np.polyval(c,v)
        sel=(np.abs(v)<=VT) if box is None else ((v>=box[0])&(v<=box[1])); dvv=abs(v[1]-v[0]); return float(y[sel].sum()*dvv),float(np.std((y-0)[fr])*np.sqrt(sel.sum())*dvv)
    res_=dict(HCOp_int=integ(nu0,sp['X961_spw0'],F0['HCOp'],box)[0],HCN_int=integ(nu0,sp['X961_spw0'],F0['HCN'],box))
    res_['HCOp_int']=integ(nu0,sp['X961_spw0'],F0['HCOp'],box)
    for k in CK:
        if k=='X961_spw0': continue
        nu=DATA[k][1]; ls=[l for l in LINES[k] if covered(nu,F0[l],150)]
        if not ls: continue
        # fix: weighted channel sum equals integral: template columns carry 1/(km/s) units -> amplitude is the integral (K km/s) after multiplying by dv?  cols*amp=K  => amp has units K*(km/s). OK.
        r,sig,nfree,A,cf,use=fit_cluster(nu,sp[k],ls,tm if box is None else None,box)
        for l,(a,e) in r.items():
            if box is not None:
                w=(vel(nu,F0[l])>=box[0])&(vel(nu,F0[l])<=box[1]); dvk=abs(vel(nu,F0[l])[1]-vel(nu,F0[l])[0]); a,e=a*w.sum()*dvk,e*w.sum()*dvk
            cov_=covered(nu,F0[l],300)
            out[(l,k)]=(a,e,sig,k,cov_)
    # HCN from spw0 using same-template box
    res_['lines']=out; return res_,sp
def go():
    # null positions: ring offsets from the midpoint of N,S,Bridge, HCO+ moment < 3% of N
    pc=REG['Bridge']; nul=[]
    kO=kern_w(vO,dv0,-350,350); m0=np.tensordot(kO,Tb0,axes=(0,0))*dv0; ny,nx=Tb0.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    pN=pxy0(REG['N']); ref=ap_spec(Tb0,pix0,pN,RAD); refN=np.trapezoid(np.clip(ref,0,None),vO[::-1]) if False else None
    cenpix=pxy0(pc); mN=m0[int(round(pN[1])),int(round(pN[0]))]
    cand=[]
    for rr in (7.0,9.0,11.0,13.0):
        for ang in np.arange(0,360,30):
            dE=rr*np.sin(np.deg2rad(ang)); dN=rr*np.cos(np.deg2rad(ang)); p=off(pc,dE,dN); x,y=pxy0(p)
            mmean=ap_spec(Tb0,pix0,(x,y),RAD); mi=np.sum(mmean[(np.abs(vO)<350)])*dv0
            cand.append((rr,ang,p,mi))
    P('null candidates HCO+ +-350 aperture integral (K km/s) vs N=%.0f:'%(np.sum(ap_spec(Tb0,pix0,pN,RAD)[np.abs(vO)<350])*dv0))
    nulls=[c for c in cand if c[3]<0.06*mN*1.0 or abs(c[3])<12]
    P('%d null apertures selected (HCO+ integral < 12 K km/s)'%len(nulls))
    return nulls
nulls=go()
regs=['N','S','Bridge','Slobe','Ridge']; tm_by={}
for rn in regs:
    pos=positions(rn); box=SPECIAL.get(rn)
    r,sp=run_region(rn,pos,'t' if box is None else 'b',box=box)
    res[rn]=r; tm_by[rn]=sp
# NULLS: same templates (N,S,Bridge) and boxes (Slobe,Ridge) applied at null positions
nullres={rn:{} for rn in regs}
for (rr,ang,p,mi) in nulls:
    spn=spectra(p)
    for rn in regs:
        box=SPECIAL.get(rn)
        tm,_=template(tm_by[rn]['X961_spw0'])
        for k in CK:
            if k=='X961_spw0': continue
            nu=DATA[k][1]; ls=[l for l in LINES[k] if covered(nu,F0[l],150)]
            if not ls: continue
            r_,sig,nf,A,cf,use=fit_cluster(nu,spn[k],ls,tm if box is None else None,box)
            for l,(a,e) in r_.items():
                if box is not None:
                    v=vel(nu,F0[l]); w=((v>=box[0])&(v<=box[1])).sum(); dvk=abs(v[1]-v[0]); a,e=a*w*dvk,e*w*dvk
                nullres[rn].setdefault((l,k),[]).append((a,e))
# combine H13CO+ from X95d spw0 (covers v>-205 km/s) and spw1 (covers v<+205 km/s): inverse-variance weighted, also in the nulls
def comb(L0,L1):
    a0,e0=np.array([x[0] for x in L0]),np.array([x[1] for x in L0]); a1,e1=np.array([x[0] for x in L1]),np.array([x[1] for x in L1])
    w0,w1=1/e0**2,1/e1**2; return list(zip((a0*w0+a1*w1)/(w0+w1),1/np.sqrt(w0+w1)))
for rn in regs:
    if ('H13COp','X95d_spw0') in nullres[rn] and ('H13COp','X95d_spw1') in nullres[rn]:
        nullres[rn][('H13COp','comb')]=comb(nullres[rn][('H13COp','X95d_spw0')],nullres[rn][('H13COp','X95d_spw1')])
        a0,e0,s0,k0,c0=res[rn]['lines'][('H13COp','X95d_spw0')]; a1,e1,s1,k1,c1=res[rn]['lines'][('H13COp','X95d_spw1')]
        w0,w1=1/e0**2,1/e1**2; res[rn]['lines'][('H13COp','comb')]=((a0*w0+a1*w1)/(w0+w1),1/np.sqrt(w0+w1),0.5*(s0+s1),'X95d_spw0+spw1',True)
P('\n== NULL apertures: scatter of fitted amplitude vs formal error (K km/s) ==')
calib={}
for rn in regs:
    for l,L in nullres[rn].items():
        a=np.array([x[0] for x in L]); e=np.array([x[1] for x in L]); z=a/e
        c=max(1.0,float(np.sqrt(np.mean(z**2)))); calib[(rn,l)]=(c,float(np.std(a)),float(np.mean(a)),len(a))
        P(f' {rn:7s} {str(l):28s} n={len(a)} mean {a.mean():+.2f} std {a.std():.2f} formal {e.mean():.2f} rms(S/N) {np.sqrt(np.mean(z**2)):.2f} max|S/N| {np.max(np.abs(z)):.1f}')
tbl={}
P('\n== WEAK LINES at 2.2" circular beam, aperture r=2.2", K km/s (RJ T_b); error = max(formal*null calibration, null std); DET>=5 sigma, marg 3-5, lim ==')
for rn in regs:
    for (l,k0),(a,e,sig,k,cov_) in res[rn]['lines'].items():
        c,sn,mn,nn=calib.get((rn,(l,k0)),(1.0,0,0,0))
        e2=max(e*c,sn); s=a/e2; cls=classify(s)
        tbl[f'{rn}|{l}|{k0}']=dict(amp=a,err_formal=e,err=e2,snr=s,cls=cls,lim3=3*e2,sig_ch=sig,cube=k,null_mean=mn,null_std=sn,calib=c)
        P(f' {rn:7s} {l:7s} {k0:10s} {a:7.2f} +- {e2:5.2f} ({s:5.1f} sig, {cls:4s}, 3sig lim {3*e2:5.1f}) [sig_ch {sig*1e3:.1f} mK, null mean {mn:+.2f} std {sn:.2f}]')
P('\n== reference integrals (X961 spw0, 2.2", |v|<500 km/s or Slobe/Ridge box): HCO+ and HCN K km/s ==')
for rn in regs: P(f' {rn}: HCO+ {res[rn]["HCOp_int"][0]:.1f}+-{res[rn]["HCOp_int"][1]:.1f}  HCN {res[rn]["HCN_int"][0]:.1f}+-{res[rn]["HCN_int"][1]:.1f}')
json.dump(dict(table=tbl,HCOp_int={rn:res[rn]['HCOp_int'] for rn in regs},HCN_int={rn:res[rn]['HCN_int'] for rn in regs},cube_list=CK,n_null=len(nulls)),open(O+'s185.json','w'),indent=1,default=float)
open(O+'s185.txt','w').write('\n'.join(log))

# ---------------- box moments (+-150 km/s as Harada+2018 Table 13; +-350 km/s) for N, S, Bridge, with the same null calibration ----------------
def boxfit(sp,box):
    out={}
    for k in CK:
        if k=='X961_spw0': continue
        nu=DATA[k][1]; ls=[l for l in LINES[k] if covered(nu,F0[l],150)]
        if not ls: continue
        r_,sig,nf,A,cf,use=fit_cluster(nu,sp[k],ls,None,box)
        for l,(a,e) in r_.items():
            v=vel(nu,F0[l]); w=((v>=box[0])&(v<=box[1])).sum(); dvk=abs(v[1]-v[0]); out[(l,k)]=(a*w*dvk,e*w*dvk)
    if ('H13COp','X95d_spw0') in out and ('H13COp','X95d_spw1') in out:
        (a0,e0),(a1,e1)=out[('H13COp','X95d_spw0')],out[('H13COp','X95d_spw1')]; w0,w1=1/e0**2,1/e1**2; out[('H13COp','comb')]=((a0*w0+a1*w1)/(w0+w1),1/np.sqrt(w0+w1))
    return out
def refint(sp0,f0,box):
    v=vel(nu0,f0); o=np.argsort(v); v=v[o]; y=sp0[o]; fr=(np.abs(v)>550)&(np.abs(v)<1500); y=y-np.polyval(np.polyfit(v[fr],y[fr],1),v)
    sel=(v>=box[0])&(v<=box[1]); dvv=abs(v[1]-v[0]); return float(y[sel].sum()*dvv)
BOXRES={}
P('\n== BOX MOMENTS (local linear baseline fitted jointly; errors = max(formal*calib, null std)); ratios to HCO+ / HCN box integral in the same window, plus 8% cross-project scale ==')
for bn,box,RAD in (('pm150',(-150,150),2.2),('pm150pk',(-150,150),0.3),('pm350',(-350,350),2.2)):
    nul={}
    for (rr,ang,p,mi) in nulls:
        for kk,(a,e) in boxfit(spectra(p),box).items(): nul.setdefault(kk,[]).append((a,e))
    for rn in ('N','S','Bridge'):
        sp=spectra(positions(rn)); R=boxfit(sp,box); IO=refint(sp['X961_spw0'],F0['HCOp'],box); IH=refint(sp['X961_spw0'],F0['HCN'],box)
        for kk,(a,e) in R.items():
            L=np.array(nul[kk]); c=max(1.0,float(np.sqrt(np.mean((L[:,0]/L[:,1])**2)))); e2=max(e*c,float(np.std(L[:,0]))); l=kk[0]
            ref=IH if l in ('H13CN','HNC','HC3N') else IO; refn='HCN' if l in ('H13CN','HNC','HC3N') else 'HCO+'
            if kk[1] in ('X95d_spw0','X95d_spw1') and l=='H13COp': continue
            rat=a/ref; er=np.hypot(e2/ref, 0.08*abs(rat))
            BOXRES[f'{bn}|{rn}|{l}']=dict(amp=a,err=e2,snr=a/e2,cls=classify(a/e2),ref=ref,refname=refn,ratio=rat,ratio_err=er,lim3_ratio=3*e2/ref+max(0,0),HCOp_int=IO,HCN_int=IH)
            P(f' {bn} {rn:6s} {l:7s} {a:6.2f}+-{e2:4.2f} ({a/e2:4.1f} sig {classify(a/e2)}) / {refn} {ref:6.1f} = {rat:.4f}+-{er:.4f}  (3sig lim {3*e2/ref:.4f})')
json.dump(BOXRES,open(O+'s185_box.json','w'),indent=1,default=float)
open(O+'s185.txt','w').write('\n'.join(log))
