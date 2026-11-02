# Q3 round 2, fix 2: Soutflow HCN/HCO+ 3-2 with the selection bias removed.  Fixed a-priori apertures (2" diameter on the 1.05" Xb4d cubes): centroid of the LINE-INTEGRATED emission (not the wing), Michiyama S-OF-Blue, and N, S, Bridge;
# velocity bins relative to each aperture's local centroid; Monte-Carlo ratio errors; injection test of the peak-selection bias with real (correlated) noise; independent-array check in Xb4b smoothed to 1.05".
# Outputs: out/q3b/ratio_*.csv, ratio.log
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
rng=np.random.default_rng(3256)
RAD=1.0; BIN=100.; NMC=20000
SYS=0.10    # position-dependent systematic floor on ln(ratio) per aperture (Xb4d vs Xb4b differ by 7-25% in core ratios; baseline/continuum/imaging); common calibration cancels in position differences
cT=(F0['HCN']/F0['HCO+'])**-2*(267.557626/265.886434)**0   # Tb(HCN)/Tb(HCO+) = (I_HCN/I_HCO+)*(nu_O/nu_H)^2 for the same beam (Rayleigh-Jeans, Jy/beam km/s -> K km/s)
cT=(F0['HCO+']/F0['HCN'])**2
log=open(OUT+'ratio.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
CUBES={}
def get(ous,sp):
    k=(ous,sp)
    if k not in CUBES: CUBES[k]=vcube(sp,ous,1.05,crop=9. if ous=='Xb4d' else 7.)
    return CUBES[k]
def winsum_map(c,lo,hi):
    m=(c['v']>=lo)&(c['v']<hi); return c['a'][m].sum(0)*c['dv']
# ---- a-priori aperture centres
def centroid_line_integrated(ous='Xb4d',box=(-2.5,2.5,-6.5,-3.6),win=(2325.,2900.)):
    """centroid (arcsec offsets) of the normalised HCN+HCO+ line-integrated map (whole local line, not the wing) of the S-extension region"""
    tot=0
    for sp in ('HCN','HCO+'):
        c=get(ous,sp); m=winsum_map(c,*win); m=np.where((c['dx']>box[0])&(c['dx']<box[1])&(c['dy']>box[2])&(c['dy']<box[3]),m,0); m=np.clip(m,0,None); tot=tot+m/m.sum()
    # use only pixels above 50% of the peak to avoid the S-nucleus tail dominating
    t=np.where(tot>0.5*tot.max(),tot,0); return float((t*c['dx']).sum()/t.sum()),float((t*c['dy']).sum()/t.sum())
def peak_wing(ous='Xb4d',box=(-2.5,2.5,-6.5,-3.6),win=(2425.,2575.)):
    c=get(ous,'HCN'); m=winsum_map(c,*win); m=np.where((c['dx']>box[0])&(c['dx']<box[1])&(c['dy']>box[2])&(c['dy']<box[3]),m,-1e9); j=np.unravel_index(np.argmax(m),m.shape); return float(c['dx'][j]),float(c['dy'][j])
APS={'N':POS['N'],'S':POS['S'],'Bridge':POS['Bridge'],'Soutflow (wing peak, selected)':peak_wing(),'Soutflow (line-integrated centroid)':centroid_line_integrated(),'M18 South-OF-Blue':POS['M18 South-OF-Blue'],'mirror of Soutflow N of S':POS['mirror of Soutflow N of S']}
pr('apertures (arcsec offsets from phase centre):',{k:tuple(round(float(a),2) for a in v) for k,v in APS.items()})
def spectra(ous,x,y,rad=RAD):
    out={}
    for sp in ('HCN','HCO+'):
        c=get(ous,sp); s=aper(c,x,y,rad)*1e3; s,e=baseline(c['v'],s,FREE[sp],1); out[sp]=(c['v'],s,e,c['dv'])
    return out
def local_centroid(sp_,vs0=2775.):
    """flux-weighted centroid of the 30 km/s-smoothed, peak-normalised mean of HCN and HCO+ above 25% of peak within vs0 +-500"""
    acc=None
    for sp in ('HCN','HCO+'):
        v,s,e,dv=sp_[sp]; n=int(round(30/dv)); sm=np.convolve(s,np.ones(n)/n,'same'); sm=sm/np.nanmax(sm[(v>vs0-500)&(v<vs0+500)])
        g=np.interp(np.arange(vs0-500,vs0+501,5.),v[::-1],sm[::-1]); acc=g if acc is None else acc+g
    vv=np.arange(vs0-500,vs0+501,5.); acc/=2; m=acc>0.25
    return float((vv[m]*acc[m]).sum()/acc[m].sum())
def win_int(v,s,dv,lo,hi): k=(v>=lo)&(v<hi); return s[k].sum()*dv/1e3,int(k.sum())      # Jy/beam km/s
def win_noise(v,s,dv,free,n):
    m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    sm=np.convolve(s,np.ones(n),'valid'); mv=np.convolve(m.astype(float),np.ones(n),'valid')>n-0.5
    return robust_sigma(sm[mv])*dv/1e3
def ratio_mc(a,ea,b,eb,n=NMC):
    A=rng.normal(a,ea,n); B=rng.normal(b,eb,n); R=A/B*cT
    ok=B>0
    return R,ok
def ratio_stat(a,ea,b,eb):
    """value +- from MC (median, 16-84%) if both >3 sigma; else limit (lower: HCO+ at 3-sigma lower bound...)"""
    R,ok=ratio_mc(a,ea,b,eb)
    if a/ea>=3 and b/eb>=3:
        q=np.percentile(R[ok],[16,50,84]); return dict(R=q[1],lo=q[0],hi=q[2],kind='value',lnsig=0.5*(np.log(max(q[2],1e-3))-np.log(max(q[0],1e-3))))
    if a/ea>=3:  # HCO+ weak: lower limit on ratio (84% of the MC mass above)
        return dict(R=np.nan,lo=np.percentile(R[ok|True][(R>0)],16) if (R>0).any() else np.nan,hi=np.nan,kind='lower limit (HCO+<3sig)',lnsig=np.nan)
    return dict(R=np.nan,lo=np.nan,hi=np.nan,kind='HCN<3sig',lnsig=np.nan)
# ---- main table: ratio vs velocity relative to local centroid, same aperture for all
def table(ous):
    rows=[];loc={}
    for an,(x,y) in APS.items():
        sp_=spectra(ous,x,y); vl=local_centroid(sp_); loc[an]=vl
        for lo in np.arange(-500,500,BIN):
            r=dict(ous=ous,ap=an,vloc=vl,dv_lo=lo,dv_hi=lo+BIN)
            for sp in ('HCN','HCO+'):
                v,s,e,dv=sp_[sp]; I,n=win_int(v,s,dv,vl+lo,vl+lo+BIN); sg=win_noise(v,s,dv,FREE[sp],n); r['I_'+sp]=I; r['e_'+sp]=sg
            st=ratio_stat(r['I_HCN'],r['e_HCN'],r['I_HCO+'],r['e_HCO+']); r.update(st); rows.append(r)
    return pd.DataFrame(rows),loc
def windows(ous,spec_win=((-400,-100,'blue half'),(-100,100,'core'),(100,400,'red half'),(-300,-100,'blue -300..-100'),(-500,500,'whole line'))):
    rows=[]
    for an,(x,y) in APS.items():
        sp_=spectra(ous,x,y); vl=local_centroid(sp_)
        for lo,hi,name in spec_win:
            r=dict(ous=ous,ap=an,vloc=vl,window=name,lo=lo,hi=hi)
            for sp in ('HCN','HCO+'):
                v,s,e,dv=sp_[sp]; I,n=win_int(v,s,dv,vl+lo,vl+hi); r['I_'+sp]=I; r['e_'+sp]=win_noise(v,s,dv,FREE[sp],n)
            st=ratio_stat(r['I_HCN'],r['e_HCN'],r['I_HCO+'],r['e_HCO+']); r.update(st); rows.append(r)
    return pd.DataFrame(rows)
if __name__=='__main__':
    pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
    for ous in ('Xb4d','Xb4b'):
        T,loc=table(ous); pr(f'== velocity-resolved ratio ({ous} smoothed to 1.05", 2" diameter aperture, 100 km/s bins relative to the local centroid; Tb ratio; MC 68% interval)'); pr('local centroids (km/s):',{k:round(v) for k,v in loc.items()})
        T.round(3).to_csv(OUT+f'ratio_vs_v_{ous}.csv',index=False)
        for an in APS:
            t=T[T.ap==an]; pr(f'-- {an} (centroid {loc[an]:.0f})'); pr(t[['dv_lo','dv_hi','I_HCN','e_HCN','I_HCO+','e_HCO+','R','lo','hi','kind']].round(3).to_string(index=False))
        W=windows(ous); W.round(3).to_csv(OUT+f'ratio_windows_{ous}.csv',index=False); pr(f'== windows ({ous})'); pr(W[['ap','window','vloc','I_HCN','e_HCN','I_HCO+','e_HCO+','R','lo','hi','kind']].round(3).to_string(index=False))
        # significance of differences: ln ratio
        pr(f'== significance of ratio differences (ln R, MC sigma) {ous}')
        for wn in ('blue half','core','red half','blue -300..-100','whole line'):
            w=W[W.window==wn].set_index('ap')
            for ref in ('Soutflow (line-integrated centroid)','M18 South-OF-Blue','Soutflow (wing peak, selected)'):
                for oth in ('N','S','Bridge'):
                    a=w.loc[ref]; b=w.loc[oth]
                    if a.kind=='value' and b.kind=='value':
                        z=(np.log(a.R)-np.log(b.R))/np.hypot(np.hypot(a.lnsig,b.lnsig),np.sqrt(2)*SYS); pr(f'{wn:16s} {ref:38s} R={a.R:.2f} vs {oth:7s} R={b.R:.2f}: ln-difference {np.log(a.R/b.R):+.2f}, z={z:+.1f}')
                    else: pr(f'{wn:16s} {ref:38s} vs {oth}: {a.kind} / {b.kind}  (R_ref lo={a.lo:.2f})')

    # ---- selection-bias injection test (real correlated noise from line-free windows at the same pixels)
    def null_maps(sp,lo,hi,step=5):
        c=get('Xb4d',sp); n=int(((c['v']>=lo)&(c['v']<hi)).sum()); v=c['v']; cs=np.cumsum(np.concatenate([np.zeros((1,)+c['a'].shape[1:]),c['a']]),0); out=[]
        for i in range(0,len(v)-n,step):
            if any(v[i]>=a and v[i+n-1]<=b for a,b in FREE[sp]): out.append((cs[i+n]-cs[i])*c['dv'])
        return np.array(out)
    from scipy.ndimage import gaussian_filter
    WIN=(2425.,2575.); box=(-2.5,2.5,-6.5,-3.6)
    cH=get('Xb4d','HCN'); NH=null_maps('HCN',*WIN); NO=null_maps('HCO+',*WIN); pr(f'== injection test: {len(NH)} HCN and {len(NO)} HCO+ line-free window maps (150 km/s)')
    bm=(cH['dx']>box[0])&(cH['dx']<box[1])&(cH['dy']>box[2])&(cH['dy']<box[3])
    SH=gaussian_filter(winsum_map(cH,*WIN),0.5/cH['pix']); SH=np.where(bm,np.clip(SH,0,None),0)
    pix=cH['pix']; yy,xx=np.mgrid[:SH.shape[0],:SH.shape[1]]
    def ap_mean(M,j): return M[np.hypot(yy-j[0],xx-j[1])<=RAD/pix].mean()
    rows=[]
    for Rt in (0.6,0.9,1.2,1.7,2.5):
        meas=[];meas_ap=[];dist=[]; SO=SH*cT/Rt
        for it in range(3000):
            nh=NH[rng.integers(len(NH))]*rng.choice([-1,1]); no=NO[rng.integers(len(NO))]*rng.choice([-1,1])
            MH=np.where(bm,SH+nh,-1e9); j=np.unravel_index(np.argmax(MH),MH.shape); MO=SO+no
            meas.append(MH[j]/MO[j]*cT/cT if MO[j]>0 else np.nan); 
            a=ap_mean(SH+nh,j); b=ap_mean(MO,j); meas_ap.append(a/b*cT/cT*(1) if b>0 else np.nan)
        # ratio of Tb: I_H/I_O * cT_true: SO=SH*cT/Rt -> I_H/I_O = Rt/cT ; measured Tb ratio = (MH/MO)*cT
        meas=np.array(meas)*cT; meas_ap=np.array(meas_ap)*cT
        # true ratio at the signal peak is Rt by construction
        ok=np.isfinite(meas); oka=np.isfinite(meas_ap)
        q=np.percentile(meas[ok],[16,50,84]); qa=np.percentile(meas_ap[oka],[16,50,84])
        rows.append(dict(R_true=Rt,pix_median=q[1],pix_16=q[0],pix_84=q[2],ap2_median=qa[1],ap2_16=qa[0],ap2_84=qa[2],frac_HCO_nonpos=1-ok.mean(),peak_snr_med=np.median([SH.max()/np.std(NH[:,bm])]) ))
    I=pd.DataFrame(rows); I.round(3).to_csv(OUT+'ratio_injection.csv',index=False); pr(I.round(3).to_string(index=False))
    pr('injection reading: bias = median(measured at the selected peak pixel)/R_true; 2" aperture = same but averaged in a 2" diameter aperture about the selected peak')
