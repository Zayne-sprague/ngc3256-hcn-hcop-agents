# Q5h2: real-data sizes (range over fit radius, floor flagged, injection-calibrated 95% intervals) and line-minus-continuum offsets
# with a null-tested budget. Needs out/q5h/q5h_size_inj.csv from q5h_inject.py.
# python3 -I scripts/q5h2_real.py -> out/q5h/q5h2.log, q5h_size_real.csv, q5h_size_summary.csv, q5h_offsets.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import q5h_inject as Q
from q5h_inject import *
log=open(OH+'q5h2.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); log.write(s+'\n'); log.flush()
inj=pd.read_csv(OH+'q5h_size_inj.csv')
CL={'N core':(0.0,2.63),'S core':(-0.15,-2.60),'N-E':(1.24,2.13),'Br-E':(1.12,1.29),'S-E':(1.20,-2.52),'S-W2':(0.60,-2.75),'N-SW':(-0.55,1.9)}
RADS=(0.3,0.35,0.4,0.5,0.6)
def fitsrc(args):
    nm,cn,x,y,rad=args; img=IMG[nm]; pk=np.percentile(img[np.hypot(dx-x,dy-y)<rad],99)
    try: xf,yf,fa,fb,fg=fitone(img,x,y,rad,BC[nm],BM[nm],max(pk,1e-3))
    except Exception: return dict(src=cn,line=nm,rad=rad,ok=False)
    # peak of the fit model
    return dict(src=cn,line=nm,rad=rad,ok=True,x=xf,y=yf,fwhm_maj=fa,fwhm_min=fb,fwhm_geo=fg,peak_img=pk)
def invert(line,amp,rad,obs):
    """95% interval of the true circular FWHM consistent with the recovered geometric FWHM obs (injection quantiles, nearest amp / radius)"""
    t=inj[(inj.line==line)&(inj.rad==min(inj.rad.unique(),key=lambda r:abs(r-rad)))]
    a=min(t.amp.unique(),key=lambda a_:abs(a_-amp)); t=t[t.amp==a].sort_values('fwhm_true')
    ok=(t.g16<=obs)&(t.g84>=obs)
    return (t.fwhm_true[ok].min(),t.fwhm_true[ok].max()) if ok.any() else (np.nan,np.nan)
if __name__=='__main__':
    with Pool(16) as pool:
        jobs=[(nm,cn,x,y,r) for cn,(x,y) in CL.items() for nm in ('HCN','HCO+') for r in RADS]
        rows=pool.map(fitsrc,jobs)
    R=pd.DataFrame([r for r in rows if r['ok']]); R.to_csv(OH+'q5h_size_real.csv',index=False)
    pd.set_option('display.width',200); pr(R.round(3).to_string())
    # summary: stable = fits not at the floor (min axis >0.03") and geometric FWHM range over radius <= +-25% of median
    S=[]
    for (cn,nm),g in R.groupby(['src','line']):
        fg=g.fwhm_geo.values; fl=(g.fwhm_min<0.03).values; med=np.median(fg); rng=(fg.min(),fg.max())
        lo,hi=[],[]
        for _,r in g.iterrows():
            a,b=invert(nm,max(r.peak_img*0.9,0.3),r.rad,r.fwhm_geo); lo.append(a); hi.append(b)
        stable=(not fl.any()) and (rng[1]/rng[0]<1.5)
        S.append(dict(src=cn,line=nm,n_fit=len(g),n_floor=int(fl.sum()),geo_min=rng[0],geo_med=med,geo_max=rng[1],ratio_maxmin=rng[1]/rng[0],
                      inj_lo=np.nanmin(lo) if np.isfinite(lo).any() else np.nan,inj_hi=np.nanmax(hi) if np.isfinite(hi).any() else np.nan,stable=stable,
                      minor_min=g.fwhm_min.min(),minor_max=g.fwhm_min.max(),major_min=g.fwhm_maj.min(),major_max=g.fwhm_maj.max()))
    S=pd.DataFrame(S); S.to_csv(OH+'q5h_size_summary.csv',index=False); pr('\n# SIZE SUMMARY (geo = sqrt(maj*min) deconvolved FWHM, arcsec; 1"=194 pc)'); pr(S.round(3).to_string())
    # ---- offsets: same estimator on line m0 and cont map, radii 0.3-0.6, null-injection scatter
    ctr={'N core':(0.0,2.63),'S core':(-0.15,-2.60),'Bridge':(0.3,0.15)}
    OFF=[]
    def posfit(nm,x,y,rad):
        img=IMG[nm]; pk=max(np.percentile(img[np.hypot(dx-x,dy-y)<rad],99),1e-5)
        xf,yf,fa,fb,fg=fitone(img,x,y,rad,BC[nm],BM[nm],pk); return xf,yf,fg,pk
    for cn,(x,y) in ctr.items():
        # continuum peak position first to centre the windows
        for rad in (0.3,0.4,0.5,0.6):
            try: xc,yc,gc,pc=posfit('cont',x,y,rad)
            except Exception: continue
            for nm in ('HCN','HCO+'):
                try: xl,yl,gl,pl=posfit(nm,x,y,rad)
                except Exception: continue
                OFF.append(dict(src=cn,rad=rad,line=nm,d_east=xl-xc,d_north=yl-yc,sep=np.hypot(xl-xc,yl-yc),cont_peak_mJy=pc*1e3,cont_fwhm=gc,line_fwhm=gl,line_peak=pl))
    O_=pd.DataFrame(OFF); O_.to_csv(OH+'q5h_offsets.csv',index=False); pr('\n# line-minus-continuum offsets (arcsec; east, north), same estimator, radii 0.3-0.6'); pr(O_.round(3).to_string())
    # null scatter: cont and lines
    def nullsc(line,amp,rad,n=200):
        with Pool(16) as pool: res=np.array(pool.map(Q.inject,[(line,0.25,amp,rad,90000+i) for i in range(n)]))
        ok=np.isfinite(res[:,0]); return res[ok,0].std(),res[ok,1].std(),np.median(np.hypot(res[ok,0],res[ok,1])),np.percentile(np.hypot(res[ok,0],res[ok,1]),95)
    pr('\n# null (zero-offset injection) centroid scatter, intrinsic FWHM 0.25" (sd east, sd north, median r, 95% r)')
    NS=[]
    for cn,camp in (('N core',0.0006),('S core',0.0019),('Bridge',0.0006)):
        for rad in (0.35,0.5):
            sc=nullsc('cont',camp,rad); pr('cont',cn,camp,rad,np.round(sc,4)); NS.append(dict(kind='cont',src=cn,amp=camp,rad=rad,sd_e=sc[0],sd_n=sc[1],r50=sc[2],r95=sc[3]))
    for nm,amps in (('HCN',(0.85,0.32,0.30)),('HCO+',(0.8,0.9,0.7))):
        for cn,amp in zip(('N core','S core','Bridge'),amps):
            for rad in (0.35,0.5):
                sc=nullsc(nm,amp,rad); pr(nm,cn,amp,rad,np.round(sc,4)); NS.append(dict(kind=nm,src=cn,amp=amp,rad=rad,sd_e=sc[0],sd_n=sc[1],r50=sc[2],r95=sc[3]))
    pd.DataFrame(NS).to_csv(OH+'q5h_offset_null.csv',index=False)
