# Q14e (fix 1,2): matched-aperture moment-0 CH3OH 254.015/HCO+ 3-2 with local baseline, windows matched to HCO+ (v0,FWHM) and Harada's +-150 km/s;
# empirical nulls: (a) same estimator in shifted line-free windows of the same spw at the same pixel, (b) same window at off-source positions (ring 9-16" from the N-S midpoint, HCO+ < 3 sigma).
# Output out/q14e/m0_table.csv, null_summary.csv, ring.csv.  python3 -I q14e_m0.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
F0={'HCO+':267.557626,'CH3OH':254.015377}
SPW={'Xb4d':{'HCO+':'spw25','CH3OH':'spw29'},'Xb4f':{'HCO+':'spw16','CH3OH':'spw20'}}
BEAMS=[('Xb4d','1.7',(1.7,1.7,0.)),('Xb4d','3.0',(3.0,3.0,0.)),('Xb4d','4.45x2.95',(4.45,2.95,87.)),('Xb4d','6.3x3.8',(6.3,3.8,82.)),('Xb4f','6.3x3.8',(6.3,3.8,82.))]
RG=['N','S','Bridge','Wclump','Soutflow']
rows=[];nullr=[];ringr=[]
rng=np.random.default_rng(1)
for ous,bn,tb in BEAMS:
    D={}
    for ln in F0:
        fr,a,cw,pix,beam,nb=cut(ous,SPW[ous][ln],tb,half=20.); D[ln]=(vel(fr,F0[ln]),a,cw)
    # ring positions (null b)
    ex=(CEN[0]-0,CEN[1])
    for reg in RG:
        ra,de=REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']
        vh,ah,cw=D['HCO+']; sh=sample(ah,cw,ra,de); p=hco_fit(vh,sh); v0,fw=p[1],p[2]
        vc,ac,cwc=D['CH3OH']; sc=sample(ac,cwc,ra,de)
        for wn,(c,W) in {'HCO+FWHM':(v0,fw),'Harada150':(2775.,300.)}.items():
            Ih,eh=m0(vh,sh,c,W,excl=(c-W/2,c+W/2)); Ic,ec=m0(vc,sc,c,W,excl=(c-W/2,c+W/2))
            # null (a): shifted windows in CH3OH cube, stride W/2, excluding |off|<W (+ real line wing) and any with <8 baseline chans
            vals=[]
            for off in np.arange(-1500,1501,W/2):
                if abs(off)<1.3*W: continue
                cc=c+off
                if cc-W/2<vc.min()+1.5*W/2 or cc+W/2>vc.max()-1.5*W/2: continue
                x,_=m0(vc,sc,cc,W,excl=(c-1.0*W,c+1.0*W))
                if np.isfinite(x): vals.append(x)
            vals=np.array(vals)
            rows.append(dict(ous=ous,beam=bn,region=reg,window=wn,v0=v0,fwhm=fw,W=W,I_HCO=Ih,I_CH3OH=Ic,e_stat=ec,null_n=len(vals),null_rms=np.std(vals) if len(vals)>2 else np.nan,null_mad=robust_sigma(vals) if len(vals)>2 else np.nan,null_mean=np.mean(vals) if len(vals)>2 else np.nan))
        # peak of CH3OH m0 map within 3.5" of region (HCO+FWHM window)
        c,W=v0,fw; vc,ac,cwc=D['CH3OH']; inw=np.abs(vc-c)<W/2; dv=abs(np.median(np.diff(vc)))
        mp=ac[inw].sum(0)*dv; ny,nx=mp.shape; yy,xx=np.mgrid[0:ny,0:nx]
        wc=cwc.pixel_to_world_values(xx,yy); dra=(wc[0]-ra)*np.cos(np.deg2rad(de))*3600; dde=(wc[1]-de)*3600; msk=np.hypot(dra,dde)<3.5
        j=np.argmax(np.where(msk,mp,-1e9)); iy,ix=np.unravel_index(j,mp.shape)
        rows[-2]['pk_dE']=dra[iy,ix]; rows[-2]['pk_dN']=dde[iy,ix]; rows[-2]['pk_m0']=mp[iy,ix]
    # ring null: positions at 9-16 arcsec from midpoint, 24 random, same HCO+ windows from N and S (use S window)
    vh,ah,cw=D['HCO+']; vc,ac,cwc=D['CH3OH']
    sS=sample(ah,cw,REG.loc['S','RA_deg'],REG.loc['S','Dec_deg']); pS=hco_fit(vh,sS); c,W=pS[1],pS[2]
    for k in range(40):
        r=rng.uniform(9,16); th=rng.uniform(0,2*np.pi); dra=r*np.sin(th)/3600/np.cos(np.deg2rad(CEN[1])); dde=r*np.cos(th)/3600
        s=sample(ah,cw,CEN[0]+dra,CEN[1]+dde); s2=sample(ac,cwc,CEN[0]+dra,CEN[1]+dde)
        Ih,_=m0(vh,s,c,W,excl=(c-W/2,c+W/2)); Ic,_=m0(vc,s2,c,W,excl=(c-W/2,c+W/2))
        ringr.append(dict(ous=ous,beam=bn,r=r,I_HCO=Ih,I_CH3OH=Ic))
    print(ous,bn,'done',flush=True)
T=pd.DataFrame(rows); T.to_csv('/workspace/out/q14e/m0_table.csv',index=False)
R=pd.DataFrame(ringr); R.to_csv('/workspace/out/q14e/ring.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(T.round(3).to_string())
print(R.groupby(['ous','beam']).agg(mean=('I_CH3OH','mean'),sd=('I_CH3OH','std'),mad=('I_CH3OH',robust_sigma),hco_med=('I_HCO','median')).round(3))
