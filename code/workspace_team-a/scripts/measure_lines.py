# Standard line table: every target line x every region x every cube (native resolution), peak-pixel and aperture spectra, pbcor.
import sys,warnings,json,glob; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from lineset import *
import pandas as pd
from scipy.optimize import curve_fit
OUS=['Xb55','Xb47','Xb4f','Xb4d','Xb4b']
BAND={'Xb55':'B5','Xb47':'B7','Xb4f':'B6','Xb4d':'B6','Xb4b':'B6'}
SPWS={'Xb55':['spw16','spw18','spw20','spw22'],'Xb47':['spw16','spw18','spw20','spw22'],'Xb4f':['spw16','spw18','spw20','spw22'],
      'Xb4d':['spw25','spw27','spw29','spw31'],'Xb4b':['spw25','spw27','spw29','spw31']}
REF={'Xb55':'spw18','Xb47':'spw18','Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'}   # HCO+ cube used to locate peak pixel
REFF0={'Xb55':178.375056,'Xb47':356.734223,'Xb4f':267.557626,'Xb4d':267.557626,'Xb4b':267.557626}
def region_radec(dx,dy):
    return RA0+dx/3600/np.cos(np.deg2rad(DE0)), DE0+dy/3600
def vel(fr,f0): return (1-fr/f0)*C
def gauss(v,a,v0,fw): return a*np.exp(-4*np.log(2)*(v-v0)**2/fw**2)
def _fit_I(rel,s,W,order,excl,dv):
    """baseline from channels W+150..W+1150 away (excluding other lines), integrate |rel|<=W; returns I, rms, baseline-subtracted spectrum, masks"""
    bl=((np.abs(rel)>W+150)&(np.abs(rel)<W+1150))&np.isfinite(s)
    for lo,hi in excl: bl&=~((rel>lo)&(rel<hi))
    if bl.sum()<30: return None
    p=np.polyfit(rel[bl],s[bl],order if bl.sum()>60 else 0); sc=s-np.polyval(p,rel)
    r=sc[bl]
    for _ in range(4):
        sg=np.std(r); r=r[np.abs(r)<3*sg]
    w=(np.abs(rel)<=W)&np.isfinite(s)
    return np.sum(sc[w])*dv,np.std(r),sc,w,bl
def analyse(v,s,dv,f0name,W,order=1,vexcl=(),null_specs=None):
    """v: velocity rel. to rest f0 (km/s radio); s: spectrum (Jy/beam). Detection significance from an empirical null (same estimator on line-free windows)."""
    rel=v-VSYS
    r=_fit_I(rel,s,W,order,vexcl,dv)
    if r is None: return None
    I,sg,sc,w,bl=r; n=int(w.sum())
    if n<3: return None
    d=dict(rms_chan=sg,nchan_int=n,I_sum=I,eI_stat=sg*dv*np.sqrt(n),nchan_base=int(bl.sum()))
    # empirical null: same estimator applied to emission-free pixels of the same cube (spatial null)
    nulls=[]
    if null_specs is not None:
        for k in range(null_specs.shape[1]):
            rr=_fit_I(rel,null_specs[:,k],W,order,vexcl,dv)
            if rr is not None and np.isfinite(rr[0]): nulls.append(rr[0])
    nulls=np.array(nulls)
    d['n_null']=len(nulls)
    if len(nulls)>=20:
        d['eI_emp']=1.4826*np.median(np.abs(nulls-np.median(nulls))); d['null_mean']=float(np.median(nulls))
    else: d['eI_emp']=np.nan; d['null_mean']=np.nan
    e=np.nanmax([d['eI_stat'],d['eI_emp']]) if np.isfinite(d['eI_emp']) else d['eI_stat']
    d['eI_sum']=e; d['snr_I']=(I-(d['null_mean'] if np.isfinite(d['null_mean']) else 0))/e
    pk=np.max(sc[w]); ipk=np.argmax(np.where(w,sc,-1e9))
    d.update(peak=pk,snr_peak=pk/sg,v_peak=rel[ipk]+VSYS)
    if I>0:
        d['v_mom1']=float(np.sum(rel[w]*sc[w])/np.sum(sc[w])+VSYS)
        sm=np.convolve(np.where(w,sc,0),np.ones(3)/3,'same'); half=sm.max()/2
        above=np.where((sm>half)&w)[0]
        d['W50']=float((above.max()-above.min()+1)*dv) if len(above) else np.nan
    if d['snr_I']>=3:
        try:
            p0=[max(pk,sg),0.,max(120,dv*2)]
            pp,cv=curve_fit(gauss,rel[w],sc[w],p0=p0,sigma=np.full(w.sum(),sg),absolute_sigma=True,maxfev=5000,bounds=([0,-400,10],[10,400,1200]))
            e=np.sqrt(np.diag(cv)); d.update(g_amp=pp[0],g_v0=pp[1]+VSYS,g_fwhm=pp[2],e_g_amp=e[0],e_g_v0=e[1],e_g_fwhm=e[2],g_I=pp[0]*pp[2]*1.0645)
        except Exception as ex: pass
    return d

def get_cube(ous,spw):
    h=fits.open(f(ous,spw),memmap=True); return h[0].header,h[0].data[0]

def beam_pix(hd):
    pix=abs(hd['CDELT2'])*3600; return np.pi*hd['BMAJ']*3600*hd['BMIN']*3600/(4*np.log(2))/pix**2,pix

def region_pixels(ous,hdref,cuberef,wref,dx,dy):
    """peak pixel of HCO+ moment-0 near region + aperture pixel mask."""
    ra,de=region_radec(dx,dy); x0,y0=[float(a) for a in wref.world_to_pixel_values(ra,de)]
    return x0,y0

if __name__=='__main__':
    rows=[]; spec_store={}; SPEC={}
    for ous in OUS:
        band=BAND[ous]
        hdr,cr=get_cube(ous,REF[ous]); w=WCS(hdr).celestial; fr=freqaxis(hdr); vr=vel(fr/1e9,REFF0[ous]); dvr=abs(np.median(np.diff(vr)))
        beam,pix=beam_pix(hdr); bmaj=hdr['BMAJ']*3600
        # HCO+ moment-0 (vsys+-350) for peak location
        sel=np.where(np.abs(vr-VSYS)<=350)[0]
        m0=np.nansum(np.array(cr[sel.min():sel.max()+1]),axis=0)
        ny,nx=m0.shape; Y,X=np.mgrid[:ny,:nx]
        rap=1.5   # aperture radius arcsec (same for all cubes; coarse cubes are beam-limited)
        print(ous,'beam',round(bmaj,2),'pix',round(pix,3),'aperture r',round(rap,2),flush=True)
        # null pixel selection
        ra_,de_=w.pixel_to_world_values(X,Y); rr_=np.hypot((ra_-RA0)*np.cos(np.deg2rad(DE0))*3600,(de_-DE0)*3600)
        rmax={'Xb55':17.,'Xb47':9.,'Xb4f':11.,'Xb4d':11.,'Xb4b':11.}[ous]; rmin=6.
        cand=(rr_>rmin)&(rr_<rmax)&np.isfinite(m0)
        thr=np.nanmedian(np.abs(m0[cand]-np.nanmedian(m0[cand])))*1.4826*2.5+np.nanmedian(m0[cand])
        cand&=(m0<thr)
        yy,xx=np.where(cand); rng=np.random.default_rng(1)
        if len(yy)>300: sel_=rng.choice(len(yy),300,replace=False); yy,xx=yy[sel_],xx[sel_]
        o=np.argsort(yy*1000+xx); yy,xx=yy[o],xx[o]; nullpix=(xx,yy); print('  null pixels',len(xx))
        regpix={}
        for rn,(dx,dy) in REGIONS.items():
            ra,de=region_radec(dx,dy); xc,yc=[float(a) for a in w.world_to_pixel_values(ra,de)]
            rs=max(1.0,0.3*bmaj)/pix
            m=((X-xc)**2+(Y-yc)**2<=rs**2)&np.isfinite(m0)
            if m.sum()==0: continue
            idx=np.argmax(np.where(m,m0,-np.inf)); yp,xp=np.unravel_index(idx,m0.shape)
            ap=((X-xp)**2+(Y-yp)**2<=(rap/pix)**2)
            regpix[rn]=(xp,yp,ap,float(w.pixel_to_world_values(xp,yp)[0]),float(w.pixel_to_world_values(xp,yp)[1]))
        for spw in SPWS[ous]:
            hd,cube=get_cube(ous,spw); fr=freqaxis(hd)/1e9; dv=abs(np.median(np.diff(fr)))/np.mean(fr)*C
            # null pixels: emission-free (HCO+ m0 < 2 sigma_med) pixels within 6-12" (or 6-17" Xb55) of phase centre
            Xn,Yn=nullpix
            null_all=np.array(cube[:,Yn,Xn],float) if len(Xn) else None
            # targets inside this cube
            tg=[]
            for (nm,f0,Wt,bl) in TARGETS[band]:
                fs=f0*(1-VSYS/C); 
                if fr.min()+0.02<fs-Wt/C*fs*0 and fs>fr.min() and fs<fr.max(): tg.append((nm,f0,Wt,bl))
            if not tg: continue
            # velocity-masked blends: exclude all other targets' +-400 km/s windows from baseline
            for rn,(xp,yp,ap,rap_ra,rap_de) in regpix.items():
                ys,xs=np.where(ap); y0,y1,x0,x1=ys.min(),ys.max()+1,xs.min(),xs.max()+1
                sub=np.array(cube[:,y0:y1,x0:x1],float); apm=ap[y0:y1,x0:x1]
                s_pk=np.array(cube[:,yp,xp],float); s_ap=np.nanmean(sub[:,apm],axis=1)
                SPEC[f'{ous}|{spw}|{rn}|pk']=s_pk; SPEC[f'{ous}|{spw}|{rn}|ap']=s_ap; SPEC[f'{ous}|{spw}|freq']=fr
                for kind,s,fac in (('peakpix',s_pk,1.0),('aperture',s_ap,1.0)):
                    for (nm,f0,Wt,bl) in tg:
                        v=vel(fr,f0)
                        excl=[]
                        for (n2,f2,W2,b2) in tg:
                            if n2!=nm and n2.split()[0] in ('HCN','HCO+','CS','H2CS','H13CN','HOC+','NS','SO') and not n2.startswith('HCN v2'):
                                off=(vel(f2*(1-VSYS/C),f0)) # velocity of other line relative to this rest freq at vsys
                                excl.append((off-VSYS-W2-100+0,off-VSYS+W2+100))
                        # shift: rel coordinate is v-VSYS; other line centre at v_other-VSYS
                        res=analyse(v,s,dv,nm,Wt,order=(2 if 'H13CN' in nm else 1),vexcl=excl,null_specs=(null_all if kind=='peakpix' else None))
                        if res is None: continue
                        res.update(ous=ous,band=band,spw=spw,region=rn,kind=kind,line=nm,f0=f0,blends=bl,dv=dv,
                                   x=xp,y=yp,ra=rap_ra,dec=rap_de,bmaj=bmaj,bmin=hd['BMIN']*3600,aperture_arcsec=(rap if kind=='aperture' else 0),
                                   beam_pix=beam)
                        rows.append(res)
                        if kind=='peakpix': spec_store[f'{ous}|{spw}|{rn}|{nm}']=(v,s)
    np.savez('/workspace/out/region_spectra.npz',**{k.replace('|','__'):v for k,v in SPEC.items()})
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/line_table_raw.csv',index=False); print(len(df))
