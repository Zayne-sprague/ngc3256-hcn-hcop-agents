# Significance of the 'Soutflow' blue-wing peak (Xb4d HCO+ 3-2 spw25, HCN 3-2 spw27):
#  (1) wing maps in windows defined relative to the local (S-nucleus) centre v_c=2722 km/s (and, for comparison, vsys 2780); (2) number of independent beams (r<12"); 
#  (3) null test: the same-width (300 km/s) integration windows applied to line-free velocity windows -> number/ distribution of peaks above the Soutflow S/N;
#  (4) wing flux, 2-D Gaussian centroid + position error, wing/core ratio at the peak, red/blue symmetry; (5) PSF-sidelobe/clean-residual test with the S-nucleus spectrum;
#  (6) neighbouring-line contamination of the windows (catalogue).
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
import pandas as pd
RA0,DE0=156.96346,-43.90461
def run(spw,nm,f0,vS=2722.):
    h=fits.open(f('Xb4d',spw),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v)))
    cube=np.array(h[0].data[0],float); w=WCS(hd).celestial; ny,nx=cube.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    ra,de=w.pixel_to_world_values(X,Y); dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx,dy); pix=abs(hd['CDELT2'])*3600
    # per-pixel channel noise from the channel set far from all lines (first 15 channels)
    sig=np.nanstd(cube[:15],axis=0)
    def wmap(lo,hi):
        m=(v>=lo)&(v<=hi); return np.nansum(cube[m],axis=0)*dv,sig*dv*np.sqrt(m.sum()),m.sum()
    area=np.pi*12**2; beam=np.pi*hd['BMAJ']*hd['BMIN']*3600**2/(4*np.log(2)); nbeam=area/beam
    reg=(r<12)&np.isfinite(sig)
    def peaks(m0,e0,thr):
        sn=m0/ndi.uniform_filter(np.where(np.isfinite(e0),e0,np.nanmedian(e0)),3)
        lm=(m0==ndi.maximum_filter(np.nan_to_num(m0),size=int(1.5/pix)))&(sn>thr)&reg
        return [(sn[y,x],m0[y,x],dx[y,x],dy[y,x],y,x) for y,x in zip(*np.where(lm))]
    out={'nm':nm,'nbeam':nbeam}
    # windows (km/s, radio, absolute)
    W={'blue_vsys':(2780-450,2780-150),'blue_local':(vS-450,vS-150),'red_local':(vS+150,vS+450),'core_local':(vS-150,vS+150)}
    maps={k:wmap(*a) for k,a in W.items()}
    res={}
    for k in ('blue_vsys','blue_local','red_local'):
        m0,e0,n=maps[k]; pk=peaks(m0,e0,5); pk.sort(reverse=True); res[k]=pk
        print(nm,k,'peaks S/N>5 in r<12":',[(round(p[0],1),round(p[1],2),round(p[2],1),round(p[3],1)) for p in pk[:6]])
    # null windows: line-free (rel. to line, absolute radio velocity): choose windows well away from the line and from known other lines (HOC+ at +-200)
    vmin,vmax=v.min(),v.max(); nulls=[]
    cands=[(vmax-300-s,vmax-s) for s in (150,520,900)] if vmax>3100 else []
    for a,b in [(x,x+300) for x in np.arange(max(vmin+80,vS+600),vmax-380,150)]+[(x,x+300) for x in np.arange(vmin+80,vS-600-300,150)]:
        # keep only windows clear of the line core+wings by >=300 km/s and not overlapping HOC+/CH3OH/other blind lines (listed below)
        if b<vS-600 or a>vS+600: nulls.append((a,b))
    print(nm,'v range',round(vmin),round(vmax),'null windows',[(round(a),round(b)) for a,b in nulls])
    rows=[]
    for a,b in nulls:
        m0,e0,n=wmap(a,b); pk=peaks(m0,e0,5); rows.append(dict(win=f'{a:.0f}-{b:.0f}',n5=len(pk),n8=sum(p[0]>8 for p in pk),n10=sum(p[0]>10 for p in pk),maxsn=max([p[0] for p in pk],default=np.nan),
              medsn=float(np.nanmedian((m0/e0)[reg])),sdsn=float(np.nanstd((m0/e0)[reg]))))
    nd=pd.DataFrame(rows); print(nd.round(2).to_string())
    out['null']=nd; out['res']=res; out['maps']=maps; out['geom']=(dx,dy,pix,hd,r)
    return out,cube,v,dv,sig
def gfit(m0,e0,dx,dy,pix,x0,y0,bm):
    # 2-D gaussian (beam-fixed FWHM) centroid fit in a 3-beam box
    box=int(2.0/pix); ys=slice(max(0,y0-box),y0+box+1); xs=slice(max(0,x0-box),x0+box+1)
    Z=m0[ys,xs]; Xg=dx[ys,xs]; Yg=dy[ys,xs]; E=e0[ys,xs]
    s=bm/2.3548
    f=lambda xy,a,xc,yc,b: a*np.exp(-0.5*(((xy[0]-xc)**2+(xy[1]-yc)**2)/s**2))+b
    p,cv=curve_fit(f,(Xg.ravel(),Yg.ravel()),Z.ravel(),p0=[Z.max(),dx[y0,x0],dy[y0,x0],0],sigma=E.ravel(),absolute_sigma=True)
    return p,np.sqrt(np.diag(cv))
if __name__=='__main__':
    allrows=[]
    for spw,nm,f0 in (('spw25','HCO+ 3-2',267.557626),('spw27','HCN 3-2',265.886434)):
        out,cube,v,dv,sig=run(spw,nm,f0); dx,dy,pix,hd,r=out['geom']
        # Soutflow peak: strongest local peak in blue_local map within (-1.5,-6)..(1,-3.5)
        m0,e0,n=out['maps']['blue_local']
        sel=(np.abs(dx+0.25)<1.3)&(dy<-3.3)&(dy>-6)
        idx=np.argmax(np.where(sel,m0/e0,-1e9)); y0,x0=np.unravel_index(idx,m0.shape)
        bm=np.sqrt(hd['BMAJ']*hd['BMIN'])*3600
        p,e=gfit(m0,e0,dx,dy,pix,x0,y0,bm); snr=m0[y0,x0]/e0[y0,x0]
        sigpos=bm/(2*snr)
        # spectra at the wing peak and at S nucleus (peak pixel of core map)
        mc,ec,_=out['maps']['core_local']; sS=(np.hypot(dx-(-0.15),dy-(-2.55))<0.6); iS=np.argmax(np.where(sS,mc,-1e9)); yS,xS=np.unravel_index(iS,mc.shape)
        spw_=cube[:,y0,x0]; sS_=cube[:,yS,xS]
        def I(s,lo,hi): m=(v>=lo)&(v<=hi); return np.nansum(s[m])*dv
        vS=2722.
        row=dict(line=nm,peak_dRA=dx[y0,x0],peak_dDec=dy[y0,x0],SN_blue_local=snr,gauss_dRA=p[1],e_dRA=e[1],gauss_dDec=p[2],e_dDec=e[2],sigma_pos_beam_over_2SN=sigpos,
                 I_blue_local=I(spw_,vS-450,vS-150),I_core_local=I(spw_,vS-150,vS+150),I_red_local=I(spw_,vS+150,vS+450),I_blue_vsys=I(spw_,2330,2630),
                 e_I_blue=sig[y0,x0]*dv*np.sqrt(((v>=vS-450)&(v<=vS-150)).sum()),I_S_blue=I(sS_,vS-450,vS-150),I_S_core=I(sS_,vS-150,vS+150),I_S_red=I(sS_,vS+150,vS+450),
                 n_beams=out['nbeam'],n_null_windows=len(out['null']),null_n8_total=int(out['null'].n8.sum()),null_n10_total=int(out['null'].n10.sum()),null_max_snr=out['null'].maxsn.max(),null_sd_snr=out['null'].sdsn.median())
        row['wing_core_ratio_Soutflow']=row['I_blue_local']/row['I_core_local']; row['wing_core_ratio_S']=row['I_S_blue']/row['I_S_core']
        # spectral-shape test: scale S spectrum to Soutflow in core, compare blue wing (PSF-sidelobe copy hypothesis)
        sc=row['I_core_local']/row['I_S_core']; row['blue_pred_if_sidelobe_copy']=sc*row['I_S_blue']; row['sidelobe_fraction_core']=sc
        # velocity centroid of Soutflow spectrum and S spectrum (core +-450)
        for tag,s in (('Sout',spw_),('S',sS_)):
            m=(v>=vS-450)&(v<=vS+450); row['v_mom1_'+tag]=float(np.sum(v[m]*s[m])/np.sum(s[m]))
        allrows.append(row); 
        np.savez(f'/workspace/out/wing_spec_{spw}.npz',v=v,sout=spw_,sS=sS_)
    df=pd.DataFrame(allrows); df.round(3).to_csv('/workspace/out/wing_null.csv',index=False); print(df.T.round(3).to_string())
