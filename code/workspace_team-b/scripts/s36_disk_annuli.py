# Q1 (round 2): disk HCN/HCO+ and radial trend at a common 10" circular beam, with an explicit PB-noise model and zero-level (bowl) test.
# Arrays: Xb55 (J=2-1), Xb4f/Xb4d (3-2), Xb47 (4-3). Annuli centred on the N-S midpoint. PB(r) is estimated from the line-free per-pixel rms (pbcor cubes: rms ~ 1/PB),
# normalised to the central value, and compared with a Gaussian 12-m PB (FWHM=1.13 lam/D). Zero level: mean moment-0 (|v|<350 km/s, same window/baseline) of the far
# field (r_far..max, PB>0.15) -> bowl flux inside an annulus = mean_far x area/beam area. Go/no-go for an annulus: Gaussian-model PB>=0.5 (the line-free rms of these pbcor cubes does NOT rise as 1/PB, so PBrms is only reported), annulus flux >3 sigma in both lines
# (sigma includes the bowl uncertainty), and |bowl correction| <30% of the annulus flux. Output out/s36_disk.json/.txt
import sys, json; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import *
S.HALF=45.0
ANN=[(0,4),(4,8),(8,12),(12,17),(17,24)]
NU={'Xb55':177.3,'Xb4f':265.9,'Xb4d':265.9,'Xb47':354.5}
FAR={'Xb55':(30,45),'Xb4f':(26,40),'Xb4d':(17,19),'Xb47':(17,24)}
txt=['## HCN/HCO+ (Jy km/s ratio; x1.0126 for K) in annuli about the N-S midpoint, 10" beam. Flux in Jy km/s (sum/beam area). PBrms = 1/(rms/rms_centre) from line-free channels; PBg = Gaussian 12-m model']
out={}
for ous,J in (('Xb55',2),('Xb4f',3),('Xb4d',3),('Xb47',4)):
    o=S.run(ous,J,10.0,windows=(350,)); H,O=o['HCN'],o['HCO+']; tag=f'{ous}_J{J}'
    fw=1.13*299792458/(NU[ous]*1e9)/12*206265
    rr=H['rr']; rms=H['rms_chan_map']; c0=np.nanmedian(rms[rr<3])
    txt.append(f'## {tag}: map half-width {rr.max()/1.0:.0f}" ; Gaussian PB FWHM {fw:.1f}"')
    rows={}
    far_lo,far_hi=FAR[ous]
    zero={}
    for mol,R in (('HCN',H),('HCO+',O)):
        mk=(R['rr']>=far_lo)&(R['rr']<far_hi)&(R['rms_chan_map']>0)
        v=R['m0_350'][mk]; zero[mol]=(np.mean(v),np.std(v)/np.sqrt(max(mk.sum()/R['beamarea'],1)),np.std(v),mk.sum())
    for lo,hi in ANN:
        r_mid=0.5*(lo+hi); mk=(H['rr']>=lo)&(H['rr']<hi)&(H['rms_chan_map']>0)
        if mk.sum()<5: continue
        pbr=c0/np.nanmedian(H['rms_chan_map'][mk]); pbg=np.exp(-4*np.log(2)*(r_mid/fw)**2)
        S_={}
        for mol,R in (('HCN',H),('HCO+',O)):
            f=R['m0_350'][mk].sum()/R['beamarea']; npix=mk.sum()
            # noise: line-free-channel pixel-sum scatter (correlated pixels) from the cube itself
            lf=R['cube_lf'][:,mk].sum(1)/R['beamarea']; sig=lf.std()*R['dv']*np.sqrt((abs(R['v'])<350).sum())
            bowl=zero[mol][0]*npix/R['beamarea']; ebowl=zero[mol][1]*npix/R['beamarea']
            S_[mol]=(f,sig,bowl,ebowl)
        (fh,sh,bh,ebh),(fo,so,bo,ebo)=S_['HCN'],S_['HCO+']
        r,e,_=ratio(fh,sh,fo,so,0.0)
        fhc,foc=fh-bh,fo-bo; ehc=np.hypot(sh,ebh); eoc=np.hypot(so,ebo)
        rc,erc,_=ratio(fhc,ehc,foc,eoc,0.0) if fhc>0 and foc>0 else (np.nan,np.nan,np.nan)
        sig_ok=(fhc>3*ehc) and (foc>3*eoc); bowl_ok=abs(bh)<0.3*abs(fh) and abs(bo)<0.3*abs(fo)
        go=(pbg>=0.5) and sig_ok and bowl_ok
        rows[f'{lo}-{hi}']=dict(PBrms=pbr,PBg=pbg,HCN=fh,eHCN=sh,HCOp=fo,eHCOp=so,bowlHCN=bh,bowlHCOp=bo,ebowlHCN=ebh,ebowlHCOp=ebo,R=r,eR=e,R_bowlcorr=rc,eR_bowlcorr=erc,go=bool(go))
        txt.append(f'  {lo:2d}-{hi:2d}" PBrms {pbr:.2f} PBg {pbg:.2f}  HCN {fh:8.2f}+-{sh:.2f} HCO+ {fo:8.2f}+-{so:.2f}  bowl(HCN,HCO+) {bh:6.2f},{bo:6.2f}+-{ebh:.2f},{ebo:.2f}  R={r:.3f}+-{e:.3f}  R(bowl-corr)={rc:.3f}+-{erc:.3f}  {"GO" if go else "NO-GO"}')
    txt.append(f'  far-field zero level ({far_lo}-{far_hi}") mean m0 [Jy/beam km/s]: HCN {zero["HCN"][0]:.3f}+-{zero["HCN"][1]:.3f} (pix rms {zero["HCN"][2]:.2f}, npix {zero["HCN"][3]}), HCO+ {zero["HCO+"][0]:.3f}+-{zero["HCO+"][1]:.3f} (pix rms {zero["HCO+"][2]:.2f})')
    out[tag]=dict(rows=rows,zero={k:list(map(float,v)) for k,v in zero.items()},pb_fwhm=fw)
    # sectors in 8-14" annulus (4 quadrants about the midpoint; +x = RA pixel axis (E to the left))
    n=H['rr'].shape[0]; Y,X=np.mgrid[:n,:n]; cx,cy=H['pxy'](*CEN); ang=np.degrees(np.arctan2(-(X-cx),(Y-cy)))  # 0=N, +90=E
    for nm,(a0,a1) in {'N':(-45,45),'E':(45,135),'S':(135,225),'W':(-135,-45)}.items():
        aa=ang%360; a0_,a1_=a0%360,a1%360
        sec=(aa>=a0_)&(aa<a1_) if a0_<a1_ else ((aa>=a0_)|(aa<a1_))
        mk=sec&(H['rr']>=6)&(H['rr']<14)&(H['rms_chan_map']>0)
        if mk.sum()<5: continue
        fh=H['m0_350'][mk].sum()/H['beamarea']; fo=O['m0_350'][mk].sum()/O['beamarea']
        txt.append(f'  sector {nm} 6-14": HCN {fh:7.2f} HCO+ {fo:7.2f} R={fh/fo:.3f}')
    print('\n'.join(txt[-14:]),flush=True)
# continuum check: fraction of the r<8" flux that sits in each annulus (10" beam), per array, from the line-free continuum image
txt.append('## continuum (spw*.cont pbcor, 10" beam): annulus flux / (r<8" flux)  [line: HCO+ same ratio printed above via fluxes]')
for ous in ('Xb55','Xb4f','Xb4d','Xb47'):
    fn=files(f'*{ous}.ngc_3256_sci.spw*_*.cont.I.pbcor.fits')[0]; d,hd=load2d(fn); ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    b0=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); cx,cy=wc.world_to_pixel_values(*CEN); n=d.shape[0]
    P=2*n; T=S.transfer((P,P),ps,b0,(10.,10.,0.))[:,:P//2+1]; a=np.zeros((P,P),np.float32); a[:n,:n]=np.nan_to_num(d)
    cc=np.fft.irfft2(np.fft.rfft2(a)*T,s=(P,P))[:n,:n]; Y,X=np.mgrid[:n,:n]; rr=np.hypot(X-cx,Y-cy)*ps; ba=1.1331*100/ps**2
    ref=cc[rr<8].sum()/ba*(100/(b0[0]*b0[1]))
    fr=[ (cc[(rr>=lo)&(rr<hi)].sum()/ba*(100/(b0[0]*b0[1])))/ref for lo,hi in ANN[:4]]
    txt.append(f'  {ous} cont: '+'  '.join(f'{lo}-{hi}": {f:.3f}' for (lo,hi),f in zip(ANN[:4],fr)))
    out[f'cont_{ous}']=dict(zip([f'{lo}-{hi}' for lo,hi in ANN[:4]],fr))
print('\n'.join(txt[-5:]))
json.dump(out,open('/workspace/out/s36_disk.json','w'),indent=1,default=float)
open('/workspace/out/s36_disk.txt','w').write('\n'.join(txt))
