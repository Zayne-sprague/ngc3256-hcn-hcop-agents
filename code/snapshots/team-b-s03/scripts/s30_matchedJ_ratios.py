# Q1: HCN/HCO+ at matched J. Each pbcor cube is cropped, per-pixel baseline-subtracted (line-free channels), convolved (analytic Fourier Gaussian, flux
# rescaled by beam-area ratio as in s26) to a common circular beam, integrated over |v-VSYS|<W km/s, and read at region positions / apertures.
# Output: out/s30_ratios.json, out/s30_ratios.txt, out/s30_maps_<tag>.npz  (maps for the regression script s31)
import sys, json; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy.ndimage import map_coordinates
HALF=24.0
CEN=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]
F0={'HCN':{2:177.26111e9,3:265.88618e9,4:354.50547e9},'HCO+':{2:178.37506e9,3:267.55763e9,4:356.73424e9}}
SPW={ # (ous, J) -> (HCN spw, HCO+ spw)
 ('Xb55',2):('spw16','spw18'),('Xb4f',3):('spw18','spw16'),('Xb4d',3):('spw27','spw25'),('Xb47',4):('spw16','spw18')}
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def transfer(shape,ps,b0,tgt):
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-8,None))@U.T
    fy=np.fft.fftfreq(shape[0],ps); fx=np.fft.fftfreq(shape[1],ps); FX,FY=np.meshgrid(fx,fy)
    q=np.stack([-FX.ravel(),FY.ravel()]); return np.exp(-2*np.pi**2*np.einsum('ij,ij->j',q,C2@q)).reshape(shape)
def load(ous,spw,f0):
    fn=files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]; h_,hd,w,nu=cube_info(fn)
    v=vrad(nu,f0)-VSYS; ps=abs(hd['CDELT1'])*3600; wc=WCS(hd).celestial
    cx,cy=wc.world_to_pixel_values(*CEN); r=int(min(HALF,(min(hd["NAXIS1"],hd["NAXIS2"])//2-2)*ps)/ps); x0,y0=int(round(cx-r)),int(round(cy-r)); n=2*r
    ok=np.where((v>-1000)&(v<700))[0]; i0,i1=ok[0],ok[-1]+1
    cube=np.nan_to_num(np.array(h_[0].data[0,i0:i1,max(y0,0):y0+n,max(x0,0):x0+n],dtype=np.float32))
    assert cube.shape[1:]==(n,n),(ous,spw,cube.shape)
    return cube,v[i0:i1],ps,wc,x0,y0,(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA'])
def conv(cube,ps,b0,tgt):
    n=cube.shape[1]; P=2*n; T=transfer((P,P),ps,b0,tgt)[:,:P//2+1]; fac=tgt[0]*tgt[1]/(b0[0]*b0[1]); out=np.empty_like(cube)
    for i in range(cube.shape[0]):
        a=np.zeros((P,P),np.float32); a[:n,:n]=cube[i]; out[i]=np.fft.irfft2(np.fft.rfft2(a)*T,s=(P,P))[:n,:n]*fac
    return out
def products(cube,v,W):
    ref=(v<-600)|((v>450)&(v<600)); core=abs(v)<W; dv=abs(np.median(np.diff(v)))
    c=cube-np.median(cube[ref],axis=0)
    return c,core,dv,ref
def run(ous,J,fwhm,windows=(250,350,500)):
    out={}
    for mol,spw in zip(('HCN','HCO+'),SPW[(ous,J)]):
        cube,v,ps,wc,x0,y0,b0=load(ous,spw,F0[mol][J])
        tgt=(fwhm,fwhm,0.0) if fwhm else b0
        if fwhm and (b0[0]>fwhm or b0[1]>fwhm): raise SystemExit(f'target beam smaller than native {b0}')
        cc=conv(cube,ps,b0,tgt) if fwhm else cube
        c,core,dv,ref=products(cc,v,350)
        n=c.shape[1]; Y,X=np.mgrid[:n,:n]
        def pxy(ra,dec): px,py=wc.world_to_pixel_values(ra,dec); return px-x0,py-y0
        pc=pxy(*CEN); rr=np.hypot(X-pc[0],Y-pc[1])*ps
        beamarea=1.1331*tgt[0]*tgt[1]/ps**2
        res={'beam':tgt,'native':b0,'dv':dv,'ps':ps}
        for W in windows:
            m=abs(v)<W; m0=c[m].sum(0)*dv; sig=c[ref].std(0)   # per-pixel per-channel rms from line-free channels
            res[f'm0_{W}']=m0; res[f'err_{W}']=sig*np.sqrt(m.sum())*dv
            if W==350: res['rms_chan_map']=sig
        res['cube_lf']=c[ref]; res['rr']=rr; res['beamarea']=beamarea; res['pxy']=pxy; res['c']=c; res['v']=v
        out[mol]=res
    return out
def val_at(m,xy): return float(map_coordinates(m,[[xy[1]],[xy[0]]],order=1)[0])
def peak_near(m,xy,rad_pix):
    n=m.shape[0]; Y,X=np.mgrid[:n,:n]; mk=np.hypot(X-xy[0],Y-xy[1])<=rad_pix
    idx=np.unravel_index(np.argmax(np.where(mk,m,-1e9)),m.shape); return idx[1],idx[0]
def ratio(a,ea,b,eb,cal=0.05):
    r=a/b; e=abs(r)*np.hypot(ea/a,eb/b) if a>0 and b>0 else np.nan
    return r,e,abs(r)*np.hypot(e/abs(r),cal) if np.isfinite(e) else np.nan
KFAC=lambda J:(F0['HCO+'][J]/F0['HCN'][J])**2   # R_Tb = R_Jy * (nu_HCO+/nu_HCN)^2
if __name__=='__main__':
    JOBS=[('Xb55',2,10.0),('Xb4f',3,10.0),('Xb4d',3,10.0),('Xb47',4,10.0),('Xb4f',3,6.5),('Xb4d',3,6.5),('Xb47',4,6.5),('Xb4d',3,None),('Xb4f',3,None),('Xb47',4,None),('Xb55',2,None)]
    allres=[]; txt=[]
    for ous,J,fwhm in JOBS:
        o=run(ous,J,fwhm); H,O=o['HCN'],o['HCO+']
        tag=f'{ous}_J{J}_'+(f'{fwhm:g}as' if fwhm else 'native'); b=H['beam']
        txt.append(f'## {tag}: target beam {b[0]:.2f}x{b[1]:.2f}  native HCN {H["native"][0]:.2f}x{H["native"][1]:.2f} HCO+ {O["native"][0]:.2f}x{O["native"][1]:.2f} dv {H["dv"]:.1f} km/s')
        rows={}
        locs={k:REGIONS[k] for k in ('N','S','Bridge','E1','E2','W','Sblue')}; locs['mid']=(CEN[0],CEN[1],'N-S midpoint')
        for name,(ra,dec,_) in locs.items():
            for kind in ('pos','peak'):
                vals={}
                for W in (250,350,500):
                    s={}
                    for mol,R in (('HCN',H),('HCO+',O)):
                        xy=R['pxy'](ra,dec)
                        if kind=='peak': xy=peak_near(O['m0_350'],xy,max(1.0,1.0/R['ps']))  # location of the HCO+ peak within 1" (same pixel used for both lines)
                        s[mol]=(val_at(R[f'm0_{W}'],xy),val_at(R[f'err_{W}'],xy),xy)
                    vals[W]=s
                if kind=='peak':   # same pixel for both: recompute HCN at the HCO+ peak pixel
                    pass
                a,ea,_=vals[350]['HCN']; bb,eb,_=vals[350]['HCO+']
                r,e,et=ratio(a,ea,bb,eb)
                rows[f'{name}_{kind}']=dict(HCN=a,eHCN=ea,HCOp=bb,eHCOp=eb,R_Jy=r,eR_stat=e,eR_tot=et,R_K=r*KFAC(J),eR_K=et*KFAC(J),
                    R_Jy_w250=vals[250]['HCN'][0]/vals[250]['HCO+'][0],R_Jy_w500=vals[500]['HCN'][0]/vals[500]['HCO+'][0],
                    snr_min=min(a/ea,bb/eb) if ea>0 and eb>0 else np.nan)
                txt.append(f'  {name:7s}{kind:5s} HCN {a:7.3f}+-{ea:.3f}  HCO+ {bb:7.3f}+-{eb:.3f} Jy/b km/s  R_Jy={r:.3f}+-{e:.3f}(stat) +-{et:.3f}(tot)  R_K={r*KFAC(J):.3f}  [w250 {rows[f"{name}_{kind}"]["R_Jy_w250"]:.3f} w500 {rows[f"{name}_{kind}"]["R_Jy_w500"]:.3f}]  S/N {rows[f"{name}_{kind}"]["snr_min"]:.0f}')
        # apertures about the midpoint and annulus
        for lab,(lo,hi) in {'r<4':(0,4),'r<8':(0,8),'r<12':(0,12),'8<r<20':(8,20)}.items():
            S={}
            for mol,R in (('HCN',H),('HCO+',O)):
                mk=(R['rr']>=lo)&(R['rr']<hi); tot=R['m0_350'][mk].sum()/R['beamarea']
                lf=R['cube_lf'][:,mk].sum(1)/R['beamarea']; dv=R['dv']; ncore=int((abs(R['v'])<350).sum())
                S[mol]=(tot,lf.std()*dv*np.sqrt(ncore))
            r,e,et=ratio(S['HCN'][0],S['HCN'][1],S['HCO+'][0],S['HCO+'][1])
            rows[f'ap_{lab}']=dict(HCN=S['HCN'][0],eHCN=S['HCN'][1],HCOp=S['HCO+'][0],eHCOp=S['HCO+'][1],R_Jy=r,eR_stat=e,eR_tot=et,R_K=r*KFAC(J),eR_K=et*KFAC(J))
            txt.append(f'  ap {lab:7s} HCN {S["HCN"][0]:8.2f}+-{S["HCN"][1]:.2f}  HCO+ {S["HCO+"][0]:8.2f}+-{S["HCO+"][1]:.2f} Jy km/s  R_Jy={r:.3f}+-{e:.3f} (tot {et:.3f}) R_K={r*KFAC(J):.3f}')
        allres.append(dict(tag=tag,ous=ous,J=J,fwhm=fwhm,beam=list(b[:2]),rows=rows))
        np.savez(f'/workspace/out/s30_maps_{tag}.npz',HCN=H['m0_350'],HCOp=O['m0_350'],eHCN=H['err_350'],eHCOp=O['err_350'],rr=H['rr'],ps=H['ps'],
                 cx=H['pxy'](*CEN)[0],cy=H['pxy'](*CEN)[1],beam=np.array(b[:2]))
        print('\n'.join(txt[-1:]),flush=True)
    json.dump(allres,open('/workspace/out/s30_ratios.json','w'),indent=1,default=float)
    open('/workspace/out/s30_ratios.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
