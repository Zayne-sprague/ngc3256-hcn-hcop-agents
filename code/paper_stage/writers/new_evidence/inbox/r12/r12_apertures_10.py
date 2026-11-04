"""r12 (i): HCN/HCO+ J=1-0 N/S contrast on the short-spacing-complete Band 3 cubes (2015.1.00412 X145, 2015.1.00993 X961) at circular beams 1.5, 2.2, 2.5".
Same estimators as r11_apertures.py: centre = HCO+ moment-0 peak within 0.4" of nominal (window +-350 km/s about 2780, also used for the +-150 window),
peak pixel; circular apertures r=0.25,0.5,0.75" and r=beam (spline on 0.05" grid, as r11); RJ T_b (flux-conserving convolution to the circular beam, s180 helpers);
baseline = median of channels >600 km/s from HCN and HCO+; errors from line-free channels of the aperture-mean spectrum x dv x sqrt(sum k^2).
A beam is skipped if it is smaller than the native beam (X961 1.5": native 2.10x1.67).  Usage: python3 -I r12_apertures_10.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import load, to_tb_conv
from regions import REGIONS
from scipy.ndimage import map_coordinates
C=299792.458; V0=2780.; F0={'HCN':88.6318,'HCOp':89.1885}
DS={'X145':'/archive/files/2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',
    'X961':'/archive/files/2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits'}
def vel(nu,f0): return C*(1-nu/f0)-V0
def kw(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
SUB=0.05; res={}; lines=[]
for dn,path in DS.items():
    cut,nu,pix,bm,pxy=load(path,half=14)
    for th in (1.5,2.2,2.5):
        if th<max(bm[:2])*1.001: continue
        Tb=to_tb_conv(cut,nu,pix,bm,th,np.mean(nu))
        vH=vel(nu,F0['HCN']); vO=vel(nu,F0['HCOp']); free=(np.abs(vH)>600)&(np.abs(vO)>600)
        Tb=Tb-np.median(Tb[free],axis=0); ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]
        dvO=abs(vO[1]-vO[0]); k350=kw(vO,dvO,-350,350); m0=np.tensordot(k350,Tb,axes=(0,0))*dvO
        g=np.arange(-th*1.2,th*1.2+1e-9,SUB)/pix; gy,gx=np.meshgrid(g,g,indexing='ij'); rr=np.hypot(gy,gx)*pix
        for rn in ('N','S'):
            px,py=pxy(REGIONS[rn][:2]); mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); cx,cy=float(X.flat[j]),float(Y.flat[j])
            o={}
            for ln,v in (('HCN',vH),('HCOp',vO)):
                dv=abs(v[1]-v[0])
                cube=np.stack([map_coordinates(Tb[i],[gy+cy,gx+cx],order=3) for i in range(Tb.shape[0])])
                sps=[('peak',Tb[:,int(cy),int(cx)])]+[(f'r{r}',cube[:,rr<=r].mean(axis=1)) for r in (0.25,0.5,0.75,th)]
                for lab,sp in sps:
                    sig=np.std(sp[free])
                    for wn,w in (('w350',350),('w150',150)):
                        k=kw(v,dv,-w,w); o[f'{ln}|{lab}|{wn}']=[float((k*sp).sum()*dv),float(sig*np.sqrt((k**2).sum())*dv)]
            res[f'{dn}|{th}|{rn}']=o
        for wn in ('w350','w150'):
            for lab in ('peak','r0.25','r0.5','r0.75','rbeam'):
                L=f'r{th}' if lab=='rbeam' else lab; row={}
                for rn in ('N','S'):
                    h=res[f'{dn}|{th}|{rn}'][f'HCN|{L}|{wn}']; q=res[f'{dn}|{th}|{rn}'][f'HCOp|{L}|{wn}']
                    R=h[0]/q[0]; row[rn]=(R,abs(R)*np.hypot(h[1]/h[0],q[1]/q[0]))
                c=row['N'][0]/row['S'][0]; ce=c*np.hypot(row['N'][1]/row['N'][0],row['S'][1]/row['S'][0])
                lines.append(f'{dn} beam {th}" {wn} {lab:6s} N {row["N"][0]:.3f}+-{row["N"][1]:.3f} S {row["S"][0]:.3f}+-{row["S"][1]:.3f} N/S {c:.3f}+-{ce:.3f}')
                res[f'NS|{dn}|{th}|{wn}|{lab}']=[c,ce,row['N'][0],row['S'][0]]
open('r12_apertures_10.txt','w').write('\n'.join(lines)); print('\n'.join(lines)); json.dump(res,open('r12_apertures_10.json','w'))
