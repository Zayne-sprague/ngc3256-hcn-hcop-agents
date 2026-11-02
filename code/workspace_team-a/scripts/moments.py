# Moment maps (0: integrated over vsys±VW, pbcor) for all strong-line cubes; saves FITS-like npz + peak lists + png
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage as ndi
VSYS=2780.; VW=350.
LINES=[('Xb4f','spw16','HCO+ 3-2',267.557626),('Xb4f','spw18','HCN 3-2',265.886434),('Xb47','spw16','HCN 4-3',354.505477),('Xb47','spw18','HCO+ 4-3',356.734223),
 ('Xb47','spw22','CS 7-6',342.882850),('Xb55','spw16','HCN 2-1',177.261111),('Xb55','spw18','HCO+ 2-1',178.375065),
 ('Xb4d','spw25','HCO+ 3-2',267.557626),('Xb4d','spw27','HCN 3-2',265.886434),('Xb4b','spw25','HCO+ 3-2',267.557626),('Xb4b','spw27','HCN 3-2',265.886434)]
def mom(ous,spw,name,f0,vlo=-VW,vhi=VW,tag='core'):
    h=fits.open(f(ous,spw),memmap=True); hd=h[0].header; fr=freqaxis(hd)
    v=(1-fr/1e9/f0)*C; dv=abs(np.median(np.diff(v)))
    sel=np.where((v>=VSYS+vlo)&(v<=VSYS+vhi))[0]
    cube=np.array(h[0].data[0,sel.min():sel.max()+1],float)
    off=np.array(h[0].data[0,0:15],float)  # first 15 channels line-free (all lines at 2780 within band interior)
    sig=np.nanstd(off,axis=0)  # per-pixel per-channel noise from line-free channels
    m0=np.nansum(cube,axis=0)*dv
    n=len(sel); e0=sig*dv*np.sqrt(n)
    return hd,m0,e0,dv,n
if __name__=='__main__':
    fig,axs=plt.subplots(3,4,figsize=(20,15)); axs=axs.ravel()
    for a,(ous,spw,name,f0) in zip(axs,LINES):
        if ous=='Xb4b': continue
        hd,m0,e0,dv,n=mom(ous,spw,name,f0); w=WCS(hd).celestial
        np.savez(f'/workspace/out/m0_{ous}_{spw}.npz',m0=m0,e0=e0)
        pix=abs(hd['CDELT2'])*3600; ny,nx=m0.shape
        ext=[-(nx/2)*pix,(nx/2)*pix,-(ny/2)*pix,(ny/2)*pix]
        a.imshow(m0,origin='lower',cmap='magma'); a.set_title(f'{ous} {spw} {name}')
        # local maxima
        lm=(m0==ndi.maximum_filter(np.nan_to_num(m0),size=max(3,int(3./pix))))&(m0>5*np.nanmedian(e0))
        ys,xs=np.where(lm); pk=sorted([(m0[y,x],x,y) for x,y in zip(xs,ys)],reverse=True)[:6]
        print(ous,spw,name,'dv',round(dv,2),'nch',n,'median e0',round(float(np.nanmedian(e0)),3))
        for val,x,y in pk:
            ra,de=[float(q) for q in w.pixel_to_world_values(x,y)]
            print(f'    peak {val:.2f} Jy/b km/s (S/N {val/np.nanmedian(e0):.0f}) at {ra:.5f} {de:.5f}   dRA={((ra-156.96346)*np.cos(np.deg2rad(-43.9046)))*3600:.2f}" dDec={(de+43.90461)*3600:.2f}"')
            a.plot(x,y,'c+')
    plt.tight_layout(); plt.savefig('/workspace/out/m0_maps.png',dpi=50)
