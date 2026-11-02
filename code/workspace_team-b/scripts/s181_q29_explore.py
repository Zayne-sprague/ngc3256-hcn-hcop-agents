import sys; sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from regions import REGIONS
import json
out={}
for key in CUBES:
    cut,nu,pix,bm,pxy=load(key,half=12)
    print(key,cut.shape,'pix %.3f beam %.2f %.2f %.0f'%(pix,*bm),'nu %.3f-%.3f dnu %.2f MHz'%(nu[0],nu[-1],(nu[1]-nu[0])*1e3))
    Tb=to_tb_conv(cut,nu,pix,bm,0.95,np.mean(nu))
    ny,nx=Tb.shape[1:]
    for rn in ('N','S'):
        x,y=pxy(REGIONS[rn][:2]); ix,iy=int(round(x)),int(round(y))
        sp=Tb[:,iy-1:iy+2,ix-1:ix+2].mean(axis=(1,2))
        np.savetxt('/workspace/out/q29/spec_%s_%s.txt'%(key,rn),np.c_[nu,sp])
        # rest freq of strong peaks assuming vsys
        srt=np.argsort(sp)[::-1][:6]
        print(rn,[(round(nu[i]/(1-V0/C),3),round(float(sp[i]),3)) for i in sorted(srt)])
        print('  robust rms (K):',1.4826*np.median(np.abs(sp-np.median(sp))))
