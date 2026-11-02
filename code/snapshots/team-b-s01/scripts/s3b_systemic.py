# Systemic velocity from the integrated (r=8") spectra of the strong HCN/HCO+ lines (radio def., LSRK)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
L={'HCN43':(354.50547,'Xb47.spw16.cube.pbcor'),'HCO+43':(356.73424,'Xb47.spw18.cube.pbcor'),'HCN32':(265.88618,'Xb4f.spw18.cube.pbcor'),'HCO+32':(267.55763,'Xb4f.spw16.cube.pbcor'),'HCN21':(177.26111,'Xb55.spw16.cube.pbcor'),'HCO+21':(178.37506,'Xb55.spw18.cube.pbcor')}
for k,(f0,nm) in L.items():
    z=np.load('/workspace/out/spectra/%s.npz'%nm); nu=z['nu']; v=vrad(nu,f0*1e9)
    for key in ['tot_r8','N_pix','S_pix']:
        y=z[key]; m=(abs(v-VSYS)<600)
        off=(abs(v-VSYS)>650)&(abs(v-VSYS)<3000)
        if off.sum()<10: print(k,key,'no baseline'); continue
        yy=y-np.median(y[off]); ymax=yy[m].max()
        sel=m&(yy>0.5*ymax); idx=np.where(sel)[0]
        vm=np.sum(v[m]*yy[m])/np.sum(yy[m])
        print(k,key,'moment1 %.0f  half-max-mid %.0f  width(0.5max) %.0f'%(vm,0.5*(v[idx[0]]+v[idx[-1]]),abs(v[idx[0]]-v[idx[-1]])))
