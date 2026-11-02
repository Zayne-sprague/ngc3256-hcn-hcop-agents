# Measure HCO+ 3-2, HCN 3-2 and CH3OH 2(0)-1(-1)E (Xb4d, 0.93", pbcor) at the CH3OH/HCN-only peak W of S ("Wclump", from regions2.py) and at the Bridge, with the same estimators.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
from ch3oh_matched import fit_one
from ch3oh import toK
RA0,DE0=156.96346,-43.90461
POS={'Wclump':(-1.80,-2.40),'Bridge':(0.30,0.15),'S':(-0.15,-2.55),'N':(0.0,2.55)}
rows=[]
for nm,(dx,dy) in POS.items():
    ra=RA0+dx/3600/np.cos(np.deg2rad(DE0)); de=DE0+dy/3600
    for spw,line,f0 in (('spw25','HCO+ 3-2',267.557626),('spw27','HCN 3-2',265.886434),('spw29','CH3OH 2(0)-1(-1)E',254.015377)):
        h=fits.open(f('Xb4d',spw),memmap=True); hd=h[0].header; w=WCS(hd).celestial; x,y=[int(round(float(a))) for a in w.world_to_pixel_values(ra,de)]
        # peak pixel of the +-350 map within 0.45" 
        cube=h[0].data[0]; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v)))
        m=(np.abs(v-2780)<350); sub=np.array(cube[m][:,y-3:y+4,x-3:x+4],float); mm=np.nansum(sub,axis=0); j,i=np.unravel_index(np.nanargmax(mm),mm.shape); s=np.array(cube[:,y-3+j,x-3+i],float)
        bm=(hd['BMAJ']*3600,hd['BMIN']*3600); yK=toK(s,fr,bm)
        # velocity centroid / width from HCO+/HCN: use the +-350 window moment; line fit for CH3OH with v0=2780 (Bridge) / 2760
        rel=v-2780; bl=(np.abs(rel)>500)&(np.abs(rel)<1500)&np.isfinite(s); p=np.polyfit(rel[bl],s[bl],1); sc=s-np.polyval(p,rel); w_=np.abs(rel)<=350
        I=np.nansum(sc[w_])*dv; sg=np.std(sc[bl]) ; vc=np.sum(v[w_]*sc[w_])/np.sum(sc[w_])
        rec=dict(pos=nm,line=line,pix_dRA=float((WCS(hd).celestial.pixel_to_world_values(x-3+i,y-3+j)[0]-RA0)*np.cos(np.deg2rad(DE0))*3600),I_Jykms_pm350=I,e_stat=sg*dv*np.sqrt(w_.sum()),v_mom1=vc)
        if 'CH3OH' in line:
            A,eA,sg2=fit_one(fr,yK,(254.015377,),v0=vc if False else 2780.,fw=130.); k2j=bm[0]*bm[1]*(254.015377*(1-2780/C))**2/1.222e6
            rec.update(I_fit_Jykms=A*k2j,eI_fit=eA*k2j,snr_fit=A/eA)
        rows.append(rec)
df=pd.DataFrame(rows); df.round(3).to_csv('/workspace/out/wclump.csv',index=False); print(df.round(3).to_string())
