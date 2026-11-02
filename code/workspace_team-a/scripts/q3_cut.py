# Q3 step 4: N-S cut through S nucleus -> Soutflow -> Sext (x=-0.2"): line centroid, width, blue-wing fraction and HCN/HCO+ vs distance; Xb4b (0.25") morphology of the Soutflow blue lobe.
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_pos import *
import pandas as pd
from scipy.optimize import curve_fit
c1={sp:cube(sp,fw=1.05) for sp in ('HCN','HCO+')}
cT=Tb_factor(265.886434/(1+VS/C),1.05)/Tb_factor(267.557626/(1+VS/C),1.05)
rows=[]
for y in np.arange(-2.55,-7.6,-0.45):
    r=dict(dDec=y,dist_from_S=-2.55-y)
    for sp,c in c1.items():
        s=spec(c,-0.2,y); v=c['v']; fm=chan_noise(c); e=1.4826*np.median(np.abs(s[fm]-np.median(s[fm])))
        m=(v>2250)&(v<2950); I=s[m].sum()*c['dv']; r['I_'+sp]=I; r['vmean_'+sp]=(s[m]*v[m]).sum()/s[m].sum() if I>0 else np.nan
        r['sd_'+sp]=np.sqrt(np.sum(s[m]*(v[m]-r['vmean_'+sp])**2)/s[m].sum()) if I>0 else np.nan
        for k,(lo,hi) in {'bl':(2300,2575),'co':(2575,2900)}.items():
            mm=(v>=lo)&(v<hi); r[f'{k}_{sp}']=s[mm].sum()*c['dv']; r[f'e{k}_{sp}']=e*c['dv']*np.sqrt(mm.sum())
    r['R_blue']=r['bl_HCN']/r['bl_HCO+']*cT if r['bl_HCN']>3*r['ebl_HCN'] and r['bl_HCO+']>3*r['ebl_HCO+'] else np.nan
    r['R_core']=r['co_HCN']/r['co_HCO+']*cT if r['co_HCN']>3*r['eco_HCN'] and r['co_HCO+']>3*r['eco_HCO+'] else np.nan
    rows.append(r)
T=pd.DataFrame(rows);T.round(3).to_csv('/workspace/out/q3/q3_cut_S_Soutflow.csv',index=False)
pd.set_option('display.width',250);print(T[['dDec','dist_from_S','I_HCN','I_HCO+','vmean_HCN','vmean_HCO+','sd_HCN','sd_HCO+','bl_HCN','bl_HCO+','R_blue','R_core']].round(2).to_string())
# --- Xb4b morphology: blue window 2300-2575 (HCN, HCO+) and core-local window, native beam, crop +-8" around the S nucleus
for ous in ('Xb4b',):
    for sp in ('HCN','HCO+'):
        c=cube(sp,ous,crop=8.); v=c['v']
        for wn,(lo,hi) in {'blue 2300-2575':(2300,2575),'S-lobe 2575-2850':(2575,2850)}.items():
            m0,n=wsum(c,lo,hi); sg=wnoise(c,lo,hi)[0]; sn=m0/sg
            # box around Soutflow +-1.5"
            box=(np.abs(c['dx']+0.2)<1.5)&(np.abs(c['dy']+4.5)<1.5)
            from scipy import ndimage as ndi
            lab,nl=ndi.label((sn>3)&box)
            tot=m0[(sn>3)&box].sum()/ (OM(c['beam'])/c['pix']**2)    # Jy km/s: sum Jy/beam km/s over pixels / pixels per beam
            print(ous,sp,wn,'peak S/N in box %.1f; pixels>3sig %d (%.2f arcsec2); connected regions %d; flux(>3sig px) %.3f Jy km/s; peak %.3f Jy/beam km/s at (%.2f,%.2f)'%(sn[box].max(),((sn>3)&box).sum(),((sn>3)&box).sum()*c['pix']**2,nl,tot,m0[box].max(),c['dx'][box][np.argmax(m0[box])],c['dy'][box][np.argmax(m0[box])]))
            # resolution: smooth to 1.05" and compare with Xb4d value at the same position
            np.savez(f'/workspace/out/q3/xb4b_{sp}_{lo}.npz',m0=m0,sn=sn,dx=c['dx'],dy=c['dy'])
