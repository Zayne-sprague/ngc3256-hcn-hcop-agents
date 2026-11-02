# Blend models (beam apertures, Xb47): (1) CS 7-6 + H2CS 10-9 + 29SiO 8-7 with the same-region HCO+/HCN 4-3 profile as template (NNLS + bootstrap)
# (2) H13CN 4-3 sitting on the blue wing of CO 3-2 (band edge): template*A + low-order polynomial
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
from scipy.optimize import nnls
rng=np.random.default_rng(1)
z={k:np.load(f'/workspace/out/regspec/Xb47.{k}.cube.pbcor.npz') for k in ('spw16','spw18','spw20','spw22')}
def vel(nu,f0): return C*(1-nu/(f0*1e9))
def tmpl(spw,f0,reg,vg):
    nu=z[spw]['nu']; y=z[spw][reg+'_ap']; v=vel(nu,f0)-VSYS; o=np.argsort(v); v,y=v[o],y[o]
    y=y-np.median(y[(abs(v)>450)&(abs(v)<900)]); y=np.convolve(y,np.ones(3)/3,'same')
    t=np.interp(vg,v,y,left=0,right=0); t[abs(vg)>350]=0; t=np.clip(t,0,None); return t/(t.sum()*abs(vg[1]-vg[0]))
out=[]
# ---- (1) CS complex
L=[('CS 7-6',342.88285),('H2CS 10_0,10-9_0,9',342.94644),('29SiO 8-7',342.98066)]
for reg in ('N','S'):
    nu=z['spw22']['nu']; y=z['spw22'][reg+'_ap']; dv=abs(np.median(np.diff(nu)))/nu.mean()*C
    f_cs=342.88285*(1-VSYS/C); sel=abs(C*(1-nu/(f_cs*1e9)+0)-0)<2000
    vcs=vel(nu,342.88285)-VSYS; m=abs(vcs)<1000; v=vcs[m]; yy=y[m]
    fl=abs(v)>550; rms=np.std(yy[fl]); 
    for tname,(tspw,tf0) in {'HCO+4-3':('spw18',356.73424),'HCN4-3':('spw16',354.50547)}.items():
        vg=np.arange(-1500,1500,0.5); 
        def design(nl):
            cols=[]
            for nm,f in L[:nl]:
                sh=vel(nu[m],f)-vel(nu[m],342.88285)    # relative shift in velocity
                t=tmpl(tspw,tf0,reg,vg); cols.append(np.interp(v-(-sh)*0+ (vel(nu[m],f)-VSYS),vg,t,left=0,right=0)*dv*0+np.interp(vel(nu[m],f)-VSYS,vg,t,left=0,right=0)*1.0)
            cols.append(np.ones_like(v)); return np.array(cols).T
        for nl in (1,3):
            A=design(nl)[:, :]; 
            # constant allowed to be +-: handle by NNLS on [A, +1, -1]
            Af=np.hstack([A,-A[:,-1:]]); x,rn=nnls(Af/rms,yy/rms); 
            res=yy-Af@x; chi2=np.sum((res/rms)**2); dof=len(v)-(nl+1)
            bs=[]
            for _ in range(300):
                yb=Af@x+rng.permutation(res); xb,_=nnls(Af/rms,yb/rms); bs.append(xb[:nl])
            bs=np.array(bs)
            for i in range(nl): out.append(dict(model=f'{nl}-line',template=tname,region=reg,line=L[i][0],W_Jykms=x[i],err_boot=bs[:,i].std(),chi2r=chi2/dof))
# ---- (2) H13CN 4-3 on the blue wing of CO 3-2 (band edge at rest 345.60 GHz): the red half of H13CN overlaps the steeply rising CO wing, so use only the
# LOW-FREQUENCY (velocity-redshifted) half of the H13CN profile (v_H13CN in [+30,+300] km/s), the CO 3-2 line lies on the high-frequency side, so the CO wing is small there; any CO wing there is >=0, so the matched-filter amplitude is an UPPER bound on H13CN.
for reg in ('N','S'):
    nu=z['spw20']['nu']; y=z['spw20'][reg+'_ap']; fr=rest_of(nu)/1e9; v=vel(nu,345.33977)-VSYS
    vg=np.arange(-1500,1500,0.5); t=np.interp(v,vg,tmpl('spw16',354.50547,reg,vg),left=0,right=0)
    fl=(fr>343.9)&(fr<345.05); base=np.median(y[fl]); rms=np.std(y[fl])
    d=y[fl]-y[fl].mean(); rho=np.sum(d[1:]*d[:-1])/np.sum(d*d)
    for lo,hi,tag in [(30,300,'low-freq half'),(120,300,'far low-freq part')]:
        m=(v>=lo)&(v<=hi)&np.isfinite(y); W=np.sum(t[m]*(y[m]-base))/np.sum(t[m]**2); sW=rms/np.sqrt(np.sum(t[m]**2))*np.sqrt((1+max(rho,0))/(1-max(rho,0)))
        out.append(dict(model='H13CN matched-filter',template='HCN4-3',region=reg,line=f'H13CN 4-3 {tag} (upper bound; CO-wing contaminated)',W_Jykms=W,err_boot=sW,chi2r=np.nan))
df=pd.DataFrame(out); df.to_csv('/workspace/out/s7d_blends.csv',index=False)
pd.set_option('display.width',200); print(df.round(3).to_string())
