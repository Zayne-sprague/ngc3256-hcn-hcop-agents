# Matched-filter (template = same-region HCN/HCO+ profile of the same array/beam) amplitudes / 3-sigma limits for weaker and expected lines
# + blend models: CS 7-6 / H2CS / 29SiO; H13CN 4-3 on the CO 3-2 blue wing.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
from scipy.optimize import nnls, curve_fit
Z={}
def spec(ous,spw,reg): 
    k=(ous,spw)
    if k not in Z: Z[k]=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz')
    return Z[k]['nu'],Z[k][reg+'_ap']
def vel(nu,f0): return C*(1-nu/(f0*1e9))
def template(ous,spw,f0,reg,dvgrid):
    nu,y=spec(ous,spw,reg); v=vel(nu,f0)-VSYS; o=np.argsort(v); v=v[o]; y=y[o]
    base=np.median(y[(abs(v)>500)&(abs(v)<900)]); y=y-base
    t=np.interp(dvgrid,v,np.convolve(y,np.ones(3)/3,'same')); t[abs(dvgrid)>350]=0; t=np.clip(t,0,None)
    return t/(t.sum()*abs(dvgrid[1]-dvgrid[0]))     # unit integral [1/(km/s)]
# (ous, spw, label, rest GHz, template line (spw,f0), flag)
W=[('Xb47','spw20','HC15N 4-3',344.20016,('spw16',354.50547),''),
   ('Xb47','spw20','SO 8_8-7_7',344.31061,('spw16',354.50547),''),
   ('Xb47','spw22','H2CS 10_0,10-9_0,9',342.94644,('spw16',354.50547),'blended with CS 7-6 / 29SiO 8-7'),
   ('Xb47','spw16','HCN v2=1f? c-HCCCH 4(3,2)-3(0,3) 354.240',354.2402,('spw16',354.50547),'on HCN blue wing'),
   ('Xb47','spw18','HCN v2=1f 4-3 (l=1f)',356.2559,('spw16',354.50547),'-400 km/s from HCO+ 4-3'),
   ('Xb4f','spw22','SO 5_6-4_5',251.82577,('spw18',265.88618),''),
   ('Xb4f','spw20','NS 11/2-9/2 Om=1/2 f (253.968)',253.96793,('spw18',265.88618),'blended: SO2 253.957, c-HCCCH 254.013'),
   ('Xb4f','spw20','NS 11/2-9/2 Om=1/2 e (253.571)',253.57052,('spw18',265.88618),'partner of 253.968; should be comparable'),
   ('Xb4f','spw20','SiS 14-13',254.10321,('spw18',265.88618),''),
   ('Xb4f','spw20','HC3N 28-27',254.69970,('spw18',265.88618),''),
   ('Xb4f','spw16','HCN v2=1f 3-2 (l=1f)',267.1993,('spw18',265.88618),'-400 km/s from HCO+ 3-2'),
   ('Xb4f','spw20','CH3OH?',254.0,('spw18',265.88618),'')]
W+=[('Xb55','spw18','HOC+ 2-1',178.97214,('spw18',178.37506),'-50 km/s from H2CS/other? check; ~600 km/s from HCO+ 2-1 (no overlap)'),
    ('Xb55','spw16','H41beta RRL',177.72284,('spw16',177.26111),''),
    ('Xb55','spw22','CH3CN 9-8 K=0,1 (165.568)',165.56825,('spw16',177.26111),'')]
W=[w for w in W if w[2]!='CH3OH?']
# same Band-6 lines in the 1" (Xb4d) and 0.25" (Xb4b) data: spw mapping Xb4f->Xb4d/Xb4b: 16->25, 18->27, 20->29, 22->31
MAP={'spw16':'spw25','spw18':'spw27','spw20':'spw29','spw22':'spw31'}
W+= [(o2,MAP[sp],lab,f0,(MAP[ts],tf0),fl) for o2 in ('Xb4d','Xb4b') for (o,sp,lab,f0,(ts,tf0),fl) in W if o=='Xb4f']
dvg=np.arange(-600,600.1,3.3)
rows=[]
for ous,spw,lab,f0,(tspw,tf0),flag in W:
    for reg in ('N','S'):
        nu,y=spec(ous,spw,reg); v=vel(nu,f0)-VSYS; o=np.argsort(v); v,y=v[o],y[o]
        ok=np.isfinite(y)&(abs(v)<900)
        if (abs(v[ok])>0).sum()<60 or v[ok].max()<250 or v[ok].min()>-250: continue
        t=template(ous,tspw,tf0,reg,np.arange(-1000,1000,abs(np.median(np.diff(v)))) )  # on a grid
        dv=abs(np.median(np.diff(v))); tg=np.arange(-1000,1000,dv)
        T=np.interp(v,tg,t)
        # baseline: line-free, |v| in 450..900 (clipped) -> subtract median
        fl=ok&(abs(v)>450)
        yb=y[fl]
        for _ in range(4):
            med=np.median(yb); s=1.4826*np.median(np.abs(yb-med)); yb=yb[abs(yb-med)<3.5*s]
        base=np.median(yb); rms=np.std(yb)
        m=ok&(abs(v)<=600)
        yy=y[m]-base; TT=T[m]
        Wf=np.sum(TT*yy)/np.sum(TT**2); sW=rms/np.sqrt(np.sum(TT**2))
        # rho correction
        d=yb-yb.mean(); rho=np.sum(d[1:]*d[:-1])/np.sum(d*d); sW*=np.sqrt((1+max(rho,0))/(1-max(rho,0)))
        status='detected' if Wf/sW>=5 else ('marginal (3-5 sigma)' if Wf/sW>=3 else 'not detected')
        rows.append(dict(ous=ous,line=lab,rest_GHz=f0,region=reg,W_Jykms=Wf,sigW=sW,snr=Wf/sW,UL3_Jykms=3*sW,rms_chan_mJy=rms*1e3,status=status,flag=flag))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s7c_weak_lines.csv',index=False)
pd.set_option('display.width',250,'display.max_colwidth',60); print(df.round(3).to_string())
