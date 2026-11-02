# Line measurements for every (cube, region, spectrum-type): peak, 1- and 2-Gaussian fits, integrated intensity, per-channel rms, T_b
import sys,warnings; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import pandas as pd
from scipy.optimize import curve_fit
LINES={ # (cube tag, line label, rest GHz, half window for W [km/s])
 ('Xb55','spw16'):[('HCN 2-1',177.26111,350)],('Xb55','spw18'):[('HCO+ 2-1',178.37506,350)],
 ('Xb55','spw22'):[('CH3OH Jk-1-Jk0 E blend (2,2 @165.061)',165.06113,250)],
 ('Xb47','spw16'):[('HCN 4-3',354.50547,350)],('Xb47','spw18'):[('HCO+ 4-3',356.73424,350)],('Xb47','spw22'):[('CS 7-6 (+H2CS 10-9, 29SiO 8-7 blend)',342.88285,350)],
 ('Xb4f','spw18'):[('HCN 3-2',265.88618,350)],('Xb4f','spw16'):[('HCO+ 3-2',267.55763,350)],
 ('Xb4d','spw27'):[('HCN 3-2',265.88618,350)],('Xb4d','spw25'):[('HCO+ 3-2',267.55763,350)],
 ('Xb4b','spw27'):[('HCN 3-2',265.88618,350)],('Xb4b','spw25'):[('HCO+ 3-2',267.55763,350)],
}
def g(v,a,v0,fw): return a*np.exp(-4*np.log(2)*((v-v0)/fw)**2)
def m1(v,a,v0,fw,c): return g(v,a,v0,fw)+c
def m2(v,a1,v1,w1,a2,v2,w2,c): return g(v,a1,v1,w1)+g(v,a2,v2,w2)+c
def measure(nu,y,f0,half,bmaj,bmin):
    v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
    win=abs(v-VSYS)<=900; v,y=v[win],y[win]; ok=np.isfinite(y); v,y=v[ok],y[ok]
    if len(v)<50: return None
    fl=abs(v-VSYS)>450
    if fl.sum()<30: return None
    yb=y[fl]
    for _ in range(5):
        med=np.median(yb); s=1.4826*np.median(np.abs(yb-med)); yb=yb[abs(yb-med)<3.5*s]
    base=np.median(yb); rms=np.std(yb)
    d=yb-np.mean(yb); rho=np.sum(d[1:]*d[:-1])/np.sum(d*d)    # lag-1 autocorrelation of line-free noise (channel-ordered, approx)
    yy=y-base; core=abs(v-VSYS)<=half; n=core.sum()
    W=np.sum(yy[core])*dv
    # empirical sigma of n-channel box sum in line-free channels
    ff=np.where(fl)[0]; sums=[yy[i:i+n].sum() for i in range(0,len(yy)-n) if fl[i:i+n].all()]
    sigW=(np.std(sums) if len(sums)>4 else rms*np.sqrt(n))*dv
    sm=np.convolve(yy,np.ones(3)/3,'same'); ip=np.argmax(np.where(core,sm,-np.inf))
    out=dict(rms_chan=rms,rho1=rho,dv=dv,peak_raw=yy[core].max(),peak_sm3=sm[ip],v_peak=v[ip],W=W,sigW=sigW,snr_W=W/sigW,n_core=n)
    f=np.sqrt((1+max(rho,0))/(1-max(rho,0)))   # error inflation for correlated noise (approximate)
    # 1-Gaussian
    try:
        p0=[max(sm[ip],3*rms),v[ip] if abs(v[ip]-VSYS)<150 else VSYS,200,0]
        p,cv=curve_fit(m1,v,y-base,p0=p0,sigma=np.full(len(v),rms),absolute_sigma=True,maxfev=20000,bounds=([0,VSYS-300,15,-5*rms],[np.inf,VSYS+300,800,5*rms]))
        e=np.sqrt(np.diag(cv))*f; res=y-base-m1(v,*p); chi=np.sum((res/rms)**2); k=4
        out.update(g1_A=p[0],g1_A_e=e[0],g1_v0=p[1],g1_v0_e=e[1],g1_fwhm=p[2],g1_fwhm_e=e[2],g1_chi2r=chi/(len(v)-k),g1_I=p[0]*p[2]*1.0645,g1_bic=chi+k*np.log(len(v)))
    except Exception as ex: pass
    # 2-Gaussian (double-horned profiles): init from two halves
    try:
        best=None
        for dvv in (60,100,140):
            p0=[max(sm[ip],3*rms)*0.8,VSYS-dvv,120,max(sm[ip],3*rms)*0.8,VSYS+dvv,120,0]
            try:
                p,cv=curve_fit(m2,v,y-base,p0=p0,sigma=np.full(len(v),rms),absolute_sigma=True,maxfev=30000,bounds=([0,VSYS-400,10,0,VSYS-400,10,-5*rms],[np.inf,VSYS+400,600,np.inf,VSYS+400,600,5*rms]))
                chi=np.sum(((y-base-m2(v,*p))/rms)**2)
                if best is None or chi<best[0]: best=(chi,p,cv)
            except Exception: pass
        chi,p,cv=best; e=np.sqrt(np.diag(cv))*f; k=7
        o=np.argsort([p[1],p[4]]); P=[p[0:3],p[3:6]]; E=[e[0:3],e[3:6]]
        out.update(g2_Ablue=P[o[0]][0],g2_vblue=P[o[0]][1],g2_wblue=P[o[0]][2],g2_Ared=P[o[1]][0],g2_vred=P[o[1]][1],g2_wred=P[o[1]][2],
                   g2_chi2r=chi/(len(v)-k),g2_bic=chi+k*np.log(len(v)),g2_Ablue_e=E[o[0]][0],g2_vblue_e=E[o[0]][1],g2_Ared_e=E[o[1]][0],g2_vred_e=E[o[1]][1])
    except Exception: pass
    return out
if __name__=='__main__':
    rows=[]
    for (ous,spw),lines in LINES.items():
        z=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz'); nu=z['nu']; bmaj,bmin=float(z['bmaj']),float(z['bmin'])
        for lab,f0,half in lines:
            for reg in REGIONS:
                if ous in ('Xb4f','Xb47','Xb55') and reg not in ('N','S'): continue   # beam >> region separation
                for st in ('pix','ap'):
                    key=f'{reg}_{st}'
                    if key not in z: continue
                    o=measure(nu,z[key],f0,half,bmaj,bmin)
                    if o is None: continue
                    nu_obs=f0*(1-VSYS/C)
                    o['Tb_peak_K']=1.222e3*(o['peak_sm3']*1e3)/(nu_obs**2*bmaj*bmin)   # Rayleigh-Jeans: T[K]=1.222e3 S[mJy/beam]/(nu[GHz]^2 bmaj*bmin[arcsec^2]), observed-frame nu
                    o['W_Kkms']=1.222e3*(o['W']*1e3)/(nu_obs**2*bmaj*bmin)
                    ra,dec=z[reg+'_xy']; o.update(ous=ous,spw=spw,line=lab,rest_GHz=f0,region=reg,spec=st,beam=f'{bmaj:.2f}x{bmin:.2f}',ra=float(ra),dec=float(dec))
                    rows.append(o)
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/s7_line_fits.csv',index=False)
    pd.set_option('display.width',280,'display.max_columns',40,'display.max_rows',600)
    cols=['ous','line','region','spec','rms_chan','peak_sm3','W','sigW','snr_W','g1_v0','g1_fwhm','g1_A','g1_chi2r','g2_chi2r','Tb_peak_K']
    print(df[cols].round(3).to_string())
