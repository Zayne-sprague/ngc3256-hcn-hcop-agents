"""Q6: beam-filling implied by the 0.25" peak brightness temperatures: RADEX (LVG, LAMDA, H2, Tbg 2.73 K) T_R(3-2) vs (Tk, n, N/dv).
eta_min = Tb_obs / T_R(model); thermalised thick limit J(Tk)-J(Tbg). python3 -I s153_q6_radex_filling.py"""
import subprocess,os,json,numpy as np
os.makedirs('/workspace/out/q6/radex',exist_ok=True); os.chdir('/workspace/out/q6/radex')
NU={'hco+':267.5576,'hcn':265.8864}
def radex(mol,Tk,n,N,dV):
    inp=f"/opt/lamda/{mol}.dat\nout_{mol}.txt\n170 370\n{Tk}\n1\nH2\n{n}\n2.73\n{N}\n{dV}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    for l in open(f'out_{mol}.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--' and abs(float(p[4])-NU[mol]/1.0)<0.3 and int(p[0])==3:
            return float(p[8]),float(p[7]),float(p[6])
    return np.nan,np.nan,np.nan
def J(T,nu): x=4.799e-11*nu*1e9; return x/np.expm1(x/T)
obs={'HCO+ S':('hco+',2.81,100.),'HCO+ N':('hco+',2.75,130.),'HCN N':('hcn',2.46,160.),'HCN S':('hcn',1.18,100.),'HCO+ E1':('hco+',2.05,60.),'HCO+ E2':('hco+',2.04,55.)}
rows=[]
for name,(mol,tb,dv) in obs.items():
    for Tk in (30,60,120):
        Jmax=J(Tk,NU[mol])-J(2.73,NU[mol])
        for n in (1e4,1e5,1e6,1e7):
            for N in (1e13,1e14,1e15,1e16,1e17):
                tr,tau,tex=radex(mol,Tk,n,N,dv)
                rows.append(dict(line=name,Tk=Tk,n=n,N=N,dv=dv,TR=tr,tau=tau,Tex=tex,eta=tb/tr if tr>0 else np.nan,Tb_obs=tb,JmaxLTE=Jmax,eta_min_LTE=tb/Jmax))
import pandas as pd
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q6/s153_radex_filling.csv',index=False)
pd.set_option('display.width',200,'display.max_rows',500)
print(df.groupby(['line','Tk']).agg(Tb_obs=('Tb_obs','first'),Jmax=('JmaxLTE','first'),eta_min_LTE=('eta_min_LTE','first'),TRmax=('TR',lambda x: np.nanmax(x))).round(2).to_string())
sub=df[(df.line=='HCO+ S')&(df.Tk==60)]; print(sub.pivot(index='N',columns='n',values='TR').round(2).to_string()); print(sub.pivot(index='N',columns='n',values='tau').round(2).to_string())
# eta range over grid cells that are allowed (TR>=Tb_obs) for the fitted regime n 1e4-1e6, Tk 30-120
for name in obs:
    s=df[(df.line==name)&(df.n<=1e6)&(df.TR>=df.Tb_obs)]
    print(name,'eta in cells with TR>=Tb (n<=1e6):','min %.3f  median %.3f'%(s.eta.min(),s.eta.median()),'N cells',len(s))
