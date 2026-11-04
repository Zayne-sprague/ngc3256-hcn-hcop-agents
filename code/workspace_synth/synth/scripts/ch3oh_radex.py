# Non-LTE (RADEX, LAMDA e-ch3oh.dat) grid for the matched-beam CH3OH intensities (out/ch3oh_matched.csv): n(H2), Tkin, N(E-CH3OH) at fixed dv=150 km/s, free common beam-filling eta (<=1) -> chi2.
import sys,subprocess,os,tempfile,itertools,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
D=pd.read_csv('/workspace/out/ch3oh_matched.csv')
# observed set (K km/s): Q J=1-3 sum, J4, J5, J6, J7, J8 (Xb55), 2(0)-1(-1) and 5(-2,4)-4(-1,3) (Xb4f matched). 10% calibration added.
OBS=[('Q123',(165.050175,165.061130,165.099240),D[(D.line=='Q J=1-3')].iloc[0]),('J4',(165.190475,),D[D.line=='Q J=4'].iloc[0]),('J5',(165.369341,),D[D.line=='Q J=5'].iloc[0]),
     ('J6',(165.678649,),D[D.line=='Q J=6'].iloc[0]),('J7',(166.169098,),D[D.line=='Q J=7'].iloc[0]),('J8',(166.898566,),D[D.line=='Q J=8'].iloc[0]),
     ('2_0-1_-1',(254.015377,),D[(D.line=='2(0)-1(-1) E')&(D.ous=='Xb4f')].iloc[0]),('5_-2-4_-1',(266.838148,),D[(D.line=='5(-2,4)-4(-1,3) E')&(D.ous=='Xb4f')].iloc[0])]
def radex(T,n,N,dv=150.,wd='/tmp/rx'):
    os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/e-ch3oh.dat\n{wd}/o.txt\n150 300\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n{dv}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True)
    rows=[]
    for l in open(f'{wd}/o.txt').read().split('\n')[11:]:
        p=l.split()
        if len(p)>=13 and p[1]=='--':
            rows.append((float(p[4]),float(p[8]),float(p[-3])))   # freq, tau, K km/s ; careful columns: name,'--',name,Eup,freq,wavel,Tex,tau,TR,popu,popl,flux,flux
    return np.array(rows)
def model(T,n,N):
    r=radex(T,n,N)   # cols: freq, tau, flux  -> recompute properly below
    return r
def parse(T,n,N):
    wd='/tmp/rx'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/e-ch3oh.dat\n{wd}/o.txt\n150 300\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n150\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True)
    out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit():
            try: out.append((float(p[3]),float(p[5]),float(p[6]),float(p[10])))  # freq, Tex, tau, K km/s
            except: pass
    return np.array(out)
if __name__=='__main__':
    fo=np.array([[o[2].I_Kkms,np.hypot(o[2].eI_Kkms,0.1*max(o[2].I_Kkms,0))] for o in OBS])
    res=[]
    for T in [10,15,20,30,50,80,120,200]:
        for ln in np.arange(2,8.01,0.5):
            for lN in np.arange(12,16.01,0.5):
                m=parse(T,10**ln,10**lN)
                if len(m)==0: continue
                pr=[]
                for name,frs,_ in OBS:
                    s=0
                    for f0 in frs:
                        k=np.argmin(abs(m[:,0]-f0)); s+= m[k,3] if abs(m[k,0]-f0)<0.01 else np.nan
                    pr.append(s)
                pr=np.array(pr)
                if not np.all(np.isfinite(pr)): continue
                y,e=fo[:,0],fo[:,1]
                eta=np.clip(np.sum(pr*y/e**2)/max(np.sum(pr*pr/e**2),1e-30),0,1e3)
                chi2=np.sum(((y-eta*pr)/e)**2)
                res.append(dict(Tkin=T,log_n=ln,log_N=lN,eta=eta,chi2=chi2,**{f'p_{o[0]}':p for o,p in zip(OBS,pr)}))
    R=pd.DataFrame(res); R.to_csv('/workspace/out/ch3oh_radex_grid.csv',index=False)
    b=R.sort_values('chi2').head(15); pd.set_option('display.width',250); print(b.round(3).to_string())
    print('observed:',dict(zip([o[0] for o in OBS],fo[:,0].round(3))))
    ok=R[R.chi2<R.chi2.min()+11.8]   # ~3-parameter 68%? use delta chi2<3.5 (1 sig,3 par=3.53) and 7.8(2sig)
    for lab,dc in (('1sig(3par)',3.53),('2sig',8.0)):
        s=R[R.chi2<R.chi2.min()+dc]; print(lab,'Tkin',sorted(s.Tkin.unique()),'log n',s.log_n.min(),s.log_n.max(),'logN',s.log_N.min(),s.log_N.max(),'eta',s.eta.min().round(3),s.eta.max().round(3))
