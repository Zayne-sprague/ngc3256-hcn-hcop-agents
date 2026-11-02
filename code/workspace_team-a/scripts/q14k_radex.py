# Q14k (fix 5): non-LTE (RADEX, LAMDA e-ch3oh.dat / ch3oh_a.dat, p-H2 collider (Rabli & Flower), Tbg 2.73 K, FWHM 150 km/s) brackets for the CH3OH conventions.
# For each (Tkin, n) and N(E)=N(A)=N_tot/2 (equal ortho/para): predicted K km/s of the 254.015 E line, the Q-branch J=1-3 blend (165.05-165.10 E), the 96.74 2k-1k group (E+A) and the 241.79 5k-4k group (E+A).
# Tests: (1) R96/254 and R241/96 vs LTE; (2) the observed 254/blend ratio at S (0.31, 0.22-0.37, q14b) and Harada's S 241.79/96.74 = 0.70+-0.1 -> allowed (Tkin,n); (3) N_tot needed to give I(254)=I_obs per solution, vs LTE at Trot=11 K.
# Output out/q14e/radex_ch3oh_grid.csv, radex_ch3oh_solutions.csv. python3 -I q14k_radex.py
import sys,subprocess,os,itertools; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
DV=150.
def run(a):
    mol,T,n,N=a; wd=f'/tmp/rq14_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n80 300\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit():
            try: out.append((float(p[3]),float(p[5]),float(p[6]),float(p[10])))   # freq, Tex, tau, K km/s
            except Exception: pass
    return np.array(out)
def grp(m,sel):
    k=(m[:,0]>sel[0])&(m[:,0]<sel[1]); return m[k,3].sum() if k.any() else np.nan
def one(a):
    T,n,N=a; me=run(('e-ch3oh',T,n,N/2)); ma=run(('ch3oh_a',T,n,N/2))
    if len(me)==0 or len(ma)==0: return None
    i254=grp(me,(254.01,254.02)); bl=grp(me,(165.04,165.11)); g96=grp(me,(96.73,96.76))+grp(ma,(96.73,96.76)); g241=grp(me,(241.69,241.93))+grp(ma,(241.69,241.93))
    k=np.argmin(abs(me[:,0]-254.015377)); 
    return dict(Tkin=T,log_n=np.log10(n),log_Ntot=np.log10(N),I254=i254,Iblend=bl,I96=g96,I241=g241,tau254=me[k,2],Tex254=me[k,1])
if __name__=='__main__':
    jobs=list(itertools.product([20.,30.,50.,100.],10**np.arange(3.5,7.01,0.5),10**np.arange(12.5,15.51,0.25)))
    with ProcessPoolExecutor(12) as ex: r=list(ex.map(one,jobs,chunksize=10))
    G=pd.DataFrame([x for x in r if x]); G['R254_blend']=G.I254/G.Iblend; G['R96_254']=G.I96/G.I254; G['R241_96']=G.I241/G.I96
    G.to_csv('/workspace/out/q14e/radex_ch3oh_grid.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_rows',300)
    # LTE reference (same N): table is in q14f lte_ratios.csv (CMB-contrast on)
    H=pd.read_csv('/workspace/out/q14f/lte_ratios.csv'); H=H[H.cmb_contrast]
    print('LTE R96/254, R241/96:',H[H['T'].isin([8,10,15,20,30,50])][['T','R96_254','R241_96']].round(3).values.tolist())
    # for each (Tkin,n): interpolate in N to give I254=0.57 (S, Xb4f 6.3x3.8 deblended) and report predicted ratios there
    rows=[]
    for (T,ln),g in G.groupby(['Tkin','log_n']):
        g=g.sort_values('log_Ntot')
        for Iobs in (0.38,0.57):
            if g.I254.max()<Iobs or g.I254.min()>Iobs: continue
            lN=np.interp(np.log(Iobs),np.log(np.clip(g.I254.values,1e-9,None)),g.log_Ntot.values)
            f=lambda c: np.interp(lN,g.log_Ntot.values,g[c].values)
            rows.append(dict(Tkin=T,log_n=ln,I254_obs=Iobs,log_Ntot=lN,R254_blend=f('R254_blend'),R96_254=f('R96_254'),R241_96=f('R241_96'),tau254=f('tau254'),Tex254=f('Tex254')))
    S=pd.DataFrame(rows); S.to_csv('/workspace/out/q14e/radex_ch3oh_solutions.csv',index=False)
    print(S[S.I254_obs==0.57].round(3).to_string())
