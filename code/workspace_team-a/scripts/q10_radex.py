# Q10 v2: RADEX (LAMDA hcn.dat, hco+.dat, H2, Tbg 2.73 K, dv=150 km/s) response of HCN/HCO+ 3-2 (Tb ratio) to n, Tk, X=N(HCN)/N(HCO+), N(HCO+) (opacity).
# python3 -I q10_radex.py -> out/q10v2/radex.log, radex_grid.npy
import os,subprocess,itertools,re,json,numpy as np
from concurrent.futures import ProcessPoolExecutor
TK=[20,30,50,80,150.]; LN=np.arange(3,7.01,0.5); LNC=[14.0,15.0,16.0,16.5,17.0]; XS=[0.3,1,3,10,30]; DV=150.
def fl(x):
    x=re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',x)
    try: return float(x)
    except Exception: return float('nan')
def run(a):
    mol,T,ln,lN=a; wd=f'/tmp/rx10_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n60 400\n{T}\n1\nH2\n{10**ln:.4e}\n2.73\n{10**lN:.4e}\n{DV}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); d={}
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0].isdigit(): d[int(p[0])]=(fl(p[7]),fl(p[6]))   # upper J -> (T_R, tau); J=2..4 means 2-1 is upper 2
    return [d.get(2,(np.nan,np.nan)),d.get(3,(np.nan,np.nan)),d.get(4,(np.nan,np.nan))] if False else [d.get(j,(np.nan,np.nan)) for j in (1,2,3,4)]
if __name__=='__main__':
    # radex rows: upper-J label 1..: "1 -- 0", "2 -- 1", "3 -- 2", "4 -- 3"
    LNall=sorted(set(list(LNC)+[x+np.log10(1) for x in LNC]+[lc+np.log10(x) for lc in LNC for x in XS]))
    path='/workspace/out/q10v2/radex_grid.npy'
    keys=[(m,T,ln,lN) for m in ('hcn','hco+') for T in TK for ln in LN for lN in LNall]
    if os.path.exists(path): arr=np.load(path)
    else:
        with ProcessPoolExecutor(4) as ex: res=list(ex.map(run,keys,chunksize=40))
        arr=np.array(res); np.save(path,arr)
    idx={k:i for i,k in enumerate(keys)}
    log=[]
    def P(*a):
        s=' '.join(str(x) for x in a); print(s); log.append(s)
    def get(m,T,ln,lN): return arr[idx[(m,T,ln,lN)]]   # [(TR,tau) for J=1..4 upper index]
    def Rcalc(T,ln,lNo,X):
        lNh=lNo+np.log10(X); h=get('hcn',T,ln,min(LNall,key=lambda z:abs(z-lNh))); o=get('hco+',T,ln,lNo)
        return h[2][0]/o[2][0],h[2][1],o[2][1],h[3][0]/h[2][0],o[3][0]/o[2][0],h[2][0]/h[1][0],o[2][0]/o[1][0]
    P('RADEX Tb-ratio R32=T_R(HCN 3-2)/T_R(HCO+ 3-2) [tauHCN32, tauHCO+32 | 4-3/3-2 HCN,HCO+ | 3-2/2-1 HCN,HCO+]; Tbg=2.73 K, dv=150 km/s, H2 collider. Observed Tb-ratio: N core 0.93-0.97, S core 0.40-0.46, outer N (2-3") ~0.30-0.4; ladders (10" beam / deblended): HCN 4-3/3-2 N 0.39 S 0.27; HCO+ 4-3/3-2 N 0.48 S 0.44; HCN 3-2/2-1 N 0.35 S 0.33; HCO+ 3-2/2-1 N 0.43 S 0.40 (Q1b)')
    for T in (30.,80.):
        for lNo in (15.0,16.0,16.5):
            P(f'\n Tk={T:.0f} K, log N(HCO+)={lNo}:   rows log n; columns X=N(HCN)/N(HCO+) = {XS}')
            for ln in LN:
                row=[]
                for X in XS:
                    R,th,to,_,_,_,_=Rcalc(T,ln,lNo,X); row.append(f'{R:5.2f}')
                R,th,to,h43,o43,h32,o32=Rcalc(T,ln,lNo,3)
                P(f'   log n={ln:.1f}: R32 '+' '.join(row)+f'   | X=3: tau(HCN32)={th:.1f}, tau(HCO+32)={to:.1f}, 4-3/3-2 {h43:.2f}/{o43:.2f}, 3-2/2-1 {h32:.2f}/{o32:.2f}')
    # solutions for R32=1.0+-0.1 vs 0.40+-0.05 and ladder consistency
    P('\n Grid points with R32 = 1.0+-0.1 (N core) / 0.42+-0.05 (S core / N outskirts) that ALSO satisfy HCN 4-3/3-2 and HCO+ 4-3/3-2 within 0.08 of N: (0.39,0.48) resp. S: (0.27,0.44); tau32 given')
    for lab,Rt,eR,(l1,l2) in (('N',1.0,0.1,(0.39,0.48)),('S',0.42,0.05,(0.27,0.44))):
        sol=[]
        for T in TK:
            for ln in LN:
                for lNo in LNC:
                    for X in XS:
                        R,th,to,h43,o43,_,_=Rcalc(T,ln,lNo,X)
                        if abs(R-Rt)<eR and abs(h43-l1)<0.08 and abs(o43-l2)<0.08: sol.append((T,ln,lNo,X,th,to))
        s=np.array(sol)
        if len(s):
            P(f' {lab}: {len(s)} solutions; Tk {sorted(set(s[:,0]))}; log n {s[:,1].min():.1f}-{s[:,1].max():.1f}; log N(HCO+) {s[:,2].min():.1f}-{s[:,2].max():.1f}; X {sorted(set(s[:,3]))}; tau32(HCN) {s[:,4].min():.1f}-{s[:,4].max():.1f}')
        else: P(f' {lab}: none')
    # one-parameter trajectories: start from the S-like state and move ONE parameter to reach R=1
    P('\n One-parameter trajectories (state S-like: Tk=50..80, log n=4, N(HCO+)=1e16, X=3). Which single change gives R32: 0.4 -> 1.0, and the other observables it predicts')
    for T in (30.,80.):
        for lbl,vary in (('X',[(T,4.0,16.0,x) for x in XS]),('n',[(T,ln,16.0,3) for ln in LN]),('N(HCO+) at X=3',[(T,4.0,lN,3) for lN in LNC])):
            P(f'  Tk={T:.0f}: vary {lbl}: '+' '.join(f'[{v[-1] if lbl=="X" else (v[1] if lbl=="n" else v[2])}: R={Rcalc(*v)[0]:.2f} tau={Rcalc(*v)[1]:.0f} 43/32={Rcalc(*v)[3]:.2f}/{Rcalc(*v)[4]:.2f}]' for v in vary))
    open('/workspace/out/q10v2/radex.log','w').write('\n'.join(log))
