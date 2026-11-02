# Q1 physical read-out of the HCN/HCO+ ladder: RADEX (LVG, LAMDA hcn.dat / hco+.dat, Tbg 2.73 K, dv=150 km/s FWHM) grid; the same Tkin, n(H2), beam filling for both molecules;
# free column N(HCO+) and abundance ratio A=N(HCN)/N(HCO+). Observables: R_K(J)=HCN/HCO+ in K km/s for J=2-1,3-2,4-3 at the matched 10" beam, N-S midpoint (Xb4d-3-2 variant and Xb4f variant averaged),
# and 6.5" N and S (J=3,4). 3 ratios vs 4 parameters: report the (n,T) region with chi2 min over (N,A) < 4 and the A range there. Not a full Q3 analysis.
import subprocess,os,json,itertools,numpy as np
os.makedirs('/workspace/out/radex_q1',exist_ok=True); os.chdir('/workspace/out/radex_q1')
LO={'hcn':{2:177.2611,3:265.8864,4:354.5055},'hco+':{2:178.3751,3:267.5576,4:356.7342}}
def radex(mol,Tk,nH,N,dV=150.):
    inp=f"/opt/lamda/{mol}.dat\nout_{mol}.txt\n170 370\n{Tk}\n1\nH2\n{nH}\n2.73\n{N}\n{dV}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    I={};tau={}
    for l in open(f'out_{mol}.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--':
            for J,f in LO[mol].items():
                if abs(float(p[4])-f)<2e-3: I[J]=float(p[11]); tau[J]=float(p[7])
    return I,tau
OBS={ # label: {J:(R_K, err)} ; err = stat+5% cal in quadrature + 5% model/beam-matching allowance
 'mid10':{2:(0.588,0.04),3:(0.445,0.04),4:(0.335,0.03)},
 'N6.5':{3:(0.55,0.04),4:(0.43,0.04)}, 'S6.5':{3:(0.40,0.03),4:(0.32,0.03)}}
Ts=[15,20,30,50,80,120,200]; ns=[1e4,1e5,3e5,1e6,3e6,1e7,3e7,1e8]
NO=10**np.arange(12,16.01,0.25); NH=10**np.arange(11,17.01,0.1)   # N(HCO+) and N(HCN) grids; A=NH/NO
res={}
for lab in OBS: res[lab]={}
for T in Ts:
    for n in ns:
        tO={N:radex('hco+',T,n,N)[0] for N in NO}; tH={N:radex('hcn',T,n,N)[0] for N in NH}
        for lab,ob in OBS.items():
            bc=(1e9,None)
            for N in NO:
                Io=tO[N]
                if len(Io)<3: continue
                for M in NH:
                    Ih=tH[M]
                    if len(Ih)<3: continue
                    c=sum(((Ih[J]/Io[J]-r)/e)**2 for J,(r,e) in ob.items())
                    if c<bc[0]: bc=(c,(N,M/N))
            res[lab][(T,n)]=bc
txt=['RADEX LVG, dv=150 km/s; chi2 minimum over (N(HCO+),A) per (Tk,n); allowed = chi2<4 (2 d.o.f. for mid10, 0-1 for the others, so only a consistency screen)']
out={}
for lab,best in res.items():
    ok=[(k,v) for k,v in best.items() if v[0]<4]
    txt.append(f'## {lab}: {len(ok)}/{len(best)} (Tk,n) cells allowed')
    if ok:
        Aall=[v[1][1] for k,v in ok]; Tk=sorted({k[0] for k,v in ok}); nn=sorted({k[1] for k,v in ok})
        txt.append(f'   A=N(HCN)/N(HCO+) range {min(Aall):.2f}-{max(Aall):.2f}; Tk {Tk[0]}-{Tk[-1]} K; n {nn[0]:.0e}-{nn[-1]:.0e}')
        out[lab]=dict(Arange=[min(Aall),max(Aall)],Tk=[Tk[0],Tk[-1]],n=[nn[0],nn[-1]],cells=[(k,v[0],v[1][1],v[1][0]) for k,v in ok])
    for (T,n),(c,(N,A)) in sorted(best.items()):
        if c<4: txt.append(f'   Tk={T:4d} n={n:.0e} chi2={c:.2f} N(HCO+)={N:.0e} A={A:.2f}')
json.dump(out,open('/workspace/out/s33_radex_q1.json','w'),indent=1,default=float)
open('/workspace/out/s33_radex_q1.txt','w').write('\n'.join(txt)); print('\n'.join(txt[:60]))
