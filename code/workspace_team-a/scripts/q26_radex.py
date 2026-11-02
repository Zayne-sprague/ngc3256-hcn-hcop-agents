# Q26: expected R32/R21 (HCN/HCO+ ratio of Tb, same beam => same as flux-ratio-of-ratios) from RADEX, X=N(HCN)/N(HCO+) tuned so R21=0.54*(nuHCO/nuHCN)^2 inverse.
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts'); import q25_radex as Q, os,subprocess
def run(mol,T,n,N,dv=150.):
    wd='/tmp/rx26'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n50 400\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n{dv}\n0\n'
    subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): out.append((int(p[1]),Q.fl(p[7])))
    return [o[1] for o in sorted(out)[:3]]   # TR for 1-0,2-1,3-2
f=(89.188525/88.631602)**2  # flux ratio = TR ratio * (nuHCN/nuHCO+)^2 for same J roughly (nu ratios ~1.006 per J; cancels in ratio of ratios)
print('T n  NHCO  X  R21  R32  R32/R21')
for T in (20,50,150):
  for ln in (4,5,6):
    for lN in (14.5,15.5,16.5):
        o=run('hco+',T,10**ln,10**lN); best=None
        for X in np.logspace(-1.5,1,26):
            h=run('hcn',T,10**ln,X*10**lN)
            r21=h[1]/o[1]; r32=h[2]/o[2]
            if best is None or abs(np.log(r21/0.54))<abs(np.log(best[1]/0.54)): best=(X,r21,r32)
        if abs(np.log(best[1]/0.54))<0.1: print(T,ln,lN,f'X={best[0]:.2f} R21={best[1]:.2f} R32={best[2]:.2f} R32/R21={best[2]/best[1]:.2f}')
