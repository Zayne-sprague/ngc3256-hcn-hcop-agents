# Q35: CS J=1-0..7-6 RADEX LVG grid, same conventions as s62_radex_grid.py (dv 150 km/s, Tbg 2.73). Output out/q35/grid_cs.npz
# Q3 round 4: RADEX LVG grid extended to Tkin 500 K and including J=1-0 (needed for the 1-0 test): HCN, HCO+ J=1-0..4-3.
# Grid: Tk 10..500 (13), log n(H2) 3..8 (0.25), log N 11.5..19 (0.25, total column in cm-2 for dv=150 km/s), Tbg 2.73 K, LAMDA /opt/lamda/{hcn,hco+}.dat
# Output out/q35/grid_cs.npz: arrays [mol: T,n,N,J(1..4),(TR,Tex,tau)]
# usage: python3 -I s62_radex_grid.py
import subprocess,os,sys,tempfile,numpy as np
from multiprocessing import Pool
LO={'cs':{1:48.9910,2:97.9810,3:146.9690,4:195.9542,5:244.9356,6:293.9122,7:342.8829}}
TK=np.array([10,14,20,28,40,56,80,112,160,224,300.,400.,500.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(11.5,19.01,0.25); DV=150.
def one(args):
    mol,Tk,ln=args; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d)
    out=np.full((len(LNC),7,3),np.nan)
    for i,lN in enumerate(LNC):
        inp=f"/opt/lamda/{mol}.dat\no.txt\n80 370\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{DV}\n0\n"
        subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,(J,f) in enumerate(LO[mol].items()):
                        if abs(float(p[4])-f)<2e-3: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    return out
if __name__=='__main__':
    jobs=[(m,T,n) for m in ('cs',) for T in TK for n in LN]
    with Pool(8) as P: res=P.map(one,jobs,chunksize=4)
    G={m:np.array([r for (mm,T,n),r in zip(jobs,res) if mm==m]).reshape(len(TK),len(LN),len(LNC),7,3) for m in ('cs',)}
    np.savez('/workspace/out/q35/grid_cs.npz',TK=TK,LN=LN,LNC=LNC,cs=G['cs'])
    print('done; nan fraction',{m:float(np.isnan(G[m]).mean()) for m in G})
