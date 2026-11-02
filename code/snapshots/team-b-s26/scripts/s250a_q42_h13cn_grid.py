# Q42: RADEX grid for H13CN with its own LAMDA file (h13cn@xpol.dat), same conventions as s220a (dv 150 km/s, Tbg 2.73, H2, Tk 30pts 10-300 K, log n 3-8, log N(H13CN) 11.5-19).
# Lines J=1-0 (86.34), 3-2 (259.01), 4-3 (345.34 GHz). Output out/q42/grid_h13cn.npz arr[T,n,N,line,(TR,Tex,tau)]
# usage: python3 -I s250a_q42_h13cn_grid.py
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, subprocess, tempfile, os
from multiprocessing import Pool
import s220a_q39_fine_grid as g
LO=[86.3402,259.0118,345.3398]
def one(a):
    Tk,ln=a; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d); out=np.full((len(g.LNC),3,3),np.nan)
    for i,lN in enumerate(g.LNC):
        inp=f"/opt/lamda/h13cn@xpol.dat\no.txt\n80 370\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{g.DV}\n0\n"
        subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,f in enumerate(LO):
                        if abs(float(p[4])-f)<2e-2: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    subprocess.run(['rm','-rf',d]); return out
if __name__=='__main__':
    jobs=[(T,n) for T in g.TK for n in g.LN]
    with Pool(16) as P: res=P.map(one,jobs,chunksize=2)
    G=np.array(res).reshape(len(g.TK),len(g.LN),len(g.LNC),3,3)
    np.savez('/workspace/out/q42/grid_h13cn.npz',TK=g.TK,LN=g.LN,LNC=g.LNC,h13cn=G); print('nan frac',np.isnan(G).mean())
