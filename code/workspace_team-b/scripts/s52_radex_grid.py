# Q3: RADEX grid for HCN and HCO+ (LAMDA hcn.dat, hco+.dat; H2 collisions; Tbg 2.73 K) over Tkin, n(H2), N per dv=150 km/s.
# Output: out/q3/grid_<geom>.npz with T_R, Tex, tau for J=2-1,3-2,4-3 on the grid (NaN where RADEX fails).
# usage: python3 -I s52_radex_grid.py lvg|slab
import subprocess,os,sys,tempfile,numpy as np
from multiprocessing import Pool
geom=sys.argv[1]; exe={'lvg':'radex_lvg','slab':'radex_slab','sphere':'radex_sphere'}[geom]
LO={'hcn':{2:177.2611,3:265.8864,4:354.5055},'hco+':{2:178.3751,3:267.5576,4:356.7342}}
TK=np.array([10,14,20,28,40,56,80,112,160,224,300.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(11.5,19.01,0.25); DV=150.
def one(args):
    mol,Tk,ln=args; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d)
    out=np.full((len(LNC),3,3),np.nan)   # N, J(2,3,4), (TR,Tex,tau)
    for i,lN in enumerate(LNC):
        inp=f"/opt/lamda/{mol}.dat\no.txt\n170 370\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{DV}\n0\n"
        subprocess.run([exe],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,(J,f) in enumerate(LO[mol].items()):
                        if abs(float(p[4])-f)<2e-3: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    return out
if __name__=='__main__':
    jobs=[(m,T,n) for m in ('hcn','hco+') for T in TK for n in LN]
    with Pool(4) as P: res=P.map(one,jobs,chunksize=4)
    G={m:np.array([r for (mm,T,n),r in zip(jobs,res) if mm==m]).reshape(len(TK),len(LN),len(LNC),3,3) for m in ('hcn','hco+')}
    np.savez(f'/workspace/out/q3/grid_{geom}.npz',TK=TK,LN=LN,LNC=LNC,hcn=G['hcn'],hcop=G['hco+'])
    print(geom,'done; nan fraction',{m:float(np.isnan(G[m]).mean()) for m in G})
