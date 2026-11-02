"""Q62: RADEX (LVG, dv=150 km/s, Tbg 2.73) grid for 12CO and 13CO J=1-0,2-1,3-2: Tk x log n(H2) x log N per 150 km/s.
usage: python3 -I s372_q62_co_grid.py [lvg|slab|sphere] -> out/q62/co_grid_<geom>.npz   (T_R [K], Tex, tau)"""
import subprocess,os,sys,tempfile,numpy as np
from multiprocessing import Pool
geom=sys.argv[1] if len(sys.argv)>1 else 'lvg'; exe={'lvg':'radex_lvg','slab':'radex_slab','sphere':'radex_sphere'}[geom]
FR={'co':[115.2712,230.538,345.796],'13co':[110.2014,220.3987,330.5880]}
TK=np.array([10,15,20,30,40,60,80,100,133,200,300.]); LN=np.arange(2.5,6.01,0.25); LNC=np.arange(14.,19.51,0.25); DV=150.
def one(a):
    mol,Tk,ln=a; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d)
    out=np.full((len(LNC),3,3),np.nan)
    for i,lN in enumerate(LNC):
        inp=f"/opt/lamda/{mol}.dat\no.txt\n100 400\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{DV}\n0\n"
        subprocess.run([exe],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,f in enumerate(FR[mol]):
                        if abs(float(p[4])-f)<2e-3: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    return out
if __name__=='__main__':
    jobs=[(m,T,n) for m in ('co','13co') for T in TK for n in LN]
    with Pool(14) as P: res=P.map(one,jobs,chunksize=2)
    G={m:np.array([r for (mm,T,n),r in zip(jobs,res) if mm==m]).reshape(len(TK),len(LN),len(LNC),3,3) for m in ('co','13co')}
    np.savez(f'/workspace/out/q62/co_grid_{geom}.npz',TK=TK,LN=LN,LNC=LNC,co=G['co'],c13=G['13co'])
    print(geom,'done; nan frac',{m:float(np.isnan(G[m]).mean()) for m in G})
