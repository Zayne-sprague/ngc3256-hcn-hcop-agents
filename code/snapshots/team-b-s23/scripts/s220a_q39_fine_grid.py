# Q39: fine RADEX LVG grid for HCN, HCO+, CS (same conventions as s62_radex_grid.py / s210a: dv 150 km/s, Tbg 2.73 K, LAMDA hcn/hco+/cs, H2 collider).
# Grid: Tk 30 pts log-spaced 10-300 K (0.051 dex), log n 3-8 step 0.1, log N 11.5-19 step 0.1 (interpolated cubically to 0.05 dex downstream).
# Output out/q39/grid_fine.npz: arr[mol][T,n,N,line,(TR,Tex,tau)]
# usage: python3 -I s220a_q39_fine_grid.py
import subprocess,os,tempfile,numpy as np,sys,time
from multiprocessing import Pool
LO={'hcn':[88.6316,265.8864,354.5055],'hco+':[89.1885,267.5576,356.7342],'cs':[97.9810,244.9356,342.8829]}   # J=1(2)-0, 3-2, 4-3 | CS 2-1, 5-4, 7-6
TK=np.logspace(1,np.log10(300),30); LN=np.round(np.arange(3,8.001,0.1),2); LNC=np.round(np.arange(11.5,19.001,0.1),2); DV=150.
def one(a):
    mol,Tk,ln=a; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d); out=np.full((len(LNC),3,3),np.nan)
    for i,lN in enumerate(LNC):
        inp=f"/opt/lamda/{mol}.dat\no.txt\n80 370\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{DV}\n0\n"
        subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,f in enumerate(LO[mol]):
                        if abs(float(p[4])-f)<2e-3: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    subprocess.run(['rm','-rf',d]); return out
if __name__=='__main__':
    mols=sys.argv[1:] or ['hcn','hco+','cs']
    if mols==['test']:
        t=time.time(); r=one(('cs',50.,5.0)); print(time.time()-t); print(r[::10]); sys.exit()
    jobs=[(m,T,n) for m in mols for T in TK for n in LN]
    with Pool(16) as P: res=P.map(one,jobs,chunksize=2)
    G={m:np.array([r for (mm,T,n),r in zip(jobs,res) if mm==m]).reshape(len(TK),len(LN),len(LNC),3,3) for m in mols}
    np.savez('/workspace/out/q39/grid_fine.npz',TK=TK,LN=LN,LNC=LNC,hcn=G['hcn'],hcop=G['hco+'],cs=G['cs'])
    print('done; nan fraction',{m:float(np.isnan(G[m]).mean()) for m in G})
