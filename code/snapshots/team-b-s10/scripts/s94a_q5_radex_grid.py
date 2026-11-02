# Q5: RADEX LVG grids (dv=150 km/s, Tbg 2.73 K, H2 collisions) for CS 7-6, SO 8_8-7_7 & 6_5-5_4 on the Q3 (Tkin,n,N) grid (out/q3/grid_lvg.npz axes). usage: python3 -I s94a_q5_radex_grid.py
import subprocess,os,sys,tempfile,numpy as np
from multiprocessing import Pool
MOL={'cs':{'file':'cs.dat','lines':{'CS7-6':342.8829},'rng':'330 360'},'so':{'file':'so@lique.dat','lines':{'SO6_5-5_4':251.8258,'SO8_8-7_7':344.3106},'rng':'240 360'}}
g=np.load('/workspace/out/q3/grid_lvg.npz'); TK,LN,LNC=g['TK'],g['LN'],g['LNC']; DV=150.
def one(args):
    mol,Tk,ln=args; d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d); M=MOL[mol]
    out=np.full((len(LNC),len(M['lines']),3),np.nan)
    for i,lN in enumerate(LNC):
        inp=f"/opt/lamda/{M['file']}\no.txt\n{M['rng']}\n{Tk}\n1\nH2\n{10**ln}\n2.73\n{10**lN}\n{DV}\n0\n"
        subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
        try:
            for l in open('o.txt'):
                p=l.split()
                if len(p)>=12 and p[1]=='--':
                    for k,f in enumerate(M['lines'].values()):
                        if abs(float(p[4])-f)<3e-3: out[i,k]=(float(p[8]),float(p[6]),float(p[7]))
        except Exception: pass
    return out
if __name__=='__main__':
    res={}
    for mol in MOL:
        jobs=[(mol,T,n) for T in TK for n in LN]
        with Pool(4) as P: r=P.map(one,jobs,chunksize=4)
        res[mol]=np.array(r).reshape(len(TK),len(LN),len(LNC),len(MOL[mol]['lines']),3); print(mol,'nan frac',np.isnan(res[mol]).mean(),flush=True)
    np.savez('/workspace/out/q5/radex_cs_so.npz',TK=TK,LN=LN,LNC=LNC,cs=res['cs'],so=res['so'])
