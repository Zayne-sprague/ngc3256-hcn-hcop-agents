# Q39: explicit-nuisance version of the fine-grid fit (ln eta grid 0.05 on [-3.6,0] and 7 Gauss-Hermite nodes per band scale), so that one-gas and two-gas
# models have comparable evidences / chi2 and N/S contrasts can share the band scale z (same calibration at N and S).
# Models: 'g3' one gas (HCN,HCO+,CS; Tk,n,eta common);  'g1' gas 1 = HCN+HCO+ ; 'g2' gas 2 = CS (own Tk,n,eta). Two-gas = g1 x g2 (given z): free, or with n2<=n1 (CS from the lower-density gas).
# usage: python3 -I s220d_q39_twogas.py  -> out/q39/tg_<pos>_<model>.npz ; s220e summarises
import sys,numpy as np,json,time
sys.argv=sys.argv[:1]+['none']
sys.path.insert(0,'/workspace/scripts')
import importlib; F=importlib.import_module('s220b_q39_fine_fit')
from multiprocessing import Pool
from numpy.polynomial.hermite_e import hermegauss
O='/workspace/out/q39/'
SB=(0.10,0.15); E=np.round(np.arange(-3.6,0.0001,0.05),3); nE=len(E)
x,wq=hermegauss(7); wq=wq/wq.sum()
Z3,Z6=np.meshgrid(x*SB[0],x*SB[1],indexing='ij'); ZW=np.outer(wq,wq).ravel(); Z3=Z3.ravel(); Z6=Z6.ravel(); nZ=len(ZW)
ZPEN=(Z3/SB[0])**2+(Z6/SB[1])**2
O3=(E[:,None]+Z3[None,:]).ravel(); O6=(E[:,None]+Z6[None,:]).ravel()      # (E,Z) flattened, index e*nZ+z
ZI=np.tile(np.arange(nZ),nE)
SPEC={'g3':['H','O','C'],'g1':['H','O'],'g2':['C']}
G=F.X
def worker_setup(pos,model):
    F.setup(pos,*SB,SPEC[model]); G.O3x=O3; G.O6x=O6
def cell_min(c):
    it,jn=c; m=np.zeros(len(O3))
    for s in G.sp:
        f=F.chi_s(it,jn,s,O3,O6); m=m+f.min(1)
    return m.min()
def cell(args):
    c,cb,thr=args; it,jn=c; sp=G.sp
    f={s:F.chi_s(it,jn,s,O3,O6) for s in sp}
    w={s:np.exp(-0.5*(f[s]-cb)) for s in sp}; S={s:w[s].sum(1) for s in sp}; mins={s:f[s].min(1) for s in sp}
    L=np.prod([S[s] for s in sp],0); tot=sum(mins.values())
    out=dict(L=L.reshape(nE,nZ),cm=tot.reshape(nE,nZ))
    # per-z marginals on N for each species (unnormalised), and profiles
    M={};P={}
    for s in sp:
        oth=np.prod([S[t] for t in sp if t!=s],0) if len(sp)>1 else np.ones(len(O3))
        ofm=sum(mins[t] for t in sp if t!=s) if len(sp)>1 else np.zeros(len(O3))
        M[s]=np.zeros((nZ,F.nN)); P[s]=np.full((nZ,F.nN),np.inf)
        for z in range(nZ):
            sel=np.arange(z,len(O3),nZ)
            M[s][z]=(w[s][sel]*oth[sel,None]).sum(0); P[s][z]=(f[s][sel]+ofm[sel,None]).min(0)+ZPEN[z]
    out['M']=M;out['P']=P
    k=tot<cb+thr; H={}
    if len(sp)>1 and k.any():
        specs=[('A','H','O')]+([('CO','C','O'),('CH','C','H')] if 'C' in sp else [])
        for key,a,b in specs:
            rest=[s for s in sp if s not in (a,b)]
            h=np.zeros((nZ,len(F.LAGS[key]))); p=np.full((nZ,len(F.LAGS[key])),np.inf)
            for z in range(nZ):
                kk=np.where(k&(ZI==z))[0]
                if len(kk)==0: continue
                rw=np.prod([S[s][kk] for s in rest],0) if rest else np.ones(len(kk))
                rf=sum(mins[s][kk] for s in rest) if rest else np.zeros(len(kk))
                h[z],p[z]=F.pair(f[a][kk],f[b][kk],w[a][kk],w[b][kk],rw,rf,F.LAGS[key]); p[z]+=ZPEN[z]
            H[key]=(h,p)
    out['H']=H
    ib=np.argmin(tot); out['best']=(tot[ib]+ZPEN[ZI[ib]],{s:float(F.NF[np.argmin(f[s][ib])]) for s in sp},float(E[ib//nZ]),int(ib%nZ))
    return out
def run(pos,model):
    t=time.time(); worker_setup(pos,model); cells=[(i,j) for i in range(len(F.TK)) for j in range(len(F.LN))]
    with Pool(16) as P_:
        m=np.array(P_.map(cell_min,cells,chunksize=20)); cb=float(m.min()); res=P_.map(cell,[(c,cb,45.) for c in cells],chunksize=10)
    nT,nn=len(F.TK),len(F.LN); L=np.array([r['L'] for r in res]).reshape(nT,nn,nE,nZ)       # exp(-chi/2 + cb/2)
    CM=np.array([r['cm'] for r in res]).reshape(nT,nn,nE,nZ)+ZPEN
    d=dict(cb=cb,L=L.sum(2),CMmin=CM.min(2))   # (T,n,Z): summed over eta
    d['Lz']=(L*ZW).sum((0,1,2))
    for s in G.sp:
        d['M_'+s]=np.sum([r['M'][s] for r in res],0); d['P_'+s]=np.min([r['P'][s] for r in res],0)
    for key in ('A','CO','CH'):
        hs=[r['H'][key] for r in res if key in r['H']]
        if hs: d['h_'+key]=np.sum([h[0] for h in hs],0); d['p_'+key]=np.min([h[1] for h in hs],0)
    bs=min((r['best'] for r in res),key=lambda t:t[0]); d['best']=json.dumps(dict(chi2=float(bs[0]),N=bs[1],eta=bs[2],z=[float(Z3[bs[3]]),float(Z6[bs[3]])]))
    np.savez(O+f'tg_{pos}_{model}.npz',**d); print(pos,model,d['best'],f'{time.time()-t:.0f}s',flush=True)
if __name__=='__main__':
    for pos in 'NSC':
        for model in ('g3','g1','g2'): run(pos,model)
