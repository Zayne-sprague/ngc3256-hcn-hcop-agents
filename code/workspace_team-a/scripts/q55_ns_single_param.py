# Q55: N vs S pair test with the Q40i machinery: which ONE parameter (T, n, N(HCN), X) must differ between N and S to reproduce R10,R21,R32?
# Reuses q40h_radex tables and out/q40b/q40f_table.json (same inputs as Q40e); adds the N-S pair that Q40i did not run, plus tau and
# vibrational/IR-pumping-free notes. Output: out/q55/ns_single_param.{log,json}
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, os, warnings; warnings.filterwarnings('ignore')
import q40h_radex as R
os.makedirs('/workspace/out/q55',exist_ok=True)
T=json.load(open('/workspace/out/q40b/q40f_table.json'))
H=R.table('hcn',R.LNH); O=R.table('hco+',R.LNO)
LX=np.arange(-1,1.01,0.25); nT,nn,nH=H.shape[:3]
RJ=np.full((4,nT,nn,nH,len(LX)),np.nan); TAU=np.full((nT,nn,nH,len(LX)),np.nan)
for ix,lx in enumerate(LX):
    for ih,lh in enumerate(R.LNH):
        io=np.where(np.abs(R.LNO-(lh-lx))<1e-6)[0]
        if len(io)==0: continue
        for J in range(4): RJ[J,:,:,ih,ix]=H[:,:,ih,J,0]/O[:,:,io[0],J,0]
        TAU[:,:,ih,ix]=H[:,:,ih,2,1]
def obs(reg):
    d={}
    for J,nm in ((0,'Michiyama'),(1,'Xb55'),(2,'Xb4f')):
        c,e,s=T[f'R|{reg}|{nm}']; d[J]=(c,np.sqrt(e*e+s*s+(0.05*c)**2))
    return d
def chi2(reg):
    c=np.zeros(RJ.shape[1:])
    for J,(v,e) in obs(reg).items(): c+=((RJ[J]-v)/e)**2
    return np.where(np.isfinite(c),c,1e9)
CN,CS=chi2('N'),chi2('S')
names={0:'T',1:'n',2:'N(HCN)',3:'X'}
log=[]
def P(s): print(s); log.append(s)
for r in ('N','S'): P(f'{r}: '+'; '.join(f'{v:.3f}+-{e:.3f}' for v,e in obs(r).values())+'  (R10,R21,R32 = HCN/HCO+ Tb)')
res={'all shared':(CN+CS).min()}
for ax in range(4): res[f'{names[ax]} differs']=(CN.min(axis=ax)+CS.min(axis=ax)).min()
P('min joint chi2 (6 data): '+'; '.join(f'{k}: {v:.1f}' for k,v in res.items()))
LT,LNr,LNHa=R.TK,R.LN,R.LNH
vals={0:LT,1:LNr,2:LNHa,3:LX}
off={}
for ax in (0,1,2,3):
    tot=CN.min(axis=ax)+CS.min(axis=ax); ok=np.argwhere(tot<=tot.min()+4); d=[]
    for i in ok:
        idx=list(i); sl=[]; k=0
        for dim in range(4):
            if dim==ax: sl.append(slice(None))
            else: sl.append(idx[k]); k+=1
        a=CN[tuple(sl)]; b=CS[tuple(sl)]
        va=vals[ax][a.argmin()]; vb=vals[ax][b.argmin()]; d.append((va-vb))
    d=np.array(d); off[names[ax]]=[float(np.min(d)),float(np.median(d)),float(np.max(d))]
    P(f'  only {names[ax]} differs: value_N - value_S over joint family (chi2<=min+4): min {d.min():+.2f} median {np.median(d):+.2f} max {d.max():+.2f} ({"log10" if ax else "K"} units: T in K, others dex)')
# opacity of HCN 3-2 in the allowed set per region; and the fraction of allowed cells that are thin
for r,c in (('N',CN),('S',CS)):
    ok=(c<=c.min()+4)&(c<=7.8); t=TAU[ok]
    P(f' {r}: chi2min {c.min():.2f}, allowed cells {ok.sum()}; tau32(HCN)>1 {(t>1).mean():.2f} >3 {(t>3).mean():.2f} >10 {(t>10).mean():.2f}; thin(<0.3) chi2min {np.where(TAU<0.3,c,1e9).min():.1f}')
# does equal excitation (same T,n,N) plus X-only difference suffice? number of dex of X needed in the family
json.dump({'chi2':{k:float(v) for k,v in res.items()},'offsets':off},open('/workspace/out/q55/ns_single_param.json','w'),indent=1)
open('/workspace/out/q55/ns_single_param.log','w').write('\n'.join(log))
