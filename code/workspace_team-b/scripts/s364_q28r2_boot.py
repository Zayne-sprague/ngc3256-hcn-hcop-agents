# Q28 r2 parametric bootstrap per nucleus: simulate data from the best two-phase (GEN=2p) or best one-state-two-size (GEN=1s) model with the assumed errors,
# refit both, record chi2_2p, chi2_1s. usage: python3 -I s364_q28r2_boot.py POS GEN NREP
import sys,pickle,numpy as np
sys.path.insert(0,'/workspace/scripts')
from s362_q28r2_lib import *
import s362_q28r2_fit as F
pos,gen,nrep=sys.argv[1],sys.argv[2],int(sys.argv[3])
o=pickle.load(open('/workspace/out/q28r2/fit_base.pkl','rb')); y,sg=o['D'][pos]; O=o['O'][pos]; sge=np.sqrt(sg**2+(0.09*O**2*0)+(0.3*O)**2)
gc,trC,taC,okCm=F.buildC(); gE,trE,taE,okEm=F.buildE(3.); PC=F.config(trC,FC,THS); PE=W*np.repeat(trE[:,SEL],3,0)*np.tile(FC,len(gE))[:,None]
R=o['R1'][pos]; c=R['minC'].argmin(); e=R['argE'][c]; best2=PC[c]+PE[e]+O
g1,tr1,ta1=models(np.r_[np.arange(1.0,2.41,0.2),2.48],np.arange(3.5,7.51,0.5),np.arange(12.5,17.51,0.5),np.array([-0.5,0,0.5,1.0])); ok1=ta1[:,3]<50*0.128
S1=cscale(THS,FC); c2=np.array([0.03,0.1,0.3,1.0])
P1=W*tr1[:,SEL][:,None,None,None,:]*(S1[None,:,:,None,:]+c2[None,None,None,:,None])+O   # [M,th,fc,c2,14]
chi1=lambda yy:np.where(ok1[:,None,None,None],(((P1-yy)/sge)**2).sum(-1),1e6).min()
b1=o[f'best1_{pos}']; j=np.where((g1==b1[0]).all(1))[0][0]; best1=P1[j,list(THS).index(b1[1]),list(FC).index(b1[2]),list(c2).index(b1[3])]
mu=best2 if gen=='2p' else best1
rng=np.random.default_rng(7 if gen=='2p' else 8); Y=mu+sge*rng.standard_normal((nrep,14)); 
okC=o['okC']; okE=o['okE']; ie=np.where(okE)[0]
uc=(PC/sge).astype(np.float32); ue=((PE+O)[ie]/sge).astype(np.float32); V=(Y/sge).astype(np.float32)
qc=(uc**2).sum(1); qe=(ue**2).sum(1); lc=uc@V.T; le=ue@V.T; v2=(V**2).sum(1)
mins=np.full(nrep,1e9)
for s in range(0,len(uc),1500):
    sl=slice(s,s+1500); G=2*(uc[sl]@ue.T)+qc[sl,None]+qe[None,:]; G[~okC[sl]]=np.inf
    for r in range(nrep):
        mins[r]=min(mins[r],(G-2*lc[sl,r][:,None]-2*le[:,r][None,:]).min()+v2[r])
c1=np.array([chi1(Y[r]) for r in range(nrep)])
obs2=R['minE'].min(); obs1=o[f'chi1_{pos}']
res=dict(pos=pos,gen=gen,chi2p=mins,chi1s=c1,obs2p=obs2,obs1s=obs1)
pickle.dump(res,open(f'/workspace/out/q28r2/boot_{pos}_{gen}.pkl','wb'))
d=c1-mins
print(f'{pos} gen={gen} nrep={nrep}: chi2_2p mean {mins.mean():.1f} [{np.percentile(mins,5):.1f},{np.percentile(mins,95):.1f}] (obs {obs2:.1f}; frac sims with chi2_2p>=obs {np.mean(mins>=obs2):.2f}); chi2_1s mean {c1.mean():.1f} (obs {obs1:.1f}); Dchi2=1s-2p mean {d.mean():.1f} 95% {np.percentile(d,95):.1f} max {d.max():.1f} (obs {obs1-obs2:.1f}; P(sim>=obs) {np.mean(d>=obs1-obs2):.3f})')
