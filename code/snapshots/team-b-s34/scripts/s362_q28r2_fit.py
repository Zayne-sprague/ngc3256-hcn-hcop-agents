# Q28 round 2: joint N+S fit. usage: python3 -I s362_q28r2_fit.py VARIANT   (variants listed in VAR); writes out/q28r2/fit_<VARIANT>.pkl/.txt
import sys,pickle,time,os,numpy as np
sys.path.insert(0,'/workspace/scripts')
from s362_q28r2_lib import *
XCN=50.
VAR={ # name: (taumax, dropped data idx, sig(10,2.2,1.29,0.24), cap (X,Sigma_max Msun/pc2) or None, E-separate flag)
 'base':(3.,[],(0.10,0.15,0.15,0.15),None),
 'tau1':(1.,[],(0.10,0.15,0.15,0.15),None),
 'sig15':(3.,[],(0.15,0.20,0.20,0.25),None),
 'no10':(3.,[6,8],(0.10,0.15,0.15,0.15),None),        # drop 2.2" 1-0 (HCN, HCO+)
 'no129':(3.,[10,11],(0.10,0.15,0.15,0.15),None),
 'no024':(3.,[12,13],(0.10,0.15,0.15,0.15),None),
 'no22':(3.,[6,7,8,9],(0.10,0.15,0.15,0.15),None),
 'no10arc':(3.,[0,1,2,3,4,5],(0.10,0.15,0.15,0.15),None),  # no 10" data at all (N/S separated resolutions only)
 'capX1e-9':(3.,[],(0.10,0.15,0.15,0.15),(1e-9,77.)),
 'capX3e-9':(3.,[],(0.10,0.15,0.15,0.15),(3e-9,77.)),
 'capX1e-8':(3.,[],(0.10,0.15,0.15,0.15),(1e-8,77.)),
 'capX1e-8_14':(3.,[],(0.10,0.15,0.15,0.15),(1e-8,14.)),
 'capX1e-9_14':(3.,[],(0.10,0.15,0.15,0.15),(1e-9,14.)),
}
def buildC():
    g,tr,ta=models(np.r_[np.arange(1.0,2.41,0.2),2.48],np.arange(3.5,7.51,0.5),np.arange(13.,17.51,0.5),np.array([-0.5,0,0.5,1.0]))
    ok=ta[:,3]<XCN*0.128; return g,tr,ta,ok
def buildE(taumax):
    g,tr,ta=models(np.round(np.arange(1.0,2.501,0.075),4),np.arange(3.5,6.001,0.125),np.arange(12.5,15.501,0.25),np.array([-0.5,0,0.5,1.0]))
    return g,tr,ta,ta.max(1)<=taumax
def config(tr,fc,th,off=0.):
    S=cscale(th,fc,off)   # [nth,nfc,14]
    P=W*tr[:,None,None,SEL]*S[None]   # [M,nth,nfc,14]
    return P.reshape(-1,14)
def stage(y,sg,wt,PC,okC,PE,okE,O,rho=None,chunk=1500):
    sge=np.sqrt(sg**2+(0.3*O)**2)
    v=(y/sge*wt).astype(np.float32); uc=(PC/sge*wt).astype(np.float32); ue=((PE+O)/sge*wt).astype(np.float32)
    ie=np.where(okE)[0]; ue=ue[ie]; qc=(uc**2).sum(1); qe=(ue**2).sum(1)+(v**2).sum(); lc=uc@v; le=ue@v
    NC=len(PC); ZE=np.zeros(len(ie)); minE=np.full(len(ie),1e9); minC=np.full(NC,1e9); argE=np.zeros(NC,int); wC=np.zeros(NC) if rho is not None else None
    for s in range(0,NC,chunk):
        sl=slice(s,min(s+chunk,NC)); oc=okC[sl]
        chi=(qc[sl,None]-2*lc[sl,None])+(qe-2*le)[None,:]+2*(uc[sl]@ue.T)
        chi[~oc]=1e6
        ZE+=np.exp(-0.5*chi.astype(np.float32)).sum(0,dtype=np.float64); minE=np.minimum(minE,chi.min(0)); 
        minC[sl]=chi.min(1); argE[sl]=chi.argmin(1)
        if rho is not None: wC[sl]=np.exp(-0.5*chi.astype(np.float32))@rho.astype(np.float32)
    return dict(ZE=ZE,minE=minE,minC=minC,argE=ie[argE],wC=wC,ie=ie)
def run(name):
    taumax,drop,sig,cap=VAR[name]; t0=time.time(); log=[]
    def P(s=''): print(s,flush=True); log.append(s)
    th=THS; gC,trC,taC,okCm=buildC(); gE,trE,taE,okEm=buildE(taumax)
    nC=len(THS)*len(FC); NCtot=len(gC)*nC
    okC=np.repeat(okCm,nC); PC=config(trC,FC,THS)
    PCoff=config(trC,FC,THS,SEP)
    NEfull=len(gE)*len(FC); okE=np.repeat(okEm,len(FC)); fE=np.tile(FC,len(gE)); 
    PE=W*np.repeat(trE[:,SEL],len(FC),0)*fE[:,None]
    if cap is not None:
        X,Smax=cap
        SE=2.244e-20*np.repeat(10**gE[:,2],len(FC))*fE/X; okE=okE&(SE<=Smax)
        SC=2.244e-20*np.repeat(10**gC[:,2],nC)*np.tile(np.repeat(FC[None],len(THS),0).ravel(),len(gC))*np.tile(np.repeat(THS**2/(THS**2+100),len(FC)),len(gC))/X
        okC=okC&(SC<=Smax)
    P(f'# variant {name}: taumax {taumax} drop {drop} sig {sig} cap {cap}; compact configs {nC*len(gC)} (ok {okC.sum()}), E models {NEfull} (ok {okE.sum()})')
    wt=np.ones(14); wt[drop]=0
    offf='/workspace/out/q28r2/offsets_base.npy'
    D={p:data(p,sig) for p in 'NS'}; O={'N':np.zeros(14),'S':np.zeros(14)}
    if os.path.exists(offf) and name!='base': O=dict(zip('NS',np.load(offf)))
    npass=2 if name=='base' else 1
    for it in range(npass):
        R1={p:stage(*D[p],wt,PC,okC,PE,okE,O[p]) for p in 'NS'}
        lj=np.log(R1['N']['ZE']/NCtot)+np.log(R1['S']['ZE']/NCtot)
        rho=np.exp(lj-lj.max()); rho/=rho.sum()
        R2={}
        for p,q in (('N','S'),('S','N')):
            R2[p]=stage(*D[p],wt,PC,okC,PE,okE,O[p],rho=rho/(R1[p]['ZE']+1e-300))
        if name=='base':
            oN=(R2['S']['wC']/R2['S']['wC'].sum())@PCoff*(wt*0+1); oS=(R2['N']['wC']/R2['N']['wC'].sum())@PCoff
            oN=oN.copy(); oS=oS.copy(); oN[6:]=0; oS[6:]=0; O={'N':oN,'S':oS}; np.save(offf,np.array([oN,oS]))
            P(f' iteration {it}: offsets (S compact in N 10", N compact in S 10") K km/s: N {np.round(oN[:6],3).tolist()}  S {np.round(oS[:6],3).tolist()}')
    if name=='base':   # final pass with converged offsets
        R1={p:stage(*D[p],wt,PC,okC,PE,okE,O[p]) for p in 'NS'}
        lj=np.log(R1['N']['ZE']/NCtot)+np.log(R1['S']['ZE']/NCtot); rho=np.exp(lj-lj.max()); rho/=rho.sum()
        R2={p:stage(*D[p],wt,PC,okC,PE,okE,O[p],rho=rho/(R1[p]['ZE']+1e-300)) for p in 'NS'}
    out=dict(name=name,gC=gC,gE=gE,okC=okC,okE=okE,fE=fE,R1=R1,R2=R2,rho=rho,O=O,NCtot=NCtot,NEfull=NEfull,wt=wt,D=D,sig=sig,cap=cap,taumax=taumax)
    # evidences
    lZE={p:np.log(R1[p]['ZE'].sum()/(NCtot*NEfull)) for p in 'NS'}
    lZsh=np.log(np.sum(R1['N']['ZE']/NCtot*R1['S']['ZE']/NCtot)/NEfull)
    out.update(lZ_sep={p:float(lZE[p]) for p in 'NS'},lZ_shared=float(lZsh))
    for p in 'NS': P(f' {p}: min chi2 two-phase (own E) {R1[p]["minE"].min():.2f}; lnZ {lZE[p]:.2f}')
    cj=R1['N']['minE']+R1['S']['minE']; ij=cj.argmin()
    P(f' joint shared-E: min chi2 {cj.min():.2f} (N {R1["N"]["minE"][ij]:.2f}, S {R1["S"]["minE"][ij]:.2f}); lnZ shared {lZsh:.2f} vs separate-E {lZE["N"]+lZE["S"]:.2f}')
    # one-state two-size null (per nucleus)
    g1,tr1,ta1=models(np.r_[np.arange(1.0,2.41,0.2),2.48],np.arange(3.5,7.51,0.5),np.arange(12.5,17.51,0.5),np.array([-0.5,0,0.5,1.0])); ok1=ta1[:,3]<XCN*0.128
    S1=cscale(THS,FC); c2=np.array([0.03,0.1,0.3,1.0])
    for p in 'NS':
        y,sg=D[p]; sge=np.sqrt(sg**2+(0.3*O[p])**2)
        pred=W*tr1[:,SEL][:,None,None,None,:]*(S1[None,:,:,None,:]+c2[None,None,None,:,None])+O[p]
        chi=(((pred-y)/sge*wt)**2).sum(-1); chi[~ok1]=1e6
        l=np.log(np.exp(-0.5*chi).mean()); out[f'lZ1_{p}']=float(l); out[f'chi1_{p}']=float(chi.min())
        i=np.unravel_index(chi.argmin(),chi.shape)
        P(f' {p}: one-state two-size: min chi2 {chi.min():.2f} lnZ {l:.2f}; best Tk {10**g1[i[0],0]:.0f} n {g1[i[0],1]} N {g1[i[0],2]} A {g1[i[0],3]} theta {THS[i[1]]} fC {FC[i[2]]} cov_ext {c2[i[3]]}; ln BF(two-phase/one-state) {lZE[p]-l:.2f}')
        out[f'best1_{p}']=(g1[i[0]],THS[i[1]],FC[i[2]],c2[i[3]])
    P(f' joint ln BF shared two-phase / (one-state N x one-state S) = {lsh(out):.2f}')
    P(f' time {time.time()-t0:.0f}s')
    pickle.dump(out,open(f'/workspace/out/q28r2/fit_{name}.pkl','wb')); open(f'/workspace/out/q28r2/fit_{name}.txt','w').write('\n'.join(log))
def lsh(o): return o['lZ_shared']-o['lZ1_N']-o['lZ1_S']
if __name__=='__main__': run(sys.argv[1])
