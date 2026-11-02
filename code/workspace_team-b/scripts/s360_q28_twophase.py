# Q28: two-phase (compact + extended thin) RADEX model of N and S vs single component, multi-beam data.
# Data per nucleus (K km/s, RJ, rest frame; all from existing team products):
#   10" Gaussian beam : HCN, HCO+ 2-1, 3-2, 4-3  (s53_ladder_fit.DATA, Xb55/Xb4f/Xb47 smoothed)           rel err sqrt(stat^2+SIG10^2)
#   2.2" aperture     : HCN, HCO+ 1-0 (Q29 s183, 2015.1.00993) and 3-2 (Q34 s192 12m+7m), +-150 km/s box x 1.00927  rel err SIG22
#   1.29" (250 pc)    : HCN, HCO+ 3-2 peak T_b (Q6r2 s160.json) -> 159.7*Tb (w=150 km/s)                          rel err SIG13
#   0.24" (49 pc)     : floor: model T_R(3-2)(0.24") >= 0.85 T_b(s160) for both lines (one-sided gaussian penalty, sigma 10% of Tb)
#   H13CN 4-3 limit (s108 convention): tau43(HCN, compact) < X*0.128, X=12C/13C = 50 (30,100 variants)
# Model: RADEX LVG grid out/q3d/grid_lvg500.npz (dv=150 km/s for every component = shared line width; I=1.0645*150*T_R).
#   compact C : Tk, n, N(HCO+), A=N(HCN)/N(HCO+), Gaussian source FWHM theta; coupling to a beam b: theta^2/(theta^2+b^2); eta_c(10")=theta^2/(theta^2+100)
#   extended E: Tk, n, N(HCO+), A, fills every beam (coupling 1), optically thin: tau<=TAUMAX for all HCN/HCO+ J<=4
#   single    : one component of size theta (couples theta^2/(theta^2+b^2)); theta up to 30"
# Evidence: mean of exp(-chi2/2) over flat grids (log Tk 10-300, log n 3-8, log N(HCO+) 13-18, log A -0.5..1, log theta); Bayes factor 2C/1C includes Occam penalty.
# usage: python3 -I s360_q28_twophase.py [SIG10 SIG22 SIG13 X TAUMAX]  -> out/q28/s360_<tag>.txt/.json
import sys,json,numpy as np,itertools
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L
SIG10,SIG22,SIG13,XCN,TAUMAX=[float(x) for x in sys.argv[1:6]] if len(sys.argv)>5 else (0.10,0.15,0.15,50.,1.0)
TAG=f's{SIG10}_{SIG22}_{SIG13}_X{XCN:.0f}_t{TAUMAX}'; K=1.00927; W=159.7
import os; os.makedirs('/workspace/out/q28',exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
d=np.load('/workspace/out/q3d/grid_lvg500.npz'); TK=d['TK']; LN=d['LN']; LNC=d['LNC']
TR={m:np.nan_to_num(d[m][...,0],nan=0.0) for m in ('hcn','hcop')}; TAU={m:np.nan_to_num(d[m][...,2],nan=0.0) for m in ('hcn','hcop')}
NOg=np.arange(12.0,18.0001,0.5); Ag=np.array([-0.5,0.0,0.5,1.0])
ic=lambda x:int(round((x-LNC[0])/0.25))
def table(Tmax=300.):
    rows=[];par=[];tmax=[];tau43=[]
    for iT,t in enumerate(TK):
        if t>Tmax+1e-6: continue
        for iN,n in enumerate(LN):
            for no in NOg:
                for a in Ag:
                    if no+a>LNC[-1]+1e-6: continue
                    io=ic(no); ih=ic(no+a)
                    h=TR['hcn'][iT,iN,ih,:]; o=TR['hcop'][iT,iN,io,:]
                    th=TAU['hcn'][iT,iN,ih,:]; to=TAU['hcop'][iT,iN,io,:]
                    rows.append(np.concatenate([h,o])); par.append((np.log10(t),n,no,a)); tmax.append(max(th.max(),to.max())); tau43.append(th[3])
    return np.array(rows),np.array(par),np.array(tmax),np.array(tau43)
R,PAR,TMX,T43=table(); P(f'# grid models: {len(R)}')
# observable ordering: [HCN2,3,4,HCO2,3,4 @10"], [HCN1,HCN3,HCO1,HCO3 @2.2"], [HCN3,HCO3 @1.29"]; floors at 0.24"
J10=[1,2,3,5,6,7]; J22=[0,2,4,6]; J13=[2,6]
BE=dict(b10=10.,b22=2.2,b13=1.29,b02=0.24)
r183=json.load(open('/workspace/out/q29/s183.json')); r192=json.load(open('/workspace/out/q34/s192.json'))['res']; s160=json.load(open('/workspace/out/q6r2/s160.json'))
def data(pos):
    I10,e10=L.DATA[pos]; a=r183['2015.1.00993 X961|2.2|'+pos]; c=r192['2p2|comb|'+pos]
    I22=np.array([a['HCN']['h150_ap'][0],c['HCN32|ap|h150'][0],a['HCOp']['h150_ap'][0],c['HCOp32|ap|h150'][0]])*K
    I13=np.array([s160['HCN3-2'][pos]['250pc(1.29")']['Tb'],s160['HCO+3-2'][pos]['250pc(1.29")']['Tb']])*W
    fl=np.array([s160['HCN3-2'][pos]['native 0.24"']['Tb'],s160['HCO+3-2'][pos]['native 0.24"']['Tb']])
    y=np.concatenate([I10,I22,I13]); sg=np.concatenate([y[:6]*np.sqrt(e10**2+SIG10**2),I22*SIG22,I13*SIG13])
    return y,sg,fl
BEAMS=np.array([10]*6+[2.2]*4+[1.29]*2); SEL=J10+J22+J13
def chi_pred(y,sg,pred,fl,predfl):
    chi=(((pred-y)/sg)**2).sum(-1)
    chi=chi+((np.maximum(0.85*fl-predfl,0)/(0.1*fl))**2).sum(-1)
    return chi
THS=np.array([0.15,0.2,0.3,0.5,0.8,1.2,1.8,2.7,4.,6.]); THS1=np.concatenate([THS,[12.,30.]]); CG=np.array([0.03,0.1,0.3,1.0])
cp=lambda th,b:th**2/(th**2+b**2)
def single(pos):
    y,sg,fl=data(pos); out=[]
    for th in THS1:
        f=cp(th,BEAMS); pred=W*R[:,SEL]*f[None,:]; pf=R[:,[2,6]]*cp(th,0.24)
        c=chi_pred(y,sg,pred,fl,pf); c=np.where(T43<XCN*0.128,c,1e6); out.append(c)
    return np.array(out)    # [theta,model]
def single2(pos):
    # one gas state (Tk,n,N,A) with two source sizes: compact gaussian theta (peak covering c1<=1) + extended covering c2<=1 (coupling 1)
    y,sg,fl=data(pos); best=[]
    for th in THS:
        f=cp(th,BEAMS)
        for c1 in CG:
            for c2 in CG:
                pred=W*R[:,SEL]*(c1*f+c2)[None,:]; pf=R[:,[2,6]]*(c1*cp(th,0.24)+c2)
                c=chi_pred(y,sg,pred,fl,pf); best.append(np.where(T43<XCN*0.128,c,1e6))
    return np.array(best)
def double(pos,chunk=40):
    y,sg,fl=data(pos); ie=np.where(TMX<=TAUMAX)[0]; ic_=np.where(T43<XCN*0.128)[0]
    P(f'  {pos}: compact models {len(ic_)}, extended thin models {len(ie)}')
    Rc=R[ic_]; Re=R[ie]
    best=(1e9,None); Zacc=0.; n=0; marg={}
    chis=np.empty((len(THS),len(ic_),len(ie)),np.float32)
    for k,th in enumerate(THS):
        f=cp(th,BEAMS)
        Pc=W*Rc[:,SEL]*f[None,:]; Pe=W*Re[:,SEL]
        fc02=cp(th,0.24)*Rc[:,[2,6]]; fe=Re[:,[2,6]]
        for s in range(0,len(ic_),chunk):
            pred=Pc[s:s+chunk,None,:]+Pe[None,:,:]; pf=fc02[s:s+chunk,None,:]+fe[None,:,:]
            chis[k,s:s+chunk]=chi_pred(y,sg,pred,fl,pf)
    return chis,ic_,ie
def q(x,p,qs=(16,50,84)):
    o=np.argsort(x); c=np.cumsum(p[o]); c=c/c[-1]; return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in qs]
res={}
for pos in ('N','S'):
    P(f'\n######## {pos}'); y,sg,fl=data(pos); P(' data (10": HCN2,3,4 HCO2,3,4 | 2.2": HCN1,HCN3,HCO1,HCO3 | 1.29": HCN3,HCO3 K km/s):  '+str(np.round(y,2).tolist())+f'  floors 0.24" T_b {fl.round(2).tolist()}')
    c1=single(pos); Z1=np.exp(-0.5*c1).mean(); i=np.unravel_index(c1.argmin(),c1.shape)
    P(f' single: min chi2 {c1.min():.1f} (theta {THS1[i[0]]}, Tk {10**PAR[i[1],0]:.0f}, log n {PAR[i[1],1]}, logNO {PAR[i[1],2]}, A {PAR[i[1],3]}); lnZ1 {np.log(Z1):.2f}')
    p1=np.exp(-0.5*(c1-c1.min())); p1/=p1.sum(); P('  single posterior theta: '+' '.join(f'{t}:{p1[k].sum():.2f}' for k,t in enumerate(THS1)))
    c1b=single2(pos); Z1b=np.exp(-0.5*c1b).sum()/(len(c1b)*len(R)); i1b=np.unravel_index(c1b.argmin(),c1b.shape)
    P(f' single state, two sizes (c1,c2 in 0.03-1): min chi2 {c1b.min():.1f} (Tk {10**PAR[i1b[1],0]:.0f}, log n {PAR[i1b[1],1]}, logNO {PAR[i1b[1],2]}, A {PAR[i1b[1],3]}); lnZ1b {np.log(Z1b):.2f}')
    chis,ic_,ie=double(pos); c2=np.minimum(chis,1e4)
    # prior: all pairs of grid models x theta; masked pairs have likelihood 0
    Z2=np.exp(-0.5*c2,dtype=np.float32).sum(dtype=np.float64)/(len(THS)*len(R)*len(R))
    Z1=np.exp(-0.5*c1).sum()/(len(THS1)*len(R))
    i=np.unravel_index(c2.argmin(),c2.shape)
    cmin2=float(c2.min()); P(f' single lnZ1 {np.log(Z1):.2f}; two-phase lnZ2 {np.log(Z2):.2f}; ln BF(2/1) = {np.log(Z2/Z1):.2f}; ln BF(2/1b) = {np.log(Z2/Z1b):.2f}; min chi2 two-phase {cmin2:.1f} (dof ~{len(y)+2-7}-... ) vs single {c1.min():.1f}')
    kc,ke=ic_[i[1]],ie[i[2]]
    P(f' best two-phase: theta {THS[i[0]]}; C: Tk {10**PAR[kc,0]:.0f} log n {PAR[kc,1]:.2f} logNO {PAR[kc,2]} A {PAR[kc,3]} ; E: Tk {10**PAR[ke,0]:.0f} log n {PAR[ke,1]:.2f} logNO {PAR[ke,2]} A {PAR[ke,3]}')
    w=c2; w-=c2.min(); w*=-0.5; np.exp(w,out=w); w/=w.sum(dtype=np.float64); wt=w.sum((1,2),dtype=np.float64); wc=w.sum((0,2),dtype=np.float64); we=w.sum((0,1),dtype=np.float64); R_={}
    P('  posterior theta: '+' '.join(f'{t}:{wt[k]:.2f}' for k,t in enumerate(THS)))
    for nm,ww,idx in (('C',wc,ic_),('E',we,ie)):
        for j,lab,tr in ((0,'Tk',lambda v:10**v),(1,'log n',lambda v:v),(2,'log N(HCO+)',lambda v:v),(3,'log A',lambda v:v)):
            v=PAR[idx,j]; qq=q(tr(v),ww); R_[f'{nm}_{lab}']=qq; P(f'  {nm} {lab}: {qq[1]:.2f} [{qq[0]:.2f},{qq[2]:.2f}]')
    pc30=float(wc[10**PAR[ic_,0]<=30.01].sum()); pe50=float(we[10**PAR[ie,0]>50].sum()); pc50=float(wc[10**PAR[ic_,0]>50].sum())
    P(f'  P(Tk_C<=30)={pc30:.3f}  P(Tk_C>50)={pc50:.3f}  P(Tk_E>50)={pe50:.3f}; prior P(Tk_C<=30)={float((10**PAR[ic_,0]<=30.01).mean()):.2f}')
    # core flux fraction at 10" in 3-2 and eta_c, gas surface density
    wf=[];Sg=[]
    th_idx=np.arange(len(THS))[:,None,None]
    ii=np.argwhere(w>1e-6*w.max()); ww=w[tuple(ii.T)]
    for X in (1e-9,3e-9,1e-8):
        pass
    th=THS[ii[:,0]]; eta=th**2/(th**2+100); NOc=10**PAR[ic_[ii[:,1]],2]; NOe=10**PAR[ie[ii[:,2]],2]
    fr=eta*R[ic_[ii[:,1]],2]/(eta*R[ic_[ii[:,1]],2]+R[ie[ii[:,2]],2])
    P(f'  compact fraction of HCN 3-2 flux in 10" beam: {q(fr,ww)[1]:.2f} [{q(fr,ww)[0]:.2f},{q(fr,ww)[2]:.2f}]; eta_c(10") {q(eta,ww)[1]:.4f}')
    for X in (1e-9,3e-9,1e-8):
        Sig=2.244e-20*(eta*NOc+NOe)/X; qq=q(Sig,ww); Sc=2.244e-20*eta*NOc/X; P(f'  X(HCO+)={X:.0e}: Sigma_gas(10")={qq[1]:.0f} [{qq[0]:.0f},{qq[2]:.0f}] Msun/pc2 (compact part {q(Sc,ww)[1]:.0f}); P(<=77)={float(ww[Sig<=77].sum()):.2f} P(<=14)={float(ww[Sig<=14].sum()):.2f}')
    Xmin=q(2.244e-20*(eta*NOc+NOe)/77.,ww); P(f'  X(HCO+) needed for Sigma<=77 Msun/pc2 (dust bound, Td=30): median {Xmin[1]:.1e} [{Xmin[0]:.1e},{Xmin[2]:.1e}]; for Sigma<=14 (Td=100): {np.median(q(2.244e-20*(eta*NOc+NOe)/14.,ww)):.1e}')
    R_.update(lnBF1b=float(np.log(Z2/Z1b)),chi2_1b=float(c1b.min()),lnZ1=float(np.log(Z1)),lnZ2=float(np.log(Z2)),lnBF=float(np.log(Z2/Z1)),chi2_1=float(c1.min()),chi2_2=cmin2,P_TcLe30=pc30,P_TcGt50=pc50,P_TeGt50=pe50,theta_post=wt.tolist(),Cfrac=q(fr,ww))
    res[pos]=R_
json.dump(res,open(f'/workspace/out/q28/s360_{TAG}.json','w'),indent=1); open(f'/workspace/out/q28/s360_{TAG}.txt','w').write('\n'.join(log))
