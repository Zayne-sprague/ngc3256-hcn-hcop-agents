# Q3 round 4, fixes 4-5: (A) J=1-0 test against the 10" ladder posterior and as a datum; compact-component fit at 1.7"; (B) derived quantities: filling factor and size, HCN/HCO+ filling ratio q, N(H2) for X(HCO+), Sigma, pressure, critical densities, cold-branch table, line-width consistency.
# 1-0 inputs (ASSUMPTIONS, labelled): Harada+2018 Table 13 1.7" HCN/HCO+ 1-0 intensities (K km/s) for N and S as parsed by team-a (/snapshots/team-a-s04/out/q4c/harada_table13_rows.csv); our own Xb4d 3-2 convolved to 1.7" (team-a ours_1p7arcsec.csv).
# usage: python3 -I s66_j10_derived.py -> out/q3d/s66_j10_derived.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L, s62_lib as S
from scipy.special import logsumexp,log_ndtr
OUT='/workspace/out/q3d/'; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(); T=10**F.lT; nm=F.nm
W0=S.prior_w(F,'logT300')
def pos(k,cls='centre'):
    I,e=L.DATA[k]; return dict(y=np.log(I),rel=e,cls=cls)
G1=np.isfinite(F.M1).all(1)
def post(lnp,w=W0,o=None): return S.Post(F,np.where(np.isfinite(lnp),lnp,-1e30),w,o)
def q16(Ps,x,qs=(16,50,84)): return Ps.wq(x,qs)
R10m=np.where(G1,np.exp(F.M1[:,0]-F.M1[:,1]),1.0); R32m=np.exp(F.M[:,1]-F.M[:,4]); R21m=np.exp(F.M[:,0]-F.M[:,3]); R43m=np.exp(F.M[:,2]-F.M[:,5])
tH=10**F.tauH[:,1]; A_=10**F.lA
BR={'thin HCN (tau32<1)':tH<1,'intermediate (1-5)':(tH>=1)&(tH<5),'thick HCN (tau32>=5)':tH>=5}
BR2={'A<=1':A_<=1.01,'1<A<=3':(A_>1.01)&(A_<=3.2),'A>3':A_>3.2}
def branch_line(Ps,name):
    s=f'   {name:34s} '+' | '.join(f'{k}: {Ps.prob(m):.2f}' for k,m in BR.items())+'   ;  '+' | '.join(f'{k}: {Ps.prob(m):.2f}' for k,m in BR2.items())
    return s
def brief(Ps):
    d=Ps.summary(); return (f"Tk {d['Tk'][2]:.0f} [{d['Tk'][1]:.0f},{d['Tk'][3]:.0f}] (5%: {Ps.wq(T,(5,))[0]:.0f}); log n {d['log n'][2]:.2f} [{d['log n'][1]:.2f},{d['log n'][3]:.2f}]; log A {d['log A'][2]:.1f} [{d['log A'][1]:.1f},{d['log A'][3]:.1f}]; tau(HCN3-2) {10**d['tau(HCN3-2)'][2]:.1f} [{10**d['tau(HCN3-2)'][1]:.1f},{10**d['tau(HCN3-2)'][3]:.1f}]; P(Tk>50)={d['P']['Tk>50'][0]:.2f}; P(n>1e5)={d['P']['n>1e5'][0]:.2f}")
# ------------------------------------------------------------------ (A) 1-0 test
H10={'N':(197.23,205.25),'S':(56.68,87.87)}; O32={'N':(63.61,85.06),'S':(20.28,53.72)}      # 1-0 (Harada), 3-2 (ours) at 1.7": (HCN, HCO+) K km/s
R10o={k:v[0]/v[1] for k,v in H10.items()}; R32o={k:v[0]/v[1] for k,v in O32.items()}
P('# (A) J=1-0 test. Observed at 1.7" (assumption-labelled literature 1-0): HCN/HCO+ 1-0 N %.3f, S %.3f (Michiyama+2018, independent: 0.94+-0.13, 0.60+-0.08); 3-2 (ours, 1.7") N %.3f, S %.3f; HCN 3-2/1-0 N %.2f S %.2f; HCO+ 3-2/1-0 N %.2f S %.2f'%(R10o['N'],R10o['S'],R32o['N'],R32o['S'],O32['N'][0]/H10['N'][0],O32['S'][0]/H10['S'][0],O32['N'][1]/H10['N'][1],O32['S'][1]/H10['S'][1]))
P('  10" ladder data: HCN/HCO+ 3-2 at 10" (Xb4f): N %.3f, S %.3f, Bridge %.3f -> the 10" aperture is not the 1.7" core; the 1-0 test below is therefore a consistency/branch test of the single zone, not a like-for-like prediction.'%tuple(L.DATA[k][0][1]/L.DATA[k][0][4] for k in ('N','S','Bridge')))
base={}
for k in ('N','S'):
    o=S.gfit(F,[pos(k)],keep_samples=True); base[k]=o
    Ps=post(o['lnpost'],W0,o); res[f'10arcsec|{k}']=Ps.summary()
    P(f'\n## 10" ladder posterior at {k} (band 10%, Q12 recovery; flat log Tk<=300): '+brief(Ps))
    P('   posterior-predictive HCN/HCO+ 1-0 (T_R ratio; filling cancels if equal): median %.2f [16-84: %.2f-%.2f; 2.5-97.5: %.2f-%.2f]; observed %s %.3f. HCN/HCO+ 3-2 predicted %.2f [%.2f-%.2f] (10" observed %.3f)'%(*Ps.wq(R10m,(50,16,84,2.5,97.5))[:3],*Ps.wq(R10m,(2.5,97.5)),k,R10o[k],*Ps.wq(R32m,(50,16,84)),L.DATA[k][0][1]/L.DATA[k][0][4]))
    for sg in (0.10,0.25):
        # extra datum: ln R10 model vs observed (equal-filling, same-band ratio; sigma 10% or 25%)
        lnp=np.where(G1,o['lnpost']-0.5*((np.log(R10m)-np.log(R10o[k]))/sg)**2,-1e30)
        Pr=post(lnp); lnZ_ratio=float(logsumexp(lnp+np.log(np.where(W0>0,W0,1e-300)))-logsumexp(o['lnpost']+np.log(np.where(W0>0,W0,1e-300))))
        res[f'10arcsec+R10|{k}|sig{sg}']=Pr.summary(); res[f'10arcsec+R10|{k}|sig{sg}']['lnBayes']=lnZ_ratio
        P(f'   + datum R10(obs {R10o[k]:.2f}, sigma {sg:.0%}): '+brief(Pr)+f'   [ln(Z_with/Z_without)={lnZ_ratio:+.2f}]')
        P(branch_line(Pr,f'   branch masses with R10 (sig {sg:.0%})'))
    P(branch_line(Ps,'   branch masses 10" ladder alone'))
    P(branch_line(S.Post(F,np.zeros(nm),W0),'   prior'))
# ------------------------------------------------------------------ compact fit at 1.7": HCN,HCO+ 1-0 and 3-2 (4 data), shared band scales for 1-0 and 3-2 (10% each), 5% per-line, eta<=1 flat in ln
def compact(k,sb=(0.10,0.10),sl=0.05,extra=None):
    cols=np.stack([F.M1[:,0],F.M[:,1],F.M1[:,1],F.M[:,4]],1)   # HCN1,HCN3,HCO1,HCO3 (ln of 1.0645*150*TR)
    y=np.log(np.array([H10[k][0],O32[k][0],H10[k][1],O32[k][1]]))
    Cm=np.diag([sl**2]*4)+np.outer([1,0,1,0],[1,0,1,0])*sb[0]**2+np.outer([0,1,0,1],[0,1,0,1])*sb[1]**2
    Ci=np.linalg.inv(Cm); a=float(np.ones(4)@Ci@np.ones(4)); u=Ci@np.ones(4)
    r=y[None,:]-cols; g=r@Ci; q0=np.einsum('ij,ij->i',g,r); bb=g.sum(1); eh=bb/a; chi=q0-bb*bb/a
    l=-0.5*chi+log_ndtr(np.sqrt(a)*(0-eh)); l=np.where(np.isfinite(l),l,-1e30); chi=np.where(np.isfinite(chi),chi,1e30)
    return l,eh,chi
P('\n## compact component at 1.7": 4 data (HCN, HCO+ at 1-0 [Harada, assumption] and 3-2 [ours]), 1-0 and 3-2 scales free to 10% each, 5% per line, one common eta<=1 -> 4 data vs 4 params + eta: weakly constraining, shown against the prior')
for k in ('N','S'):
    l,eh,chi=compact(k); Pc=post(l); res[f'compact|{k}']=Pc.summary()
    P(f'[{k}] min chi2 {chi.min():.2f}; '+brief(Pc)); P(branch_line(Pc,f'   {k} compact'))
    # overlap with the 10" posterior: product of the two (same gas, ignoring different beams) 
    lj=l+base[k]['lnpost']; Pj=post(lj); P(f'[{k}] compact x 10" ladder (assumes the 10" gas also emits the 1.7" 1-0 and 3-2 with the same Tk,n,N,A): '+brief(Pj)+f' ; best joint chi2(10") {base[k]["chi"][np.argmax(lj)]:.2f} chi2(1.7") {chi[np.argmax(lj)]:.2f}'); P(branch_line(Pj,f'   {k} joint'))
    res[f'joint|{k}']=Pj.summary()
# ------------------------------------------------------------------ (B) derived quantities (N alone representative; S for contrast)
KPC=40e3*np.pi/180/3600   # Mpc? pc per arcsec at 40 Mpc
PC_AS=40e6*np.pi/180/3600
P(f'\n# (B) derived quantities. Scale: 1" = {PC_AS:.0f} pc at 40 Mpc; 10" beam omega 113.3 arcsec^2 (equivalent circle of diameter 12.0").')
for k in ('N','S'):
    o=base[k]; Ps=post(o['lnpost'],W0,o); seta,wt=Ps.sample_eta(); sid=o['sidx']
    le=seta[:,0]/np.log(10); eta=10**le
    def wqs(x,qs=(2.5,16,50,84,97.5)):
        oo=np.argsort(x); cs=np.cumsum(wt[oo]); return [float(x[oo][min(np.searchsorted(cs,q/100),len(x)-1)]) for q in qs]
    bm_diam=np.sqrt(113.31*4/np.pi)  # arcsec of equivalent circle
    size=bm_diam*np.sqrt(eta)*PC_AS   # pc, equivalent-circle diameter of the emitting area
    lNO=F.lNO[sid]; N_H=F.lNH[sid]
    P(f'\n## {k}: filling factor and size (single zone; eta = emitting solid angle / 10" beam; flat in ln eta <=1)')
    P('   log10 eta '+' '.join(f'{v:.2f}' for v in wqs(le))+' (2.5,16,50,84,97.5 pct);  eta median %.3f [%.3f,%.3f]'%tuple(10**v for v in wqs(le,(50,16,84))))
    P('   equivalent diameter of the emitting area %.0f pc [%.0f,%.0f] (16-84%%), 2.5-97.5%%: %.0f-%.0f pc; (at 1.7" resolution the N core is unresolved (p0007: HCN core FWHM <~0.3")) '%(*wqs(size,(50,16,84)),*wqs(size,(2.5,97.5))))
    res[f'eta|{k}']=dict(log10eta=wqs(le),size_pc=wqs(size))
    # N(H2)
    for X in (1e-9,1e-8):
        lNH2=lNO-np.log10(X); beam=lNH2+le
        Sig=np.log10(2.244e-20)+beam    # Msun/pc2 (mu=2.8 per H2 incl He)
        P(f'   X(HCO+)={X:.0e}: log N(H2) of the emitting gas {wqs(lNH2,(16,50,84))[1]:.2f} [{wqs(lNH2,(16,50,84))[0]:.2f},{wqs(lNH2,(16,50,84))[2]:.2f}] (95%: {wqs(lNH2)[0]:.1f}-{wqs(lNH2)[4]:.1f}); beam-averaged log <N(H2)> {wqs(beam,(16,50,84))[1]:.2f} [{wqs(beam,(16,50,84))[0]:.2f},{wqs(beam,(16,50,84))[2]:.2f}]; Sigma_mol(10" beam) log {wqs(Sig,(16,50,84))[1]:.2f} [{wqs(Sig,(16,50,84))[0]:.2f},{wqs(Sig,(16,50,84))[2]:.2f}] Msun/pc2')
        res[f'NH2|{k}|{X}']=dict(logN=wqs(lNH2),beam=wqs(beam),logSigma=wqs(Sig))
    lNdv=lNO-np.log10(150)
    P(f'   N(HCO+)/dv: log {wqs(lNdv,(16,50,84))[1]:.2f} [{wqs(lNdv,(16,50,84))[0]:.2f},{wqs(lNdv,(16,50,84))[2]:.2f}] (cm-2 per km/s; 95%: {wqs(lNdv)[0]:.1f}-{wqs(lNdv)[4]:.1f}); N(HCN)/dv log {wqs(F.lNH[sid]-np.log10(150),(16,50,84))[1]:.2f}; the grid uses dv=150 km/s but only N/dv matters: for single clouds of FWHM 54 km/s (Brunetti+2021 sigma 23 km/s, as cited by the review, not re-checked) the cloud column is N/dv x 54 = {wqs(lNdv+np.log10(54),(16,50,84))[1]:.2f} dex (same tau and excitation).')
    lP=F.ln[sid]+F.lT[sid]
    P(f'   P/k = n Tk: log {wqs(lP,(16,50,84))[1]:.2f} [{wqs(lP,(16,50,84))[0]:.2f},{wqs(lP,(16,50,84))[2]:.2f}] K cm-3 (95%: {wqs(lP)[0]:.1f}-{wqs(lP)[4]:.1f}); the review quotes Brunetti+2021 turbulent pressures 1e5-1e10 K cm-3 (not re-checked here): the whole posterior lies inside this range, so pressure does not discriminate.')
# q: HCN filling relative to HCO+ (ln q flat in [ln0.1, ln10], eta_HCN=q eta_HCO+ <=1)
P('\n## HCN/HCO+ filling ratio q (HCN flux = q x eta x model) for N and S: marginalised over ln q flat in [0.1,10] with eta_HCO+<=1 and q eta<=1; in thin lines only A*q is fixed by the data')
qg=np.exp(np.linspace(np.log(0.1),np.log(10),13))
for k in ('N','S'):
    lnps=[];Zq=[]
    for q in qg:
        p=pos(k); p['y']=p['y'].copy(); p['y'][:3]-=np.log(q)
        o=S.gfit(F,[p],keep_samples=False,cap=min(0.0,-np.log(q))); lnps.append(o['lnpost']); Zq.append(logsumexp(o['lnpost']+np.log(np.where(W0>0,W0,1e-300))))
    Zq=np.array(Zq); wq_=np.exp(Zq-Zq.max()); wq_/=wq_.sum()
    lnp=logsumexp(np.stack(lnps),axis=0,b=1.0/len(qg))
    Ps=post(lnp)
    cq=np.cumsum(wq_); qq=lambda f:float(qg[min(np.searchsorted(cq,f),len(qg)-1)])
    P(f'[{k}] posterior of q: median {qq(.5):.2f} [{qq(.16):.2f},{qq(.84):.2f}] ; P(q<1)={wq_[qg<1].sum():.2f} P(q>1)={wq_[qg>1].sum():.2f} (prior flat in ln q: P(q<1)=0.5); with free q: '+brief(Ps)); P(branch_line(Ps,f'   {k} free q'))
    res[f'q|{k}']=dict(q=qg.tolist(),w=wq_.tolist(),summary=Ps.summary())
# critical densities from LAMDA
def lamda(fn):
    L_=[l.strip() for l in open(fn) if not l.strip().startswith('!')]   # LAMDA: value lines only
    nl=int(L_[2]); lev=[(float(L_[3+j].split()[1]),float(L_[3+j].split()[2])) for j in range(nl)]
    j=3+nl; nt=int(L_[j]); tr=[L_[j+1+t].split() for t in range(nt)]; A={(int(t[1])-1,int(t[2])-1):float(t[3]) for t in tr}
    j=j+1+nt; j+=1                       # number of partners
    j+=1                                 # partner name line
    ncol=int(L_[j]); ntmp=int(L_[j+1]); temps=[float(x) for x in L_[j+2].split()]; rows=[L_[j+3+r].split() for r in range(ncol) if True]
    # transition rows follow after the 'number of coll trans' / 'number of temps' / temps header lines
    rows=[L_[j+3+r].split() for r in range(ncol)]
    C={(int(r[1])-1,int(r[2])-1):np.array([float(x) for x in r[3:3+ntmp]]) for r in rows}
    return lev,A,temps,C
P('\n## critical densities n_cr(J->J-1) = A / sum_l C(J->l) (H2 collisions, LAMDA, optically thin, downward rates only); with line trapping the effective values are lower by ~1/tau')
cr={}
for mol,fn in (('HCN','hcn'),('HCO+','hco+')):
    lev,A,temps,C=lamda(f'/opt/lamda/{fn}.dat')
    for Tk in (30,100,200):
        row=[]
        for J in (1,2,3,4):
            u=J; tot=0
            for (a,b),v in C.items():
                if a==u and b<u: tot+=np.interp(Tk,temps,v)
            row.append(A[(u,J-1)]/tot)
        P(f'   {mol:5s} Tk={Tk:3d} K: n_cr(1-0,2-1,3-2,4-3) = '+', '.join(f'{v:.1e}' for v in row)); cr[f'{mol}|{Tk}']=row
res['ncrit']=cr
# cold, thick branch table (ratios-only: band scales free) and the baseline posterior split by Tk
P('\n## cold/thick branch: ratios-only fit (band scales free, no recovery term; eta free) at N, and baseline fit split by Tk')
o=S.gfit(F,[pos('N')],keep_samples=False,sig_b=(10,10,10),use_rec=False); Pr=post(o['lnpost'])
def tab(Ps,mask,name):
    pm=Ps.p*mask; m=pm.sum()
    if m<=0: P(f'   {name}: no mass'); return
    pp=pm/m
    def qs(x): 
        oo=np.argsort(x); cs=np.cumsum(pp[oo]); return [float(x[oo][min(np.searchsorted(cs,q),len(x)-1)]) for q in (.16,.5,.84)]
    P(f'   {name}: mass {m:.2f}; Tk {qs(T)[1]:.0f} [{qs(T)[0]:.0f},{qs(T)[2]:.0f}]; log n {qs(F.ln)[1]:.2f} [{qs(F.ln)[0]:.2f},{qs(F.ln)[2]:.2f}]; log N(HCO+)/dv {qs(F.lNO)[1]-np.log10(150):.2f} [{qs(F.lNO)[0]-np.log10(150):.2f},{qs(F.lNO)[2]-np.log10(150):.2f}]; log A {qs(F.lA)[1]:.1f}; tau(HCN3-2) {10**qs(F.tauH[:,1])[1]:.1f} [{10**qs(F.tauH[:,1])[0]:.1f},{10**qs(F.tauH[:,1])[2]:.1f}]; tau(HCO+3-2) {10**qs(F.tauO[:,1])[1]:.1f}; P/k log {qs(F.ln+F.lT)[1]:.2f}')
    res[f'branch|{name}']=dict(mass=float(m),Tk=qs(T),logn=qs(F.ln),tauH=qs(F.tauH[:,1]),logNO=qs(F.lNO))
Pb=post(base['N']['lnpost'])
for nm_,Ps in (('baseline N (10%)',Pb),('ratios-only N',Pr)):
    P(' '+nm_)
    for lab,m in (('Tk<20 K',T<20),('20-50 K',(T>=20)&(T<50)),('Tk>=50 K',T>=50),('Tk>=100 K',T>=100),('log n>=6',F.ln>=6)): tab(Ps,m,f'{nm_} {lab}')
P('   best chi2 (profile over eta) of cold models (Tk<20 K) vs all, baseline N fit:')
for lab,m in (('Tk<20 K',T<20),('Tk<20 K & log n>=6',(T<20)&(F.ln>=6)),('all',np.ones(nm,bool))): P(f'      {lab}: min chi2 = {base["N"]["chi"][m].min():.2f}')
json.dump(res,open(OUT+'s66_j10_derived.json','w'),indent=1,default=float); open(OUT+'s66_j10_derived.txt','w').write('\n'.join(log))
