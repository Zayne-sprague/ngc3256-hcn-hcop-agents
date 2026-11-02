# Q28 r2 analysis of fit_<variant>.pkl: profiles, posteriors, branch tables. usage: python3 -I s363_q28r2_analyse.py base [more variants]
import sys,pickle,numpy as np
sys.path.insert(0,'/workspace/scripts')
from s362_q28r2_lib import *
import s362_q28r2_fit as F
def wq(x,w,qs=(16,50,84)):
    o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in qs]
def analyse(name):
    o=pickle.load(open(f'/workspace/out/q28r2/fit_{name}.pkl','rb')); log=[]
    def P(s=''): print(s); log.append(s)
    P(f'######## variant {name}: taumax {o["taumax"]} cap {o["cap"]} dropped {np.where(o["wt"]==0)[0].tolist()} sig {o["sig"]}')
    for l in open(f'/workspace/out/q28r2/fit_{name}.txt').read().split('\n')[1:]: P(l)
    gE,fE,gC=o['gE'],o['fE'],o['gC']; ie=o['R1']['N']['ie']
    pE=np.c_[gE[ie//3],fE[ie]]  # logTk,logn,logN,logA,fE
    rho=o['rho']; chiJ=o['R1']['N']['minE']+o['R1']['S']['minE']; c0=chiJ.min()
    P(f' shared-E posterior (joint N+S), 16/50/84%: '+'; '.join(f'{l} {wq(10**pE[:,0] if k==0 else pE[:,k],rho)[1]:.3g} [{wq(10**pE[:,0] if k==0 else pE[:,k],rho)[0]:.3g},{wq(10**pE[:,0] if k==0 else pE[:,k],rho)[2]:.3g}]' for k,l in enumerate(['Tk_E','log n_E','log N(HCO+)_E','log A_E','f_E'])))
    P(f'  P(Tk_E>50)={rho[10**pE[:,0]>50].sum():.3f}  P(Tk_E>100)={rho[10**pE[:,0]>100].sum():.3f}  P(n_E<1e5)={rho[pE[:,1]<5].sum():.3f}  best-fit E: Tk {10**pE[chiJ.argmin(),0]:.0f} log n {pE[chiJ.argmin(),1]:.3f} logN {pE[chiJ.argmin(),2]} A {pE[chiJ.argmin(),3]} fE {pE[chiJ.argmin(),4]}')
    # profile likelihood: chi2_joint(Tk_E) minimised over n,N,A,fE and over both compact phases
    P(' profile chi2_joint(Tk_E)-min (min over n_E,N_E,A_E,f_E and both compact phases):')
    ut=np.unique(pE[:,0]); pr=np.array([chiJ[pE[:,0]==t].min()-c0 for t in ut]); P('   '+' '.join(f'{10**t:.0f}:{p:.1f}' for t,p in zip(ut,pr)))
    ok=ut[pr<4]; P(f'  Tk_E allowed at dchi2<4: {10**ok.min():.0f}-{10**ok.max():.0f} K; excluded dchi2>4 below {10**ok.min():.0f} K')
    un=np.unique(pE[:,1]); pn=np.array([chiJ[pE[:,1]==t].min()-c0 for t in un]); P(' profile chi2_joint(log n_E)-min: '+' '.join(f'{t:.2f}:{p:.1f}' for t,p in zip(un[::2],pn[::2])))
    ok=un[pn<4]; P(f'  log n_E allowed at dchi2<4: {ok.min():.2f}-{ok.max():.2f}')
    for nfix in (4.0,4.5,5.0,5.5):
        m=np.abs(pE[:,1]-nfix)<1e-6
        if m.any(): P(f'  chi2(Tk_E) at fixed log n_E={nfix}: '+' '.join(f'{10**t:.0f}:{chiJ[m&(pE[:,0]==t)].min()-c0:.1f}' for t in ut[::2] if (m&(pE[:,0]==t)).any()))
    # compact-phase posteriors (joint-E weighted)
    gc,trC,taC,okCm=F.buildC(); nC=18
    PCfull=F.config(trC,FC,THS); idx=np.arange(len(PCfull)); mi=idx//nC; ith=(idx//3)%len(THS); jf=idx%3
    for p in 'NS':
        w=o['R2'][p]['wC'].copy(); w[~o['okC']]=0; w/=w.sum()
        Tk=10**gc[mi,0]; P(f' {p} compact posterior: '+'; '.join(f'{l} {wq(v,w)[1]:.3g} [{wq(v,w)[0]:.3g},{wq(v,w)[2]:.3g}]' for l,v in (('Tk_C',Tk),('log n_C',gc[mi,1]),('log N(HCO+)_C',gc[mi,2]),('log A_C',gc[mi,3]),('theta"',THS[ith]),('f_C',FC[jf]))))
        pc=lambda m:w[m].sum()
        prior_cold=np.mean(Tk[o['okC']]<=30.01)
        P(f'   P(Tk_C<=30)={pc(Tk<=30.01):.3f} (prior {prior_cold:.3f})  P(30<Tk<=100)={pc((Tk>30.01)&(Tk<=100.5)):.3f}  P(Tk_C>100)={pc(Tk>100.5):.3f}')
        P(f'   P(theta<=0.3")={pc(THS[ith]<=0.3):.2f}  P(f_C<1)={pc(FC[jf]<1):.2f}')
        A=10**gc[mi,3]; P(f'   N(HCN)/N(HCO+) compact: {wq(A,w)[1]:.2f} [{wq(A,w)[0]:.2f},{wq(A,w)[2]:.2f}]')
    return o,log
def branches(o,fulllog):
    """predicted observables of the best solutions in the cold (Tk_C<=30) and warm (>100) branches, per nucleus"""
    gc,trC,taC,okCm=F.buildC(); gE,trE,taE,okE_=F.buildE(o['taumax']); nC=18; X=50.
    P=lambda s: (print(s),fulllog.append(s))
    for p in 'NS':
        w=o['R2'][p]['wC'].copy(); w[~o['okC']]=0; w/=w.sum(); minC=o['R1'][p]['minC']; argE=o['R1'][p]['argE']
        idx=np.arange(len(minC)); mi=idx//nC; ith=(idx//3)%len(THS); jf=idx%3; Tk=10**gc[mi,0]
        y,sg=o['D'][p]; wt=o['wt']
        for lab,mask in (('cold Tk_C<=30',Tk<=30.01),('mid 30-100',(Tk>30.01)&(Tk<=100.5)),('warm Tk_C>100',Tk>100.5)):
            m=mask&o['okC']; 
            if not m.any(): continue
            ib=np.where(m)[0][minC[m].argmin()]; eb=argE[ib]; me=eb//3; fe=FC[eb%3]; th=THS[ith[ib]]; fc=FC[jf[ib]]
            pred=(F.config(trC[mi[ib]:mi[ib]+1],FC,THS)[(ith[ib]*3+jf[ib])]+W*fe*trE[me][SEL]+o['O'][p])
            P(f' {p} {lab}: P(branch)={w[m].sum():.3f}; best chi2 {minC[ib]:.1f}; C: Tk {Tk[ib]:.0f} log n {gc[mi[ib],1]:.2f} logN {gc[mi[ib],2]} A {gc[mi[ib],3]} theta {th}" fC {fc}; E: Tk {10**gE[me,0]:.0f} log n {gE[me,1]:.2f} logN {gE[me,2]:.2f} A {gE[me,3]} fE {fe}')
            P(f'    predicted 0.24" T_b HCN/HCO+ = {pred[12]/W:.2f}/{pred[13]/W:.2f} K (observed {y[12]/W:.2f}/{y[13]/W:.2f}); 1.29" {pred[10]/W:.2f}/{pred[11]/W:.2f} (obs {y[10]/W:.2f}/{y[11]/W:.2f}); chi2 contributions of the 0.24" pair {(((pred[12:]-y[12:])/sg[12:])**2).sum():.1f}')
            # posterior-weighted observables within the branch (each C paired with its own best E)
            ww=w[m]/w[m].sum(); ii=np.where(m)[0]; sel=np.random.default_rng(1).choice(len(ii),size=min(4000,len(ii)*3),p=ww); ii=ii[sel]
            e=argE[ii]; me_=e//3; fe_=FC[e%3]; th_=THS[ith[ii]]; fc_=FC[jf[ii]]; mc=mi[ii]
            def T(line,beam,jc): return fc_*cp(th_,beam)*trC[mc,jc]+fe_*trE[me_,jc]
            r32_10=T('',0.5,2)/T('',0.5,0); T43=T('',0.5,3); T32=T('',0.5,2); r43_32=T43/T32
            r13c=(1-np.exp(-taC[mc,3]/X))/(1-np.exp(-taC[mc,3])+1e-12); r13e=(1-np.exp(-taE[me_,3]/X))/(1-np.exp(-taE[me_,3])+1e-12)
            fl_c=fc_*cp(th_,2.7)*trC[mc,3]; fl_e=fe_*trE[me_,3]; alpha=(fl_c*r13c+fl_e*r13e)/(fl_c+fl_e)
            Xs=1e-9; mass={X_:2.244e-20*fc_*10**gc[mc,2]/X_*1.1331*(th_*PCAS)**2 for X_ in (1e-9,3e-9,1e-8)}
            q=lambda v:'%.3g [%.3g,%.3g]'%tuple(np.percentile(v,[50,16,84]))
            P(f'    HCN 3-2/1-0 T_b ratio at 0.5": {q(r32_10)}; HCN 4-3 T_b at 0.5": {q(T43)} K; 4-3/3-2: {q(r43_32)}; H13CN/HCN 4-3 flux ratio in 2.7" (X=50): {q(alpha)} (limit 0.12); HCN v2=1 LTE-collisional exp(-1050K/Tk) {np.exp(-1050/Tk[ii]).mean():.1e} (limit 0.02-0.03)')
            P(f'    compact gas mass M(X=1e-9/3e-9/1e-8) = '+' / '.join('%.2g'%np.median(mass[k]) +' Msun' for k in mass)+f'; n_H2 {np.median(10**gc[mc,1]):.2g} cm-3, size {np.median(th_)*PCAS:.0f} pc FWHM; n x size (cm-3 pc) {np.median(10**gc[mc,1]*th_*PCAS):.2g}')
if __name__=='__main__':
    for v in sys.argv[1:]:
        o,log=analyse(v)
        if v in('base','tau1','sig15'): branches(o,log)
        open(f'/workspace/out/q28r2/analysis_{v}.txt','w').write('\n'.join(log))
