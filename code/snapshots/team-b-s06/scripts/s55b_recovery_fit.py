# Q3 round 2, fixes 1,2,4,6: single-zone RADEX fits of the 10" ladder per position (N,S,Bridge,E1,E2,W + joint disk) on the fine (interpolated) grid, with a
# one-sided, region-class-dependent flux-recovery nuisance in place of the symmetric common band scale.
# Recovery model (calibrated on the nested Xb4b/Xb4d/Xb4f set, s26 + out/s30_ratios.json: same line, same 10" beam, Xb4d/Xb4f = 0.69-0.87 at N/S/Bridge peaks, 0.32-0.78 at E1/E2/W):
#   true flux = observed / f_band,  f = 1 - L_b,  L_Xb4f = Lf (unknown, >=0, U[0,Lmax]); L_Xb47 = Lf*k47, L_Xb55 = Lf*k55, k in {0.5,1,1.5} (flat; Xb55/Xb47 have different baselines/beams);
#   HCN extra relative loss delta*k_b (HCN more extended at the disk: Xb4d/Xb4f HCN vs HCO+ = 0.65/0.85 E1, 0.34/0.62 E2, 0.6/0.77 W; ~equal at N,S,Bridge).
#   centre (N,S,Bridge): Lf in {0,.1,.2,.3,.4}, delta in {0,.1};  disk (E1,E2,W): Lf in {0,.15,.3,.45,.6}, delta in {0,.15,.3};  'stress': Lf up to 0.7 (centre 0.5), k in {0.5,1,2}.
#   Combos are marginalised with equal prior weight; the band-scale calibration nuisance is independent per band (15/10/15%) in addition.
# Disk errors include an extra PB-model term 0.3*(1-PB) per band (PB from s55a). W HCN 4-3 (2.7 sigma) is dropped from the likelihood and tested as an upper limit afterwards.
# usage: python3 -I s55b_recovery_fit.py  -> out/q3b/s55b_recovery_fit.txt/json, s55b_posteriors.npz
import sys,json,itertools,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S
from scipy.special import logsumexp
L=S.L; F=S.Fine(); PBJ=json.load(open(S.OUT+'s55a_pb_snr.json'))['pb']
log=[]
def P(s=''): print(s,flush=True); log.append(s)
CENTRE=('N','S','Bridge'); DISK=('E1','E2','W')
GRID={'centre':dict(Lf=[0,.1,.2,.3,.4],de=[0,.1],k=[.5,1,1.5]),'disk':dict(Lf=[0,.15,.3,.45,.6],de=[0,.15,.3],k=[.5,1,1.5]),
      'centre_stress':dict(Lf=[0,.1,.2,.3,.4,.5],de=[0,.1,.2],k=[.5,1,2]),'disk_stress':dict(Lf=[0,.15,.3,.45,.6,.7],de=[0,.15,.3,.45],k=[.5,1,2]),'none':dict(Lf=[0],de=[0],k=[1])}
def combos(g): return [(a,b,c,d) for a in g['Lf'] for b in g['de'] for c in g['k'] for d in g['k']]
def extra_for(pos):
    if pos not in DISK: return None
    pb=np.array(PBJ[pos]); s=0.3*(1-pb); return np.r_[s,s]       # HCN 2,3,4 ; HCO+ 2,3,4 (bands 55,4f,47)
def fit(pos_list,cls,sig_b=(0.15,0.10,0.15),drop=None,grid=None,keep=False):
    """pos_list: positions sharing (n,Tk,N,A); each has its own eta; recovery combo shared across positions. returns Post, summary, combo posterior"""
    g=GRID[grid or cls]; cb=combos(g); engines=[];ys=[]
    for pos in pos_list:
        I,e=L.DATA[pos]; extra=extra_for(pos); rel=e.copy()
        if pos=='W' and drop: rel=rel.copy(); rel[drop]=5.0
        engines.append(S.Engine(F,S.covmat(rel,sig_b,extra=extra))); ys.append(np.log(I))
    tot=None; lnZc=[]; musum=None
    for (Lf,de,k47,k55) in cb:
        f=S.recovery_vec(Lf,de,k47,k55); lnL=0; mus=[]
        for En,y in zip(engines,ys):
            l,mu,chi=En.lnlike(y-np.log(f)); lnL=lnL+l; mus.append(mu)
        lnZc.append(logsumexp(lnL)-np.log(len(lnL)))
        w=np.exp(lnL-0)  # unnormalised (may under/overflow -> use shift)
        tot=lnL if tot is None else np.logaddexp(tot,lnL)
    tot=tot-np.log(len(cb)); lnZc=np.array(lnZc)
    # eta posterior (last position's eta at the best-fit combo is not needed): recompute mu for the first position at weighted combos
    Pst=S.Post(F,tot); d=Pst.summary(); 
    pc=np.exp(lnZc-lnZc.max()); pc/=pc.sum(); d['combo_post']={'Lf':{str(v):float(sum(p for p,c in zip(pc,cb) if c[0]==v)) for v in g['Lf']},'delta':{str(v):float(sum(p for p,c in zip(pc,cb) if c[1]==v)) for v in g['de']},
        'k47':{str(v):float(sum(p for p,c in zip(pc,cb) if c[2]==v)) for v in g['k']},'k55':{str(v):float(sum(p for p,c in zip(pc,cb) if c[3]==v)) for v in g['k']}}
    d['lnZ']=float(logsumexp(lnZc)-np.log(len(cb))); d['chi2min_combo_best']=None
    # predicted HCN 4-3 (W test): flux K km/s (observed-frame, i.e. before recovery) at posterior models with eta ~ eta_hat of the best combo
    return Pst,d,(cb,pc)
def line(d,keys=('log n','Tk','log A','log N(HCO+)/dv','log P/k')):
    out=[]
    for k in keys:
        v=d[k]; out.append(f"{k} {v[2]:.2f} [68%: {v[1]:.2f},{v[3]:.2f}; 95%: {v[0]:.2f},{v[4]:.2f}]")
    return '; '.join(out)
def pline(d): p=d['P']; return f"P(n>1e5)={p['n>1e5']:.2f} P(n<1e4)={p['n<1e4']:.2f} P(Tk>50)={p['Tk>50']:.2f} P(Tk<30)={p['Tk<30']:.2f} P(Tk<20)={p['Tk<20']:.2f} P(A>1)={p['A>1']:.2f} P(tau_HCN32>1)={p['tauHCN3-2>1']:.2f} | prior-edge mass: "+', '.join(f'{k}={v:.2f}' for k,v in d['edge'].items() if v>0.02)
if __name__=='__main__':
    res={}; maps={}
    P('# Observed per-position 10" ladders (K km/s) [HCN 2-1,3-2,4-3 | HCO+ 2-1,3-2,4-3]:')
    for pos in CENTRE+DISK: P(f'  {pos:6s} '+' '.join(f'{v:7.3f}' for v in L.DATA[pos][0]))
    scen=[('nominal (no recovery; cal 15/10/15%)',lambda pos: dict(grid='none')),
          ('one-sided recovery (centre Lf<=0.4; disk Lf<=0.6)',lambda pos: dict(grid=None)),
          ('stress recovery (centre Lf<=0.5, disk Lf<=0.7, k 0.5-2)',lambda pos: dict(grid=('centre_stress' if pos in CENTRE else 'disk_stress'))),
          ('ratios-only (all band scales free, no recovery)',lambda pos: dict(grid='none',sig_b=(10,10,10)))]
    for sn,fn in scen:
        P(f'\n## {sn}')
        for pos in CENTRE+DISK:
            cls='centre' if pos in CENTRE else 'disk'; kw=fn(pos)
            Pst,d,_=fit([pos],cls,drop=[2] if pos=='W' else None,**kw)
            res[f'{sn}|{pos}']=d; maps[f'{sn}|{pos}']=Pst.map2d()
            P(f'[{pos}] '+line(d)); P('      '+pline(d)); 
            if sn.startswith('one-sided'): P('      recovery posterior: '+'; '.join(f'{k}: '+','.join(f'{a}:{b:.2f}' for a,b in v.items()) for k,v in d['combo_post'].items()))
        kw=fn('E1'); Pst,d,_=fit(list(DISK),'disk',drop=[2],**kw); res[f'{sn}|disk joint']=d; maps[f'{sn}|disk joint']=Pst.map2d()
        P('[disk joint E1+E2+W, shared (n,Tk,N,A) and shared recovery combo] '+line(d)); P('      '+pline(d))
        Pst,d,_=fit(list(CENTRE),'centre',**fn('N')); res[f'{sn}|centre joint']=d; maps[f'{sn}|centre joint']=Pst.map2d()
        P('[centre joint N+S+Bridge] '+line(d)); P('      '+pline(d))
    json.dump(res,open(S.OUT+'s55b_recovery_fit.json','w'),indent=1,default=float)
    np.savez(S.OUT+'s55b_posteriors.npz',LT=S.LT,LN=S.LN,**{k.replace(' ','_').replace('|','__').replace('(','').replace(')','').replace(',','').replace('<=','le').replace('=','').replace('%','').replace('/','-').replace(';','')[:120]:v for k,v in maps.items()})
    open(S.OUT+'s55b_recovery_fit.txt','w').write('\n'.join(log))
