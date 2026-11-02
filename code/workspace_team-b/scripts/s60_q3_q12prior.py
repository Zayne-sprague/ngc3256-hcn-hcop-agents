# Q3 round 3: single-zone RADEX ladder fits at 10" with the Q12 (p0008) recovery information as a relative-recovery budget, not a correction.
# Recovery enters as multipliers (true/observed) on the J=2-1 and J=4-3 fluxes relative to J=3-2 (Xb4f); a common loss of all bands is degenerate with the filling eta (<=1) and drops out.
#   centre N,S,Bridge (Q12 R55/R4f 1.05-1.18, R47/R4f 0.89-0.97):  f2 = true/obs(2-1) in U[0.85,0.95], f4 in U[1.03,1.12]
#   E1 (compact, Xb4d/Xb4f=1.01, R55,R47 ~ 1):                      f2,f4 in U[0.95,1.05]
#   E2, W (unconstrained in Q12; W HCN not usable):                 wide  f2 in U[0.6,1.0], f4 in U[1.0,1.5]  (reported as unconstrained)
# Calibration: 10%/band independent (Xb55,Xb4f,Xb47), +7% per line; disk PB term per line = max(0.3(1-PB), 0.5|1-A/B|) (A/B = smooth(pbcor) vs smooth(pbcor*PB)/PB, s55a). W HCN 4-3 (2.7 sigma) dropped.
# usage: python3 -I s60_q3_q12prior.py -> out/q3c/s60_q12prior.txt/json, s60_posteriors.npz
import sys,json,numpy as np,os
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S
from scipy.special import logsumexp
L=S.L; F=S.Fine(); J=json.load(open('/workspace/out/q3b/s55a_pb_snr.json')); PBJ=J['pb']; PT=J['pbtest']
OUT='/workspace/out/q3c/'; os.makedirs(OUT,exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
CENTRE=('N','S','Bridge'); DISK=('E1','E2','W')
def g3(a,b): return np.linspace(a,b,4)
PRI={'q12':{'centre':(g3(.85,.95),g3(1.03,1.12)),'E1':(g3(.95,1.05),g3(.95,1.05)),'E2':(g3(.6,1),g3(1,1.5)),'W':(g3(.6,1),g3(1,1.5))},
     'none':{k:(np.array([1.]),np.array([1.])) for k in ('centre','E1','E2','W')},
     'wide':{k:(g3(.6,1),g3(1,1.5)) for k in ('centre','E1','E2','W')}}
def cls(pos): return 'centre' if pos in CENTRE else pos
def extra(pos,pbmult=1.):
    if pos not in DISK: return None
    pb=np.array(PBJ[pos]); s=0.3*(1-pb)
    ab=[PT[f'{b}|{pos}'] for b in ('Xb55 (B5, J=2)','Xb4f (B6, J=3)','Xb47 (B7, J=4)')]   # [HCN,HCO+] per band
    sh=np.maximum(s,0.5*np.abs(1-np.array([a[0] for a in ab]))); so=np.maximum(s,0.5*np.abs(1-np.array([a[1] for a in ab])))
    return pbmult*np.r_[sh,so]
def fit(pos_list,prior='q12',sig_b=(0.10,0.10,0.10),pbmult=1.,drop=None):
    en=[];ys=[]
    for pos in pos_list:
        I,e=L.DATA[pos]; rel=e.copy()
        if pos=='W': rel[2]=5.0       # drop HCN 4-3
        en.append(S.Engine(F,S.covmat(rel,sig_b,extra=extra(pos,pbmult)))); ys.append(np.log(I))
    # shared correction across the positions in the list (use class of first)
    f2g,f4g=PRI[prior][cls(pos_list[0])]; tot=None
    for f2 in f2g:
        for f4 in f4g:
            f=np.array([f2,1,f4,f2,1,f4]); lnL=0
            for En,y in zip(en,ys): l,mu,chi=En.lnlike(y+np.log(f)); lnL=lnL+l
            tot=lnL if tot is None else np.logaddexp(tot,lnL)
    Pst=S.Post(F,tot-np.log(len(f2g)*len(f4g))); return Pst,Pst.summary()
def line(d):
    return '; '.join(f"{k} {d[k][2]:.2f} [{d[k][1]:.2f},{d[k][3]:.2f}; 95%: {d[k][0]:.2f},{d[k][4]:.2f}]" for k in ('log n','Tk','log A','log N(HCO+)/dv','log P/k'))
def pline(d):
    p=d['P']; return f"P(n>1e5)={p['n>1e5']:.2f} P(n<1e4)={p['n<1e4']:.2f} P(Tk>50)={p['Tk>50']:.2f} P(Tk<30)={p['Tk<30']:.2f} P(A>1)={p['A>1']:.2f} P(A>3)={p['A>3']:.2f} P(tau_HCN32>1)={p['tauHCN3-2>1']:.2f} tau(HCN3-2) med {10**d['tau(HCN3-2)'][2]:.1f} [{10**d['tau(HCN3-2)'][1]:.1f},{10**d['tau(HCN3-2)'][3]:.1f}]"
def pgt(a,b,x):    # P(a>b) for two independent marginals on grid x
    return float(sum(a[i]*b[:i].sum() for i in range(len(x)))+0.5*np.sum(a*b))
if __name__=='__main__':
    res={};maps={}
    P('# Ladder ratios (K km/s, RJ, 10" beam) and their budgeted systematic: stat is <1% at the centre; band cal 10%/band (indep) + Q12 recovery factor on T21/T32 [0.85-0.95] and T43/T32 [1.03-1.12] (centre)')
    for pos in CENTRE+DISK:
        I=L.DATA[pos][0]; P(f'  {pos:6s} T21/T32 HCN {I[0]/I[1]:.2f} HCO+ {I[3]/I[4]:.2f} | T43/T32 HCN {I[2]/I[1]:.3f} HCO+ {I[5]/I[4]:.3f} | R_K(2,3,4)={I[0]/I[3]:.3f},{I[1]/I[4]:.3f},{I[2]/I[5]:.3f}')
    scen=[('Q12-informed recovery budget, cal 10/10/10%','q12',(.1,.1,.1),1.),('no recovery term, cal 10/10/10%','none',(.1,.1,.1),1.),
          ('Q12-informed, cal 15/10/15%','q12',(.15,.1,.15),1.),('wide recovery for all (f2 .6-1, f4 1-1.5), cal 10%','wide',(.1,.1,.1),1.),
          ('Q12-informed, cal 10%, PB term x2','q12',(.1,.1,.1),2.),('ratios-only (band scales free)','none',(10,10,10),1.),
          ('Q12-informed, Xb55 (2-1) dropped (free scale)','q12',(10,.1,.1),1.),('Q12-informed, Xb47 (4-3) dropped (free scale)','q12',(.1,.1,10),1.)]
    for sn,pr,sb,pm in scen:
        P(f'\n## {sn}')
        for pos in CENTRE+DISK:
            Pst,d=fit([pos],pr,sb,pm); res[f'{sn}|{pos}']=d; maps[f'{sn}|{pos}']=Pst.map2d()
            P(f'[{pos}] '+line(d)); P('      '+pline(d))
        for nm,lst in (('centre joint N+S+Bridge',CENTRE),('E1 alone vs ',None)):
            if lst is None: continue
            Pst,d=fit(list(lst),pr,sb,pm); res[f'{sn}|{nm}']=d; maps[f'{sn}|{nm}']=Pst.map2d(); P(f'[{nm}] '+line(d)); P('      '+pline(d))
        # contrasts (independent marginals; calibration/recovery common to positions is NOT cancelled here: conservative)
        LT,LN=S.LT,S.LN; mc=maps[f'{sn}|centre joint N+S+Bridge']
        for pos in ('E1','E2','W'):
            m=maps[f'{sn}|{pos}']
            P(f'   {pos} vs centre: P(Tk_{pos}<Tk_centre)={pgt(mc.sum(1),m.sum(1),LT):.3f}  P(n_{pos}>n_centre)={pgt(m.sum(0),mc.sum(0),LN):.3f}')
        for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
            ma,mb=maps[f'{sn}|{a}'],maps[f'{sn}|{b}']
            P(f'   {a} vs {b}: P(Tk_{a}>Tk_{b})={pgt(ma.sum(1),mb.sum(1),LT):.3f}  P(n_{a}>n_{b})={pgt(ma.sum(0),mb.sum(0),LN):.3f}  (posteriors are near-identical: the data share the same prior-dominated n,Tk degeneracy)')
    json.dump(res,open(OUT+'s60_q12prior.json','w'),indent=1,default=float)
    np.savez(OUT+'s60_posteriors.npz',LT=S.LT,LN=S.LN,**{('k%03d'%i):v for i,v in enumerate(maps.values())},keys=np.array(list(maps.keys())))
    open(OUT+'s60_q12prior.txt','w').write('\n'.join(log))
