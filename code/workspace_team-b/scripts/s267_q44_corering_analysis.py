# Q44 round 2: analysis of core/ring fits (s266), same as s262.
# Q44 analysis of s261 fits (fully matched 2.2" fit); adapted from s250c (Q42). Original header: Q42 analysis of s250b fits (+ s235 baseline 'H' for N): per variant/position/band-prior:
#  lnZ(X); posterior of X=12C/13C (flat in log X on 22-100, and LN(45,0.2dex)); branch posteriors P(Tk<=30), P(30-100), P(>=100); BF hot:cold (prior-mass corrected);
#  N(HCN)/N(HCO+) marginal; dchi2(cold-hot) at fixed X; and the X contrast S/N from the 'T' variants.
# usage: python3 -I s250c_q42_analysis.py -> out/q42/analysis_q44_corering.txt, analysis_q44_corering.json
import sys,os,json,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s235_q40r2_Xmarg as S
m=S.s.m; TK=m.TK; LOGX=S.LOGX; SIGS=S.SIGS
OQ='/workspace/out/q44/'; O40='/workspace/out/q44/'
L=[]
def P(x=''): print(x); L.append(x)
def path(pos,var,k,lx):
    if var=='H0': return O40+f'fit_{pos}_H_s{k}_lX{int(round(lx*100))}.npz'
    if var=='A2': return O40+f'fit_{pos}_A2_s{k}_lX{int(round(lx*100))}.npz'
    return OQ+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
LOGXE=np.round(np.arange(1.35,2.2501,0.05),2)
def loadset(pos,var,k,xmax=2.0001):
    R={}
    for lx in LOGXE[LOGXE<=xmax]:
        f=path(pos,var,k,lx)
        if os.path.exists(f): R[lx]=np.load(f)
    return R
def lnZ(r): return -0.5*float(r['nsp'])*r['chibest']+np.log(r['Wc'].sum())
ST=0.05; lag=np.arange(-20,29)
def quant(h,qs=(0.16,0.5,0.84),lg=lag):
    x=np.concatenate([[lg[0]-0.5],lg+0.5])*ST; c=np.concatenate([[0],np.cumsum(h)/h.sum()]); return [10**float(np.interp(q,c,x)) for q in qs]
PRI={'flat logX':None,'LN(45,0.2dex)':(np.log10(45),0.2)}
OUT={}
def analyse(pos,var,k,xmax=2.0001):
    R=loadset(pos,var,k,xmax)
    if len(R)<7: return None
    lxs=np.array(sorted(R)); lz=np.array([lnZ(R[a]) for a in lxs]); res={'lnZ':dict(zip([f'{10**a:.0f}' for a in lxs],lz.round(2).tolist()))}
    for nm,pr in PRI.items():
        dx=np.gradient(lxs)                      # trapezoid weights: the X grid is 0.05 or 0.1 dex and not always uniform
        lp=np.log(dx)+(0 if pr is None else -0.5*((lxs-pr[0])/pr[1])**2)
        w=np.exp(lz+lp-(lz+lp).max()); w/=w.sum(); c=np.cumsum(w)
        xq=[10**float(np.interp(q,c,lxs)) for q in (0.16,0.5,0.84)]
        Zb=np.array([R[a]['Wc'].sum(1)*np.exp(-1.5*R[a]['chibest']-lz.max()) for a in lxs]); pw=np.exp(lp)[:,None]*Zb; tot=pw.sum()
        br={'cold<=30':pw[:,TK<=30].sum()/tot,'30-100':pw[:,(TK>30)&(TK<100)].sum()/tot,'>=100':pw[:,TK>=100].sum()/tot,'<=150':pw[:,TK<=150].sum()/tot}
        bf=br['>=100']/br['cold<=30']*(np.sum(TK<=30)/np.sum(TK>=100))
        hs=np.array([R[a]['h_A']/R[a]['h_A'].sum() for a in lxs]); hA=(w[:,None]*hs).sum(0); q=quant(hA)
        WT=(pw.sum(0) if False else None)
        Wtk=sum(np.exp(lp[i])*R[a]['Wc'].sum(1)*np.exp(-0.5*float(R[a]['nsp'])*R[a]['chibest']-lz.max()+0) for i,a in enumerate(lxs)); Wtk=Wtk/Wtk.sum(); ctk=np.cumsum(Wtk)
        Wn=sum(np.exp(lp[i])*R[a]['Wc'].sum(0)*np.exp(-0.5*float(R[a]['nsp'])*R[a]['chibest']-lz.max()) for i,a in enumerate(lxs)); Wn=Wn/Wn.sum(); cn=np.cumsum(Wn)
        tkq=[float(np.interp(q,ctk,TK)) for q in (0.16,0.5,0.84)]; nq=[float(np.interp(q,cn,m.LN)) for q in (0.16,0.5,0.84)]
        res[nm]=dict(Tk=tkq,logn=nq,X=xq,branch=br,BF_hot_cold=bf,A=q,w=w.tolist(),lxs=lxs.tolist())
    # fixed-X dchi2
    res['dchi2_cold_hot']={f'{10**a:.0f}':float(R[a]['chimin'][TK<=30].min()-R[a]['chimin'][TK>=100].min()) for a in lxs[::3]}
    res['chi2min_bestX']=float(min(R[a]['chibest'] for a in lxs))
    return res
def show(pos,var,k,r):
    for nm in PRI:
        d=r[nm]; b=d['branch']
        P(f'     Tk {d["Tk"][1]:.0f} [{d["Tk"][0]:.0f},{d["Tk"][2]:.0f}] log n {d["logn"][1]:.2f} [{d["logn"][0]:.2f},{d["logn"][2]:.2f}]'); P(f'  {pos:2s} {var:5s} bandprior {SIGS[k]} Xprior {nm:14s}: X={d["X"][1]:.0f} [{d["X"][0]:.0f},{d["X"][2]:.0f}]  P(Tk<=30)={b["cold<=30"]:.4f} P(30-100)={b["30-100"]:.3f} P(>=100)={b["≥100" if False else ">=100"]:.3f} P(<=150)={b["<=150"]:.2f}  BF(hot:cold)={d["BF_hot_cold"]:.0f}  N(HCN)/N(HCO+)={d["A"][1]:.2f} [{d["A"][0]:.2f},{d["A"][2]:.2f}]')
    P('      dchi2(cold-hot) at X: '+' '.join(f'{x}:{v:.1f}' for x,v in r['dchi2_cold_hot'].items())+f'; chi2min(best X) {r["chi2min_bestX"]:.2f}')
jobs=[(p,v) for v in ('M','M32') for p in ('Nc','Nr','Sc','Sr','Cc','Cr')]
for pos,var in jobs:
    ks=[0,1]
    for k in ks:
        r=analyse(pos,var,k)
        if r is None: continue
        OUT[f'{pos}_{var}_{k}']=r; show(pos,var,k,r)
        if True:
            re=analyse(pos,var,k,2.2501)
            if re is not None and len(re['lnZ'])>len(r['lnZ']):
                OUT[f'{pos}_{var}_{k}_ext']=re; P('   [X grid extended to 178]'); show(pos,var,k,re)
open(OQ+'analysis_q44_corering.txt','w').write('\n'.join(L)); json.dump(OUT,open(OQ+'analysis_q44_corering.json','w'),default=float)
