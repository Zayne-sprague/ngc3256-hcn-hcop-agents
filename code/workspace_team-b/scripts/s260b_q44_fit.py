# Q44: fully matched one-zone RADEX fit: ALL data at ONE 2.2" circular beam, r=2.2" aperture means, +-150 km/s box (SiO 2-1 outside), no Harada peak-pixel inputs.
#  12C: HCN, HCO+ 1-0 (Q29 s185 'pm150' reference integrals; X961, short-spacing complete) and 3-2 (Q34 s192 '2p2|comb|*|h150' 12m+7m uv-matched taper); all x(1+z)=1.0093 (rest frame) as s193.
#  13C: model R = I(13C)/I(12C) 1-0 (H13CO+/HCO+, H13CN/HCN) from Q29 s185_box 'pm150' (s250b 'T150' terms, 8% cross-project scale included); N(13C)=N(12C)/X, X nuisance (grid 0.05 dex, 22..100 or ..178).
#  Machinery = s250b (s235/s230/s220b, one-zone Tk,n,N_HCN,N_HCO+, filling factor eta<=1, band scales z3,z6 with priors) with m.g monkeypatched; error per 12C line = 1% + 5% floor (s220b).
# Variants (VAR):
#  P   : 4 x 12C (1-0, 3-2 of HCN, HCO+) + 13C 1-0 ratios                      [primary, no CS]
#  P10 : 1-0 only (3-2 lines masked) + 13C 1-0 ratios
#  PCS : P + Harada 1.7" CS 2-1, 5-4 (peak pixel; NOT matched in aperture; shown for comparison only)
#  P10C: P10 + Harada CS 2-1 only (same caveat)
#  Pw  : P with +-350 km/s 12C integrals (s193) and the s185 template-aperture 13C ratios (|v|<500)
#  Pp  : P with the 1-0 peak pixel (r=0.3", s185 pm150pk) 12C and 13C but 3-2 aperture 12C (3-2 peak pixel not stored) -> mixed: reported as a check only
# usage: python3 -I s260b_q44_fit.py POS[,POS] VAR [k list] ; env XSTEP=1 (X every 0.1 dex), XEXT=1 (X=112,141,178) -> out/q44/fit_<pos>_<var>_s<k>_lX<..>.npz
import sys,os,json,time,numpy as np
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s250b_q42_ratio_fit as B
sys.argv=_a
S235=B.S235; s=B.s; m=B.m; X=m.X
O='/workspace/out/q44/'; os.makedirs(O,exist_ok=True)
K=1.0093
# 1-0 +-150 aperture (s185 box refs), template aperture +-350 (s193: includes K already)  ; 3-2 s192
I10={'N':(86.7,103.9),'S':(50.0,86.1),'C':(70.8,99.2)}                # HCN, HCO+  r=2.2" +-150 km/s (K km/s, s185.txt 'pm150 ... / HCN ref')
I10pk={'N':(145.5,160.6),'S':(63.2,108.6),'C':(81.2,117.6)}
M32=json.load(open('/workspace/out/q34/s192.json'))['res']
def i32(pos,key):
    r=M32['2p2|comb|'+{'N':'N','S':'S','C':'Bridge'}[pos]]; return r['HCN32|ap|'+key][0],r['HCOp32|ap|'+key][0]
I10w={'N':(91.7/K,109.6/K),'S':(52.9/K,92.8/K),'C':(75.6/K,106.3/K)}   # s193 printout /K  (w350)
VAR={'P':dict(sp=['H','O'],lines=[0,1,2,3],ap='h150',r10=I10,tm='T150'),
     'P10':dict(sp=['H','O'],lines=[0,2],ap='h150',r10=I10,tm='T150'),
     'PCS':dict(sp=['H','O','C'],lines=[0,1,2,3,4,5],ap='h150',r10=I10,tm='T150'),
     'P10C':dict(sp=['H','O','C'],lines=[0,2,4],ap='h150',r10=I10,tm='T150'),
     'Pw':dict(sp=['H','O'],lines=[0,1,2,3],ap='w350',r10=I10w,tm='T'),
     'Pp':dict(sp=['H','O'],lines=[0,1,2,3],ap='h150',r10=I10pk,tm='Tp')}
# 13C ratio terms: reuse s250b DATA keys (T150 / T / Tp) -> B.terms
CFG={'var':'P','pos':'N'}
g0=m.g
def g(p,l):
    v=VAR[CFG['var']]; i=m.LINES.index(l)
    if l.startswith('HCN') or l.startswith('HCO'):
        h10,o10=v['r10'][p]; h32,o32=i32(p,v['ap'])
        val={'HCN1':h10,'HCN3':h32,'HCO+1':o10,'HCO+3':o32}[l]*K
        return val,0.01*val
    return g0(p,l)
m.g=g
setup0=m.setup
def setup2(pos,s3,s6,species):
    v=VAR[CFG['var']]; setup0(pos,s3,s6,v['sp'])
    for i in range(len(m.LINES)):
        if i not in v['lines']: X.v[i]=1e8
m.setup=setup2
from multiprocessing import Pool
def run_one(pos,logX,variant,sig):
    s.X_iso={'H':[],'O':[],'C':[]}
    s.MOD13['H']=s.lnI13(s.d['hcn'],logX); s.MOD13['O']=s.lnI13(s.d['hcop'],logX)
    cells=[(i,j) for i in range(len(m.TK)) for j in range(len(m.LN))]
    m.setup(CFG['pos'],sig[0],sig[1],['H','O','C'])
    with Pool(min(16,os.cpu_count())) as P:
        mm=np.array(P.map(m.cell_min,cells,chunksize=20)); cb=mm.min()
        res=P.map(m.cell_pass2,[(c,cb,45.) for c in cells],chunksize=20)
    Wc=np.array([r['Wc'] if r else 0 for r in res]).reshape(len(m.TK),len(m.LN))
    out=dict(chibest=cb,Wc=Wc,chimin=mm.reshape(len(m.TK),len(m.LN)))
    for key in m.LAGS:
        hs=[r['hist'][key] for r in res if r and key in r['hist']]
        if hs: out['h_'+key]=np.sum(hs,0); out['p_'+key]=np.min([r['prof'][key] for r in res if r and key in r['prof']],0)
    return out
def run(pos,var,k,lx):
    CFG['var']=var; CFG['pos']=pos
    B.CUR['terms']=B.terms(pos,VAR[var]['tm']); B.CUR['logX']=lx; B.RAT.clear()
    for (spc,li,Ro,sg) in B.CUR['terms']: B.RAT[(spc,li)]=B.ratio_model(spc,li,lx)
    return run_one(pos,lx,'N',S235.SIGS[k])
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]; ks=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [0,1]
    lxs=[2.05,2.15,2.25] if os.environ.get('XEXT') else (list(S235.LOGX[::2])+[S235.LOGX[-1]] if os.environ.get('XSTEP') else list(S235.LOGX))
    for pos in poss:
      for k in ks:
        for lx in lxs:
            fn=O+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
            if os.path.exists(fn): continue
            t=time.time(); o=run(pos,var,k,lx); np.savez(fn,**o); cm=o['chimin']; Tk=m.TK
            print(pos,var,S235.SIGS[k],lx,f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
