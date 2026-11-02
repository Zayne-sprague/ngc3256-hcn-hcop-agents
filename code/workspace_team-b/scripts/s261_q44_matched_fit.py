# Q44: FULLY MATCHED one-zone RADEX fit: every datum at ONE 2.2" circular beam, r<2.2" aperture mean, +-150 km/s (Harada-like window; avoids SiO 2-1 at +320 km/s from H13CO+ 1-0), NO Harada peak-pixel input.
# Machinery = s250b/s235/s230/s220b (imported, monkeypatched): same RADEX grid (Tk 10-300, log n 3-8, 0.05 dex N grids, eta<=1 flat, band scales s3 (1-0) s6 (3-2) shared by all lines), 13C species at N/X, same Tk,n,eta,dv=150.
# Data vector (K km/s, x(1+z)=1.00927 as s193):  HCN 1-0, HCO+ 1-0: Q29 s183 X961 2.2" h150_ap;  HCN 3-2, HCO+ 3-2: Q34 s192 12m+7m 'comb' 2.2" h150 ap;
#   H13CN 1-0, H13CO+ 1-0: Q29 s185_box 'pm150' (absolute amplitudes, null-calibrated errors);  CS 2-1, 5-4: only Harada 1.7" peak exists -> scaled by the Q44 conversion factor (2.2" aperture / Harada) measured on HCN and HCO+
#   at the same J (mean, +-25% systematic added in quadrature), flagged as the one non-matched input (variants with 'C').
# Variants (VAR):  M12 : 12C only (4 lines);  M : + H13CO+ 1-0 and H13CN 1-0 (PRIMARY);  M10 : 1-0 12C lines + 13C 1-0 only (no 3-2);  MC : M + scaled CS;  M10C : M10 + scaled CS 2-1 only;
#   M4 : M + H13CN 4-3 censored ratio at N (alpha<0.12, 4.3x2.8" Xb47; ratio term, beam-mismatched, s250b convention);  M32 : 3-2 12C lines + 13C 1-0 (no 12C 1-0 lines).
# H13CN uses its own LAMDA grid (out/q42/grid_h13cn.npz) unless H13CN_FILE=hcn.  usage: python3 -I s261_q44_matched_fit.py POS VAR [k list] ; env XSTEP=1 (X every 0.1 dex), XEXT=1 (X=112,141,178)
import sys,os,json,time,numpy as np
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s250b_q42_ratio_fit as R
sys.argv=_a
S235=R.S235; s=R.s; m=R.m; d=R.d; NF=R.NF
from multiprocessing import Pool
from scipy.interpolate import CubicSpline
O='/workspace/out/q44/'; os.makedirs(O,exist_ok=True)
K=1.00927
r183=json.load(open('/workspace/out/q29/s183.json')); r192=json.load(open('/workspace/out/q34/s192.json'))['res']; box=json.load(open('/workspace/out/q29/s185_box.json'))
HAR={'N':dict(HCN1=(197.23,.96),HCO1=(205.25,.95),HCN3=(71.38,.63),HCO3=(91.56,.62),CS2=(48.99,.64),CS5=(15.45,.63)),
     'S':dict(HCN1=(56.68,.96),HCO1=(87.87,.95),HCN3=(17.72,.63),HCO3=(43.69,.62),CS2=(22.92,.64),CS5=(4.62,.63)),
     'C':dict(HCN1=(93.15,.96),HCO1=(139.46,.94),HCN3=(19.24,.61),HCO3=(44.61,.60),CS2=(29.89,.64),CS5=(5.37,.62))}
RG={'N':'N','S':'S','C':'Bridge'}
def vec(pos):
    g=RG[pos]; a=r183['2015.1.00993 X961|2.2|'+g]; c=r192['2p2|comb|'+g]
    I=np.array([a['HCN']['h150_ap'][0],c['HCN32|ap|h150'][0],a['HCOp']['h150_ap'][0],c['HCOp32|ap|h150'][0]])*K
    e=np.array([a['HCN']['h150_ap'][1],c['HCN32|ap|h150'][1],a['HCOp']['h150_ap'][1],c['HCOp32|ap|h150'][1]])*K
    h=HAR[pos]
    f1=np.mean([I[0]/h['HCN1'][0],I[2]/h['HCO1'][0]]); f3=np.mean([I[1]/h['HCN3'][0],I[3]/h['HCO3'][0]])
    cs=[h['CS2'][0]*f1,h['CS5'][0]*f3]; ce=[np.hypot(h['CS2'][1]*f1,0.25*cs[0]),np.hypot(h['CS5'][1]*f3,0.25*cs[1])]
    return I,e,cs,ce,(f1,f3)
def iso13(pos):
    g=RG[pos]; hc=box[f'pm150|{g}|H13COp']; hn=box[f'pm150|{g}|H13CN']
    return dict(O=[(0,hc['amp']*K,hc['err']*K)],H=[(0,hn['amp']*K,hn['err']*K)])
VARS={'M12':dict(lines=(0,1,2,3),iso=False),'M':dict(lines=(0,1,2,3),iso=True),'M10':dict(lines=(0,2),iso=True),'MC':dict(lines=(0,1,2,3,4,5),iso=True),
      'M10C':dict(lines=(0,2,4),iso=True),'M4':dict(lines=(0,1,2,3),iso=True,r43=True),'M32':dict(lines=(1,3),iso=True),'M12C':dict(lines=(0,1,2,3,4,5),iso=False),'M3':dict(lines=(1,3),iso=False)}
def setup(pos,var,s3,s6):
    I,e,cs,ce,_=vec(pos); Iall=np.concatenate([I,cs]); eall=np.concatenate([e,ce])
    use=np.zeros(6,bool); use[list(VARS[var]['lines'])]=True
    m.X.y=np.log(Iall); v=(eall/Iall)**2+0.05**2; v[~use]=1e8; m.X.v=v
    m.X.o3,m.X.o6,m.X.pi=m.nuis(s3,s6)
    m.X.sp=['H','O']+(['C'] if use[4] or use[5] else [])
    q=iso13(pos) if VARS[var]['iso'] else dict(O=[],H=[])
    s.X_iso={'H':q['H'],'O':q['O'],'C':[]}
def lnI13(a,logX):
    TR=np.nan_to_num(a[...,:2,0],nan=1e-3); l=np.log(1.0645*150*np.maximum(TR,1e-3))
    return CubicSpline(m.LNC,l,axis=2,extrapolate=True)(np.clip(NF-logX,11.5,19))
def run(pos,var,k,lx):
    sg=S235.SIGS[k]; setup(pos,var,sg[0],sg[1])
    s.MOD13['H']=lnI13(R.G13['h13cn'] if R.H13FILE=='xpol' else d['hcn'],lx); s.MOD13['O']=lnI13(d['hcop'],lx)
    R.CUR['terms']=[]; R.RAT.clear()
    if VARS[var].get('r43') and pos=='N':
        R.CUR['terms']=[('H',2,0.0,0.04)]; R.RAT[('H',2)]=R.ratio_model('H',2,lx)
    cells=[(i,j) for i in range(len(m.TK)) for j in range(len(m.LN))]
    with Pool(min(16,os.cpu_count())) as P:
        mm=np.array(P.map(m.cell_min,cells,chunksize=20)); cb=mm.min()
        res=P.map(m.cell_pass2,[(c,cb,45.) for c in cells],chunksize=20)
    Wc=np.array([r['Wc'] if r else 0 for r in res]).reshape(len(m.TK),len(m.LN))
    out=dict(chibest=cb,Wc=Wc,chimin=mm.reshape(len(m.TK),len(m.LN)),nsp=len(m.X.sp))
    for key in m.LAGS:
        hs=[r['hist'][key] for r in res if r and key in r['hist']]
        if hs: out['h_'+key]=np.sum(hs,0); out['p_'+key]=np.min([r['prof'][key] for r in res if r and key in r['prof']],0)
    return out
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]; ks=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [0,1]
    if os.environ.get('SHOW'):
        for p in poss:
            I,e,cs,ce,f=vec(p); print(p,'12C',I.round(2),e.round(2),'CS',np.round(cs,2),np.round(ce,2),'conv f1,f3',np.round(f,3),'13C',iso13(p))
        sys.exit()
    tag=os.environ.get('TAGX','')
    for pos in poss:
      for k in ks:
        LX=([1.6] if os.environ.get('LXONLY') else [2.05,2.15,2.25] if os.environ.get('XEXT') else (S235.LOGX if not os.environ.get('XSTEP') else list(S235.LOGX[::2])+[S235.LOGX[-1]]))
        for lx in LX:
            fn=O+f'fit_{pos}_{var}{tag}_s{k}_lX{int(round(lx*100))}.npz'
            if os.path.exists(fn): continue
            t=time.time(); o=run(pos,var,k,lx); np.savez(fn,**o); cm=o['chimin']; Tk=m.TK
            print(pos,var,S235.SIGS[k],lx,f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
