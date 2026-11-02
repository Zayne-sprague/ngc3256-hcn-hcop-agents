# Q42: add 13C/12C RATIO constraints to the Q40r2 one-zone RADEX fit (s235/s230/s220b machinery, monkeypatched; HCN, HCO+, CS 1.7" Harada lines unchanged).
# Ratio terms depend only on (Tk, n, N_12C, X): the filling factor and band scales cancel in a same-line 13C/12C ratio. N(13C)=N(12C)/X, same Tk, n, width; X = nuisance on the 0.05 dex grid of s235.
# Model ratio R = I(13C line)/I(12C line) in K km/s (RJ T_b). H13CN uses its own LAMDA file (h13cn@xpol.dat; grid out/q42/grid_h13cn.npz, s250a) unless H13CN_FILE=hcn (HCN file at N/X, the Q40 convention);
# H13CO+ uses the HCO+ file at N/X (no H13CO+ LAMDA file installed). 4-3: observed alpha = I(H13CN 4-3)/I(HCN 4-3) in Jy km/s -> T_b ratio = alpha*(354.5055/345.3398)^2.
# Variants (VAR): data terms (obs R, sigma, line index 0=1-0,1=3-2,2=4-3); error^2 = sigma^2 + (0.10 R_mod)^2 [LAMDA] + (beamsys R_obs)^2:
#   'H4c'  : Harada H13CO+ (s235 'H') + H13CN 4-3 CENSORED: R=0 +- 0.04 in alpha (the 3-sigma limit alpha<0.12 at N from Q4/p0007)
#   'H4g'  : same with the fitted (not censored) alpha = 0.045 +- 0.04 (envelope of the accepted wing models of s41)
#   'H4d'  : H4c with a dilution factor f=0.35: the 4.3x2.8" limit constrains the 1.7" core only through f*alpha_core <= alpha_beam, sigma_alpha/f
#   'T'    : common-taper (2.2" aperture, Q29 s185) 1-0 ratios H13CO+/HCO+ and H13CN/HCN at N, S, C (no absolute 13C intensities)
#   'Tp'   : same at the 2.2" PEAK PIXEL (+-150 km/s; s185 'pm150pk'); 'T150': 2.2" aperture +-150 km/s box ('pm150')
#   'T32'  : T + Harada 1.7" H13CO+/HCO+ 3-2 ratio at N (0.0284 +- 0.0055)
#   'T432' : T32 + H13CN 4-3 censored at N (S, C: fitted alpha with their wing-model envelope)
# env: XSTEP=1 -> X grid every 0.1 dex (1.35..1.95)+2.0; XEXT=1 -> extension X=112,141,178 (log X 2.05, 2.15, 2.25) beyond the s235 grid
# usage: python3 -I s250b_q42_ratio_fit.py POS VAR [k list e.g. 0,1]   -> out/q42/fit_<pos>_<var>_s<k>_lX<..>.npz
import sys,os,numpy as np,time
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s235_q40r2_Xmarg as S235
sys.argv=_a
from scipy.interpolate import CubicSpline
s=S235.s; m=s.m; d=s.d; NF=m.NF; LNC=m.LNC
O='/workspace/out/q42/'; os.makedirs(O,exist_ok=True)
G13=np.load(O+'grid_h13cn.npz') if os.path.exists(O+'grid_h13cn.npz') else None
H13FILE=os.environ.get('H13CN_FILE','xpol')
NU13,NU12=345.3398,354.5055; F43=(NU12/NU13)**2
def lnTR(a,li): return np.log(np.maximum(np.nan_to_num(a[...,li,0],nan=1e-3),1e-3))
SPL={}
def spl(a,li,key):
    if key not in SPL: SPL[key]=CubicSpline(LNC,lnTR(a,li),axis=2,extrapolate=True)
    return SPL[key]
def ratio_model(sp,li,logX):
    a12=d['hcn'] if sp=='H' else d['hcop']
    a13=a12
    if sp=='H' and H13FILE=='xpol': a13=G13['h13cn']
    l12=spl(a12,li,(sp,li,'12'))(NF); l13=spl(a13,li,(sp,li,'13'+('x' if a13 is not a12 else 'c')))(np.clip(NF-logX,11.5,19))
    R=np.clip(np.exp(l13-l12),1e-8,1e3)
    return R*(F43 if (sp=='H' and li==2) else 1.0)
# ---- data: sp, li, Robs, sigma, beamsys (fractional extra systematic), censored flag handled by Robs=0
DATA={
 'N':dict(T=[('O',0,0.0231,0.0031),('H',0,0.0207,0.0034)],
          T32=[('O',1,2.6/91.56,0.0284*0.19)],
          H4c=[('H',2,0.0,0.04)],H4g=[('H',2,0.045,0.04)],H4d=[('H',2,0.0,0.04/0.35)],
          T4=[('H',2,0.0,0.04)],
          Tp=[('O',0,0.0372,0.0066),('H',0,0.0302,0.0058)],T150=[('O',0,0.0271,0.0049),('H',0,0.0220,0.0046)]),
 'S':dict(T=[('O',0,0.0184,0.0040),('H',0,0.0057,0.0059)],T4=[('H',2,0.01,0.05)],
          Tp=[('O',0,0.0141,0.0088),('H',0,-0.0006,0.0130)],T150=[('O',0,0.0150,0.0047),('H',0,0.0048,0.0069)]),
 'C':dict(T=[('O',0,0.0185,0.0033),('H',0,0.0237,0.0044)],T4=[('H',2,0.03,0.07)],
          Tp=[('O',0,0.0222,0.0087),('H',0,0.0265,0.0117)],T150=[('O',0,0.0182,0.0045),('H',0,0.0234,0.0058)]),
}
def terms(pos,var):
    D=DATA[pos]; t=[]
    if var in('H4c','H4g','H4d'): return D[var]
    if var in('T','Tp','T150'): return D[var]
    if var=='T32': return D['T']+D.get('T32',[])
    if var=='T432': return D['T']+D.get('T32',[])+D['T4']
    raise KeyError(var)
CUR={'terms':[],'logX':1.6,'bs':0.0}
RAT={}
base=m.chi_s
def wrapped(it,jn,sp,o3=None,o6=None):
    f=base(it,jn,sp,o3,o6)
    if f is None: return f
    for (spc,li,Ro,sg) in CUR['terms']:
        if spc!=sp: continue
        R=RAT[(spc,li)][it,jn]
        s2=sg**2+(0.10*R)**2
        f=f+((R-Ro)**2/s2)[None,:]
    return f
m.chi_s=wrapped
def run(pos,var,k,lx):
    CUR['terms']=terms(pos,var); CUR['logX']=lx; RAT.clear()
    for (spc,li,Ro,sg) in CUR['terms']: RAT[(spc,li)]=ratio_model(spc,li,lx)
    sg=S235.SIGS[k]; basevar='H' if var.startswith('H4') else 'N'
    if basevar=='H': 
        q=s.iso('H')[pos]; s.X_iso={'H':q['H'],'O':q['O'],'C':[]}
    return s.run_one(pos,lx,basevar,sg)
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]; ks=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [0,1]
    for pos in poss:
      for k in ks:
        for lx in ([2.05,2.15,2.25] if os.environ.get('XEXT') else (S235.LOGX if not os.environ.get('XSTEP') else list(S235.LOGX[::2])+[S235.LOGX[-1]])):
            fn=O+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
            if os.path.exists(fn): continue
            t=time.time(); o=run(pos,var,k,lx); np.savez(fn,**o); cm=o['chimin']; Tk=m.TK
            print(pos,var,S235.SIGS[k],lx,f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
