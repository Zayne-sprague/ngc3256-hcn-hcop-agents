# Q27 round 2 (reviewer fixes). Same data/grid/calibration as s104 (s62_lib Grid + band 10%, 7% per line, Q12 centre recovery), positions N, S, Bridge, 10" beam.
# CORRECTION found while writing this: F.tauH/F.tauO have J index 0=1-0 (grid file axis J=1..4), so s62/s64/s66/s104 'tau(HCN 3-2)' = F.tauH[:,1] is tau(2-1), and
#   the s104 H13CN limit used F.tauH[:,2] (=3-2) instead of the 4-3 line (F.tauH[:,3]). Here: tau32=tauH[:,2], tau43=tauH[:,3].
# Parts: A floor sensitivity (line width w 100-250 km/s: T_R,peak=exp(M)/(1.0645 w); factor f 0.7-1.0; 49 pc / 250 pc T_b from p0026 s160.json)
#        B separate eta_HCN, eta_HCO+ (flat ln eta [1e-3,1] each, optionally |ln q|<=ln10), X(HCN/H13CN) 30/50/100 (limit tau_H13CN(4-3)<0.128 -> tau_HCN(4-3)< X*0.128)
#        C Tk<=30 K / Tk=25 K (p0012 S 1" core) in the 10" aperture: profile chi2 and posterior odds
#        D N(H2), beam-averaged Sigma and 10" gas mass for X(HCO+) 1e-9..1e-8, vs dust Sigma from the Xb47 10" continuum (upper bound) and the CO ceiling
# Output out/q27/s108.txt/json
import sys,json,os,numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.stats import norm
import s53_ladder_fit as L, s62_lib as S
os.makedirs('/workspace/out/q27',exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(); T=10**F.lT; LO=np.log(1e-3)
s160=json.load(open('/workspace/out/q6r2/s160.json'))
tau32=F.tauH[:,2]; tau43=F.tauH[:,3]
def floors(pos,key): return s160['HCN3-2'][pos][key]['Tb'],s160['HCO+3-2'][pos][key]['Tb']
def q(x,p,qs=(16,50,84)):
    o=np.argsort(x); c=np.cumsum(p[o]); return [float(x[o][min(np.searchsorted(c,a/100),len(x)-1)]) for a in qs]
def setup(pos):
    I,e=L.DATA[pos]; return S.GEngine(F,[dict(y=np.log(I),rel=e,cls='centre')],sig_b=(.1,.1,.1))
def ll_shared(E,lo=LO,mask=None):
    l,eh,sd,chi=E.lnlike(); eh=eh[:,0]; sd=sd[0]
    ll=-0.5*chi+np.log(np.maximum(norm.cdf((0-eh)/sd)-norm.cdf((lo-eh)/sd),1e-300))-np.log(0-lo)+0.5*np.log(2*np.pi)+np.log(sd)
    if mask is not None: ll=np.where(mask,ll,-np.inf)
    return ll,eh,sd,chi
def ll_sep(E,lo=LO,lq=None,mask=None,ng=36):
    Ci=E.Ci; rM=E.y[None,:]-F.M; W=rM@Ci
    bH=np.array([1,1,1,0,0,0.]); bO=1-bH
    Q0=np.einsum('ij,ij->i',W,rM); uO=W@bO; uH=W@bH; cOO=bO@Ci@bO; cHO=bH@Ci@bO; cHH=bH@Ci@bH
    sdH=1/np.sqrt(cHH); gs=np.linspace(lo,0,ng); acc=None
    for g in gs:
        Qg=Q0-2*g*uO+g*g*cOO; sH=uH-g*cHO; mu=sH/cHH; chi=Qg-sH*sH/cHH
        a=lo if lq is None else max(lo,g-lq); b=0 if lq is None else min(0,g+lq)
        pr=np.log(np.maximum(norm.cdf((b-mu)/sdH)-norm.cdf((a-mu)/sdH),1e-300))
        term=-0.5*chi+pr; acc=term if acc is None else np.logaddexp(acc,term)
    ll=acc+np.log(gs[1]-gs[0])-2*np.log(0-lo)+0.5*np.log(2*np.pi/cHH)
    if mask is not None: ll=np.where(mask,ll,-np.inf)
    return ll
def post(ll,pr):
    w=S.prior_w(F,pr); lp=ll+np.log(np.where(w>0,w,1e-300)); lp=np.where(w>0,lp,-np.inf); p=np.exp(lp-lp.max()); lnZ=lp.max()+np.log(p.sum())-np.log((w>0).sum()); return p/p.sum(),lnZ
def summ(name,ll,pr='logT300',extra=''):
    p,lnZ=post(ll,pr)
    r=dict(Tk=q(T,p),Tk5=q(T,p,(5,))[0],logn=q(F.ln,p),PTk50=float(p[T>50].sum()),PTk30=float(p[T<=30.01].sum()),PTk100=float(p[T>100].sum()),
           tau32=[10**x for x in q(tau32,p)],tau43=[10**x for x in q(tau43,p)],lnZ=float(lnZ))
    P(f'  [{pr}] {name}: Tk {r["Tk"][1]:.0f} [{r["Tk"][0]:.0f},{r["Tk"][2]:.0f}] (5%: {r["Tk5"]:.0f}) P(Tk>50)={r["PTk50"]:.2f} P(Tk<=30)={r["PTk30"]:.3f} log n {r["logn"][1]:.2f} [{r["logn"][0]:.2f},{r["logn"][2]:.2f}] tau32 {r["tau32"][1]:.1f} [{r["tau32"][0]:.1f},{r["tau32"][2]:.1f}] tau43 {r["tau43"][1]:.1f} lnZ {lnZ:.1f} {extra}')
    return r,p
res={}; AREA=1.1331*(10*194.0)**2; rng=np.random.default_rng(3)
for pos in ('N','S','Bridge'):
    E=setup(pos); I,e=L.DATA[pos]; res[pos]=R={}
    P(f'\n########## {pos}: 10" ladder HCN/HCO+ 2-1,3-2,4-3 (K km/s) = {np.round(I,2).tolist()}')
    P('--- A: floor sensitivity (shared eta flat ln [1e-3,1]); T_R,peak=I/(1.0645 w); require model T_R(3-2, both) >= f x T_b')
    ll0,eh,sd,chi=ll_shared(E); r,_=summ('no floor',ll0); R['A_nofloor']=r
    for key,lab in (('native 0.24"','49pc'),('250pc(1.29")','250pc')):
        fH,fO=floors(pos,key); P(f' floor {lab}: T_b HCN/HCO+ 3-2 = {fH:.2f}/{fO:.2f} K')
        for w in (100,150,200,250):
            for f in (0.7,0.85,1.0):
                TRp=np.exp(F.M)/(1.0645*w); mk=(TRp[:,1]>=f*fH)&(TRp[:,4]>=f*fO)
                if mk.sum()==0: P(f'   w={w} f={f}: no model'); continue
                ll,eh_,sd_,chi_=ll_shared(E,mask=mk)
                p,_=post(ll,'logT300')
                # eta posterior: truncated normal mixture
                idx=rng.choice(len(p),20000,p=p); a=norm.cdf((LO-eh_[idx])/sd_); b=norm.cdf((0-eh_[idx])/sd_); le=eh_[idx]+sd_*norm.ppf(np.clip(a+rng.random(len(idx))*(b-a),1e-300,1-1e-16)); et=np.exp(np.percentile(le,[16,50,84]))
                rr,_=summ(f'{lab} w={w} f={f}',ll,extra=f'eta {et[1]:.3f} [{et[0]:.3f},{et[2]:.3f}]'); rr['eta']=et.tolist(); R[f'A_{lab}_w{w}_f{f}']=rr
    P('--- B: separate eta_HCN, eta_HCO+ (flat ln each [1e-3,1]) vs shared; X(HCN/H13CN) -> tau43(HCN) < X*0.128')
    for X in (None,30,50,100):
        mk=None if X is None else (tau43<np.log10(X*0.128)); lab='' if X is None else f' + H13CN X={X} (tau43<{X*0.128:.1f})'
        l1=ll_shared(E,mask=mk)[0]; r,_=summ('shared eta'+lab,l1); R[f'B_shared_X{X}']=r
        r,_=summ('separate eta, flat box'+lab,ll_sep(E,mask=mk)); R[f'B_sep_X{X}']=r
        r,_=summ('separate eta, |ln q|<=ln10'+lab,ll_sep(E,lq=np.log(10),mask=mk)); R[f'B_sepq10_X{X}']=r
    fH,fO=floors(pos,'250pc(1.29")'); TRp=np.exp(F.M)/(1.0645*150); m250=(TRp[:,1]>=0.85*fH)&(TRp[:,4]>=0.85*fO)
    r,_=summ('separate eta + 250pc floor + H13CN X=50',ll_sep(E,mask=m250&(tau43<np.log10(6.4)))); R['B_sep_floor250_H13']=r
    r,_=summ('shared eta + 250pc floor + H13CN X=50',ll_shared(E,mask=m250&(tau43<np.log10(6.4)))[0]); R['B_shared_floor250_H13']=r
    P("--- C: cold (p0012 S-like 20-30 K) gas in the 10\" aperture")
    c25=(np.abs(F.lT-np.log10(25))<0.03); w_=S.prior_w(F,'logT300'); pr_=float(w_[T<=30.01].sum()/w_.sum())
    lg,eh1,sd1,chi1=E.lnlike(); lfix=-0.5*(chi1+(eh1[:,0]/sd1[0])**2)+0.5*0
    for lab,ll in (('shared eta [1e-3,1]',ll_shared(E)[0]),('separate eta',ll_sep(E)),('eta=1 fixed',lfix)):
        p,lnZ=post(ll,'logT300'); pm=float(p[T<=30.01].sum()); odds=(pm/(1-pm+1e-300))/(pr_/(1-pr_)); R[f'C_{lab}']=dict(P_le30=pm,prior=pr_,BF=float(odds))
        P(f'  {lab}: P(Tk<=30)={pm:.4f} (prior {pr_:.3f}); Bayes factor cold/warm = {odds:.3f}')
    for lo in (np.log(1e-3),np.log(0.05),np.log(0.3)):
        l_,eh_,sd_,chi_=ll_shared(E,lo=lo); cc=chi_+((np.clip(eh_,lo,0)-eh_)/sd_)**2; k=np.argmin(np.where(c25,cc,1e9))
        P(f'  chi2(best eta in [{np.exp(lo):.3g},1]): all-T min {cc.min():.2f} (Tk {T[np.argmin(cc)]:.0f}); Tk=25 K min {cc[c25].min():.2f} (n {F.ln[k]:.2f}, tau32 {10**tau32[k]:.1f}, eta_best {np.exp(np.clip(eh_[k],lo,0)):.3f}); Tk<=30 min {cc[T<=30.01].min():.2f}')
        R[f'C_chi_lo{np.exp(lo):.3g}']=dict(all=float(cc.min()),T25=float(cc[c25].min()),T30=float(cc[T<=30.01].min()))
    P(f'--- D: column / mass of posterior branches (X(HCO+) 1e-9..1e-8; mu=2.8 per H2; D=40 Mpc; 10" beam area {AREA:.2e} pc2)')
    TR49=np.exp(F.M)/(1.0645*150); f49=floors(pos,'native 0.24"'); m49=(TR49[:,1]>=0.85*f49[0])&(TR49[:,4]>=0.85*f49[1])
    for lab,mk in (('baseline',None),('+H13CN X=50',tau43<np.log10(6.4)),('250pc floor',m250),('49pc floor',m49),('eta=1',False)):
        if mk is False:
            ll=lfix; eta_f=lambda idx:np.ones(len(idx))
        else:
            ll,eh_,sd_,_=ll_shared(E,mask=mk)
            def eta_f(idx,eh_=eh_,sd_=sd_):
                a=norm.cdf((LO-eh_[idx])/sd_); b=norm.cdf((0-eh_[idx])/sd_); return np.exp(eh_[idx]+sd_*norm.ppf(np.clip(a+rng.random(len(idx))*(b-a),1e-300,1-1e-16)))
        p,_=post(ll,'logT300'); idx=rng.choice(len(p),40000,p=p); eta=eta_f(idx); NO=10**F.lNO[idx]; out={}
        for X in (1e-9,3e-9,1e-8):
            Sig=2.244e-20*eta*NO/X; M=Sig*AREA; pc=lambda v:[float(x) for x in np.percentile(v,[16,50,84])]
            out[f'{X:.0e}']=dict(Sigma=pc(Sig),Mgas=pc(M),P_gt3e3=float((Sig>3e3).mean()),P_gt1e3=float((Sig>1e3).mean()))
            P(f'   {lab} X={X:.0e}: eta {np.median(eta):.3f}; Sigma_beam {np.median(Sig):.0f} [{pc(Sig)[0]:.0f},{pc(Sig)[2]:.0f}] Msun/pc2; M(10") {np.median(M):.1e} [{pc(M)[0]:.1e},{pc(M)[2]:.1e}] Msun; P(Sigma>3e3)={out[f"{X:.0e}"]["P_gt3e3"]:.2f}; P(>1e3)={out[f"{X:.0e}"]["P_gt1e3"]:.2f}')
        R[f'D_{lab}']=out
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter,map_coordinates
fn='/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'; hd=fits.open(fn)[0]; h=hd.header; img=np.squeeze(hd.data).astype(float)
nu=346.5; bmaj,bmin=h['BMAJ']*3600,h['BMIN']*3600; ps=abs(h['CDELT2'])*3600
Tb=1.222e3*img*1e3/(nu**2*bmaj*bmin)
sk=np.sqrt(10.0**2-bmaj*bmin)/2.3548/ps; Tb10=gaussian_filter(np.nan_to_num(Tb),sk); w=WCS(h).celestial
from regions import REGIONS
Jf=lambda Td:(4.799e-11*nu*1e9)/np.expm1(4.799e-11*nu*1e9/Td); kap=0.1*(nu/1000.)**1.5
dres={}; P('\n########## dust: Xb47 346.5 GHz continuum smoothed to 10" (all continuum = dust => upper bound; kappa_gas=0.1(nu/1THz)^1.5 cm2/g, factor 3)')
for pos in ('N','S','Bridge'):
    x,y=w.world_to_pixel_values(REGIONS[pos][0],REGIONS[pos][1]); t=float(map_coordinates(Tb10,[[y],[x]],order=1)[0])
    for Td in (30,50,100):
        tau=-np.log(1-t/(Jf(Td)-Jf(2.73))); Sig=tau/kap*4785.; dres[f'{pos}_{Td}']=dict(Tb10=t,tau=tau,Sigma=Sig,M10=Sig*AREA)
        P(f'  {pos} Td={Td}: Tb(10")={t:.3f} K tau={tau:.4f} Sigma_dust<={Sig:.0f} Msun/pc2 M(10")<={Sig*AREA:.1e} Msun')
res['dust']=dres
json.dump(res,open('/workspace/out/q27/s108.json','w'),indent=1,default=float); open('/workspace/out/q27/s108.txt','w').write('\n'.join(log))
