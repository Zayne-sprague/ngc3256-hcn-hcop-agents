# Q3 round 4, fix 3: disk (E1) vs centre excitation contrast, tested (i) directly on calibration-independent ratio-of-ratios and (ii) on the RADEX posteriors with the calibration treated as COMMON to both positions.
# (i) RR = [I_J/I_3](E1) / [I_J/I_3](centre N), per molecule, J=2 and 4.  Error terms (ln): stat; PB term at E1 per line = max(0.3(1-PB), 0.5|1-A/B|) from s55a (A/B: pbcor map smoothed vs pbcor*PB smoothed then /PB);
#     Q12 recovery difference (centre f2 U[.85,.95], f4 U[1.03,1.12]; E1 f U[.95,1.05]); the band scales (10%/band) and the per-line 7% term are COMMON to both positions in the same cubes and cancel in RR; shown both ways.
# (ii) per draw of the shared nuisances (band scales, optional shared line terms) the two positions have independent theta posteriors; draws weighted by the product of their evidences.
# usage: python3 -I s65_disk_contrast.py -> out/q3d/s65_disk_contrast.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L, s62_lib as S
from scipy.special import logsumexp
OUT='/workspace/out/q3d/'; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
PJ=json.load(open('/workspace/out/q3b/s55a_pb_snr.json')); PBJ=PJ['pb']; PT=PJ['pbtest']; SNR=PJ['snr']
def pbextra(pos):
    pb=np.array(PBJ[pos]); s=0.3*(1-pb); ab=[PT[f'{b}|{pos}'] for b in ('Xb55 (B5, J=2)','Xb4f (B6, J=3)','Xb47 (B7, J=4)')]
    sh=np.maximum(s,0.5*np.abs(1-np.array([a[0] for a in ab]))); so=np.maximum(s,0.5*np.abs(1-np.array([a[1] for a in ab]))); return np.r_[sh,so]
ex=pbextra('E1'); P('# E1 PB/beam extra sigma per datum (HCN 2,3,4; HCO+ 2,3,4): '+' '.join(f'{v:.3f}' for v in ex))
P('# E1 stat S/N of the 10" |v|<350 km/s fluxes (s55a): HCN 2-1,3-2,4-3 = '+', '.join(f'{v:.1f}' for v in SNR['E1']['HCN'])+'; HCO+ = '+', '.join(f'{v:.1f}' for v in SNR['E1']['HCO'])+'  (stat only; the HCN 4-3 PB/A-B systematic is 34% at E1)')
centre=np.mean([L.DATA[k][0] for k in ('N','S','Bridge')],0); IN,eN=L.DATA['N']; IE,eE=L.DATA['E1']
# ---------------- (i) ratio of ratios
rng=np.random.default_rng(5); NS=20000
rec_c=(np.exp(np.log(rng.uniform(.85,.95,NS))),rng.uniform(1.03,1.12,NS)); rec_e=(rng.uniform(.95,1.05,NS),rng.uniform(.95,1.05,NS))
P('\n# (i) ratio of ratios E1/centre(N). Columns: RR observed; ln-sigma from stat; from PB/A-B term; from recovery; total cal-cancelled; total if band(10%) + 7% per-line terms were NOT common (conservative); significance of |ln RR|')
tab={}
for mol,o,nm in (('H',0,'HCN'),('O',3,'HCO+')):
    for jn,j,kk in (('T21/T32',0,0),('T43/T32',2,1)):
        rr=(IE[o+j]/IE[o+1])/(IN[o+j]/IN[o+1]); lrr=np.log(rr)
        st=np.hypot(np.hypot(eE[o+j],eE[o+1]),np.hypot(eN[o+j],eN[o+1]))
        pb=np.hypot(ex[o+j],ex[o+1])
        # recovery: true ratio = obs*f_J/f_3 ; RR_true/RR_obs = (fE_J/fc_J) ; sample
        fe=rec_e[kk]; fc=rec_c[kk]; corr=np.log(fe/fc); rec_sd=corr.std(); rec_mean=corr.mean()
        tot=np.sqrt(st**2+pb**2+rec_sd**2); unc=np.sqrt(tot**2+2*(0.1**2+0.1**2)*0+ (0.1**2+0.1**2)*2*0.0 + 2*(2*0.07**2) + 2*(0.1**2+0.1**2))
        # conservative: independent band scales per position (2 positions x (J band + 3-2 band)) + independent 7% per line (2 positions x 2 lines)
        lrr_c=lrr+rec_mean
        sig=abs(lrr_c)/tot
        tab[(nm,jn)]=dict(RR=rr,RRcorr=float(np.exp(lrr_c)),st=st,pb=pb,rec=rec_sd,tot=tot,unc=unc,sig=sig)
        P(f'   {nm:5s} {jn}: RR={rr:.2f} (recovery-corrected {np.exp(lrr_c):.2f});  stat {st:.3f}  PB/AB {pb:.3f}  recovery {rec_sd:.3f}  -> total(cancelled) {tot:.3f} => |ln RR|/sigma={sig:.1f};  if calibration NOT common: sigma {unc:.3f} => {abs(lrr_c)/unc:.1f} sigma')
res['RR']={f'{a}|{b}':v for (a,b),v in tab.items()}
P('   Reading: RR<1 for T43/T32 means E1 has weaker 4-3 relative to 3-2 than the centre; RR>1 for T21/T32 means stronger 2-1. T43/T32 HCN: factor ~2 deficit; it survives band calibration (cancels) but the E1 PB/A-B term at Xb47 for HCN is 17% (HCO+ 10%).')
# ---------------- (ii) posterior contrast with shared calibration
# Per-position Gaussian terms (stat, PB(E1), 7% per line, recovery with Q12 mean shift and sd) are analytic (s62_lib.GEngine, n=1, band blocks removed); the band scales (and, optionally, a shared per-line term) are Monte-Carlo draws applied identically to both positions (Sobol, 256 draws).
import os
os.environ.setdefault('OMP_NUM_THREADS','2')
from multiprocessing import Pool
from scipy.special import log_ndtr
def cols_mask(which):
    v=np.zeros(6)
    if which=='HCN': v[3:]=30.
    if which=='HCO+': v[:3]=30.
    return v
def contrast(args):
    nm,sig_b,shared_line,which,nd=args
    F=S.Grid(); si=0.0 if shared_line>0 else 0.07; dm=cols_mask(which)
    pe=dict(y=np.log(IE),rel=eE,cls='E1',extra=np.hypot(ex,dm)); pc=dict(y=np.log(IN),rel=eN,cls='centre',extra=dm)
    En=[S.GEngine(F,[p],sig_b=(0,0,0),s_ind=si) for p in (pe,pc)]
    ws={k:S.prior_w(F,k) for k in ('logT300','linT500')}
    lw={k:np.log(np.where(w>0,w,1)) for k,w in ws.items()}; mk={k:w>0 for k,w in ws.items()}
    from scipy.stats import qmc,norm
    u=qmc.Sobol(d=3+6,scramble=True,seed=11).random(nd); z=norm.ppf(np.clip(u,1e-9,1-1e-9))
    off=(z[:,:3]*np.array(sig_b))[:,[0,1,2,0,1,2]]+z[:,3:]*shared_line
    tm=lambda p:np.bincount(F.iT,weights=p,minlength=len(S.LT)); nmg=lambda p:np.bincount(F.iN,weights=p,minlength=len(S.LN))
    Z={k:[] for k in ws}; PTm={k:[] for k in ws}; PNm={k:[] for k in ws}
    for d in range(nd):
        L2=[]
        for E in En:
            y=E.y-off[d]; g=(E.Ci@y)[None,:]-E.MC[(0,0)]; r=y[None,:]-F.M; q0=np.einsum('ij,ij->i',r,g); bb=g.sum(1); a_=E.Ai[0,0]; eh=bb*a_; chi=q0-bb*bb*a_
            L2.append(-0.5*chi+log_ndtr((E.cap-eh)/np.sqrt(a_)))
        for k in ws:
            zz=[];tt=[];nn=[]
            for l in L2:
                l=np.where(mk[k],l+lw[k],-np.inf); m=l.max(); p=np.exp(l-m); sp=p.sum(); zz.append(m+np.log(sp)); tt.append(tm(p/sp)); nn.append(nmg(p/sp))
            Z[k].append(zz[0]+zz[1]); PTm[k].append(tt); PNm[k].append(nn)
    out={}
    for k in ws:
        W=np.array(Z[k]); W=np.exp(W-W.max()); W/=W.sum(); pT=0;pn=0;ptE=0;ptC=0
        for d in range(nd):
            a,b=PTm[k][d]; pT+=W[d]*float(np.sum(b*(np.cumsum(a)-a+0.5*a))); ptE=ptE+W[d]*a; ptC=ptC+W[d]*b
            a,b=PNm[k][d]; pn+=W[d]*float(np.sum(a*(np.cumsum(b)-b+0.5*b)))
        qT=lambda m,q:float(10**S.LT[min(np.searchsorted(np.cumsum(m)/np.sum(m),q),len(S.LT)-1)])
        out[k]=dict(P_Tk_E1_lt_centre=pT,P_n_E1_gt_centre=pn,Tk_E1=[qT(ptE,q) for q in (.16,.5,.84)],Tk_c=[qT(ptC,q) for q in (.16,.5,.84)])
    return nm,which,out
if __name__=='__main__':
    P('\n# (ii) posterior contrast E1 vs centre(N): same band-scale draws applied to both positions; stat+PB(E1)+recovery (Q12 per position) always included; prior flat log Tk<=300 (and flat Tk<=500).')
    jobs=[(nm,sb,sl,w,256) for nm,sb,sl in (('shared band 10/10/10%',(.1,.1,.1),0.0),('shared band 15/10/15%',(.15,.1,.15),0.0),('shared band 10% + shared 7% per-line term (s_ind=0)',(.1,.1,.1),0.07),('no calibration error (stat+PB+recovery only)',(1e-3,1e-3,1e-3),0.0)) for w in ('both','HCN','HCO+')]
    with Pool(6) as pool: R=pool.map(contrast,jobs,chunksize=1)
    for nm,which,r in R:
        res[f'{nm}|{which}']=r
        for k,v in r.items():
            P(f'   [{nm}] [{which:5s}] prior {k}: P(Tk_E1<Tk_centre)={v["P_Tk_E1_lt_centre"]:.2f}  P(n_E1>n_centre)={v["P_n_E1_gt_centre"]:.2f}   Tk_E1 {v["Tk_E1"][1]:.0f} [{v["Tk_E1"][0]:.0f},{v["Tk_E1"][2]:.0f}]  Tk_centre(N) {v["Tk_c"][1]:.0f} [{v["Tk_c"][0]:.0f},{v["Tk_c"][2]:.0f}]')
    open(OUT+'s65_disk_contrast.txt','w').write('\n'.join(log)); json.dump(res,open(OUT+'s65_disk_contrast.json','w'),indent=1,default=float)
