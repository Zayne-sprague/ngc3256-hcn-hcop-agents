# Q8c: source-model uv-coverage forward model (reviewer fix 1).
# Each position (N, S, Bridge) = compact + extended circular Gaussian (peak Tb A_c, A_e K; FWHM s, s_e "), seen through array a with a short-baseline loss:
#   visibility high-pass H(q)=1-exp(-(q/qmin)^2) [VARIANT taper; L = 2 sqrt(ln2)/(pi qmin) = "largest scale kept", FWHM "], or a sharp cut [VARIANT sharp].
#   Peak of the image smoothed to beam B=(B1,B2) (exact beam-difference smoothing of the cleaned image, native beam b_a):
#   taper: A s^2 [1/sqrt(S1 S2) - 1/sqrt((S1+L^2)(S2+L^2))]/n_a, S_i=s^2+B_i^2, n_a=1-b1b2/sqrt((b1^2+L^2)(b2^2+L^2)) (dirty-beam normalisation: a point source is fully recovered)
#   taper_raw: same without n_a;   sharp: A s^2/sqrt(S1S2) exp(-(sqrt(S1S2)-b1b2)/L^2)
# Array parameters: L_f (Xb4f), L_d (Xb4d), L_b (Xb4b) free; Xb47 and Xb55 observe the SAME (compact) configuration as Xb4f (beam ratio = wavelength ratio, checked below)
#   -> same minimum baseline Bmin up to a lognormal(0,0.25) factor rB: L_47 = L_f (lam47/lamf)/rB47, L_55 = L_f (lam55/lamf)/rB55.
# Data (all calibrate L's and geometry JOINTLY): 3-2 native-beam peak Tb of Xb4b/Xb4d/Xb4f (+Xb4d smoothed to 1.7") at N,S,Bridge for HCN and HCO+;
#   Xb4d/Xb4f and Xb4b/Xb4f 3-2 flux ratios at 3 matched beams (6.2x3.7, 6.5 circ., 9.5x6.4); continuum Xb4d/Xb4f, Xb4b/Xb4f ratios at 2 matched beams, N,S,Bridge (own geometry per position);
#   flux scale nuisance g_d,g_b ~ lognormal(0,0.14) relative to Xb4f (10% per OUS).  emcee posterior.
# Prediction: Xb47(4-3)/Xb4f(3-2) recovery and the correction to the N/S T43/T32 contrast, over the posterior (geometry + array-parameter uncertainty).
import sys,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, emcee, time
VARFULL=sys.argv[1] if len(sys.argv)>1 else 'taper'; PHYS=VARFULL.endswith('_phys'); VAR=VARFULL.replace('_phys',''); NSTEP=int(sys.argv[2]) if len(sys.argv)>2 else 12000
O='/workspace/out/q8c/'; import os; os.makedirs(O,exist_ok=True)
rng=np.random.default_rng(int(sys.argv[3]) if len(sys.argv)>3 else 1)
REG=['N','S','Bridge']; LIN=['HCN','HCO+']
LAM={'f':1/265.9,'d':1/265.9,'b':1/265.9,'c':1/257.7,'47':1/354.5,'55':1/177.3}   # relative wavelengths (3-2 line, continuum, 4-3, 2-1)
NATB={'b':(0.25,0.235),'d':(0.95,0.75),'f':(5.9,3.6),'47':(4.2,2.75),'55':(9.0,6.0)}
def peak(A,s,B,L,nb,var=VAR):
    s2=s**2; S1=s2+B[0]**2; S2=s2+B[1]**2; L2=L**2
    if var=='sharp':
        S=np.sqrt(S1*S2); return A*s2/S*np.exp(-(S-nb[0]*nb[1])/L2)
    v=A*s2*(1/np.sqrt(S1*S2)-1/np.sqrt((S1+L2)*(S2+L2)))
    if var=='taper': v=v/(1-nb[0]*nb[1]/np.sqrt((nb[0]**2+L2)*(nb[1]**2+L2)))
    return v
def unf(A,s,B): return A*s**2/np.sqrt((s**2+B[0]**2)*(s**2+B[1]**2))
# ---------------- data
T=pd.read_csv('/workspace/out/q4/peak_tb.csv'); O17=pd.read_csv('/workspace/out/q4c/ours_1p7arcsec.csv'); RC=pd.read_csv('/workspace/out/q8/recovery_32_arrays.csv'); CT=pd.read_csv('/workspace/out/q8/continuum_matched.csv')
BM={'Xb4f(6.2x3.7,82)':(6.2,3.7),'circ6.5':(6.5,6.5),'Xb55(9.5x6.4,-88)':(9.5,6.4)}
NAT=[]  # (reg,line,arr,beam,nbeam_key,Tb,err)
for _,r in T.iterrows():
    if r.ous in ('Xb4b','Xb4d','Xb4f') and r.line in ('HCN 3-2','HCO+ 3-2'):
        b=tuple(float(x) for x in r.beam.split('x')); NAT.append((r.region,r.line.split()[0],r.ous[-1],b,r.Tb_peak_K,np.hypot(r.eTb_K,0.04*r.Tb_peak_K)))
for _,r in O17.iterrows(): NAT.append((r.region,r.line.split()[0],'d',(1.7,1.7),r.Tpeak_K,np.hypot(r.rms_K,0.04*r.Tpeak_K)))
RAT=[]  # line ratios array/Xb4f at matched beam
for _,r in RC.iterrows(): RAT.append((r.region,r.line,r['array'][-1],BM[r.beam],r.ratio_to_Xb4f,np.hypot(r.err,0.03*r.ratio_to_Xb4f)))
cpv=CT.pivot_table(index=['beam','array'],columns='region',values='mJy_beam')
CRAT=[]
for bn in ('Xb4f(6.2x3.7,82)','Xb55(9.5x6.4,-88)'):
    for a in ('Xb4d','Xb4b'):
        for reg in REG: CRAT.append((reg,a[-1],BM[bn],cpv.loc[(bn,a),reg]/cpv.loc[(bn,'Xb4f'),reg]))
# ---------------- parameter layout
names=[]
for reg in REG:
    for l in LIN: names+= [f'{reg}_{l}_{k}' for k in ('lns','lnse','lnAc','lnAe')]
for reg in REG: names+=[f'cont_{reg}_{k}' for k in ('lns','lnse','lnrE')]
names+=['lnLf','lnLd','lnLb','lngd','lngb','lnrB47','lnrB55']; ix={n:i for i,n in enumerate(names)}; ND=len(names)
lo=np.full(ND,-np.inf); hi=np.full(ND,np.inf)
for n,i in ix.items():
    if n.endswith('lns'): lo[i],hi[i]=np.log(0.08),np.log(1.5)
    if n.endswith('lnse'): lo[i],hi[i]=np.log(1.5),np.log(30)
    if n.endswith('lnAc') or n.endswith('lnAe'): lo[i],hi[i]=np.log(1e-2),np.log(1e2)
    if n.endswith('lnrE'): lo[i],hi[i]=np.log(1e-3),np.log(1e2)
    if n in ('lnLf','lnLd','lnLb'): lo[i],hi[i]=np.log(0.3),np.log(8.2 if PHYS else 60)   # PHYS: shortest baseline >=15 m -> L=0.53 lam/Bmin <= 8.2" at 266 GHz
def arrL(th,a,cont=False):
    Lf=np.exp(th[:,ix['lnLf']]); lam=LAM['c' if cont else 'f']/LAM['f']
    if a=='f': return Lf*lam
    if a=='d': return np.exp(th[:,ix['lnLd']])*lam
    if a=='b': return np.exp(th[:,ix['lnLb']])*lam
    if a=='47': return Lf*(LAM['47']/LAM['f'])/np.exp(th[:,ix['lnrB47']])
    if a=='55': return Lf*(LAM['55']/LAM['f'])/np.exp(th[:,ix['lnrB55']])
def comp(th,reg,l): return [np.exp(th[:,ix[f'{reg}_{l}_{k}']]) for k in ('lns','lnse','lnAc','lnAe')]
def model(th,reg,l,a,B,nb=None):
    s,se,Ac,Ae=comp(th,reg,l); L=arrL(th,a); nb=nb or NATB[a if a in NATB else 'f']
    return peak(Ac,s,B,L,nb)+peak(Ae,se,B,L,nb)
def gscale(th,a): return np.exp(th[:,ix['lngd']]) if a=='d' else (np.exp(th[:,ix['lngb']]) if a=='b' else 1.0)
def lnprob(th): return lnprob_core(th)
def lnprob_core(th):
    th=np.atleast_2d(th); ok=np.all((th>lo)&(th<hi),axis=1); lp=np.full(len(th),-np.inf)
    if not ok.any(): return lp
    t=th[ok]; c2=np.zeros(len(t))
    if PHYS:   # Bmin>=15 m also for Xb47 (L<=6.2") and Xb55 (L<=12.7")
        g=(arrL(t,'47')<=6.2)&(arrL(t,'55')<=12.7)
        if not g.all():
            ok2=ok.copy(); ok2[np.where(ok)[0][~g]]=False; lp=np.full(len(th),-np.inf)
            if ok2.any(): lp[ok2]=lnprob_core(th[ok2])
            return lp

    for reg,l,a,b,Tb,e in NAT:
        m=model(t,reg,l,a,b,nb=b if a!='d' or b[0]!=1.7 else NATB['d'])*gscale(t,a); c2+=((m-Tb)/e)**2
    for reg,l,a,B,r,e in RAT:
        mf=model(t,reg,l,'f',B); ma=model(t,reg,l,a,B)*gscale(t,a); c2+=((ma/mf-r)/e)**2
    for reg,a,B,r in CRAT:
        s,se,rE=[np.exp(t[:,ix[f'cont_{reg}_{k}']]) for k in ('lns','lnse','lnrE')]
        mf=peak(1.,s,B,arrL(t,'f',True),NATB['f'])+peak(rE,se,B,arrL(t,'f',True),NATB['f']); La=arrL(t,a,True)
        ma=(peak(1.,s,B,La,NATB[a])+peak(rE,se,B,La,NATB[a]))*gscale(t,a); c2+=((ma/mf-r)/np.hypot(0.01*r,0.03*r))**2
    prior=-0.5*((t[:,ix['lngd']]/0.14)**2+(t[:,ix['lngb']]/0.14)**2+(t[:,ix['lnrB47']]/0.25)**2+(t[:,ix['lnrB55']]/0.25)**2)
    lp[ok]=-0.5*c2+prior; return np.where(np.isnan(lp),-np.inf,lp)
if __name__=='__main__':
    nw=160; t0=time.time()
    # initialise from the best of many prior draws, then scatter
    P=lo+(hi-lo)*rng.random((400000,ND)); 
    for n in ('lngd','lngb','lnrB47','lnrB55'): P[:,ix[n]]=rng.normal(0,0.14 if 'g' in n[2:3] else 0.25,len(P))
    for n in ('lnLf','lnLd','lnLb'): pass
    L=np.concatenate([lnprob(P[i:i+50000]) for i in range(0,len(P),50000)]); best=np.argsort(-L)[:nw]
    p0=P[best]; print('init best lnp',L[best[0]],'median',np.median(L[best]),'time',time.time()-t0,flush=True)
    S=emcee.EnsembleSampler(nw,ND,lnprob,vectorize=True,moves=[(emcee.moves.DEMove(),0.8),(emcee.moves.DESnookerMove(),0.2)])
    S.run_mcmc(p0,NSTEP,progress=False)
    ch=S.get_chain(); lpr=S.get_log_prob(); print('acc',S.acceptance_fraction.mean(),'time',time.time()-t0)
    try: print('tau',S.get_autocorr_time(tol=0))
    except Exception as e: print(e)
    burn=NSTEP//2; flat=ch[burn::20].reshape(-1,ND); lf=lpr[burn::20].reshape(-1)
    np.savez(O+f'chain_{VARFULL}.npz',flat=flat,lp=lf,names=names,chain=ch[::10],lpall=lpr[::10])
    k=np.argmax(lf); print('max lnp',lf[k],'chi2min',-2*lf.max(),'ndata',len(NAT)+len(RAT)+len(CRAT))
    for n in ('lnLf','lnLd','lnLb','lngd','lngb','lnrB47','lnrB55'):
        v=flat[:,ix[n]]; print(n,'median',np.exp(np.median(v)),'16-84',np.exp(np.percentile(v,[16,84])))
