# Q61 step 2: (A) two-band (98,258 GHz) free-free+sync+dust decomposition at the nuclear peaks for circular beams 1.1-6.5";
# (B) multi-band (98/170/258/346) joint N+S SED MCMC at 4.5/6.5/10.5" for beta, Td_N, Td_S, Td_N/Td_S and dust fractions.
# python3 -I q61_sed.py -> out/q61/sed.json, sed.log (needs out/q61/phot.json)
import json,numpy as np,emcee,warnings; warnings.filterwarnings('ignore')
rng=np.random.default_rng(61); ph=json.load(open('/workspace/out/q61/phot.json'))
hk=4.799e-2   # h/k in K/GHz
def Bn(nu,T): return nu**3/np.expm1(hk*nu/T)
def dust_shape(nu,beta,T,nu0=257.8): return (nu/nu0)**beta*Bn(nu,T)/Bn(nu0,T)
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s); log.append(s)
FFPRI={'A':{'N':(0.63,0.41),'S':(0.74,0.31)},'B':{'N':(0.76,0.12),'S':(0.90,0.12)},'C':{'N':(0.5,0.29),'S':(0.5,0.29)}}
P('f_ff(98 GHz) priors (mean,sd): A=own H40a Q10c (N 0.63+-0.41, S 0.74+-0.31, 1.05" r<1"), B=Michiyama+2020 SED fit (0.76/0.90, +-0.12 adopted), C=uninformative (0.5+-0.29); truncated to [0,1]')
OUT={}
NMC=40000
def trunc(mu,sd,n):
    x=rng.normal(mu,sd,6*n); x=x[(x>=0)&(x<=1)]; return x[:n]
def alpha_dust(beta,T,n1=257.8,n3=98.05): return -np.log(dust_shape(n3,beta,T,n1))/np.log(n1/n3)
P('\n(A) two-band decomposition, S258 = ff + sync + dust; ff alpha=-0.1; sync alpha U(-1.2,-0.4); dust beta U(1.2,2.2), Td U(20,50) (exact Planck slope 98->258); scale errors 5% (B3) and 5% (B6d) common to N and S; noise from phot.json')
q=257.8/98.05
for fw in (1.1,1.5,2.5,4.5,6.5):
    d={s:{b:ph[f'{fw}|{b}|{s}'] for b in ('B3','B6d')} for s in 'NS'}
    g3=1+0.05*rng.standard_normal(NMC); g6=1+0.05*rng.standard_normal(NMC)
    beta=rng.uniform(1.2,2.2,NMC); Td=rng.uniform(20,50,NMC); ad=alpha_dust(beta,Td); as_=rng.uniform(-1.2,-0.4,NMC)
    for sc in 'ABC':
        res={}
        for s in 'NS':
            s3=d[s]['B3']['mJy_beam']*g3+d[s]['B3']['rms_mJy']*rng.standard_normal(NMC); s6=d[s]['B6d']['mJy_beam']*g6+d[s]['B6d']['rms_mJy']*rng.standard_normal(NMC)
            f=rng.normal(*FFPRI[sc][s],NMC); f=np.clip(f,0,1)   # clipped (mass at 0/1)
            r=s6/s3
            fs=(f*q**-0.1+(1-f)*q**ad-r)/(q**ad-q**as_); fd=1-f-fs
            ok=(fs>=0)&(fd>=0)
            res[s]=dict(dust=np.where(ok,fd*q**ad*s3,np.nan),ff=np.where(ok,f*q**-0.1*s3,np.nan),sync=np.where(ok,fs*q**as_*s3,np.nan),tot=s6,rej=1-ok.mean())
            res[s]['fd']=res[s]['dust']/s6
        ratio=res['N']['dust']/res['S']['dust']; tot=res['N']['tot']/res['S']['tot']
        pc=lambda x:'%.2f [%.2f,%.2f]'%tuple(np.nanpercentile(x,[50,16,84]))
        P(f' {fw:4.1f}" f_ff-prior {sc}: dust frac(258) N {pc(res["N"]["fd"])} S {pc(res["S"]["fd"])} | ff frac(258) N {pc(res["N"]["ff"]/res["N"]["tot"])} S {pc(res["S"]["ff"]/res["S"]["tot"])} | sync frac N {pc(res["N"]["sync"]/res["N"]["tot"])} S {pc(res["S"]["sync"]/res["S"]["tot"])} | total N/S {pc(tot)} dust N/S {pc(ratio)} (rejected N {res["N"]["rej"]:.2f} S {res["S"]["rej"]:.2f})')
        p3=lambda x:np.nanpercentile(x,[50,16,84]).tolist()
        OUT[f'A|{fw}|{sc}']=dict(fdN=p3(res['N']['fd']),fdS=p3(res['S']['fd']),dust_ratio=p3(ratio),tot_ratio=p3(tot),ffS=p3(res['S']['ff']/res['S']['tot']),syncS=p3(res['S']['sync']/res['S']['tot']),ffN=p3(res['N']['ff']/res['N']['tot']),syncN=p3(res['N']['sync']/res['N']['tot']),rejN=res['N']['rej'],rejS=res['S']['rej'])
        np.save(f'/workspace/out/q61/dustratio_{fw}_{sc}.npy',ratio[np.isfinite(ratio)])
BANDS={4.5:['B3','B6d','B7'],6.5:['B3','B6f','B7'],10.5:['B3','B5','B6f','B7']}
SC={'B3':0.10,'B5':0.10,'B6f':0.07,'B7':0.10,'B6d':0.07}
def run(fw,sc,bands,label):
    nu=np.array([ph[f'{fw}|{b}|N']['nu'] for b in bands])
    y={s:np.array([ph[f'{fw}|{b}|{s}']['mJy_beam'] for b in bands]) for s in 'NS'}; e={s:np.array([ph[f'{fw}|{b}|{s}']['rms_mJy'] for b in bands]) for s in 'NS'}
    ref=[i for i,b in enumerate(bands) if b in ('B6f','B6d')][0]
    k=[i for i in range(len(bands)) if i!=ref]; sg=np.array([SC[bands[i]] for i in k]); ng=len(k)
    S98={s:y[s][0] for s in 'NS'}
    def model(p,i):
        return np.exp(p[4+i])*dust_shape(nu,p[2],np.exp(p[i]))+p[6+i]*(nu/98.05)**p[3]+p[8+i]*(nu/98.05)**-0.1
    def lnp(p):
        Td=np.exp(p[:2])
        if not(15<Td.min() and Td.max()<70 and 1.0<p[2]<2.4 and -1.2<p[3]<-0.3 and p[6:10].min()>=0): return -np.inf
        g=np.ones(len(bands)); g[k]=np.exp(p[10:10+ng])
        lp=-0.5*np.sum((p[10:10+ng]/sg)**2)
        for j,s in enumerate('NS'):
            lp+=-0.5*np.sum(((y[s]-model(p,j)*g)/e[s])**2)
            f=p[8+j]/S98[s]
            if f>1.0: return -np.inf
            mu,sd=FFPRI[sc][s]; lp+=-0.5*((f-mu)/sd)**2
        return lp
    nd=10+ng; nw=64
    p0=np.zeros(nd); p0[0]=p0[1]=np.log(30); p0[2]=1.8; p0[3]=-0.7
    for j,s in enumerate('NS'): p0[4+j]=np.log(y[s][ref]*0.5); p0[6+j]=0.1*S98[s]; p0[8+j]=0.6*S98[s]
    pos=p0+1e-2*rng.standard_normal((nw,nd))*np.maximum(abs(p0),0.1); pos[:,6:10]=np.abs(pos[:,6:10])
    S=emcee.EnsembleSampler(nw,nd,lnp,moves=[(emcee.moves.DEMove(),0.8),(emcee.moves.DESnookerMove(),0.2)]); S.run_mcmc(pos,12000,progress=False)
    ch=S.get_chain(discard=4000,thin=10,flat=True); ac=np.nanmean(S.get_autocorr_time(tol=0))
    TdN,TdS=np.exp(ch[:,0]),np.exp(ch[:,1]); pc=lambda x:'%.2f [%.2f,%.2f]'%tuple(np.percentile(x,[50,16,84]))
    fdN=np.exp(ch[:,4])/y['N'][ref]; fdS=np.exp(ch[:,5])/y['S'][ref]; dr=np.exp(ch[:,4]-ch[:,5])
    fsN=ch[:,6]*(257.8/98.05)**ch[:,3]/y['N'][ref]; fsS=ch[:,7]*(257.8/98.05)**ch[:,3]/y['S'][ref]
    P(f' {label} {fw}" prior {sc} bands {bands}: Td_N {pc(TdN)} Td_S {pc(TdS)} Td_N/Td_S {pc(TdN/TdS)} P(TdN>TdS)={np.mean(TdN>TdS):.2f} beta {pc(ch[:,2])} | dust frac(258) N {pc(fdN)} S {pc(fdS)} | sync frac(258) N {pc(fsN)} S {pc(fsS)} | dust(258) N/S {pc(dr)} [tau_int~{ac:.0f}, acc {np.mean(S.acceptance_fraction):.2f}]')
    p3=lambda x:np.percentile(x,[50,16,84]).tolist()
    OUT[f'B|{fw}|{sc}|{label}']=dict(TdN=p3(TdN),TdS=p3(TdS),Tratio=p3(TdN/TdS),beta=p3(ch[:,2]),fdN=p3(fdN),fdS=p3(fdS),dust_ratio=p3(dr),P_TdN_gt_TdS=float(np.mean(TdN>TdS)))
P('\n(B) joint N+S multi-band MCMC (emcee 64 walkers x 12000, burn 4000): dust (nu/258)^beta B_nu(Td)/B_nu(258) (optically thin; Tb<<Td), ff alpha -0.1, sync power law, per-band lognormal scale (B3 10%, B5 10%, B7 10% vs B6f/B6d 7% reference), Td log-uniform... flat in lnTd 15-70 K, beta U(1,2.4)')
for fw,bands in BANDS.items():
    for sc in 'ABC': run(fw,sc,bands,'all')
json.dump(OUT,open('/workspace/out/q61/sed.json','w'),indent=1); open('/workspace/out/q61/sed.log','w').write('\n'.join(log))
