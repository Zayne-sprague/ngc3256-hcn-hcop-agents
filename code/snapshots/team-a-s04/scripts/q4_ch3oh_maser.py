# Q4e: is the Band-5 CH3OH Q-branch ((J_-1)-(J_0) E, 165.05-166.9 GHz) excess over the Band-6 CH3OH lines a non-LTE inversion / weak maser, thermal sub-thermal excitation, or spatial filtering?
#  A) profile test: Q J=1-3 blend (Xb55 spw22, 9.36x6.28" native) fitted as [J=1+2 at LTE weight] + [J=3] Gaussians with common v0, FWHM; compare to HCO+ 2-1 at the same beam/position and to LTE a(J3)/a(J1+2)
#  B) LTE (thin) ladder fit with error floor from the scatter of the J>=4 values (baseline systematics), Xb4f lines divided by an assumed flux recovery f (1, 0.7, 0.5)
#  C) RADEX (LAMDA e-ch3oh.dat, H2, Tbg 2.73, dv=150) grid: chi2 with and without the Q J=1-3 point; best-fit excess factor; tau sign of the Q lines (negative = inversion)
# Outputs: out/q4/ch3oh_profile.csv, ch3oh_lte_scan.csv, ch3oh_radex_grid.csv, ch3oh_maser_summary.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_lib import *
from ch3oh import c as CAT
import subprocess,os,re,itertools
from concurrent.futures import ProcessPoolExecutor
RA0,DE0=156.9634583333,-43.90461111111
TB=(9.357,6.284,-88.45)
M=pd.read_csv('/workspace/out/ch3oh_matched.csv')
Q123=(165.050175,165.061130,165.099240)
def cat_line(f0): return CAT[(CAT.f0-f0).abs()<2e-4].iloc[0]
def lte_w(f0,T): l=cat_line(f0); return l.f0*l.sijmu2*np.exp(-l.Eu/T)
# ---------------- A) profile test
def gauss(x,a,v0,w): return a*np.exp(-0.5*((x-v0)/(w/2.355))**2)
def profile():
    from scipy.optimize import least_squares
    fr,s,_,_=spectrum('Xb55','spw22','N',None,pos=(RA0,DE0)); fr2,s2,_,_=spectrum('Xb55','spw18','N',TB,pos=(RA0,DE0))
    out=[]
    # HCO+ 2-1 core Gaussian (|v-2780|<450) and HCN 2-1
    for nm,(spw,f0) in {'HCO+ 2-1':('spw18',178.375056),'HCN 2-1':('spw16',177.261111)}.items():
        fq,sq,_,_=spectrum('Xb55',spw,'N',TB,pos=(RA0,DE0)); v=vel(fq,f0); m=np.abs(v-2780)<450
        r=least_squares(lambda p:gauss(v[m],*p[:3])+p[3]+p[4]*(v[m]-2780)-sq[m]*1e3,[sq[m].max()*1e3,2780,150,0,0]); J=r.jac; sg=np.sqrt(np.sum(r.fun**2)/(m.sum()-5)); cv=np.linalg.inv(J.T@J)*sg**2
        out.append(dict(model=nm,v0=r.x[1],e_v0=np.sqrt(cv[1,1]),fwhm=r.x[2],e_fwhm=np.sqrt(cv[2,2]),a_mJy=r.x[0]))
    # Q123 blend: components J1,J2 (fixed LTE ratio at Trot) and J3 with free amplitudes; common v0,FWHM (v0 referenced to J=1 rest freq)
    for T in (6.,9.,15.,30.):
        v=vel(fr,Q123[0]); dv12=(1-Q123[1]/Q123[0])*C; dv3=(1-Q123[2]/Q123[0])*C
        w1,w2=lte_w(Q123[0],T),lte_w(Q123[1],T); f1=w1/(w1+w2); m=np.abs(v-2777)<500
        def model(p,vv=v[m]): return p[0]*(f1*gauss(vv,1,p[2],p[3])+(1-f1)*gauss(vv-dv12,1,p[2],p[3]))+p[1]*gauss(vv-dv3,1,p[2],p[3])+p[4]+p[5]*(vv-2777)
        best=None
        for p0 in ([60,10,2777,150,0,0],[30,60,2777,150,0,0],[10,10,2780,100,0,0],[90,0,2770,160,0,0]):
            try: r=least_squares(lambda p:(model(p)-s[m]*1e3),p0)
            except Exception: continue
            if best is None or r.cost<best.cost: best=r
        r=best; sg=np.sqrt(2*r.cost/(m.sum()-6)); J=r.jac; cv=np.linalg.pinv(J.T@J)*sg**2; dvv=abs(np.median(np.diff(v))); n=lambda wd: wd/2.355*np.sqrt(2*np.pi)/dvv*dvv
        W12=r.x[0]*n(r.x[3]); W3=r.x[1]*n(r.x[3]); kq=kb(Q123[0],9.357,6.284)*1e-3
        # LTE prediction of W3/W12 (thin): ratio of nu*Smu2*exp(-Eu/T)
        pr=lte_w(Q123[2],T)/(lte_w(Q123[0],T)+lte_w(Q123[1],T))
        out.append(dict(model=f'Q J=1-3 (J1:J2 at LTE {T:.0f} K)',v0=r.x[2],e_v0=np.sqrt(cv[2,2]),fwhm=r.x[3],e_fwhm=np.sqrt(cv[3,3]),a_mJy=np.nan,W12_Kkms=W12*kq,eW12=np.sqrt(cv[0,0])*n(r.x[3])*kq,W3_Kkms=W3*kq,eW3=np.sqrt(cv[1,1])*n(r.x[3])*kq,
                        W3_over_W12=W3/W12 if W12!=0 else np.nan,W3_over_W12_LTE=pr,chi2red=2*r.cost/sg**2/(m.sum()-6),rms_mJy=sg))
    return pd.DataFrame(out)
# ---------------- B) LTE with error floor and recovery factor
def lte_scan(floor=0.12,recov=(1.0,0.7,0.5)):
    use=M[(M.ous=='Xb55')|((M.ous=='Xb4f')&(M.line.isin(['2(0)-1(-1) E','5(-2,4)-4(-1,3) E'])))].copy()
    rows=[]
    for f in recov:
        y=use.I_Kkms.values.copy(); y=np.where(use.ous=='Xb4f',y/f,y); e=np.hypot(use.eI_Kkms.values*np.where(use.ous=='Xb4f',1/f,1),0.10*np.abs(y))
        hi=use.line.isin(['Q J=4','Q J=5','Q J=6','Q J=7','Q J=8']).values; e=np.where(hi,np.hypot(e,floor),e)
        for T in np.array([3,4,5,6,7,8,9,10,12,14,17,20,25,30,40,60,100.]):
            w=np.array([sum(lte_w(f0,T) for f0 in ((Q123) if r.line=='Q J=1-3' else (r.nu_rest,))) for r in use.itertuples()])
            for excl in (False,True):
                k=np.ones(len(y),bool)
                if excl: k[(use.line=='Q J=1-3').values]=False
                A=np.sum((w*y/e**2)[k])/np.sum((w*w/e**2)[k]); chi=np.sum((((y-A*w)/e)**2)[k]); q=np.where((use.line=='Q J=1-3').values)[0][0]
                rows.append(dict(recov=f,Trot=T,exclude_Q123=excl,chi2=chi,ndata=int(k.sum()),A=A,Q123_pred=A*w[q],Q123_obs=y[q],Q123_obs_over_pred=y[q]/(A*w[q]),eQ123_rel=e[q]/(A*w[q])))
    return pd.DataFrame(rows)
# ---------------- C) RADEX grid
def rdx(args):
    T,n,N=args; wd=f'/tmp/rxm_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/e-ch3oh.dat\n{wd}/o.txt\n150 300\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n150\n0\n'; subprocess.run(['radex'],input=inp,text=True,capture_output=True); out=[]
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit():
            try: out.append((float(p[3]),float(p[5]),float(p[6]),float(p[10])))   # freq, Tex, tau, K km/s
            except Exception: pass
    return np.array(out)
FREQ={'Q1':165.050175,'Q2':165.061130,'Q3':165.099240,'J4':165.190475,'J5':165.369341,'J6':165.678649,'J7':166.169098,'J8':166.898566,'2_0-1_-1':254.015377,'5_-2-4_-1':266.838148}
def radex_grid(floor=0.12,recov=(1.0,0.7,0.5)):
    TKs=[5,8,10,15,20,30,50,80,150]; LNs=np.arange(2,8.01,0.5); LNc=np.arange(12,16.51,0.25)
    jobs=list(itertools.product(TKs,10**LNs,10**LNc))
    path=OUT+'ch3oh_radex_raw.npz'
    if os.path.exists(path): raw=np.load(path,allow_pickle=True)['raw']
    else:
        with ProcessPoolExecutor(4) as ex: raw=list(ex.map(rdx,jobs,chunksize=30))
        np.savez(path,raw=np.array(raw,dtype=object))
    rows=[]
    use=M[(M.ous=='Xb55')|((M.ous=='Xb4f')&(M.line.isin(['2(0)-1(-1) E','5(-2,4)-4(-1,3) E'])))]
    names=['Q123','J4','J5','J6','J7','J8','2_0-1_-1','5_-2-4_-1']
    for f in recov:
        y=use.I_Kkms.values.copy(); isb6=(use.ous=='Xb4f').values; y=np.where(isb6,y/f,y); e=np.hypot(use.eI_Kkms.values*np.where(isb6,1/f,1),0.10*np.abs(y)); hi=use.line.isin(['Q J=4','Q J=5','Q J=6','Q J=7','Q J=8']).values; e=np.where(hi,np.hypot(e,floor),e)
        for (T,n,N),r in zip(jobs,raw):
            if len(r)==0: continue
            def g(f0,col): k=np.argmin(abs(r[:,0]-f0)); return r[k,col] if abs(r[k,0]-f0)<0.02 else np.nan
            pr=np.array([sum(g(F,3) for F in (FREQ['Q1'],FREQ['Q2'],FREQ['Q3'])),g(FREQ['J4'],3),g(FREQ['J5'],3),g(FREQ['J6'],3),g(FREQ['J7'],3),g(FREQ['J8'],3),g(FREQ['2_0-1_-1'],3),g(FREQ['5_-2-4_-1'],3)])
            taus=[g(FREQ[k],2) for k in ('Q1','Q2','Q3')]
            if not np.all(np.isfinite(pr)): continue
            for excl in (False,True):
                k=np.ones(len(y),bool)
                if excl: k[0]=False
                eta=np.clip(np.sum((pr*y/e**2)[k])/max(np.sum((pr*pr/e**2)[k]),1e-30),0,1.); chi=np.sum((((y-eta*pr)/e)**2)[k])
                rows.append(dict(recov=f,exclude_Q123=excl,Tkin=T,log_n=np.log10(n),log_N=np.log10(N),eta=eta,chi2=chi,Q123_pred=eta*pr[0],Q123_obs=y[0],tauQ1=taus[0],tauQ2=taus[1],tauQ3=taus[2],pred_2_0=eta*pr[6]))
    return pd.DataFrame(rows)
if __name__=='__main__':
    P=profile(); P.round(4).to_csv(OUT+'ch3oh_profile.csv',index=False); pd.set_option('display.width',280,'display.max_columns',40); print(P.round(3).to_string())
    S=lte_scan(); S.round(5).to_csv(OUT+'ch3oh_lte_scan.csv',index=False)
    for f in (1.0,0.7,0.5):
        for ex in (False,True):
            s=S[(S.recov==f)&(S.exclude_Q123==ex)]; i=s.chi2.idxmin(); b=s.loc[i]
            print(f'LTE recov={f} exclude_Q123={ex}: Trot={b.Trot} chi2={b.chi2:.1f}/{b.ndata-2} Q123 obs/pred={b.Q123_obs_over_pred:.2f}+-{b.eQ123_rel*b.Q123_obs_over_pred:.2f}')
    G=radex_grid(); G.to_csv(OUT+'ch3oh_radex_grid.csv',index=False)
    for f in (1.0,0.7,0.5):
        for ex in (False,True):
            g=G[(G.recov==f)&(G.exclude_Q123==ex)]; b=g.loc[g.chi2.idxmin()]
            print(f'RADEX recov={f} excl_Q123={ex}: min chi2 {b.chi2:.1f} at T={b.Tkin} logn={b.log_n} logN={b.log_N} eta={b.eta:.3g} Q123 obs/pred={b.Q123_obs/b.Q123_pred:.2f} tauQ=({b.tauQ1:.3g},{b.tauQ2:.3g},{b.tauQ3:.3g})')
        g=G[(G.recov==f)&(~G.exclude_Q123)]; m=g.chi2.min(); s=g[g.chi2<=m+8]; print('   all-data 2sig (3par):',len(s),'models; min tau over Q lines',s[['tauQ1','tauQ2','tauQ3']].min().min())
        g=G[(G.recov==f)&(G.exclude_Q123)]; m=g.chi2.min(); s=g[g.chi2<=m+4]; r=(s.Q123_obs/s.Q123_pred.replace(0,np.nan)); print('   without Q123 (d chi2<4): predicted Q123 range %.2f-%.2f (obs %.2f); obs/pred range %.2f-%.1f; models with tau<0 on any Q line: %d of %d'%(s.Q123_pred.min(),s.Q123_pred.max(),s.Q123_obs.iloc[0],r.min(),r.max(),int(((s[['tauQ1','tauQ2','tauQ3']]<0).any(axis=1)).sum()),len(s)))
    # inversion search: does ANY (Tkin,n,N) in the grid give negative tau (or Tex>Tkin-like population inversion) on the 165 GHz Q lines or 2(0)-1(-1)?
    g=G[(G.recov==1.0)&(~G.exclude_Q123)]
    print('grid-wide tau(Q J=1,2,3) minima:',g.tauQ1.min(),g.tauQ2.min(),g.tauQ3.min(),'; number of grid cells with any negative Q tau:',int((g[['tauQ1','tauQ2','tauQ3']]<0).any(axis=1).sum()),'of',len(g))
    neg=g[(g[['tauQ1','tauQ2','tauQ3']]<0).any(axis=1)]
    if len(neg): print(neg.groupby('Tkin').size().to_string()); print(neg[['Tkin','log_n','log_N','tauQ1','tauQ2','tauQ3']].head(10).to_string())
