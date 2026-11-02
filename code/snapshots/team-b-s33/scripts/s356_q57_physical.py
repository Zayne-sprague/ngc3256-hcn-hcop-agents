# Q57 round 2: physical priors on the s353 likelihood cubes, evidence/chi2/A_N/A_S/contrast-profile per hypothesis.
# Constraints (mask W[n,NO] on the (log n, log N(HCO+)) plane): log n >= NMIN, log n <= NMAX, and path length L = N(HCO+)/(xmax n) <= LMAX
#   (xmax = [HCO+]/[H2] upper bound 1e-8; N(HCO+) is the in-clump column since the filling factor eta multiplies the model).
# usage: python3 -I s356_q57_physical.py [SYS MOD]  -> out/q57/s356_s<SYS>_m<MOD>.txt
import sys,os,numpy as np
a=sys.argv[1:]
if a: os.environ['Q57_SYS'],os.environ['Q57_MOD']=a[0],a[1]; sys.argv=['x']
sys.path.insert(0,'/workspace/scripts')
import s353_q57_joint_bracket as S
PC=3.086e18
CN0=S.build('N').astype(float); CS0=S.build('S').astype(float)   # [g,X,T,n,NO,A]
LXg=S.LXg; GW=S.GW; LNs=S.LNs; NOg=S.NOg; LAg=S.LAg; nT=len(S.TK)
out=[]
def P(t): print(t,flush=True); out.append(t)
def mask(nmin,nmax,Lpc,xmax=1e-8):
    W=np.ones((len(LNs),len(NOg)))
    for i,ln in enumerate(LNs):
        for j,no in enumerate(NOg):
            if ln<nmin-1e-6 or ln>nmax+1e-6: W[i,j]=0
            if Lpc and no>np.log10(xmax)+ln+np.log10(Lpc*PC): W[i,j]=0
    return W
def analyse(tag,xlo,xhi,W,shared_contrast=True):
    xs=(LXg>=xlo-1e-6)&(LXg<=xhi+1e-6); nX=xs.sum(); nA=len(LAg)
    CN=CN0[:,xs]*W[None,None,None,:,:,None]; CS=CS0[:,xs]*W[None,None,None,:,:,None]
    cntW=W.sum(); cntO=W.sum(0)   # allowed (n) per NO
    def Z(i):
        N_,S_=CN[i],CS[i]   # [X,T,n,NO,A]
        z={}
        z['Hexc']=((N_.sum((1,2))*S_.sum((1,2))).sum())/(nT**2*(cntO**2).sum()*nA*nX)
        z['Hcol']=(N_.sum((1,2,3))*S_.sum((1,2,3))).sum()/((nT*cntW)**2*nA*nX)
        z['HcolX']=(N_.sum((0,1,2,3))*S_.sum((0,1,2,3))).sum()/((nT*cntW*nX)**2*nA)
        z['Habun']=(N_.sum(4)*S_.sum(4)).sum()/(nT*nX*cntW*nA**2)
        z['HfullXs']=(N_.sum((1,2,3,4))*S_.sum((1,2,3,4))).sum()/((nT*cntW*nA)**2*nX)
        z['Hfull']=N_.sum()*S_.sum()/((nT*cntW*nA*nX)**2)
        return z
    zz=[Z(i) for i in range(3)]; Zf={k:sum(GW[i]*zz[i][k] for i in range(3)) for k in zz[0]}
    N_,S_=CN[1],CS[1]
    mx={'Hfull':N_.max()*S_.max(),'HfullXs':(N_.max((1,2,3,4))*S_.max((1,2,3,4))).max(),
        'Hexc':(N_.max((1,2))*S_.max((1,2))).max(),'Hcol':(N_.max((1,2,3))*S_.max((1,2,3))).max(),
        'HcolX':(N_.max((1,2,3)).max(0)*S_.max((1,2,3)).max(0)).max(),'Habun':(N_.max(4)*S_.max(4)).max()}
    P(f'## {tag}: X {10**xlo:.0f}-{10**xhi:.0f}; allowed (n,NO) cells {int(cntW)}/{W.size}')
    P('   BF vs Hfull: '+', '.join(f'{k} {v/Zf["Hfull"]:.3g}' for k,v in Zf.items())+f'  | BF(Habun/Hexc)={Zf["Habun"]/Zf["Hexc"]:.3g}, Habun/Hcol={Zf["Habun"]/Zf["Hcol"]:.3g}')
    P('   delta(-2lnL) vs Hfull: '+', '.join(f'{k} {2*np.log(mx["Hfull"]/v):.1f}' for k,v in mx.items()))
    # best-cell location under Hfull
    for nm,C in (('N',CN[1]),('S',CS[1])):
        b=np.unravel_index(C.argmax(),C.shape); P(f'   Hfull best {nm}: Tk={S.TK[b[1]]:.0f} log n={LNs[b[2]]:.1f} logN(HCO+)={NOg[b[3]]:.1f} A={10**LAg[b[4]]:.2f} X={10**LXg[xs][b[0]]:.0f}')
    pAN=CN.sum((0,1,2,3,4)); pAN/=pAN.sum(); pAS=CS.sum((0,1,2,3,4)); pAS/=pAS.sum()
    q=lambda p:10**LAg[np.minimum(np.searchsorted(np.cumsum(p),[0.16,0.5,0.84]),nA-1)]
    J=pAN[:,None]*pAS[None,:]; ratio=(LAg[:,None]-LAg[None,:]).ravel(); w=J.ravel(); o=np.argsort(ratio); cw=np.cumsum(w[o])/w.sum()
    r=10**ratio[o][np.searchsorted(cw,[0.05,0.16,0.5,0.84,0.95])]
    P(f'   Hfull A_N {q(pAN)}, A_S {q(pAS)}; A_N/A_S (5,16,50,84,95%) {np.round(r,2)}; P(>1)={w[ratio>1e-3].sum()/w.sum():.3f} P(>1.3)={w[ratio>np.log10(1.3)].sum()/w.sum():.3f}')
    # contrast profile (1-0 HCN/HCO+ model contrast, g=0): Hexc and Hcol; bins in log10 contrast
    bins=np.arange(-0.4,0.41,0.02); nb=len(bins)
    prof={'Hexc':np.zeros(nb),'Hcol':np.zeros(nb)}
    xl=LXg[xs]; ratio_c={}
    for ix,X in enumerate(xl):
        M=S.model_cube(X)     # [6,T,n,NO,A] (eta cancels in ratio)
        lc=np.log10(M[0]/M[2]); lc3=np.log10(M[1]/M[3])
        bi=np.clip(np.round((lc-bins[0])/0.02).astype(int),0,nb-1)
        for mode in ('Hexc','Hcol'):
            for iA in range(nA):
                # per nucleus: max likelihood per contrast-bin, over T,n (and NO if Hcol)
                if mode=='Hexc':
                    for iO in range(len(NOg)):
                        fN=np.zeros(nb);fS=np.zeros(nb)
                        np.maximum.at(fN,bi[:,:,iO,iA].ravel(),N_[ix,:,:,iO,iA].ravel())
                        np.maximum.at(fS,bi[:,:,iO,iA].ravel(),S_[ix,:,:,iO,iA].ravel())
                        prof[mode]=np.maximum(prof[mode],maxplus(fN,fS,bins))
                else:
                    fN=np.zeros(nb);fS=np.zeros(nb)
                    np.maximum.at(fN,bi[:,:,:,iA].ravel(),N_[ix,:,:,:,iA].ravel())
                    np.maximum.at(fS,bi[:,:,:,iA].ravel(),S_[ix,:,:,:,iA].ravel())
                    prof[mode]=np.maximum(prof[mode],maxplus(fN,fS,bins))
    # profile over contrast c=cN/cS; compare to global max (Hfull-like) -> delta(-2lnL)
    gmax=mx['Hfull']
    for mode in prof:
        d=2*np.log(gmax/np.maximum(prof[mode],1e-300)); P(f'   profile 1-0 model contrast {mode} (delta -2lnL vs Hfull max; min {d.min():.1f} at c={10**bins[d.argmin()]:.2f}): '+' '.join(f'{10**b:.2f}:{x:.1f}' for b,x in zip(bins[::3],d[::3]) if 0.95<10**b<1.8))
        # contrast at which delta relative to hypothesis min = 1,4 (upper)
        dd=d-d.min(); ok1=bins[(dd<=1)]; ok4=bins[dd<=4]
        P(f'     {mode}: contrast within dchi2<=1 of its best: {10**ok1.min():.2f}-{10**ok1.max():.2f}; <=4: {10**ok4.min():.2f}-{10**ok4.max():.2f}; observed 1-0 contrast 1.45')
def maxplus(fN,fS,bins):
    # out[b]= max over iN-iS in bin b of fN[iN]*fS[iS]
    nb=len(bins); o=np.zeros(nb)
    for i in range(nb):
        if fN[i]==0: continue
        j=np.arange(nb); k=np.round((bins[i]-bins[j]-bins[0])/0.02).astype(int); ok=(k>=0)&(k<nb)
        np.maximum.at(o,k[ok],fN[i]*fS[j[ok]])
    return o
if __name__=='__main__':
    for xlo,xhi in ((1.6,2.0),(1.3,2.3)):
        analyse('unconstrained',xlo,xhi,mask(3.0,8.0,None))
        analyse('log n>=3.8',xlo,xhi,mask(3.8,8.0,None))
        analyse('log n>=4.2',xlo,xhi,mask(4.2,8.0,None))
        analyse('log n 3-6 flat (alt prior)',xlo,xhi,mask(3.0,6.0,None))
        analyse('log n>=3.8, L<=450pc, x<=1e-8',xlo,xhi,mask(3.8,8.0,450))
        analyse('log n>=3.8, L<=100pc (N core), x<=1e-8',xlo,xhi,mask(3.8,8.0,100))
        analyse('log n>=4.2, L<=100pc, x<=1e-8',xlo,xhi,mask(4.2,8.0,100))
        analyse('log n>=3.8, L<=100pc, x<=3e-9 (tight)',xlo,xhi,mask(3.8,8.0,100,3e-9))
    tag=f"s{S.FSYS}_m{S.FMOD}"; open(S.O+f's356_{tag}.txt','w').write('\n'.join(out))
