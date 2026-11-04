"""F4 helper: full profile-likelihood curves (delta(-2lnL) vs model 1-0 HCN/HCO+ N/S contrast, equal A) from the q57 likelihood cubes.
Same code as s356_q57_physical.py:analyse (contrast-profile block), but the whole curve is saved (s356 printed only 4 points)."""
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s356_q57_physical as M
S=M.S
def curves(xlo,xhi,W):
    xs=(M.LXg>=xlo-1e-6)&(M.LXg<=xhi+1e-6); nA=len(M.LAg)
    CN=M.CN0[:,xs]*W[None,None,None,:,:,None]; CS=M.CS0[:,xs]*W[None,None,None,:,:,None]
    N_,S_=CN[1],CS[1]
    gmax=N_.max()*S_.max()
    bins=np.arange(-0.4,0.41,0.02); nb=len(bins)
    prof={'Hexc':np.zeros(nb),'Hcol':np.zeros(nb)}
    for ix,X in enumerate(M.LXg[xs]):
        Mc=S.model_cube(X); lc=np.log10(Mc[0]/Mc[2]); bi=np.clip(np.round((lc-bins[0])/0.02).astype(int),0,nb-1)
        for iA in range(nA):
            for iO in range(len(M.NOg)):
                fN=np.zeros(nb);fS=np.zeros(nb)
                np.maximum.at(fN,bi[:,:,iO,iA].ravel(),N_[ix,:,:,iO,iA].ravel()); np.maximum.at(fS,bi[:,:,iO,iA].ravel(),S_[ix,:,:,iO,iA].ravel())
                prof['Hexc']=np.maximum(prof['Hexc'],M.maxplus(fN,fS,bins))
            fN=np.zeros(nb);fS=np.zeros(nb)
            np.maximum.at(fN,bi[:,:,:,iA].ravel(),N_[ix,:,:,:,iA].ravel()); np.maximum.at(fS,bi[:,:,:,iA].ravel(),S_[ix,:,:,:,iA].ravel())
            prof['Hcol']=np.maximum(prof['Hcol'],M.maxplus(fN,fS,bins))
    return {k:(2*np.log(gmax/np.maximum(v,1e-300))).tolist() for k,v in prof.items()}|{'bins':bins.tolist()}
out={}
out['X40-100 unconstrained']=curves(1.6,2.0,M.mask(3.0,8.0,None)); print('a',flush=True)
out['X40-100 log n>=3.8']=curves(1.6,2.0,M.mask(3.8,8.0,None)); print('b',flush=True)
out['X20-200 unconstrained']=curves(1.3,2.3,M.mask(3.0,8.0,None)); print('c',flush=True)
json.dump(out,open('f4_profiles.json','w'))
