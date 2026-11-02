# Q49: recompute Q22e core mass fractions / path lengths from Q42 dynesty weighted samples (shared envelope, baseline, cmb), per seed.
# Same conventions as scripts/q22b_post.py (mu=2.8, D=40 Mpc, core area from HCN 3-2 FWHM, envelope = filling x 10" circle).
import glob,numpy as np
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
mu=2.8*1.6605e-24; pc=3.086e18; Msun=1.989e33; PC=194.; SIZE={'S':0.455,'N':0.507}
for tag in ('base_cmb','drop1_cmb'):
  fs=sorted(glob.glob(f'/workspace/out/q42/runs2/q42_{tag}_sh_s*.npz'))
  for Xab in (1e-9,1e-8):
    res={'N':[],'S':[]};Ls={'N':[],'S':[]};ME=[]
    for f in fs:
        d=np.load(f,allow_pickle=True);S,w,nm=d['S'],d['w']/d['w'].sum(),list(d['names']); P=lambda n:S[:,nm.index(n)]
        ne=10**P('ne');NHe=10**P('NHe')/Xab;we=10**P('lwe'); Me=we*NHe*mu*np.pi*(5*PC*pc)**2/Msun; ME.append(wq(Me,w,.5))
        for r in 'NS':
            NHc=10**P('NHc_'+r)/Xab; nc=10**P('nc_'+r); fc=10**P('lfc_'+r); Lc=NHc/nc/pc
            Mc=fc*NHc*mu*np.pi*(SIZE[r]/2*PC*pc)**2/Msun; fr=Mc/(Mc+Me)
            res[r].append(wq(fr,w,[.16,.5,.84])); Ls[r].append(wq(Lc,w,.5))
    for r in 'NS':
        a=np.array(res[r]); print(tag,'X=%g'%Xab,r,'nseeds',len(fs),'core frac med per seed',a[:,1].round(4),'mean[16,84] %.3f [%.3f,%.3f]'%tuple(a.mean(0)[[1,0,2]]),'Lc med pc',np.round(Ls[r],2),'Menv',['%.2g'%m for m in ME])
