# Q25: chi2 fit of the four observables per annulus [R10, R21, HCN 2-1/1-0, HCO+ 2-1/1-0] (10" Tb, Michiyama 1-0 + Xb55 2-1) on the RADEX grid (T,n,N_HCN,N_HCO+).
# errors: R10 and R21 4% (intra-band: noise+beam), ladders 14% (1-0 vs 2-1 cross-band flux scale ~10% each). Report best fit and the range of logX, log n within dchi2<2.3 (2 params) per fixed Tk.
import sys; sys.path.insert(0,'/workspace/scripts')
from q25_radex import *
H=table('hcn'); O=table('hco+'); d=json.load(open('/workspace/out/q25/q25_ratio10.json')); f=1.0126
src=sys.argv[1] if len(sys.argv)>1 else 'Michiyama'
h10=np.array(d[f'{src}|HCN 1-0']['F']); o10=np.array(d[f'{src}|HCO+ 1-0']['F']); h21=np.array(d['Xb55|HCN 2-1']['F']); o21=np.array(d['Xb55|HCO+ 2-1']['F'])
OBS=np.array([h10/o10*f,h21/o21*f,h21/h10*.25,o21/o10*.25]); ERR=np.array([.04,.04,.14,.14])
R10=H[...,None,0]/O[...,None,:,0]; R21=H[...,None,1]/O[...,None,:,1]
LH=(H[...,1]/H[...,0])[...,None]*np.ones(len(LNC)); LO=((O[...,1]/O[...,0])[...,None,:])*np.ones((len(LNC),1))
P=np.stack([R10,R21,LH,LO],-1); P=np.where(np.isfinite(P),P,1e9)
LX=LNC[:,None]-LNC[None,:]
for ia,name in ((0,'0-4'),(1,'4-8'),(3,'8-12')):
    c=(((P-OBS[:,ia])/(ERR*OBS[:,ia]))**2).sum(-1)   # [T,n,NH,NO]
    i=np.unravel_index(np.argmin(c),c.shape); print(f'{src} {name}": obs R10,R21,LH,LO={np.round(OBS[:,ia],3)}; global best chi2={c[i]:.1f} Tk={TK[i[0]]:.0f} logn={LN[i[1]]:.2f} logNH={LNC[i[2]]:.2f} logNO={LNC[i[3]]:.2f} logX={LX[i[2],i[3]]:+.2f}')
    for it in (2,3,4):
        cc=c[it]; j=np.unravel_index(np.argmin(cc),cc.shape); ok=cc<=cc.min()+2.3
        n,a,b=np.where(ok); lx=LX[a,b]
        print(f'   Tk={TK[it]:.0f}: best chi2={cc.min():.1f} logn={LN[j[0]]:.2f} logX={LX[j[1],j[2]]:+.2f}; dchi2<2.3: logn {LN[n].min():.2f}..{LN[n].max():.2f}, logX {lx.min():+.2f}..{lx.max():+.2f}, logNHCO+ {LNC[b].min():.2f}..{LNC[b].max():.2f}')
