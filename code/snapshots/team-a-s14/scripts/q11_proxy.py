# Q11 forward model B: the observed Xb55 2-1 10" HCN and HCO+ maps (nearly fully recovered, see q11_profile f~0.2) are used as the TRUE 10" sky proxy;
# an extra intrinsic HCN fall eps*clip((r-4)/6,0,1) is imposed; the sky is multiplied by a Gaussian PB, filtered with the 3-2 transfer function
# T(u)=clip((u-ua)/ua,0,1), ua=f*Bmin/lambda(265.9 GHz), divided by PB, annuli measured. The predicted 3-2 annulus ratios (R0 scaled to the observed 0-4" value) are compared with observed.
import sys,json; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
c=299792458.; tb=(10.,10.,0.); ANN=[(0,4),(4,8),(8,12)]; XM=-0.075
def mp(key):
    d=load(key,nullspec=False); k=kernel(d['hd'],tb,d['beam']); s=fftconvolve(np.nan_to_num(d['m0']),k,mode='same')*OM(tb)/OM(d['beam']); return d,s
H,sH=mp('HCN 2-1'); O,sO=mp('HCO+ 2-1'); pix=H['pix']; omp=OM(tb)/pix**2
n=sH.shape[0]; npad=512
R=np.hypot(H['dx']-XM,H['dy']); Rc=np.hypot(H['dx'],H['dy'])
def pad(a): P=np.zeros((npad,npad)); P[:n,:n]=a; return P
u=np.fft.fftfreq(npad,d=pix/206265.)/1e3; U=np.hypot(*np.meshgrid(u,u))
PBd=lambda fw: np.exp(-4*np.log(2)*Rc**2/fw**2)
def ann(s): return np.array([s[(R>=a)&(R<b)].sum()/omp for a,b in ANN])
def filt(s,pbfw,ua,pbcut=0.05):
    pb=PBd(pbfw); im=np.real(np.fft.ifft2(np.fft.fft2(pad(s*pb))*np.clip((U-ua)/ua,0,1)))[:n,:n]
    return np.where(pb>pbcut,im/np.maximum(pb,pbcut),0)
obs={}
for tag,(a,b) in {'Xb4f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),'Xb4d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d')}.items():
    h,sh=mp(a); o,so=mp(b); r=np.hypot(h['dx']-XM,h['dy']); om=OM(tb)/h['pix']**2
    f=lambda s:np.array([s[(r>=x)&(r<y)].sum()/om for x,y in ANN]); obs[tag]=f(sh)/f(so)
print('observed 3-2 R',{k:v.round(3) for k,v in obs.items()},' observed 2-1 R',(ann(sH)/ann(sO)).round(3))
res=[]
for lab,pbfw,Bm in (('Xb4f 7m-like (PB37, Bmin 8.7m)',37.,8.7),('Xb4f 12m-like (PB22, Bmin 15m)',22.,15.),('Xb4d 12m (PB22, Bmin 15m)',22.,15.)):
  tag='Xb4d' if 'Xb4d' in lab else 'Xb4f'
  for f in (0.1,0.2,0.3,0.4):
    ua=f*Bm/(c/265.9e9)/1e3
    for eps in (0.,0.12,0.25,0.4):
        h=filt(sH*(1-eps*np.clip((Rc-4)/6,0,1)),pbfw,ua); o=filt(sO,pbfw,ua)
        a=ann(h); b=ann(o); Rm=a/b
        res.append(dict(setup=lab,f=f,ua=ua,eps=eps,flux_HCN=a.tolist(),flux_HCO=b.tolist(),R=Rm.tolist(),Rnorm=(Rm*obs[tag][0]/Rm[0]).tolist(),obs=obs[tag].tolist()))
        print(f'{lab:32s} f={f:.1f} ua={ua:4.1f}kl eps={eps:.1f} HCO+ flux {b.round(1)} HCN {a.round(1)} R {Rm.round(3)} -> R scaled to obs 0-4" {np.round(Rm*obs[tag][0]/Rm[0],3)} (obs {obs[tag].round(3)})')
json.dump(res,open('/workspace/out/q11/q11_proxy.json','w'),indent=1)
