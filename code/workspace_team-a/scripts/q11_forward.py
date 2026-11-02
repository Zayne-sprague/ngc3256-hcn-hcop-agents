# Q11 forward model of interferometer recovery for the 10" annulus HCN/HCO+ ratios (2-1 Xb55 vs 3-2 Xb4f/Xb4d).
# sky -> x primary beam (Gaussian, array-dependent) -> Fourier transfer function T(u) (0 below u_a=f*Bmin/lambda, linear ramp to 1 at 2u_a)
# -> / PB -> 10" circular Gaussian -> annulus fluxes about the N-S midpoint. Sky: HCO+ = N/S Gaussian cores + exponential disk (fit to the observed
# 2-1 HCO+ annulus profile); HCN = HCO+ x R(r), R(r)=R0*(1-eps*clip((r-4)/6,0,1)) (eps = true fractional fall of the ratio between r<=4 and r>=10").
# eps is then constrained by the observed 2-1 and 3-2 annulus ratios jointly (chi2).
import sys,json,itertools; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
c=299792458.; N=384; PIX=0.5
ANN=[(0,4),(4,8),(8,12),(12,16),(16,24)]
x=(np.arange(N)-N//2)*PIX; X,Y=np.meshgrid(x,x); XM,YM=-0.075,0.0
Rm=np.hypot(X-XM,Y-YM); Rc=np.hypot(X,Y)
u=np.fft.fftfreq(N,d=PIX/206265.)/1e3; U=np.hypot(*np.meshgrid(u,u))
def gauss(fw): s=fw/2.3548; return np.exp(-0.5*(np.hypot(X-0,Y)/s)**2)
def sky_hco(h,sc,fc):
    core=lambda y0: np.exp(-0.5*((X-0)**2+(Y-y0)**2)/(sc/2.3548)**2)
    cores=core(2.55)+core(-2.55); disk=np.exp(-Rc/h)
    cores/=cores.sum(); disk/=disk.sum(); return fc*cores+(1-fc)*disk
def Rr(r,R0,eps): return R0*(1-eps*np.clip((r-4)/6.,0,1))
def Tf(freq_ghz,Bmin,f,ramp=2.):
    ua=f*Bmin/(c/(freq_ghz*1e9))/1e3
    return np.clip((U-ua)/((ramp-1)*ua),0,1)
def observe(sky,Tm,pbfw):
    pb=np.exp(-4*np.log(2)*(X**2+Y**2)/pbfw**2)
    im=np.real(np.fft.ifft2(np.fft.fft2(sky*pb)*Tm))
    ok=pb>0.05; im=np.where(ok,im/np.maximum(pb,1e-3),0)
    K=np.fft.fft2(np.fft.ifftshift(gauss(10.)))*0+np.exp(-2*(np.pi*np.hypot(*np.meshgrid(u*1e3,u*1e3))/206265.)**2*(10./2.3548)**2)
    return np.real(np.fft.ifft2(np.fft.fft2(im)*K))
def annuli(img,zero=None):
    z=np.mean(img[(Rm>=zero[0])&(Rm<zero[1])]) if zero else 0
    return np.array([(img[(Rm>=a)&(Rm<b)]-z).sum() for a,b in ANN])
def obs_annuli(key,tb=(10.,10.,0.),fl=None):
    d=load(key,nullspec=False); k=kernel(d['hd'],tb,d['beam']); s=fftconvolve(np.nan_to_num(d['m0']),k,mode='same')*OM(tb)/OM(d['beam'])
    r=np.hypot(d['dx']-XM,d['dy']-YM); omp=OM(tb)/d['pix']**2
    z=np.nanmean(s[(r>=fl[0])&(r<fl[1])]) if fl else 0
    return np.array([(s[(r>=a)&(r<b)]-z).sum()/omp for a,b in ANN]),d
if __name__=='__main__':
    obs={}
    for tag,(kh,ko) in {'21':('HCN 2-1','HCO+ 2-1'),'32f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),'32d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d')}.items():
        obs[tag]=(obs_annuli(kh)[0],obs_annuli(ko)[0]); print(tag,'HCN',obs[tag][0].round(2),'HCO+',obs[tag][1].round(2),'R',(obs[tag][0]/obs[tag][1]).round(3))
    # setups: (name, freq, Bmin(m), PB FWHM")
    SET={'21':[('Xb55 7m-like',177.3,8.7,55.),('Xb55 12m-like',177.3,15.,33.)],
         '32f':[('Xb4f 7m-like',265.9,8.7,37.),('Xb4f 12m-like',265.9,15.,22.)]}
    results=[]
    # 1) morphology from the 2-1 HCO+ profile
    best={}
    for (nm,fr,Bm,pbfw) in SET['21']:
      for f in (0.5,1.0,1.5):
        Tm=Tf(fr,Bm,f); bb=None
        for h,sc,fc in itertools.product((4,6,8,12,18,26),(3,5,8),(0.1,0.25,0.4,0.55,0.7)):
            m=annuli(observe(sky_hco(h,sc,fc),Tm,pbfw)); y=obs['21'][1]
            A=(m*y).sum()/(m*m).sum(); chi=(((A*m-y)/(0.03*y+0.5))**2).sum()
            if bb is None or chi<bb[0]: bb=(chi,h,sc,fc,A)
        best[(nm,f)]=bb; print('morph',nm,'f',f,'chi2 %.1f h=%d sc=%d fc=%.2f'%bb[:4])
    # 2) for each (2-1 setup, f21) and (3-2 setup, f32) predict ratio profiles for eps grid, using the 2-1-fitted morphology
    EPS=(0.,0.1,0.2,0.3,0.4)
    sig21=np.array([0.01,0.01,0.015,0.03,0.04]); sig32=np.array([0.015,0.015,0.1,0.15,0.2])
    out=[]
    for (nm,fr,Bm,pbfw) in SET['21']:
     for f in (0.5,1.0,1.5):
      chi,h,sc,fc,A=best[(nm,f)]; sky=sky_hco(h,sc,fc)
      for (n3,fr3,B3,pb3) in SET['32f']:
       for f3 in (0.5,1.0,1.5):
        T3=Tf(fr3,B3,f3)
        for eps in EPS:
          R0r=obs['21'][0][0]/obs['21'][1][0]
          # R0 is set per band so that the 0-4" ratio matches (J-dependent excitation offset is a free constant)
          for tag,TT,pbw,fo,nmm in (('21',Tf(fr,Bm,f),pbfw,fr,nm),('32',T3,pb3,fr3,n3)):
            pass
          def pred(TT,pbw,idx):
              hco=observe(sky,TT,pbw); hcn=observe(sky*Rr(Rc,1.,eps),TT,pbw)
              return annuli(hcn)/annuli(hco)  # relative to R0=1
          p21=pred(Tf(fr,Bm,f),pbfw,0); p32=pred(T3,pb3,1)
          o21=obs['21'][0]/obs['21'][1]; o32=obs['32f'][0]/obs['32f'][1]
          # free R0 per band (excitation/abundance offset) fitted by weighted lstsq
          def fit(p,o,s): w=1/s**2; R0=(w*p*o).sum()/(w*p*p).sum(); return R0,(((R0*p-o)/s)**2).sum(),R0*p
          R21,c21,m21=fit(p21,o21,sig21); R32,c32,m32=fit(p32,o32,sig32)
          out.append(dict(s21=nm,f21=f,s32=n3,f32=f3,eps=eps,chi21=c21,chi32=c32,m21=m21.round(3).tolist(),m32=m32.round(3).tolist(),R21=R21,R32=R32))
    json.dump(dict(obs={k:[v[0].tolist(),v[1].tolist()] for k,v in obs.items()},best={f'{k[0]}|{k[1]}':list(v) for k,v in best.items()},grid=out),open('/workspace/out/q11/q11_forward.json','w'),indent=1,default=float)
    import collections
    print('\n eps  : min chi2(2-1) / min chi2(3-2) / min joint over recovery settings; and best model ratios 8-12" (2-1, 3-2)')
    for eps in EPS:
        g=[o for o in out if o['eps']==eps]; b=min(g,key=lambda o:o['chi21']+o['chi32'])
        print(f"eps={eps:.1f} chi21 min {min(o['chi21'] for o in g):6.1f}  chi32 min {min(o['chi32'] for o in g):7.1f} joint min {b['chi21']+b['chi32']:7.1f} [{b['s21']} f{b['f21']} / {b['s32']} f{b['f32']}] m21 {b['m21']} m32 {b['m32']}")
