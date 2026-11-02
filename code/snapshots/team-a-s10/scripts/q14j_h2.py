# Q14j (fix 4): dust-based N(H2) at matched beams and abundance proxies N(X)/N(H2).
# Continuum (pbcor) of Xb55 (B5 ~170 GHz), Xb4f (B6 258 GHz), Xb47 (B7 346.5 GHz) smoothed to common beams (6.3x3.8 PA82 for B6/B7; 9.5x6.4 PA ~ for all three) and sampled at the regions.
# Free-free/synchrotron removal: S_nu = D(nu/nu0)^(2+beta) + F (nu/nu0)^alpha_f ; MC over beta 1.5-2, alpha_f -0.4..-0.1 (Harada+2018 3mm indices N -0.39, S -0.37), flux cal 5% (Harada) / 7%.
# N(H2)=S_dust/(Omega_beam B_nu(Td) kappa_nu mu m_H); kappa_nu=0.1 (nu/1THz)^beta cm2 per g of gas (gas-to-dust 100 included), mu=2.8; Td 20-50 K (assumption).
# Output out/q14e/h2_*.csv. python3 -I q14j_h2.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
import q14f_lte as q
def img(ous,name,tb):
    p=f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.{name}.cont.I.pbcor.fits'; h=fits.open(p,memmap=True); hd=h[0].header; d=np.nan_to_num(np.array(h[0].data).squeeze(),nan=0.).astype(np.float32)
    nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); k=kernel(hd,tb,nb)
    if k is None: return None
    d=fftconvolve(d,k.astype(np.float32),mode='same')*OM(tb)/OM(nb); w=WCS(hd).celestial; nu=hd.get('CRVAL3',None)
    return d,w,nb
FILES={'B5':('Xb55','spw16_18_20_22',170.0),'B6':('Xb4f','spw16_18_20_22',257.7),'B7':('Xb47','spw16_18_20_22',346.53)}
REGS=['N','S','Bridge','Wclump','Soutflow']
rng=np.random.default_rng(11)
def hB(nu,T): x=4.799e-11*nu*1e9/T; return 2*6.626e-27*(nu*1e9)**3/2.998e10**2/np.expm1(x)
res=[]
for tbn,tb in (('6.3x3.8',(6.3,3.8,82.)),('9.6x6.5',(9.6,6.5,82.))):
    S={}
    for b,(ous,nm,nu) in FILES.items():
        r=img(ous,nm,tb)
        if r is None: continue
        d,w,nb=r
        for reg in REGS:
            px,py=[float(z) for z in w.world_to_pixel_values(REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg'])]
            S[(b,reg)]=float(map_coordinates(d,[[py],[px]],order=1)[0])
        # rms from corner
        S[(b,'rms')]=float(np.std(d[:30,:30]))
    om=OM(tb)*2.3504e-11   # sr
    for reg in REGS:
        s6,s7=S[('B6',reg)],S[('B7',reg)]; s5=S.get(('B5',reg),np.nan)
        # MC two-band split (B6,B7) and with B5 (if available)
        for Td in (20.,30.,50.):
            Ns=[];fd=[]
            for k in range(4000):
                beta=rng.uniform(1.5,2.0); af=rng.uniform(-0.4,-0.1); c6,c7=rng.normal(1,0.07,2)
                a6,a7=s6*c6,s7*c7; r=(346.53/257.7)
                D6=(a7-a6*r**af)/(r**(2+beta)-r**af); D7=max(D6,0)*r**(2+beta)
                kap=0.1*(346.53/1000)**beta
                NH2=D7*1e-23/(om*hB(346.53,Td)*kap*2.8*1.6726e-24)
                Ns.append(NH2); fd.append(D7/a7)
            Ns=np.array(Ns); res.append(dict(beam=tbn,region=reg,Td=Td,S258_mJy=s6*1e3,S346_mJy=s7*1e3,S170_mJy=s5*1e3,alpha_258_346=np.log(s7/s6)/np.log(346.53/257.7),dust_frac346_med=np.median(fd),dust_frac_lo=np.percentile(fd,16),dust_frac_hi=np.percentile(fd,84),
                                             NH2_med=np.median(Ns),NH2_16=np.percentile(Ns,16),NH2_84=np.percentile(Ns,84)))
H=pd.DataFrame(res); H.to_csv('/workspace/out/q14e/h2_dust.csv',index=False); pd.set_option('display.width',250); print(H.round(3).to_string())
