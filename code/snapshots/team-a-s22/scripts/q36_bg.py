# Q36 step 3: frequency-dependent background mean intensity J_bg(nu) [K] seen by the HCN/HCO+ gas in each core: CMB + local continuum.
# Tb_RJ(258 GHz): S core intrinsic 0.80 K (deconvolved 0.31x0.26", S_int=3.5 mJy, Xb4b cont), N core 0.18 K (0.24" beam peak; N is resolved 1.6x1.1", mean 0.11 K).
# dust Tb ~ nu^beta (beta=1.8, thin), free-free Tb ~ nu^-2.1; dust fraction at 258 GHz from q36_sed (0.5 +-0.1; ~1 for a compact core if ff is extended: 'dusty' bracket).
# geometry eta: mean intensity at the centre of a uniform emitting sphere = 0.5 x line-of-sight (diameter) brightness; 1.0 = gas fully inside/behind the whole column.
import numpy as np
H=6.62607e-27; KB=1.380649e-16
def Jp(nu,T): x=H*nu*1e9/(KB*T); return (H*nu*1e9/KB)/np.expm1(x)
def Tbg_eq(nu,J): return (H*nu*1e9/KB)/np.log1p((H*nu*1e9/KB)/J)       # blackbody T giving mean intensity J at nu [GHz]
TB258={'S':0.80,'N':0.18}
ENVTB={'envmid':0.02,'envhi':0.07}   # local envelope dust Tb_RJ at 258 GHz [K]; observed 10.5" beam-mean is 0.005-0.006 K (q36_cont), so these assume sub-beam dust filling 0.3 / 0.08
SCEN={'cmb':dict(eta=0,fd=0,clump=1),'nominal':dict(eta=0.7,fd=0.5,clump=1),'dusty':dict(eta=1.0,fd=1.0,clump=1),'high':dict(eta=1.0,fd=1.0,clump=2.0)}
def Jbg(reg,scen,nu_GHz,beta=1.8):
    if reg=='E': return Jp(np.asarray(nu_GHz,float),2.725)+ENVTB[scen]*(np.asarray(nu_GHz,float)/258.)**beta
    s=SCEN[scen]; nu=np.asarray(nu_GHz,float); tb=TB258[reg]*s['clump']
    return Jp(nu,2.725)+s['eta']*tb*(s['fd']*(nu/258.)**beta+(1-s['fd'])*(nu/258.)**-2.1)
if __name__=='__main__':
    nus=[88.63,177.26,265.89,354.51]; L=['J_bg [K] and equivalent blackbody Tbg at HCN 1-0..4-3:']
    for reg in 'SN':
        for sc in SCEN:
            J=Jbg(reg,sc,nus); L.append(f'{reg} {sc:8s} J_bg='+' '.join(f'{j:.3f}' for j in J)+'  Tbg_eq='+' '.join(f'{Tbg_eq(n,j):.2f}' for n,j in zip(nus,J)))
    for sc in ENVTB:
        J=Jbg('E',sc,nus); L.append(f'ENV {sc:8s} J_bg='+' '.join(f'{j:.3f}' for j in J)+'  Tbg_eq='+' '.join(f'{Tbg_eq(n,j):.2f}' for n,j in zip(nus,J)))
    for T in (2.725,5.,10.,20.): L.append(f'uniform Tbg={T}: J_bg='+' '.join(f'{Jp(n,T):.3f}' for n in nus))
    print('\n'.join(L)); open('/workspace/out/q36/q36_bg.log','w').write('\n'.join(L))
