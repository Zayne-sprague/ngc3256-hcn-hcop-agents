# Q40c: Harada (2015.1.00412) vs Michiyama (2015.1.00993) 1-0: uv-coverage diagnostic from the image products.
# FT of the moment-0 map on a common grid (beam-deconvolved analytically), azimuthal mean |V|(u) per line; ratio H/M shows where Harada loses power (effective shortest baseline).
# Then: apply a common uv high-pass (u>umin klambda, applied to both 1-0 maps, deconvolved-beam-corrected) and a common low-pass to 10" and re-measure plain R at N,S,E1,W; also a constant-offset match.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import q25_ratio10 as Q
from scipy.ndimage import map_coordinates
import json
NG,PIXG=320,0.4   # common grid: 128"
def cov2(maj,mnr,pa):
    th=np.deg2rad(pa); u=np.array([np.sin(th),np.cos(th)]); v=np.array([np.cos(th),-np.sin(th)])
    return (maj**2/(8*np.log(2)))*np.outer(u,u)+(mnr**2/(8*np.log(2)))*np.outer(v,v)
def regrid(d,arr):
    # native image -> common grid centred on RA0,DE0 (offsets in arcsec; x = east)
    pix=d['pix']; ny,nx=arr.shape
    # native offsets: dx[j,i] approx linear in i,j; get origin via nearest to 0,0
    dx,dy=d['dx'],d['dy']; j0,i0=np.unravel_index(np.argmin(np.hypot(dx,dy)),dx.shape)
    gx=(np.arange(NG)-NG/2)*PIXG; X,Y=np.meshgrid(gx,gx)
    sx=np.sign(dx[j0,min(i0+1,nx-1)]-dx[j0,i0]); sy=np.sign(dy[min(j0+1,ny-1),i0]-dy[j0,i0])
    ci=i0+sx*X/pix; cj=j0+sy*Y/pix
    return map_coordinates(np.nan_to_num(arr),[cj,ci],order=1,mode='constant')
def ft(d,arr):
    g=regrid(d,arr); F=np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g)))*PIXG**2
    u=np.fft.fftshift(np.fft.fftfreq(NG,PIXG)); U,V=np.meshgrid(u,u)
    S=cov2(*d['beam']); q=S[0,0]*U**2+2*S[0,1]*U*V+S[1,1]*V**2; Bf=np.exp(-2*np.pi**2*q)
    return F,Bf,U,V
KL=206.265  # cycles/arcsec -> klambda
if __name__=='__main__':
    res={}; D={}
    for nm in ('Harada','Michiyama'):
        for ln,f0 in Q.REST.items():
            d=Q.load10(Q.P[nm],f0); D[(nm,ln)]=d
    # radial amplitude of |F|/Bf (beam-corrected) in klambda bins
    bins=np.array([0,2,4,6,8,10,12,15,20,25,30,40])
    lines=[]
    for ln in Q.REST:
        prof={}
        for nm in ('Harada','Michiyama'):
            F,Bf,U,V=ft(D[(nm,ln)],D[(nm,ln)]['m0']); r=np.hypot(U,V)*KL
            prof[nm]=[np.mean(np.abs(F[(r>=a)&(r<b)])/np.maximum(Bf[(r>=a)&(r<b)],1e-3)) for a,b in zip(bins[:-1],bins[1:])]
            # null-window noise level in the same bins
            nz=[]
            for n in D[(nm,ln)]['nulls']:
                Fn,_,_,_=ft(D[(nm,ln)],n); nz.append([np.mean(np.abs(Fn[(r>=a)&(r<b)])/np.maximum(Bf[(r>=a)&(r<b)],1e-3)) for a,b in zip(bins[:-1],bins[1:])])
            prof[nm+'_noise']=np.mean(nz,axis=0).tolist()
        h,m=np.array(prof['Harada']),np.array(prof['Michiyama'])
        lines.append(f'{ln}: |V| (Jy km/s, beam-corrected) per uv bin [kl] {bins.tolist()}\n  Harada   {np.round(h,2).tolist()} (noise {np.round(prof["Harada_noise"],2).tolist()})\n  Michiyama {np.round(m,2).tolist()} (noise {np.round(prof["Michiyama_noise"],2).tolist()})\n  H/M {np.round(h/m,2).tolist()}')
        res[ln]=prof
    print('\n'.join(lines)); open('/workspace/out/q40b/q40c_uv.log','w').write('\n'.join(lines))
    json.dump(res,open('/workspace/out/q40b/q40c_uv.json','w'))
