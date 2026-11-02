# Q25: HCN/HCO+ 1-0 (archive Band 3: Harada 2015.1.00412.S, Michiyama 2015.1.00993.S) vs 2-1 (Xb55) at a common 10" beam, annuli 0-4/4-8/6-9/8-12" about the N-S midpoint.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import json
tb=(10.,10.,0.); ANN=[(0,4),(4,8),(6,9),(8,12)]
P={'Harada':'/archive/files/2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',
   'Michiyama':'/archive/files/2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits'}
REST={'HCN 1-0':88.631602,'HCO+ 1-0':89.188526}
def load10(path,f0):
    h=fits.open(path,memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C
    dv=abs(np.median(np.diff(v))); pix=abs(hd['CDELT2'])*3600; ny,nx=hd['NAXIS2'],hd['NAXIS1']
    def win(lo,hi):
        m=np.where((v>=lo)&(v<=hi))[0]
        return np.nansum(np.array(h[0].data[0,m[0]:m[-1]+1],float),axis=0)*dv,len(m)
    m0,n=win(VSYS-HW,VSYS+HW); nulls=[]
    vo=[(1-REST[k]*(1-VSYS/C)/f0)*C for k in REST]   # both lines' velocities in this frame: exclude +-900 km/s around each
    hw=HW
    for lo in np.arange(v.min()+60,v.max()-2*hw-60,hw):
        if all(lo+2*hw<x-900+VSYS*0 or lo>x+900 for x in vo): nulls.append(win(lo,lo+2*hw)[0])
    dx,dy=offsets(hd,ny,nx)
    return dict(m0=m0,nulls=nulls,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),pix=pix,dx=dx,dy=dy,hd=hd,vr=(v.min(),v.max()),nch=n)
def sm(d,a=None):
    k=kernel(d['hd'],tb,d['beam']); a=np.nan_to_num(d['m0'] if a is None else a)
    return fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
def ann(d,s):
    r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb)/d['pix']**2
    return np.array([np.nansum(s[(r>=a)&(r<b)])/omp for a,b in ANN]), r, omp
if __name__=='__main__':
    res={}
    sets={}
    for nm,p in P.items():
        for ln,f0 in REST.items():
            d=load10(p,f0); s=sm(d); fl,r,omp=ann(d,s)
            nz=[ann(d,sm(d,n))[0] for n in d['nulls']]
            e=np.std(nz,axis=0) if len(nz)>2 else np.full(len(ANN),np.nan)
            print(f'{nm:9s} {ln:8s} beam {d["beam"][0]:.2f}x{d["beam"][1]:.2f} pix {d["pix"]:.2f} v {d["vr"][0]:.0f}..{d["vr"][1]:.0f} nch {d["nch"]} nnull {len(nz)}  flux Jy km/s {fl.round(2)} +- {e.round(2)}')
            res[f'{nm}|{ln}']=dict(F=fl.tolist(),e=e.tolist()); sets[(nm,ln)]=(d,s)
    D21=[load(k,nullspec=False) for k in ('HCN 2-1','HCO+ 2-1')]
    S21=[smooth(d,10.) for d in D21]
    for k,d,s in zip(('HCN 2-1','HCO+ 2-1'),D21,S21):
        fl=ann(d,s)[0]; res[f'Xb55|{k}']=dict(F=fl.tolist()); print('Xb55',k,fl.round(2))
    print('annuli',ANN)
    for nm in P:
        h=np.array(res[f'{nm}|HCN 1-0']['F']); o=np.array(res[f'{nm}|HCO+ 1-0']['F']); eh=np.array(res[f'{nm}|HCN 1-0']['e']); eo=np.array(res[f'{nm}|HCO+ 1-0']['e'])
        R=h/o; eR=abs(R)*np.hypot(eh/h,eo/o); print(f'R10 {nm}: {R.round(3)} +- {eR.round(3)}; 8-12/4-8 = {R[3]/R[1]:.3f}')
        res[f'R10|{nm}']=dict(R=R.tolist(),eR=eR.tolist())
    h21=np.array(res['Xb55|HCN 2-1']['F']); o21=np.array(res['Xb55|HCO+ 2-1']['F']); R21=h21/o21; print('R21 Xb55',R21.round(3),'8-12/4-8',round(R21[3]/R21[1],3))
    json.dump(res,open('/workspace/out/q25/q25_ratio10.json','w'),indent=1)
