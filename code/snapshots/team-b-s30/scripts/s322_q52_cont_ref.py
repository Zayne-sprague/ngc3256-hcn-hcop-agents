"""Q52: N/S of the 1.1 mm continuum (dust-column proxy) at the same circular beams as the line ratios, and N/S of each line integral itself.
usage: python3 -I s322_q52_cont_ref.py -> out/q52/s322.txt"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from regions import REGIONS
from astropy.io import fits
R=json.load(open('/workspace/out/q52/s320.json'))
f='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'
L=[]
for th in (1.4,2.5):
    cut,nu,pix,bm,pxy=load(f,half=10)
    if cut.ndim==2: cut=cut[None]
    Tb=to_tb_conv(cut,np.array([257.752]),pix,bm,th,257.752)[0]
    v={}
    for rn in 'NS':
        x,y=pxy(REGIONS[rn][:2]); v[rn]=float(Tb[int(round(y)),int(round(x))])
    L.append(f'cont 1.1mm th={th}: N {v["N"]*1e3:.1f} mK  S {v["S"]*1e3:.1f} mK  N/S {v["N"]/v["S"]:.2f}')
    J='32'
    for ln in ('HCO+','HCN','HNC','CN'):
        a,b=R[f'{J}|{th}|N|{ln}']['350'],R[f'{J}|{th}|S|{ln}']['350']; L.append(f'  3-2 {ln} N/S = {a[0]/b[0]:.2f}')
    L.append(f'  3-2 HCO+/cont: N {R[f"{J}|{th}|N|HCO+"]["350"][0]/v["N"]:.0f} S {R[f"{J}|{th}|S|HCO+"]["350"][0]/v["S"]:.0f}  (K km/s per K)')
    for ln in ('HCN','HNC','CN'):
        L.append(f'  {ln}/cont: N {R[f"{J}|{th}|N|{ln}"]["350"][0]/v["N"]:.0f} S {R[f"{J}|{th}|S|{ln}"]["350"][0]/v["S"]:.0f}')
for th in (2.5,3.5):
    for ln in ('HCO+','HCN','HNC'):
        a,b=R[f'10|{th}|N|{ln}']['350'],R[f'10|{th}|S|{ln}']['350']; L.append(f'1-0 th={th} {ln} N/S = {a[0]/b[0]:.2f}')
    cn=lambda rn:R[f'10|{th}|{rn}|CN_G1']['350'][0]+R[f'10|{th}|{rn}|CN_G2']['350'][0]
    L.append(f'1-0 th={th} CN N/S = {cn("N")/cn("S"):.2f}')
open('/workspace/out/q52/s322.txt','w').write('\n'.join(L)); print('\n'.join(L))
