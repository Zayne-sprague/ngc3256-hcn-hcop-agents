"""Q53: 13CO(2-1) (2015.1.00412 X155 spw0, 0.85") as a column reference independent of the 1.1 mm continuum; N/S of X/13CO and HCO+/13CO,
total and core/wings. Beam 2.5" (1-0 lines) and 1.4" (3-2 lines vs 13CO convolved to 1.4"). usage: python3 -I s331_q53_13co_ref.py -> out/q53/s331.txt"""
import sys, numpy as np
sys.path.insert(0,'/workspace/scripts')
import s320_q52_chem_ratios as m, s330_q53_velbins as s   # s330 re-runs on import (cached, quick)
m.LINES['21']={'13CO':(m.A+'2015.1.00412.S/member.uid___A001_X33e_X155.ari_l.ngc_3256_sci.spw0_218561MHz.12m.cube.I.pbcor.fits',220.3987)}
res=s.res; out=[]
for th in (1.4,2.5):
    for rn in 'NS': res[f'21|{th}|{rn}|13CO']=s.spec('21','13CO',th,rn)
def I(J,th,rn,ln,b): return res[f'{J}|{th}|{rn}|{ln}'][b]
for J,th,lns in (('10',2.5,['HCN','HCO+','HNC','CN_G2']),('32',1.4,['HCN','HCO+','HNC']),('32',2.5,['HCN','HCO+','HNC'])):
  for b in ('all','core','blue','red'):
    for ln in lns:
        r={}
        for rn in 'NS':
            X=I(J,th,rn,ln,b); R=I('21',th,rn,'13CO',b); q=X[0]/R[0]; r[rn]=(q,q*np.hypot(X[1]/X[0],R[1]/R[0]),R[0])
        ns=r['N'][0]/r['S'][0]; e=ns*np.hypot(r['N'][1]/r['N'][0],r['S'][1]/r['S'][0])
        out.append(f'J{J} th={th} {b:5s} {ln:6s}/13CO(2-1) N={r["N"][0]:.3f} S={r["S"][0]:.3f} N/S={ns:.2f}+-{e:.2f}  13CO I N={r["N"][2]:.1f} S={r["S"][2]:.1f}')
for th in (1.4,2.5):
    a=I('21',th,'N','13CO','all'); c=I('21',th,'S','13CO','all'); out.append(f'13CO(2-1) th={th} N={a[0]:.1f}+-{a[1]:.1f} S={c[0]:.1f}+-{c[1]:.1f} N/S={a[0]/c[0]:.2f}')
open('/workspace/out/q53/s331.txt','w').write('\n'.join(out)); print('\n'.join(out))
