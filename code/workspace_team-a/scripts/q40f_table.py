# Q40f: final tables. Central = matched-uv (u>8,10 kl) common-beam decomposition (q40e_sys.pkl); sys = RMS over variants (+umin 6/8/10); noise from filtered null windows.
# Also: plain 10" position values (q40.json) and native-beam decontaminated (q40_decontam.json) for comparison. Outputs out/q40b/q40f_table.md/.json
import sys; sys.path.insert(0,'/workspace/scripts')
import pickle, json, warnings, numpy as np
warnings.filterwarnings('ignore')
O=pickle.load(open('/workspace/out/q40b/q40e_sys.pkl','rb'))
NU={1:(88.631602,89.188526),2:(177.261111,178.375056),3:(265.886434,267.557626)}; KF={J:(b/a)**2 for J,(a,b) in NU.items()}
SETJ={'Harada':1,'Michiyama':1,'Xb55':2,'Xb4f':3}
VN=sorted({k[1] for k in O}); L=[]
def P(s=''): print(s); L.append(s)
def amp(o,reg):
    n=o['names']
    if reg=='N+S': return o['amp'][n.index('N')]+o['amp'][n.index('S')]
    return o['amp'][n.index(reg)] if reg in n else np.nan
def err(o,reg):
    n=o['names']
    if reg=='N+S': return np.hypot(o['err'][n.index('N')],o['err'][n.index('S')])
    return o['err'][n.index(reg)]
def Rv(umin,vn,nm,reg):
    h,q=O[(umin,vn,nm,'HCN')],O[(umin,vn,nm,'HCO+')]; a,b=amp(h,reg),amp(q,reg)
    return a/b*KF[SETJ[nm]] if (a>0 and b>0) else np.nan
def Rn(umin,nm,reg):
    h,q=O[(umin,'base',nm,'HCN')],O[(umin,'base',nm,'HCO+')]; a,b=amp(h,reg),amp(q,reg); ea,eb=err(h,reg),err(q,reg)
    if not (a>0 and b>0): return np.nan
    r=a/b*KF[SETJ[nm]]; return r*np.hypot(ea/a,eb/b)
UM=(8,10)
def quant(f):
    """f(umin,vn)->value. returns central (mean of base at umin 8,10), noise-free sys (rms over all umin 6/8/10 x variants about central)"""
    cen=np.nanmean([f(u,'base') for u in UM]); s=[f(u,v) for u in (6,8,10) for v in VN]; s=np.array(s)
    return cen,np.sqrt(np.nanmean((s[np.isfinite(s)]-cen)**2)),int(np.isfinite(s).sum())
def noise_of(fn): return np.nanmean([fn(u) for u in UM])
res={}
P('## Table A. HCN/HCO+ (T_B scale) per region, matched uv u>8-10 kl, common 8.64x6.0 beam; value (noise | sys)')
P('| region | 1-0 Michiyama | 1-0 Harada | 2-1 Xb55 | 3-2 Xb4f |'); P('|---|---|---|---|---|')
for reg in ('N','S','N+S','E1','E2','W'):
    row=[]
    for nm in ('Michiyama','Harada','Xb55','Xb4f'):
        c,s,n=quant(lambda u,v: Rv(u,v,nm,reg)); e=noise_of(lambda u: Rn(u,nm,reg)); res[f'R|{reg}|{nm}']=(c,e,s)
        row.append(f'{c:.2f} ({e:.2f}|{s:.2f})' if np.isfinite(c) else 'n/m')
    P(f'| {reg} | '+' | '.join(row)+' |')
P('\n## Table B. clump/reference contrast Rclump/Rref per J (noise | sys); same-J so no calibration term. * = deviates from 1 by >3 sigma (noise+sys)')
P('| clump | ref | 1-0 Michiyama | 1-0 Harada | 2-1 | 3-2 |'); P('|---|---|---|---|---|---|')
def contrast(nm,c,ref):
    f=lambda u,v: Rv(u,v,nm,c)/Rv(u,v,nm,ref); cen,s,n=quant(f)
    # noise: propagate independent noise of the two R (shared nuclei correlated ignored)
    rc=np.nanmean([Rv(u,'base',nm,c) for u in UM]); rr=np.nanmean([Rv(u,'base',nm,ref) for u in UM])
    e=cen*np.hypot(noise_of(lambda u: Rn(u,nm,c))/rc,noise_of(lambda u: Rn(u,nm,ref))/rr)
    return cen,e,s
for c in ('E1','W'):
    for ref in ('N','S','N+S'):
        row=[]
        for nm in ('Michiyama','Harada','Xb55','Xb4f'):
            cen,e,s=contrast(nm,c,ref); res[f'C|{c}|{ref}|{nm}']=(cen,e,s); sig=abs(cen-1)/np.hypot(e,s)
            row.append(f'{cen:.2f} ({e:.2f}|{s:.2f}){"*" if sig>3 else ""}')
        P(f'| {c} | {ref} | '+' | '.join(row)+' |')
P('\n## Table C. ladder/double ratios (cal 10% on the absolute ones only; clump/ref double ratios are calibration-free)')
def dbl(c,num,den):
    f=lambda u,v: Rv(u,v,num,c)/Rv(u,v,den,c); cen,s,n=quant(f)
    e=cen*np.hypot(noise_of(lambda u: Rn(u,num,c))/np.nanmean([Rv(u,'base',num,c) for u in UM]),noise_of(lambda u: Rn(u,den,c))/np.nanmean([Rv(u,'base',den,c) for u in UM]))
    return cen,e,s
P('| region | R21/R10 (M) | R32/R21 | R32/R10 (M) |'); P('|---|---|---|---|')
for reg in ('N','S','N+S','E1','W'):
    row=[]
    for num,den in (('Xb55','Michiyama'),('Xb4f','Xb55'),('Xb4f','Michiyama')):
        cen,e,s=dbl(reg,num,den); res[f'D|{reg}|{num}/{den}']=(cen,e,s); row.append(f'{cen:.2f} ({e:.2f}|{s:.2f}; cal 0.10 -> tot {np.sqrt(e*e+s*s+(0.1*cen)**2):.2f})')
    P(f'| {reg} | '+' | '.join(row)+' |')
P('\nclump/reference of the double ratios (cal-free):')
def ddc(c,ref,num,den):
    f=lambda u,v: (Rv(u,v,num,c)/Rv(u,v,den,c))/(Rv(u,v,num,ref)/Rv(u,v,den,ref)); cen,s,n=quant(f)
    nz=np.hypot.reduce([noise_of(lambda u: Rn(u,nm,r))/np.nanmean([Rv(u,'base',nm,r) for u in UM]) for nm in (num,den) for r in (c,ref)])
    return cen,cen*nz,s
for c in ('E1','W'):
    for ref in ('N','S','N+S'):
        row=[]
        for num,den in (('Xb55','Michiyama'),('Xb4f','Xb55'),('Xb4f','Michiyama')):
            cen,e,s=ddc(c,ref,num,den); res[f'DC|{c}|{ref}|{num}/{den}']=(cen,e,s); row.append(f'{num}/{den}: {cen:.2f} ({e:.2f}|{s:.2f}; {(cen-1)/np.hypot(e,s):+.1f}sig)')
        P(f'  {c}/{ref}: '+'; '.join(row))
# comparison tables: plain 10" and native decontaminated
j=json.load(open('/workspace/out/q40/q40.json'))['R']; dc=json.load(open('/workspace/out/q40/q40_decontam.json'))['res']
P('\n## Table D. other methods, E1,W contrast vs N, S, N+S per J (Michiyama | Harada | 2-1 | 3-2 Xb4f)')
for lab,get in (('plain 10" position values',lambda nm,r: j[f'{nm}|{r}'][0] if r in('N','S','E1','W') else None),('native-beam decontaminated',lambda nm,r: dc[f'{r}|{nm}'][0])):
    P(f'{lab}:')
    for c in ('E1','W'):
        for ref in ('N','S','N+S'):
            vals=[]
            for nm in ('Michiyama','Harada','Xb55','Xb4f'):
                if lab.startswith('plain') and ref=='N+S':
                    # N+S not tabulated in plain: use mean of N and S R (plain flux weighted would need fluxes) -> skip
                    vals.append('-'); continue
                vals.append(f'{get(nm,c)/get(nm,ref):.2f}')
            P(f'  {c}/{ref}: '+' | '.join(vals))
json.dump({k:list(map(float,v)) for k,v in res.items()},open('/workspace/out/q40b/q40f_table.json','w'),indent=1)
open('/workspace/out/q40b/q40f_table.md','w').write('\n'.join(L))
