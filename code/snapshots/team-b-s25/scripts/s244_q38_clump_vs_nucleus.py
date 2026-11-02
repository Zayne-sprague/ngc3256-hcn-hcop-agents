"""Q38 fix 2: same-method (s58d NNLS, native-beam moment-0, V0 + variants) HCN/HCO+ ladder for the nuclear complex (nuc_total, N, N+S) and the clumps E1/W, R32/R21 and R43/R21 per region
and the clump/nucleus double ratio with errors: stat (sign-flip), model-variant spread, and 5% per-line calibration on each R (7% on each ladder ratio, independent bands -> 10% on the ratio of ratios).
Also the 10" position values (s240 A) for N.  usage: python3 -I s244_q38_clump_vs_nucleus.py -> out/q38/s244.txt/json"""
import json,numpy as np
O='/workspace/out/'; cf=json.load(open(O+'q12b/clumpfit.json')); S=json.load(open(O+'q38/s240.json'))
KF=1.0126; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def flux(line,ous,reg):
    d=cf[f'{line}|{ous}']; v=d['V0']; 
    f=(lambda x:x[reg]) if reg not in('N+S','nuc_total') else (lambda x:x['N']+x['S'] if reg=='N+S' else x['nuc_total'])
    val=f(v); var=[f(d[k]) for k in d]
    err=v['err'].get(reg,None) if reg in v['err'] else (np.hypot(v['err']['N'],v['err']['S']) if reg=='N+S' else np.sqrt(sum(v['err'][n]**2 for n in ('N','S','Bridge','haloB')+(('haloA',) if 'haloA' in v['err'] else ()))))
    return val,float(np.hypot(err,np.std(var)))
OUS={2:'Xb55',3:'Xb4f',4:'Xb47'}; res={}
for reg in ('nuc_total','N+S','N','E1','W'):
    R={}
    for J,o in OUS.items():
        h=flux('HCN',o,reg); c=flux('HCO+',o,reg); r=h[0]/c[0]*KF; er=r*np.hypot(h[1]/max(h[0],1e-9),c[1]/c[0]); R[J]=(r,er)
    cal=0.05
    r32=R[3][0]/R[2][0]; e32=r32*np.hypot(R[3][1]/R[3][0],R[2][1]/R[2][0]); e32c=r32*np.hypot(e32/r32,np.sqrt(2)*cal*np.sqrt(2))
    r43=R[4][0]/R[2][0]; e43=r43*np.hypot(R[4][1]/max(R[4][0],1e-9),R[2][1]/R[2][0])
    res[reg]=dict(R=R,r32_21=(r32,e32c),r43_21=(r43,e43))
    P(f' {reg:9s} R21 {R[2][0]:.3f}+-{R[2][1]:.3f}  R32 {R[3][0]:.3f}+-{R[3][1]:.3f}  R43 {R[4][0]:.3f}+-{R[4][1]:.3f} | R32/R21 {r32:.3f}+-{e32c:.3f} (stat+variants {e32:.3f}) | R43/R21 {r43:.3f}+-{e43:.3f}')
P('\n clump/nucleus double ratio of R32/R21 (nucleus = nuc_total and N):')
for ref in ('nuc_total','N'):
    for c in ('E1','W'):
        a=res[c]['r32_21']; b=res[ref]['r32_21']; q=a[0]/b[0]; e=q*np.hypot(a[1]/a[0],b[1]/b[0]); z=(q-1)/e
        P(f'  {c}/{ref}: {q:.2f}+-{e:.2f}  ({z:+.1f} sigma from 1)'); res[f'{c}/{ref}']=(q,e)
P('\n R43/R21: nucleus (nuc_total) vs clump upper limits: ',)
a=res['nuc_total']['r43_21']; P(f'  nuc_total R43/R21 = {a[0]:.2f}+-{a[1]:.2f}; E1 clump R43<0.058 (2 sigma) -> R43/R21 < {0.058/res["E1"]["R"][2][0]:.2f}; W < {0.328/res["W"]["R"][2][0]:.2f}')
P('\n 10" position values (blended, s240 A): R32/R21 at each position (Xb4f/Xb55)')
for p in ('N','E1','E2','W'):
    r2=S['A'][f'Xb55|{p}']['R_K']['0']; r3=S['A'][f'Xb4f|{p}']['R_K']['0']; r4=S['A'][f'Xb47|{p}']['R_K']['0']; P(f'  {p}: R21 {r2:.3f} R32 {r3:.3f} R43 {r4:.3f}  R32/R21 {r3/r2:.2f} R43/R21 {r4/r2:.2f}')
json.dump(res,open(O+'q38/s244.json','w'),indent=1,default=float); open(O+'q38/s244.txt','w').write('\n'.join(log))
