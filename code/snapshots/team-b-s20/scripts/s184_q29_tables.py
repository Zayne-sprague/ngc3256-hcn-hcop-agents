"""Q29: summary tables from s182.json (X94a/X94e Band 3, 0.77/0.63" native) and s183.json (cross-dataset 1-0, Xb4d 3-2) -> out/q29/s184_tables.txt, s184.json"""
import json, numpy as np
O='/workspace/out/q29/'; r=json.load(open(O+'s182.json')); r3=json.load(open(O+'s183.json')); log=[]
def P(s=''): print(s); log.append(s)
def rat(a,b):
    x=a[0]/b[0]; return x, abs(x)*np.hypot(a[1]/a[0] if a[0] else np.nan,b[1]/b[0])
out={}
P('## T1: HCN/HCO+ T_b ratio (K km/s ratio at equal beam), 1-0 vs 3-2 (Xb4d), aperture r=beam FWHM at HCO+ peak, +-350 km/s (lobe/ridge: Q17 windows)')
P('beam  region  | 1-0 X94a (0.77")   | 1-0 X145 (1.38x1.14")  | 1-0 X961 (2.1x1.7") | 3-2 Xb4d | contrast N/S(1-0,X145 or X94a) ')
for th in ('1.00','1.29','1.70','2.20','1.50'):
    for rn in ('N','S','Bridge','Slobe','Ridge'):
        row=[]
        if th in r:
            o=r[th]['b3'][rn]['ap']; row.append(('X94a',rat(o['HCN']['win'],o['HCOp']['win'])))
        for nm,key in (('X145','2015.1.00412 X145'),('X961','2015.1.00993 X961')):
            k=f'{key}|{float(th)}|{rn}'
            if k in r3: o=r3[k]; row.append((nm,rat(o['HCN']['w350_ap'],o['HCOp']['w350_ap'])))
        x4=None
        if th in r: o=r[th]['x4d'][rn]['win']; x4=rat(o['HCN32'],o['HCOp32'])
        elif th=='1.50': o=r3['xb4d_32']['1.50'][rn]['win']; x4=rat(o['HCN32'],o['HCOp32'])
        elif th=='2.20': o=r3['xb4d_32']['2.20'][rn]['win']; x4=rat(o['HCN32'],o['HCOp32'])
        P(f'{th}" {rn:7s}| '+' | '.join(f'{n} {v[0]:.2f}+-{v[1]:.2f}' for n,v in row)+(f' | 3-2 {x4[0]:.2f}+-{x4[1]:.2f}' if x4 else ''))
        out[f'T1|{th}|{rn}']={'1-0':{n:v for n,v in row},'3-2':x4}
P()
P('## T2: N/S contrast of HCN/HCO+ (aperture, +-350): 1-0 vs 3-2, same beam')
for th,key in (('1.00',None),('1.29',None),('1.50','2015.1.00412 X145'),('1.70','2015.1.00412 X145'),('2.20','2015.1.00412 X145')):
    c={}
    for tag in ('X94a','X145','X961'):
        try:
            if tag=='X94a':
                if th not in r: continue
                R=[rat(r[th]['b3'][g]['ap']['HCN']['win'],r[th]['b3'][g]['ap']['HCOp']['win']) for g in ('N','S')]
            else:
                kk='2015.1.00412 X145' if tag=='X145' else '2015.1.00993 X961'
                R=[rat(r3[f'{kk}|{float(th)}|{g}']['HCN']['w350_ap'],r3[f'{kk}|{float(th)}|{g}']['HCOp']['w350_ap']) for g in ('N','S')]
        except KeyError: continue
        c[tag]=(R[0][0]/R[1][0],R[0][0]/R[1][0]*np.hypot(R[0][1]/R[0][0],R[1][1]/R[1][0]))
    if th in r: X=r[th]['x4d']
    else: X=r3['xb4d_32'][th]
    R32=[X[g]['win']['HCN32'][0]/X[g]['win']['HCOp32'][0] for g in ('N','S')]
    P(f'{th}": 1-0 '+', '.join(f'{k} {v[0]:.2f}+-{v[1]:.2f}' for k,v in c.items())+f' | 3-2 {R32[0]/R32[1]:.2f}'); out[f'T2|{th}']={'1-0':c,'3-2':R32[0]/R32[1]}
P()
P('## T3: r31 = 3-2/1-0 T_b ratio (aperture +-350, matched beam; X145 1-0 and Xb4d 3-2 at 1.5"; X94a at 1.0")')
for th,src in (('1.50','X145'),('1.00','X94a')):
    for rn in ('N','S','Bridge','Slobe','Ridge'):
        if src=='X145':
            k=f'2015.1.00412 X145|1.5|{rn}'; a=r3[k]['HCN']['w350_ap']; b=r3[k]['HCOp']['w350_ap']; x=r3['xb4d_32']['1.50'][rn]['win']
        else:
            o=r['1.00']['b3'][rn]['ap']; a=o['HCN']['win'];b=o['HCOp']['win']; x=r['1.00']['x4d'][rn]['win']
        rh=rat(x['HCN32'],a); ro=rat(x['HCOp32'],b)
        P(f'{th}" ({src}) {rn:7s}: r31 HCN {rh[0]:.3f}+-{rh[1]:.3f}  HCO+ {ro[0]:.3f}+-{ro[1]:.3f}  (1-0 HCN {a[0]:.1f}+-{a[1]:.1f}, HCO+ {b[0]:.1f}+-{b[1]:.1f}; 3-2 HCN {x["HCN32"][0]:.1f} HCO+ {x["HCOp32"][0]:.1f})')
        out[f'T3|{th}|{rn}']=dict(r31_HCN=rh,r31_HCO=ro)
P()
P('## T4: X94a/X145 recovery ratio of the 1-0 aperture flux (+-350) at matched beam (flux/short-spacing check)')
for th in ('1.7','2.2'):
    for rn in ('N','S','Bridge'):
        a=r3[f'2016.1.00965 X94a|{th}|{rn}']; b=r3[f'2015.1.00412 X145|{th}|{rn}']; c=r3.get(f'2015.1.00993 X961|{th}|{rn}')
        s=f'{th}" {rn:7s} X94a/X145: HCN {a["HCN"]["w350_ap"][0]/b["HCN"]["w350_ap"][0]:.2f} HCO+ {a["HCOp"]["w350_ap"][0]/b["HCOp"]["w350_ap"][0]:.2f}'
        if c: s+=f' ; X961/X145: HCN {c["HCN"]["w350_ap"][0]/b["HCN"]["w350_ap"][0]:.2f} HCO+ {c["HCOp"]["w350_ap"][0]/b["HCOp"]["w350_ap"][0]:.2f}'
        P(s)
P()
P('## T5: weak lines, X94a/X94e Band 3 (template fit = HCO+ 1-0 profile; amplitude = integrated K km/s; peak-pixel and aperture), with 3-sigma limits')
def lim(v,e): return f'{v:.1f}+-{e:.1f}'+('' if v>3*e else f' (<{max(v,0)+3*e:.1f}, 3sig)')
for th in ('1.00','1.70'):
  for md in ('pix','ap'):
    for rn in ('N','S','Bridge','Slobe','Ridge'):
        o=r[th]['b3'][rn][md]; f=o['fit']
        hco=o['HCOp']['win'][0]; hcn=o['HCN']['win'][0]; hnc=o['HNC']['win'][0]; hc3=o['HC3N']['win'][0]
        P(f'{th}" {md:3s} {rn:7s}: HCO+ {hco:.0f} HCN {hcn:.0f} | H13CO+ {lim(*f["H13COp"])} SiO {lim(*f["SiO"])} H13CN {lim(*f["H13CN"])} HN13C {lim(*f["HN13C"])} CCH(3/2) {lim(*f["CCHa"])} | HNC {hnc:.0f} HC3N {hc3:.1f}+-{o["HC3N"]["win"][1]:.1f}')
        out[f'T5|{th}|{md}|{rn}']=dict(fit=f,HCO=hco,HCN=hcn,HNC=hnc,HC3N=hc3)
P()
P('## T6: ratios from X94a/X94e (+-350 aperture) HNC/HCN; Harada Table 13 (1.7", +-150) for comparison: N 0.559 S 0.495 ; HC3N/HCN N 0.058 S 0.089; SiO/HCO+ N 0.031 S 0.031; H13CO+/HCO+ N 0.046')
for th in ('1.00','1.70'):
    for rn in ('N','S','Bridge'):
        o=r[th]['b3'][rn]['pix']; a=r[th]['b3'][rn]['ap']
        P(f'{th}" {rn}: HNC/HCN pix {o["HNC"]["h150"][0]/o["HCN"]["h150"][0]:.2f} ap350 {a["HNC"]["win"][0]/a["HCN"]["win"][0]:.2f} | HC3N/HCN pix150 {o["HC3N"]["h150"][0]/o["HCN"]["h150"][0]:.3f}')
open(O+'s184_tables.txt','w').write('\n'.join(log)); json.dump(out,open(O+'s184.json','w'),indent=1,default=float)
