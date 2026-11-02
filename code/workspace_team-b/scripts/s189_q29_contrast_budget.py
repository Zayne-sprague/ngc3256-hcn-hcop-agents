"""Q29 s189: error budget for the N/S contrast of HCN/HCO+ (1-0 vs 3-2) and for r31 = T_b(3-2)/T_b(1-0) per line, at 1.5" (X145 1-0) and 2.2" (X961 1-0) circular beams vs Xb4d 3-2.
Budget terms: (a) statistical (from s183 aperture errors, incl. N,S independent); (b) window/aperture: +-150 vs +-350 km/s apertures (half-difference); (c) absolute flux scale 10%/band: cancels in the N/S contrast of one line's ratio
(same pair of lines, same pointing: scale multiplies N and S equally) and in the N/S ratio of any single line; NOT in r31 (different projects) -> r31 N/S ratios are also free of it but absolute r31 carry 10%x10%;
(d) 3-2 recovery differential Xb4d vs compact Xb4f at the Xb4f 6" beam (out/q19/s102.txt: R3=4f/4d HCN N 1.158 S 1.355 Bridge 1.174; HCO+ N 1.169 S 1.296 Bridge 1.199): taken as an UPPER bound of the 3-2 loss at N,S,Bridge at 2.2" (and 0 as lower bound);
the 1-0 data (X145, X961) are short-spacing complete (X961/X145 agree 0.95-0.99).  usage: python3 -I s189_q29_contrast_budget.py -> out/q29/s189.txt/json"""
import json, numpy as np
O='/workspace/out/q29/'; r=json.load(open(O+'s183.json')); log=[]
def P(s=''): print(s,flush=True); log.append(s)
R3={'HCN':{'N':1.158,'S':1.355,'Bridge':1.174},'HCO+':{'N':1.169,'S':1.296,'Bridge':1.199}}
out={}
for th,k10 in (('1.50','2015.1.00412 X145'),('2.20','2015.1.00993 X961')):
    P(f'\n######## beam {th}" ({k10} vs Xb4d)')
    v={}
    for rn in ('N','S','Bridge'):
        a=r[f'{k10}|{th.rstrip("0") if th=="1.50" else th}|{rn}'] if f'{k10}|{th.rstrip("0") if th=="1.50" else th}|{rn}' in r else None
    key=lambda rn: f'{k10}|{float(th)}|{rn}'
    for rn in ('N','S','Bridge'):
        a=r[key(rn)]; x=r['xb4d_32'][th][rn]
        for w,w10,w32 in (('350','w350_ap','win'),('150','h150_ap','h150')):
            v[(rn,w)]=dict(H10=a['HCN'][w10][0],O10=a['HCOp'][w10][0],eH10=a['HCN'][w10][1],eO10=a['HCOp'][w10][1],H32=x[w32]['HCN32'][0],O32=x[w32]['HCOp32'][0],eH32=x[w32]['HCN32'][1],eO32=x[w32]['HCOp32'][1])
    for w in ('350','150'):
        d={rn:v[(rn,w)] for rn in ('N','S','Bridge')}
        c10=lambda a,b: (d[a]['H10']/d[a]['O10'])/(d[b]['H10']/d[b]['O10'])
        c32=lambda a,b: (d[a]['H32']/d[a]['O32'])/(d[b]['H32']/d[b]['O32'])
        def err10(a,b): return np.sqrt(sum((d[q][e]/d[q][m])**2 for q in (a,b) for e,m in (('eH10','H10'),('eO10','O10'))))*c10(a,b)
        def err32(a,b): return np.sqrt(sum((d[q][e]/d[q][m])**2 for q in (a,b) for e,m in (('eH32','H32'),('eO32','O32'))))*c32(a,b)
        for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
            # recovery-corrected 3-2: multiply each line at each position by R3
            cc=c32(a,b)*(R3['HCN'][a]/R3['HCO+'][a])/(R3['HCN'][b]/R3['HCO+'][b])
            out[f'{th}|{w}|{a}/{b}']=dict(c10=c10(a,b),e10=err10(a,b),c32=c32(a,b),e32=err32(a,b),c32_rec=cc)
            P(f' +-{w:3s} {a}/{b}: 1-0 {c10(a,b):.3f}+-{err10(a,b):.3f} | 3-2 {c32(a,b):.3f}+-{err32(a,b):.3f}, recovery-bound-corrected 3-2 {cc:.3f} | D=3-2/1-0 {c32(a,b)/c10(a,b):.3f} (corrected {cc/c10(a,b):.3f})')
    # r31 per line at +-350 and +-150
    for w in ('350','150'):
        for ln,(k1,k3) in (('HCN',('H10','H32')),('HCO+',('O10','O32'))):
            row=f' r31 {ln} +-{w}: '
            for rn in ('N','S','Bridge'):
                d=v[(rn,w)]; rr=d[k3]/d[k1]; e=rr*np.hypot(d['e'+k1[0]+'10']/d[k1],d['e'+k3[0]+'32']/d[k3])
                rc=rr*R3[ln][rn]; out[f'{th}|r31|{ln}|{w}|{rn}']=dict(r31=rr,err=e,r31_rec_bound=rc)
                row+=f'{rn} {rr:.3f}+-{e:.3f} (rec-bound {rc:.3f}) | '
            nsr=out[f'{th}|r31|{ln}|{w}|N']['r31']/out[f'{th}|r31|{ln}|{w}|S']['r31']; nsc=out[f'{th}|r31|{ln}|{w}|N']['r31_rec_bound']/out[f'{th}|r31|{ln}|{w}|S']['r31_rec_bound']
            P(row+f'N/S {nsr:.3f} (rec-bound {nsc:.3f})'); out[f'{th}|r31NS|{ln}|{w}']=dict(obs=nsr,rec_bound=nsc)
P('\nBudget summary at 2.2" (X961 vs Xb4d), N/S contrast of HCN/HCO+: see lines above. Absolute-flux-scale: cancels (see docstring).')
open(O+'s189.txt','w').write('\n'.join(log)); json.dump(out,open(O+'s189.json','w'),indent=1,default=float)
