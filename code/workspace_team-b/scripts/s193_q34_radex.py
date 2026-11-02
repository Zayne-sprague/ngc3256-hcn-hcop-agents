"""Q34 s193: the s186 RADEX evidence machinery (Q29) re-run with the uv-matched 3-2 (s192, 'comb' = Xb4d 12m + Xb4f 7m, tapered to 1.5"/2.2") replacing the Xb4d image-plane 3-2.
Same 1-0 data (X145 1.5", X961 2.2"), same priors/band-scale marginalisation. Reports posteriors per position and ln BF for N-vs-S etc. (excitation vs abundance hypotheses), next to the old Xb4d-based numbers
(out/q29/s186.json).  usage: python3 -I s193_q34_radex.py -> out/q34/s193.txt/json"""
import sys, json, numpy as np
src=open('/workspace/scripts/s186_q29_radex_evidence.py').read(); src=src[:src.index("res={}")]
src=src.replace("O='/workspace/out/q29/'; log=[]","O='/workspace/out/q29/'; log=[]")
exec(src)
M=json.load(open('/workspace/out/q34/s192.json'))['res']; old=data
def data(th,reg):
    bn={'1.50':'1p5','2.20':'2p2'}[th]; m=M[f'{bn}|comb|{reg}']
    I0,e0=old(th,reg); I=I0.copy(); e=e0*I0
    I[1]=m['HCN32|ap|w350'][0]*K; I[3]=m['HCOp32|ap|w350'][0]*K; e[1]=m['HCN32|ap|w350'][1]*K; e[3]=m['HCOp32|ap|w350'][1]*K
    return I,e/I
res={}; USE={'3-2 only':np.array([False,True,False,True]),'1-0 + 3-2':np.array([True,True,True,True])}
O2='/workspace/out/q34/'; o186=json.load(open('/workspace/out/q29/s186.json'))
for th in ('1.50','2.20'):
    P(f'\n######## beam {th}" uv-matched 3-2 (comb). Data K km/s (x{K:.4f}):')
    for g in ('N','S','Bridge'):
        y,e=data(th,g); P(f'  {g}: HCN1 {y[0]:.1f} HCN3 {y[1]:.1f} HCO1 {y[2]:.1f} HCO3 {y[3]:.1f}')
    for g in ('N','S','Bridge'):
        p=post1(th,g,USE['1-0 + 3-2']); res[f'post|{th}|{g}']=p; q=o186[f'post|{th}|1-0 + 3-2|{g}']
        P(f'  [1-0+3-2] {g}: Tk {p["Tk"][1]:.0f} [{p["Tk"][0]:.0f},{p["Tk"][2]:.0f}] P(Tk>50)={p["PTk50"]:.2f}; log n {p["logn"][1]:.2f} [{p["logn"][0]:.2f},{p["logn"][2]:.2f}]; log A {p["logA"][1]:.2f} [{p["logA"][0]:.2f},{p["logA"][2]:.2f}]; log N(HCO+)/150km/s {p["logNO"][1]:.1f}; tau(HCN 1-0) {10**p["tauHCN10"][1]:.1f}')
        P(f'     (old Xb4d-based: Tk {q["Tk"][1]:.0f} [{q["Tk"][0]:.0f},{q["Tk"][2]:.0f}]; log n {q["logn"][1]:.2f}; log A {q["logA"][1]:.2f} [{q["logA"][0]:.2f},{q["logA"][2]:.2f}])')
    for pair in (('N','S'),('N','Bridge'),('S','Bridge')):
        for un in ('3-2 only','1-0 + 3-2'):
            E=evid(th,pair,USE[un]); res[f'ev|{th}|{pair[0]}-{pair[1]}|{un}']=E; Eo=o186[f'ev|{th}|{pair[0]}-{pair[1]}|{un}']
            P(f'  ln BF vs free [{pair[0]} vs {pair[1]}, {un}]: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free')+'   | old Xb4d: '+'  '.join(f'{k} {v:+.2f}' for k,v in Eo.items() if k not in ('free',)))
open(O2+'s193.txt','w').write('\n'.join(log)); json.dump(res,open(O2+'s193.json','w'),indent=1,default=float)
