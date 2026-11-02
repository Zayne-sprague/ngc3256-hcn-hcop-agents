"""Q29 s187b: as s187 (with H13CO+ term, X 40-100) for N vs S and N vs Bridge at 2.2", but with the Xb4d 3-2 fluxes multiplied by the Xb4d->Xb4f recovery bound R3 (HCN N 1.158 S 1.355 Bridge 1.174; HCO+ N 1.169 S 1.296 Bridge 1.199; out/q19/s102.txt)
(scale 1.0 = as measured, 0.5 = half the bound, 1.0 of the bound = full).  Tests whether the common-excitation preference survives the 3-2 recovery differential.  usage: python3 -I s187b_q29_radex_recovery.py -> out/q29/s187b.txt"""
import sys, json, numpy as np
s187=open('/workspace/scripts/s187_q29_radex_h13.py').read(); s187=s187[:s187.index("res={}; ALL")]
exec(s187)
R3={'HCN':{'N':1.158,'S':1.355,'Bridge':1.174},'HCO+':{'N':1.169,'S':1.296,'Bridge':1.199}}
data0=data
def mk(f):
    def d(th,reg):
        I,e=data0(th,reg); I=I.copy(); I[1]*=1+f*(R3['HCN'][reg]-1); I[3]*=1+f*(R3['HCO+'][reg]-1); return I,e
    return d
ALL=np.array([True]*4); th='2.20'; out={}
for f in (0.0,0.5,1.0):
    data=mk(f); P(f'\n#### recovery bound fraction f={f} (3-2 x [1+f(R3-1)])')
    ex={g:lnLR(g) for g in ('N','S','Bridge')}
    for g in ('N','S','Bridge'):
        p=post2(th,g,ALL,ex[g]); out[f'post|{f}|{g}']=p
        P(f'  {g}: Tk {p["Tk"][1]:.0f} [{p["Tk"][0]:.0f},{p["Tk"][2]:.0f}]; log n {p["logn"][1]:.2f}; log A {p["logA"][1]:.2f} [{p["logA"][0]:.2f},{p["logA"][2]:.2f}]; tau(HCO+ 1-0) {10**p["tauHCO10"][1]:.1f}; log N(HCO+) {p["logNO"][1]:.2f}; log N(HCN) {p["logNH"][1]:.2f}')
    for pair in (('N','S'),('N','Bridge')):
        for sl in (0.05,0.10):
            E=evid2(th,pair,ALL,(ex[pair[0]],ex[pair[1]]),sl=sl); out[f'ev|{f}|{pair}|{sl}']=E
            P(f'  ln BF vs free [{pair[0]} vs {pair[1]}] sl={sl}: '+'  '.join(f'{k} {v:+.2f}' for k,v in E.items() if k!='free'))
open(O+'s187b.txt','w').write('\n'.join(log)); json.dump(out,open(O+'s187b.json','w'),indent=1,default=float)
