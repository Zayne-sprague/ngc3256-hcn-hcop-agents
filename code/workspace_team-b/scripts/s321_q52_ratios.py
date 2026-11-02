"""Q52: ratios/double ratios from out/q52/s320.json. usage: python3 -I s321_q52_ratios.py -> out/q52/s321.txt, s321.json"""
import json, numpy as np
R=json.load(open('/workspace/out/q52/s320.json')); out={}; L=[]
def I(J,th,rn,ln,w='350'):
    r=R.get(f'{J}|{th}|{rn}|{ln}'); return None if not r else r[w]
def cnI(J,th,rn,w='350'):
    if J=='32': return I(J,th,rn,'CN',w)
    a,b=I(J,th,rn,'CN_G1',w),I(J,th,rn,'CN_G2',w); return [a[0]+b[0],float(np.hypot(a[1],b[1]))]
def get(J,th,rn,ln,w='350'): return cnI(J,th,rn,w) if ln=='CN' else I(J,th,rn,ln,w)
def ratio(a,b):
    r=a[0]/b[0]; return r,r*np.hypot(a[1]/a[0],b[1]/b[0])
for J,ths in (('10',[2.5,3.5]),('32',[1.4,2.5])):
  for th in ths:
    for num,den in [('HCN','HCO+'),('HNC','HCO+'),('CN','HCO+'),('HNC','HCN'),('CN','HCN')]+([('HNC_412','HCN_412'),('HCN_412','HCO+_412'),('HNC_412','HCO+_412')] if J=='10' else []):
        row={}
        for rn in 'NS':
            for w in ('350','250'):
                a,b=get(J,th,rn,num,w),get(J,th,rn,den,w)
                if a and b: row[rn+w]=ratio(a,b)
        if 'N350' not in row: continue
        sysw=lambda rn: abs(row[rn+'350'][0]-row[rn+'250'][0])
        rN,eN=row['N350']; rS,eS=row['S350']
        ns=rN/rS; ens=ns*np.hypot(eN/rN,eS/rS); ens=float(np.hypot(ens,ns*np.hypot(sysw('N')/rN,sysw('S')/rS)))
        L.append(f'J={J} {th}" {num}/{den}: N {rN:.3f}+-{eN:.3f} (w250 {row["N250"][0]:.3f}) S {rS:.3f}+-{eS:.3f} (w250 {row["S250"][0]:.3f}) N/S {ns:.2f}+-{ens:.2f}')
        out[f'{J}|{th}|{num}/{den}']=dict(N=rN,eN=eN,S=rS,eS=eS,NS=ns,eNS=ens)
open('/workspace/out/q52/s321.txt','w').write('\n'.join(L)); json.dump(out,open('/workspace/out/q52/s321.json','w'),indent=1); print('\n'.join(L))
