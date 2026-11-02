"""Q53: velocity-resolved N/S of X/HCO+ at matched circular beam. Reuses s320 loading/baseline recipe.
Bins (km/s about v_sys): core |v|<100, blue -350..-100, red 100..350, wings = blue+red. usage: python3 -I s330_q53_velbins.py -> out/q53/s330.txt/json"""
import sys, json, numpy as np, os
sys.path.insert(0,'/workspace/scripts')
import s320_q52_chem_ratios as m
from s180_q29_load import to_tb_conv
os.makedirs('/workspace/out/q53',exist_ok=True)
BINS={'core':(-100,100),'blue':(-350,-100),'red':(100,350),'all':(-350,350),'core50':(-50,50)}
def spec(J,ln,th,rn):
    path,f0=m.LINES[J][ln]; pos,v0=m.POS[rn]
    win=(-350,350); masks=[]
    if J=='10' and ln.startswith('CN'):
        other='CN_G1' if ln=='CN_G2' else 'CN_G2'; masks=[(m.LINES[J][other][1],)+m.WIN[other]]
        win=m.WIN[ln]
    if J=='10' and ln in('HCN','HCO+'): masks=[(m.LINES[J]['HCO+' if ln=='HCN' else 'HCN'][1],-350,350)]
    cut,nu,pix,bm,pxy=m.getcube(path,th)
    v=m.vel(nu,f0,v0); lo,hi=win; sel=(v>lo-1400)&(v<hi+1400)
    Tb=to_tb_conv(cut[sel],nu[sel],pix,bm,th,float(np.mean(nu))); vs=v[sel]; nus=nu[sel]
    x,y=pxy(pos); sp=Tb[:,int(round(y)),int(round(x))]; dv=abs(vs[1]-vs[0])
    free=(vs<lo-250)|(vs>hi+250)
    for (g0,glo,ghi) in masks:
        vg=m.vel(nus,g0,v0); free&=~((vg>glo-250)&(vg<ghi+250))
    free&=(vs>lo-1200)&(vs<hi+1200)
    p=np.polyfit(vs[free],sp[free],1); sp=sp-np.polyval(p,vs); sig=np.std(sp[free])
    out={}
    for b,(a,c) in BINS.items():
        w=(vs>=a)&(vs<=c); out[b]=[float(sp[w].sum()*dv),float(sig*np.sqrt(w.sum())*dv)]
    return out
res={}; txt=[]
cfg=[('10',['HCN','HCO+','HNC','CN_G2','CN_G1'],[2.5,3.5]),('32',['HCN','HCO+','HNC'],[1.4,2.5])]
for J,lns,ths in cfg:
  for th in ths:
    for rn in 'NS':
      for ln in lns: res[f'{J}|{th}|{rn}|{ln}']=spec(J,ln,th,rn)
def comb(J,th,rn,ln,b):
    if b=='wings': a=res[f'{J}|{th}|{rn}|{ln}']; return [a['blue'][0]+a['red'][0],np.hypot(a['blue'][1],a['red'][1])]
    return res[f'{J}|{th}|{rn}|{ln}'][b]
for J,lns,ths in cfg:
  for th in ths:
    for b in ['core','core50','blue','red','wings','all']:
      for ln in [l for l in lns if l!='HCO+']:
        r={}
        for rn in 'NS':
            X=comb(J,th,rn,ln,b); H=comb(J,th,rn,'HCO+',b)
            q=X[0]/H[0]; e=q*np.hypot(X[1]/X[0],H[1]/H[0]); r[rn]=(q,e,H[0])
        ns=r['N'][0]/r['S'][0]; ens=ns*np.hypot(r['N'][1]/r['N'][0],r['S'][1]/r['S'][0])
        txt.append(f'J{J} th={th} {b:6s} {ln:6s}/HCO+ N={r["N"][0]:.3f}+-{r["N"][1]:.3f} S={r["S"][0]:.3f}+-{r["S"][1]:.3f} N/S={ns:.2f}+-{ens:.2f}  HCO+ I N={r["N"][2]:.1f} S={r["S"][2]:.1f}')
json.dump(res,open('/workspace/out/q53/s330.json','w'),indent=1)
open('/workspace/out/q53/s330.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
