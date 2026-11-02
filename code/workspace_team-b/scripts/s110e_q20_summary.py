# Q20 step e: collate s110b-d: matched-beam 4-3/3-2 Tb ratios of the ridge and lobe wing excess, HCN/HCO+ pattern, profile shape per line, and upper limits for CS 7-6, SO, H13CN 4-3, HC3N.
# Errors: Monte-Carlo propagation of the independent fit errors (s110c); 3-2 and 4-3 fit errors share the template/mask so ratio errors are, if anything, overstated.   usage: python3 -I s110e_q20_summary.py
import json, numpy as np
OUT='/workspace/out/q20'; r=json.load(open(OUT+'/s110c.json')); dd=json.load(open(OUT+'/s110d.json')); rng=np.random.default_rng(3); txt=[]; P=lambda s:(txt.append(s),print(s))
def draw(k,N=20000): e=r[k]; return rng.normal(e['s'],e['se'],N)*e['pk']
def q(x): return f'{np.median(x):.2f} (+{np.percentile(x,84)-np.median(x):.2f}/-{np.median(x)-np.percentile(x,16):.2f})'
NU=(354.505/265.886)**2
P('--- matched-beam (4.5"x3.0") wing-excess flux, 3-2 vs 4-3 (E1, Jy/beam km/s at the 3-2 template peak) and Tb ratio r43/32 = (F43/F32)/%.3f'%NU)
R={}
for mol in ('HCN','HCO+'):
    for pos,wn in (('ridge','red'),('lobe','blue')):
        a=draw(f'E1|{mol}32|{wn}'); b=draw(f'E1|{mol}43|{wn}'); ra=b/a/NU; ok=a>0
        R[(mol,pos)]=np.where(ok,ra,np.nan)
        P(f'{mol:4s} {pos:5s}: F32 {np.mean(a):+.3f}+-{np.std(a):.3f}  F43 {np.mean(b):+.3f}+-{np.std(b):.3f}  significance 4-3 {np.mean(b)/np.std(b):.1f} sigma; r43/32 = {q(ra[ok])}')
for mol in ('HCN','HCO+'):
    x=R[(mol,'ridge')]/R[(mol,'lobe')]; x=x[np.isfinite(x)]; P(f'{mol}: (r43/32 ridge)/(r43/32 lobe) = {q(x)} ; P(ridge<lobe)={np.mean(x<1):.2f}')
# stacked 4-3 amplitude (inverse-variance mean of s over HCN43 and HCO+43; units of own 3-2 excess)
for pos,wn in (('ridge','red'),('lobe','blue')):
    for J in ('32','43'):
        ss=np.array([r[f'E1|{m}{J}|{wn}']['s'] for m in ('HCN','HCO+')]); ee=np.array([r[f'E1|{m}{J}|{wn}']['se'] for m in ('HCN','HCO+')]); w=1/ee**2
        P(f'stack HCN+HCO+ {J} {pos}: s = {np.sum(w*ss)/w.sum():.2f} +- {1/np.sqrt(w.sum()):.2f}  ({np.sum(w*ss)/w.sum()*np.sqrt(w.sum()):.1f} sigma)')
P('--- HCN/HCO+ flux ratio of the wing excess (matched beam)')
core=np.load(OUT+'/s110b_maps.npz'); 
for J in ('32','43'):
    for pos,wn in (('ridge','red'),('lobe','blue')):
        a=draw(f'E1|HCN{J}|{wn}'); b=draw(f'E1|HCO+{J}|{wn}'); x=(a/b)[(b>0)]; P(f'  J={J} {pos}: HCN/HCO+ excess {q(x)}')
P('--- upper limits (3 sigma = max(best,0)+3 se) on the red excess at the ridge and blue at the lobe, Jy/beam km/s at 4.5"x3.0", and as ratio to HCN 3-2 excess at the same place (Jy); core ratio = line core/HCN3-2 core at the ridge position')
import sys; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
c0,v0,ps7,wc7,xx,yy,b07=load_cube('Xb47','spw16',354.50547e9,-100,100,half=24.0); px=pxy_fn(wc7,xx,yy)
ix,iy=[int(round(t)) for t in px(*RED)]; jx,jy=[int(round(t)) for t in px(*LOBE)]
for line in ('CS76','SO65','H13CN43'):
    for pos,wn,k,(cx,cy) in (('ridge','red',0,(ix,iy)),('lobe','blue',1,(jx,jy))):
        if line=='H13CN43' and wn=='blue': continue
        e=r[f'E1|{line}|{wn}']; ul=max(e['s'],0)+3*e['se']; cr=core[f'{line}|core'][cy,cx]/core['HCN32|core'][cy,cx]
        P(f'  {line:8s} {pos:5s}: s={e["s"]:+.2f}+-{e["se"]:.2f}; 3s UL on s (=F/F_HCN32 excess) < {ul:.2f}, i.e. < {ul*e["pk"]:.2f} Jy/beam km/s; core ratio at same pixel {cr:.3f} -> UL is {ul/cr:.0f}x the core ratio')
P('  (HC3N 28-27 and SO 6_5-5_4 at 1.05"; s110d): SO ridge 3s UL 0.082 Jy km/s vs HCN 3-2 ridge 0.257 -> <0.32; lobe <0.089/0.325 = <0.27; core ratio SO/HCN ~ %.3f. HC3N ridge <0.131 -> <0.51; lobe <0.117 -> <0.36'%(dd['SO65|core']['S']/dd['HCN32|core']['S']))
P('--- profile shape: 40 km/s bins (|v|), flux of the excess at the template peak, ridge (red) and lobe (blue); shape statistic = centroid |v| of bins 160-360 and fraction in 160-200')
for mol in ('HCN','HCO+'):
    for pos in ('ridge','lobe'):
        for J,k in (('3-2','f32'),('4-3','f43')):
            vs=np.arange(160,360,40)+20; F=np.array([r[f'bin|{mol}|{pos}|{a}'][k] for a in range(160,360,40)]); f=F[:,0]; e=F[:,1]
            cen=[];frac=[]
            for _ in range(5000):
                x=rng.normal(f,e); 
                if x.sum()>0.1: cen.append((x*vs).sum()/x.sum()); frac.append(x[0]/x.sum())
            P(f'  {mol} {pos} {J}: bins 160..360 '+' '.join(f'{a:+.2f}' for a in f)+f' (sd ~{np.mean(e):.2f}); total {f.sum():+.2f}+-{np.sqrt((e**2).sum()):.2f}; centroid |v| {np.median(cen):.0f} (+{np.percentile(cen,84)-np.median(cen):.0f}/-{np.median(cen)-np.percentile(cen,16):.0f}); frac in 160-200 {np.median(frac):.2f} (+{np.percentile(frac,84)-np.median(frac):.2f}/-{np.median(frac)-np.percentile(frac,16):.2f}) [N_valid={len(cen)}/5000]')
open(OUT+'/s110e.txt','w').write('\n'.join(txt))
