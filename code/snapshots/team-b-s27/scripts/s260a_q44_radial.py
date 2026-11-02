"""Q44(a): radial test of H13CO+/HCO+ and H13CN/HCN 1-0 at N, S, Bridge at ONE 2.2" beam (X961/X95d, 2015.1.00993.S short-spacing complete), +-150 km/s box
(SiO 2-1 is 320 km/s away, outside the box). Aperture means (r = 0.3" ~ peak pixel, 1.1", 2.2", 3.3" and the 2.2-4.4" annulus-like difference computed from r=3.3 and 1.1) with the s185 machinery
(same templates/baselines/null-aperture error calibration, nulls evaluated with the same radius).  exec()s s185 up to the box-moment section (no re-edit of s185).
usage: python3 -I s260a_q44_radial.py -> out/q44/s260a.json, s260a.txt"""
import json,sys
src=open('/workspace/scripts/s185_q29_weaklines.py').read()
src=src[:src.index("BOXRES={}")]
src=src.replace("open(O+'s185.json','w')","open('/tmp/_s185_dummy.json','w')").replace("open(O+'s185.txt','w')","open('/tmp/_s185_dummy.txt','w')")
exec(src)
OUTR={}; L2=[]
def Q(s=''): print(s,flush=True); L2.append(s)
box=(-150,150)
for RAD in (0.3,1.1,2.2,3.3):
    nul={}
    for (rr,ang,p,mi) in nulls:
        for kk,(a,e) in boxfit(spectra(p),box).items(): nul.setdefault(kk,[]).append((a,e))
    for rn in ('N','S','Bridge'):
        sp=spectra(positions(rn)); R=boxfit(sp,box); IO=refint(sp['X961_spw0'],F0['HCOp'],box); IH=refint(sp['X961_spw0'],F0['HCN'],box)
        for kk,(a,e) in R.items():
            l=kk[0]
            if l not in ('H13COp','H13CN'): continue
            if l=='H13COp' and kk[1]!='comb': continue
            L=np.array(nul[kk]); c=max(1.0,float(np.sqrt(np.mean((L[:,0]/L[:,1])**2)))); e2=max(e*c,float(np.std(L[:,0])))
            ref=IH if l=='H13CN' else IO
            OUTR[f'{RAD}|{rn}|{l}']=dict(amp=a,err=e2,ref=ref,HCOp=IO,HCN=IH,ratio=a/ref,ratio_err=float(np.hypot(e2/ref,0.08*abs(a/ref))) if False else e2/ref)
            Q(f' r={RAD:.1f}" {rn:6s} {l:6s} {a:6.2f}+-{e2:4.2f} K km/s  / {"HCN" if l=="H13CN" else "HCO+"} {ref:6.1f} = {a/ref:.4f}+-{e2/ref:.4f}  (HCO+ {IO:.1f} HCN {IH:.1f})')
json.dump(OUTR,open('/workspace/out/q44/s260a.json','w'),indent=1,default=float)
open('/workspace/out/q44/s260a.txt','w').write('\n'.join(L2))
