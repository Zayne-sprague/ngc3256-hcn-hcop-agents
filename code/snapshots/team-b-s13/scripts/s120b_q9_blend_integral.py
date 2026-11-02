# Q9: model-independent integral of the Xb55 CH3OH J_-1-J_0 E (J=1-3) blend at N,S,Bridge,W (peak pixel and beam aperture),
# for a grid of baseline treatments / windows -> window+baseline systematic. Compare team-a p0001 (2.23 N, reviewer 1.84) and our earlier 0.97 (N) / 1.38 (S).
import sys,json; sys.path.insert(0,'/workspace/scripts')
from common import *
z=np.load('/workspace/out/regspec/Xb55.spw22.cube.pbcor.npz'); nu=z['nu']
f0=165.061130; v=C*(1-nu/(f0*1e9)); o=np.argsort(v); v=v[o]; dv=abs(np.median(np.diff(v)))
VC={'N':2792,'S':2721,'Bridge':2757,'W':2757}   # region centres (s13 region vsys; Bridge/W: mean)
# lines known in window (J=4..8 in J=2 frame) to exclude from baseline
LINEF=[165.190475,165.369341,165.678649,166.169098]
out={}
res={}
for reg in ('N','S','Bridge','W'):
  for st in ('pix','ap'):
    y=z[f'{reg}_{st}'][o]*1.0; vc=VC[reg]; rows=[]
    for half in (300,400,500,600):
      for base in ('const<2400','lin2','lin1','poly2','const_all'):
        m=np.isfinite(y); line=abs(v-vc)<half
        excl=np.zeros_like(line)
        for f in LINEF: excl|=abs(v-C*(1-f/f0)-vc)<130
        if base=='const<2400': fm=(v<2400)&~excl&m&(v>200)   # blue side only
        elif base=='const_all': fm=m&~line&~excl&(abs(v-vc)<2500)
        else: fm=m&~line&~excl&(abs(v-vc)<{'lin2':1200,'lin1':800,'poly2':1500}[base])&(abs(v-vc)>half)
        deg={'const<2400':0,'const_all':0,'lin2':1,'lin1':1,'poly2':2}[base]
        p=np.polyfit(v[fm]-vc,y[fm],deg); b=np.polyval(p,v-vc)
        I=np.sum((y-b)[line&m])*dv; rms=np.std((y-b)[fm]); e=rms*np.sqrt(line.sum())*dv*1.3
        rows.append((half,base,I,e))
    Is=np.array([r[2] for r in rows]); e0=np.median([r[3] for r in rows])
    res[f'{reg}_{st}']=dict(rows=rows,med=float(np.median(Is)),lo=float(Is.min()),hi=float(Is.max()),stat=float(e0))
    print(f'{reg:6s} {st:3s} I[Jy/b km/s]: median {np.median(Is):.2f} range {Is.min():.2f}-{Is.max():.2f} stat +-{e0:.2f} | half=400: '+' '.join(f'{r[1]}:{r[2]:.2f}' for r in rows if r[0]==400))
json.dump(res,open('/workspace/out/q9/s120b.json','w'),indent=1,default=float)
