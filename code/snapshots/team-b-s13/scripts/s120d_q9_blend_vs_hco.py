# Q9: region-integrated (Xb55 beam-aperture 'ap' spectra) CH3OH J=1-3 blend / HCO+ 2-1 ratio in K km/s (same beam, so only nu^2 factor), blend baseline systematic from s120b.
import sys,json; sys.path.insert(0,'/workspace/scripts')
from common import *
B=json.load(open('/workspace/out/q9/s120b.json'))
z=np.load('/workspace/out/regspec/Xb55.spw18.cube.pbcor.npz'); nu=z['nu']; f0=178.37506
v=C*(1-nu/(f0*1e9)); o=np.argsort(v); v=v[o]; dv=abs(np.median(np.diff(v)))
VC={'N':2792,'S':2721,'Bridge':2757,'W':2757}
res={}
for reg in ('N','S','Bridge','W'):
    y=z[f'{reg}_ap'][o]; vc=VC[reg]; m=np.isfinite(y)
    Is=[]
    for half in (500,700):
        for lo,hi in ((900,1500),(1200,2000)):
            fm=m&(abs(v-vc)>lo)&(abs(v-vc)<hi); p=np.polyfit(v[fm]-vc,y[fm],1)
            Is.append(np.sum((y-np.polyval(p,v-vc))[m&(abs(v-vc)<half)])*dv)
    IH=np.median(Is); sH=np.ptp(Is)
    rows=[r[2] for r in B[f'{reg}_ap']['rows'] if r[0] in (400,500)]
    Ib=np.median(rows); lo,hi=np.percentile(rows,[16,84])
    k=(f0/165.06113)**2
    res[reg]=dict(IHCO21=float(IH),IHCO21_spread=float(sH),Iblend=float(Ib),blend_16_84=[float(lo),float(hi)],ratioK=float(Ib/IH*k),ratioK_16_84=[float(lo/IH*k),float(hi/IH*k)])
    print(reg,res[reg])
json.dump(res,open('/workspace/out/q9/s120d.json','w'),indent=1)
