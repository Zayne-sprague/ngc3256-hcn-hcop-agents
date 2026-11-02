# CH3OH J_-1-J_0 E Q-branch (J=1..5) in Xb55 spw22: joint template (HCO+ 2-1 profile) NNLS fit of the 5 lines for N and S beam apertures
# and a Boltzmann (optically thin LTE) rotation temperature from the amplitude ratios (a_J ~ S mu^2 exp(-Eu/T)).
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
exec(open('/workspace/scripts/s12_line_table.py').read().split("rows=[]\nfor lab,f0,lst,note in LINES:")[0])
from scipy.optimize import nnls, curve_fit
LJ={1:(165.050175,23.367,5.38627),2:(165.061130,28.011,8.94855),3:(165.099240,34.976,12.46907),4:(165.190475,44.263,15.92557),5:(165.369341,55.871,19.29095)}
out=[]
for reg in ('N','S'):
    T=make_template(ous:='Xb55',reg,'ap'); nu,y=spec('Xb55','spw22',reg,'ap')
    f0=LJ[2][0]; v=vel(nu,f0); o=np.argsort(v); v,y=v[o],y[o]; ok=np.isfinite(y); v,y=v[ok],y[ok]; dv=abs(np.median(np.diff(v)))
    tn=lambda vv: np.interp(vv,T['v'],T['t'],left=0,right=0)/T['t'].max()
    vc=0.5*(T['hi']+T['lo']); half=0.5*(T['hi']-T['lo'])
    shift={J:C*(1-f*(1-vc/C)/f0) for J,(f,E,S_) in LJ.items()}      # velocity of line J in the J=2 frame if the gas is at vc
    free=abs(v-vc)<2500
    for J in LJ: free&=abs(v-shift[J])>half+100
    base=np.median(y[free]); yy=y-base; rms=np.std(yy[free])
    cols=[tn(v-shift[J]+vc) for J in LJ]
    A=np.array(cols).T; m=abs(v-vc)<900
    amp,rn=nnls(A[m]/rms,yy[m]/rms)
    cov=np.linalg.inv((A[m].T@A[m])/rms**2); err=np.sqrt(np.diag(cov))
    I=amp*np.array([ (tn(np.arange(-1000,1000,dv)+vc)).sum()*dv ]*5)
    Ie=err*(tn(np.arange(-1000,1000,dv)+vc)).sum()*dv
    Eu=np.array([LJ[j][1] for j in LJ]); S=np.array([LJ[j][2] for j in LJ])
    def lt(E,lnN,T): return lnN-E/T
    ok=I>0
    # fit ln(I/S) = a - Eu/Trot, weights by errors
    x=Eu; yv=np.log(np.clip(I,1e-3,None)/S); sg=Ie/np.clip(I,1e-3,None)
    p,cv=curve_fit(lt,x,yv,p0=[0,30],sigma=sg,absolute_sigma=True)
    pred=np.exp(p[0]-Eu/p[1])*S; chi=np.sum(((I-pred)/Ie)**2)
    print(f'{reg}: I_J (Jy/beam km/s) '+', '.join(f'J={j}:{i:.3f}+-{e:.3f}' for j,i,e in zip(LJ,I,Ie))+f' | Trot={p[1]:.0f}+-{np.sqrt(cv[1,1]):.0f} K, chi2={chi:.1f} (3 dof)')
    out.append((reg,I.tolist(),Ie.tolist(),float(p[1]),float(np.sqrt(cv[1,1])),float(chi)))
import json; json.dump(out,open('/workspace/out/s17_ch3oh_ladder.json','w'))
