# LTE (optically thin, beam-filling-matched) test: CH3OH 2(0)-1(-1)E 254.0154 GHz (Xb4f, smoothed to Xb55 beam) vs J=1-3 J(-1)-J(0) E blend (Xb55)
import numpy as np, json
S=json.load(open('/workspace/out/s21_ch3oh254.json'))
# Splatalogue/JPL (out/linelist.csv) values: (rest GHz, Eu K, log10 Aij, gu)
L254=(254.015377,20.09,-4.72074,5)
BL=[(165.050175,23.37,-4.62911,3),(165.061130,28.01,-4.63041,5),(165.099240,34.98,-4.63216,7)]
def W(l,T): f,E,la,g=l; return g*10**la*np.exp(-E/T)/f**2
def R(T): return W(L254,T)/sum(W(b,T) for b in BL)       # W(254)/W(blend), K km/s ratio
Ts=np.array([4,5,6,8,10,12,15,20,30,50,100.])
print('LTE prediction W254/Wblend (K km/s, thin):',dict(zip(Ts,np.round(R(Ts),3))))
nu254=251.68; nub=163.53   # observed GHz of the two lines
sm=S['smoothed']; mad=S['smoothed_mad']
# Xb55 blend (J=1-3, HCO+ template, +-900 km/s) from report #66: N 0.97, S 1.38, +-0.07 stat, +-30% template systematic
bl={'N':(0.97,0.07),'S':(1.38,0.07)}; bl_mid=((0.97+1.38)/2,np.hypot(0.07,0.3*1.18))
rows={}
for reg in ('N','S'):
    a=sm[reg]['I254_55']; ea=np.hypot(mad,0.05*a)
    b,eb=bl[reg]; eb=np.hypot(eb,0.3*b)
    r_Jy=a/b; er=r_Jy*np.hypot(ea/max(abs(a),1e-9),eb/b) if a>0 else np.nan
    rK=r_Jy*(nub/nu254)**2
    rows[reg]=dict(I254_55=a,eI254=ea,Iblend=b,eIblend=eb,R_Jy=r_Jy,R_K=rK,eR_K=er*(nub/nu254)**2 if a>0 else np.nan,UL3_R_K=3*ea/b*(nub/nu254)**2)
    print(reg,rows[reg])
# Bridge/mid-point: 254 line at the Bridge pixel (between N and S) vs blend mean of the N and S pixels in the 9" beam
a=sm['Bridge']['I254_55']; ea=np.hypot(mad,0.05*a); b,eb=bl_mid
rK=a/b*(nub/nu254)**2; erK=rK*np.hypot(ea/a,eb/b); rows['Bridge']=dict(I254_55=a,eI254=ea,Iblend=b,eIblend=eb,R_K=rK,eR_K=erK)
print('Bridge/mid',rows['Bridge'])
# temperature constraint (Bridge/mid and S)
for reg in ('S','Bridge'):
    r=rows[reg]['R_K']; e=rows[reg]['eR_K']; Tg=np.linspace(3,200,2000); m=np.abs(R(Tg)-r)<e*1.0
    lo=np.abs(R(Tg)-(r+e)).argmin(); hi=np.abs(R(Tg)-(r-e)).argmin()
    print(reg,'R_K=%.3f+-%.3f -> LTE T:'%(r,e),'best',Tg[np.abs(R(Tg)-r).argmin()],' 1sig range',sorted([Tg[lo],Tg[hi]]),' max R(T)=%.3f at T->inf %.3f'%(R(Tg).max(),R(1e4)))
json.dump(rows,open('/workspace/out/s22_ch3oh_lte_ratio.json','w'),indent=1,default=float)
