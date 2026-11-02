# Q9: CH3OH excitation with the blend baseline/window systematic propagated (s120b), 254.015 GHz line at the Xb55 beam (s21),
# optically thin LTE W254/Wblend -> Trot (Monte Carlo), the J=4,5 limits and the 'Q-branch excess' (observed blend / blend predicted from Trot set by the other lines).
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
rng=np.random.default_rng(1)
B=json.load(open('/workspace/out/q9/s120b.json')); S21=json.load(open('/workspace/out/s21_ch3oh254.json'))
sm=S21['smoothed']; mad=S21['smoothed_mad']
L254=(254.015377,20.09,-4.72074,5)
BL=[(165.050175,23.37,-4.62911,3),(165.061130,28.01,-4.63041,5),(165.099240,34.98,-4.63216,7)]
J4=(165.190475,44.26,-4.6340,9); J5=(165.369341,55.87,-4.6349,11)   # gu=2J+1 (E-type (-1,J)-(0,J) series), A from out/linelist.csv-class values (approx.)
def W(l,T): f,E,la,g=l; return g*10**la*np.exp(-E/T)/f**2
def Wb(T): return sum(W(b,T) for b in BL)
Tg=np.linspace(3,60,1141)
Rg=W(L254,Tg)/Wb(Tg)           # K km/s ratio 254/blend (monotonic decreasing)
nuc=(163.53/251.68)**2         # Jy->K factor between the two lines (same beam)
out={}
for reg in ('N','S','Bridge','W'):
    key={'N':'N','S':'S','Bridge':'Bridge','W':'W'}[reg]
    rows=[r for r in B[f'{reg}_ap']['rows'] if r[0] in (400,500)]
    Ib=np.array([r[2] for r in rows]); sb=np.median([r[3] for r in rows])
    cen=np.median(Ib); lo,hi=Ib.min(),Ib.max()
    s=sm.get(key) if isinstance(sm,dict) else None
    if s is None: print(reg,'no 254 smoothed'); out[reg]=dict(Iblend=cen,range=[lo,hi]); continue
    a=s['I254_55'] if 'I254_55' in s else None
    ea=np.hypot(mad,0.05*a)
    # Monte Carlo
    N_=20000; bs=rng.choice(Ib,N_)+rng.normal(0,sb,N_); as_=a+rng.normal(0,ea,N_)
    r=as_/bs*nuc
    ok=(bs>0)&(r>0)
    T=np.interp(r[ok],Rg[::-1],Tg[::-1])
    pct=np.percentile(T,[16,50,84]) if ok.sum()>0 else None
    # fraction of draws with ratio above the max/ below min achievable
    out[reg]=dict(Iblend=float(cen),blend_range=[float(lo),float(hi)],I254=float(a),eI254=float(ea),R_K=float(np.median(r[ok])),
        Tlo_med_hi=[float(x) for x in pct],frac_R_below_T60=float((r<Rg[-1]).mean()),
        R_K_16_84=[float(x) for x in np.percentile(r,[16,84])])
    print(reg,json.dumps(out[reg]))
print('LTE W254/Wblend:',{T:round(float(W(L254,T)/Wb(T)),3) for T in (4,6,8,10,12,15,20,30)})
print('J4/blend, J5/blend LTE:',{T:(round(float(W(J4,T)/Wb(T)),3),round(float(W(J5,T)/Wb(T)),3)) for T in (6,10,15,20,30)})
# Q-branch 'excess': team-a states Q-branch is 1.6x above any model. Observed J4 (3-sigma UL) vs blend: allowed T from J4 UL
for reg,(ul4) in (('N',0.194),('S',0.262)):
    b=out[reg]['Iblend']; r=ul4/b*(165.19/165.06)**2*0   # placeholder, see note
json.dump(out,open('/workspace/out/q9/s120c.json','w'),indent=1)
