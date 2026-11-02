# Q23 round 2 (reviewer fixes). Convention: v = radio velocity; LOW frequency spw edge = ch 0 = RED (+v) side. The ramp is on the RED side.
# (1) direct test red(+420..+600) vs blue(-600..-420) levels: observed vs ramp-model prediction (in-sample fit, leave-one-region-out, and no-ramp controls)
# (2) ramp shape: triangle(c0), exponential(tau), flat-to-window-start (pessimistic), Xb55-empirical (s133/profile of Xb55 spw18) fitted on line-free ch < chfit
# (3) estimators: E0 = s30 (median of v<-600 & +450..600, +-350); E1 = blue-only median (v<-450), +-350; E2 = blue-only linear fit (-900<v<-450), +-350; E3 = E1 baseline, +-250 window (window check).
#     'templ' = every estimator recomputed after subtracting each fitted ramp shape (ramp under window + red reference channels). Integrals, HCN/HCO+ ratios, T43/T32 etc at N,S,Bridge.
import sys, json, numpy as np; sys.path.insert(0,'/workspace/scripts')
from common import vrad, VSYS
S131=json.load(open('/workspace/out/q22/s131.json'))
F0={'HCN':{2:177.26111e9,3:265.88618e9,4:354.50547e9},'HCO+':{2:178.37506e9,3:267.55763e9,4:356.73424e9}}
LINES=[('Xb55',16,'HCN',2),('Xb55',18,'HCO+',2),('Xb4f',18,'HCN',3),('Xb4f',16,'HCO+',3),('Xb47',16,'HCN',4),('Xb47',18,'HCO+',4)]
REG=['N','S','Bridge']; CAL=0.07
def load(ous,spw,mol,J,kind='pix'):
    z=np.load(f'/workspace/out/regspec/{ous}.spw{spw}.cube.pbcor.npz'); nu=z['nu']; v=vrad(nu,F0[mol][J])-VSYS
    Y={r:np.nan_to_num(z[f'{r}_{kind}'])*1e3 for r in REG}; cont={r:S131[f'{ous}.spw{spw}.{r}']['cont'] for r in REG}
    return v,np.arange(len(nu)),Y,cont
# empirical template (fraction of continuum) from Xb55 spw18 stacked low-edge profile, 10-ch bins (s133.txt), channels at 1.95 MHz
EMP_BINS=np.array([-0.147,-0.112,-0.104,-0.096,-0.138,-0.046,-0.125,-0.058,-0.075,0.045,0.083,0.099]); EMP_CEN=np.arange(12)*10+5
def shapes(ch,dvscale):
    """unit-amplitude shapes (value 1 at ch0) -> dict name -> array. dvscale: ch width relative to Xb55 1.95MHz (for empirical template in MHz)."""
    S={}
    for c0 in (50,100,200,400): S[f'tri{c0}']=np.clip(1-ch/c0,0,None)
    for tau in (15,30,60): S[f'exp{tau}']=np.exp(-ch/tau)
    e=np.interp(ch*dvscale,EMP_CEN,-EMP_BINS/0.147,right=0.); e[ch*dvscale>95]=0; S['emp']=np.clip(e,0,None)   # empirical Xb55 spw18 shape, 0 beyond ch 95 (positive bins = line wing)
    return S
def fitamp(Y,cont,fm,sh):
    num=sum(np.sum(cont[r]*sh[fm]*Y[r][fm]) for r in REG); den=sum(np.sum((cont[r]*sh[fm])**2) for r in REG)
    f=num/den; res=np.concatenate([Y[r][fm]-f*cont[r]*sh[fm] for r in REG]); s=res.std()*1.3
    return f,s/np.sqrt(den),float(np.sum(res**2)/res.var()/1.0)
def estimators(y,v,W0=350):
    dv=abs(np.median(np.diff(v))); out={}
    ref=(v<-600)|((v>450)&(v<600)); b=np.median(y[ref]); m=abs(v)<W0
    out['E0']=((y[m]-b).sum()*dv/1e3, np.std(y[ref]-b)*1.3*np.sqrt(m.sum())*dv/1e3)
    bl=v<-450; b1=np.median(y[bl]); out['E1']=((y[m]-b1).sum()*dv/1e3, np.std(y[bl]-b1)*1.3*np.sqrt(m.sum())*dv/1e3)
    bl2=(v<-450)&(v>-900); p=np.polyfit(v[bl2],y[bl2],1); base=np.polyval(p,v)      # linear baseline on blue side only, -900<v<-450, same +-350 window
    out['E2']=((y[m]-base[m]).sum()*dv/1e3, np.std(y[bl2]-np.polyval(p,v[bl2]))*1.3*np.sqrt(m.sum())*dv/1e3)
    m3=abs(v)<250                                                                    # window check (not a ramp test): E1 baseline, +-250 window
    out['E3']=((y[m3]-b1).sum()*dv/1e3, np.std(y[bl]-b1)*1.3*np.sqrt(m3.sum())*dv/1e3)
    return out
res={'direct':{},'fits':{},'est':{}}
print('=== (1) DIRECT TEST: mean level red(+420..+600) and blue(-600..-420) [mJy/ch] per region; model-predicted red level; control = spws with no ramp')
for ous,spw,mol,J in LINES:
    v,ch,Y,cont=load(ous,spw,mol,J); dv=abs(np.median(np.diff(v))); chl=int(np.argmin(abs(v)))
    red=(v>420)&(v<600); blue=(v<-420)&(v>-600); nr=red.sum(); nb=blue.sum()
    rms={r:np.std(Y[r][v<-450]) for r in REG}
    # fit range: line-free red side, v>+380 (ch< chl - 380/dv) ; the full ramp-bearing line-free channels
    chfit=chl-int(np.ceil(380/dv)); fm=ch<chfit; dvs=dv*(3.91 if ous!='Xb55' else 1.95)/dv/1.95  # channel width in units of Xb55 channels
    chw=(3.91 if ous!='Xb55' else 1.95)/1.95
    sh=shapes(ch,chw)
    key=f'{ous}.spw{spw}.{mol}{J}'; d={'chl':chl,'chfit':int(chfit),'dv':dv,'rows':{}}
    print(f'\n{key}: line ch {chl}, dv {dv:.2f} km/s/ch, red band ch {ch[red].min()}-{ch[red].max()} ({nr} ch), line-free fit range ch<{chfit}')
    fits={}
    for nm,s in sh.items():
        f,fe,c2=fitamp(Y,cont,fm,s); fits[nm]=(f,fe,c2)
    flat=sum(np.sum(Y[r][fm]**2) for r in REG); 
    for nm,(f,fe,c2) in fits.items():
        res_=sum(np.sum((Y[r][fm]-f*cont[r]*sh[nm][fm])**2) for r in REG); s2=np.mean([rms[r]**2 for r in REG])*1.3**2
        fits[nm]=(f,fe,float((flat-res_)/s2))   # chi2 improvement over zero-ramp (corr. inflated)
    d['fits']={k:[float(x) for x in v_] for k,v_ in fits.items()}
    print('  ramp amplitude f (x cont) +-err ; delta-chi2 vs no ramp : '+'; '.join(f'{k} {v_[0]:+.3f}+-{v_[1]:.3f} ({v_[2]:.0f})' for k,v_ in fits.items()))
    best=max(fits,key=lambda k:fits[k][2]); d['best']=best
    for r in REG:
        o=Y[r][red].mean(); b=Y[r][blue].mean(); e_r=rms[r]*1.3/np.sqrt(nr); e_b=rms[r]*1.3/np.sqrt(nb)
        pred=fits[best][0]*cont[r]*sh[best][red].mean()
        # leave-one-region-out: amplitude from the other two regions
        oth=[q for q in REG if q!=r]
        num=sum(np.sum(cont[q]*sh[best][fm]*Y[q][fm]) for q in oth); den=sum(np.sum((cont[q]*sh[best][fm])**2) for q in oth); floo=num/den
        predl=floo*cont[r]*sh[best][red].mean()
        d['rows'][r]=dict(red=o,blue=b,diff=b-o,err=float(np.hypot(e_r,e_b)),pred_red=float(pred),pred_red_loo=float(predl),cont=cont[r])
        print(f'  {r:6s} cont {cont[r]:5.1f}: red {o:+.2f}+-{e_r:.2f} blue {b:+.2f}+-{e_b:.2f} blue-red {b-o:+.2f}+-{np.hypot(e_r,e_b):.2f} | predicted red level ({best}) {pred:+.2f}, leave-one-region-out {predl:+.2f}')
    res['direct'][key]=d
    # ---- estimators
    res['est'][key]={}
    for r in REG:
        E=estimators(Y[r],v); row={k:list(map(float,x)) for k,x in E.items()}
        # ramp-template integral under E2/E0 windows: integral of (ramp - ramp-baseline) -> bias; use best amplitudes +-1sig for each shape
        row['templ']={}
        for nm,(f,fe,c2) in fits.items():
            ramp=f*cont[r]*sh[nm]
            e3=estimators(Y[r]-ramp,v)
            row['templ'][nm]={k:float(e3[k][0]) for k in e3}
        res['est'][key][r]=row
json.dump(res,open('/workspace/out/q23/s145.json','w'),default=float)
# ---- summary of estimators
keys={('HCN',2):'Xb55.spw16.HCN2',('HCO+',2):'Xb55.spw18.HCO+2',('HCN',3):'Xb4f.spw18.HCN3',('HCO+',3):'Xb4f.spw16.HCO+3',('HCN',4):'Xb47.spw16.HCN4',('HCO+',4):'Xb47.spw18.HCO+4'}
# nominal shape: Xb55-empirical (MHz-scaled); pessimistic bound: tri400 (ramp ~flat through the window red wing)
print('\n=== (2) INTEGRALS, Jy/beam km/s (peak pixel, native beam, +-350 km/s unless noted). E0 s30; E1 blue-only median; E2 blue-only linear; E3 = E1 with +-250 window')
print('    corrected = integral after subtracting the fitted ramp (shapes emp = nominal, exp30, tri100, tri400 = pessimistic), with E0 and E1 baselines')
for key in res['est']:
    print(key)
    for r in REG:
        row=res['est'][key][r]; t=row['templ']
        print(f'  {r:6s} E0 {row["E0"][0]:7.3f}+-{row["E0"][1]:.3f} E1 {row["E1"][0]:7.3f} E2 {row["E2"][0]:7.3f} E3 {row["E3"][0]:7.3f} | corrected E0: '+' '.join(f'{k} {t[k]["E0"]:.3f}' for k in ('emp','exp30','tri100','tri400'))+' | corrected E1: '+' '.join(f'{k} {t[k]["E1"]:.3f}' for k in ('emp','exp30','tri100','tri400')))
def val(mol,J,r,est,shape=None):
    row=res['est'][keys[(mol,J)]][r]; return (row['templ'][shape][est] if shape else row[est][0]), row[est][1]
def relerr(mol,J,r,est):
    v_,e=val(mol,J,r,est); return np.hypot(e/abs(v_),CAL)
print('\n=== (3) RATIOS: shift of each estimator relative to E0 (and in units of the total ratio error = stat + 7% cal per line)')
def rat(num,den,r,est,shape=None):
    return val(num[0],num[1],r,est,shape)[0]/val(den[0],den[1],r,est,shape)[0]
RAT=[('HCN/HCO+ 2-1',('HCN',2),('HCO+',2)),('HCN/HCO+ 3-2',('HCN',3),('HCO+',3)),('HCN/HCO+ 4-3',('HCN',4),('HCO+',4)),
     ('HCN T43/T32',('HCN',4),('HCN',3)),('HCN T21/T32',('HCN',2),('HCN',3)),('HCO+ T43/T32',('HCO+',4),('HCO+',3)),('HCO+ T21/T32',('HCO+',2),('HCO+',3))]
for r in REG:
    print(r)
    for nm,a,b in RAT:
        R0=rat(a,b,r,'E0'); sig=np.hypot(relerr(*a,r,'E0'),relerr(*b,r,'E0')); s_=f'  {nm:14s} E0 {R0:.3f}'
        for est in ('E1','E2','E3'):
            Re=rat(a,b,r,est); s_+=f' | {est} {100*(Re/R0-1):+5.1f}% ({(Re/R0-1)/sig:+.2f}s)'
        for shp in ('emp','tri400'):
            for est in ('E0','E1'):
                Re=rat(a,b,r,est,shp); s_+=f' | {est}-{shp} {100*(Re/R0-1):+5.1f}% ({(Re/R0-1)/sig:+.2f}s)'
        print(s_)
