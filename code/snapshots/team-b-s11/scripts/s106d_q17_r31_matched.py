# Q17 round 2, step d: matched-beam 3-2/1-0 ratio of the red ridge vs the blue lobe, RADEX consequences, revised mass range (reviewer fixes 4,5).
# 3-2 cubes (Xb4d, pbcor) are convolved to a circular 1.9" beam (~ geometric mean of M18's 2.1"x1.7" beam; flux rescaled) and the flux is taken in a circular aperture of diameter 2.1" (M18 'inner 2.1"') centred on the published M18 positions
# (OF-Red, South-OF-Blue), in M18's windows ([240,360] red, [-360,-180] blue), baseline A. Errors: 3-2 = sd of random control apertures (field) (+) matched-core systematic = |model core-scaled wing| scatter;
# 1-0 = M18 table 3 (25%/30% quoted). Position tolerance: shifts of +-0.5" (M18 centroid uncertainty ~ beam/SNR + our offset 0.45").  Also the ridge flux after subtracting the core-scaled wing prediction (s106a).
# RADEX (grid out/q3d/grid_lvg500.npz): allowed (n,Tk,N) for ridge and lobe r31 (HCN & HCO+, chi2<4 on both), Tk 10..50 K, log n 3.5..7; the N_HCN/I(3-2) conversion factor ratio gives the mass range at equal 3-2 flux.
# usage: python3 -I s106d_q17_r31_matched.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
from scipy.interpolate import RegularGridInterpolator as RGI
OUT='/workspace/out/q17r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
def hms(h,m,s): return (h+m/60+s/3600)*15
OFR=(hms(10,27,51.20),-(43+54/60+15.64/3600)); SOB=(hms(10,27,51.17),-(43+54/60+21.43/3600))
nu32={'HCN':265.8864,'HCO+':267.5576}; nu10={'HCN':88.6316,'HCO+':89.1885}
M18={'HCN':{'OFR':(0.08,0.02),'SOB':(0.07,0.02)},'HCO+':{'OFR':(0.08,0.02),'SOB':(0.10,0.02)}}
BF=1.9
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,BF,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
    rad=1.05; h=(np.hypot(*np.mgrid[-int(rad/ps+1):int(rad/ps+1)+1,-int(rad/ps+1):int(rad/ps+1)+1])*ps<=rad).astype(float)/R['ba']
    def Wm(a,b):
        k=np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    Mcore=fftconvolve(Wm(-150,150),h,mode='same'); Y,X=np.mgrid[:n,:n]
    mx,my=R['pxy'](*offset_pos(S0,0.1,2.55)); step=int(round(1.0/ps)); grid=np.zeros((n,n),bool); grid[::step,::step]=True
    def dist(pos): x,y=R['pxy'](*pos); return np.hypot(X-x,Y-y)*ps
    far=(dist(S0)>2.5)&(dist(N0)>2.5)&(dist(OFR)>2.5)&(dist(SOB)>2.5)
    sel=(np.hypot(X-mx,Y-my)*ps<10)&grid&far
    for tag,pos,(a,b),key in (('OF-Red',OFR,(240,360),'OFR'),('South-OF-Blue',SOB,(-360,-180),'SOB')):
        W=fftconvolve(Wm(a,b),h,mode='same'); ixy=tuple(int(round(t)) for t in R['pxy'](*pos)[::-1]); f0=W[ixy]
        sg=float(np.std(W[sel&(np.abs(Mcore)<0.3)])); cr=Mcore[ixy]
        mc=sel&(Mcore>0.5*cr)&(Mcore<2*cr); sgm=float(np.std(W[mc])) if mc.sum()>=5 else sg
        A=np.vstack([Mcore[sel],np.ones(sel.sum())]).T; co=np.linalg.lstsq(A,W[sel],rcond=None)[0]; leak=co[0]*cr+co[1]
        # shifts
        sh=[]
        for dE,dN in ((0.5,0),(-0.5,0),(0,0.5),(0,-0.5)):
            p2=offset_pos(pos,dE,dN); ix2=tuple(int(round(t)) for t in R['pxy'](*p2)[::-1]); sh.append(float(W[ix2]))
        s10,e10=M18[mol][key]; conv_=(nu10[mol]/nu32[mol])**2
        r=f0/s10*conv_; er=r*np.hypot(np.hypot(sg,0.5*abs(leak))/f0,e10/s10)
        rx=(f0-leak)/s10*conv_
        P(f'{mol} {tag} {a}..{b} @1.9" beam, d=2.1" aperture: 3-2 flux {f0:.3f} (control sd {sg:.3f}; matched-core sd {sgm:.3f}; core-scaled leakage model {leak:+.3f}; core {cr:.1f}); shifts +-0.5": {min(sh):.3f}..{max(sh):.3f}; M18 1-0 {s10}+-{e10}; r31 = {r:.2f}+-{er:.2f} (after leakage subtraction {rx:.2f}); shifted-aperture r31 range {min(sh)/s10*conv_:.2f}..{max(sh)/s10*conv_:.2f}')
        res[f'{mol}|{tag}']=dict(f32=float(f0),sg=sg,sgm=sgm,leak=float(leak),r31=float(r),er=float(er),r31_leaksub=float(rx),shift_lo=float(min(sh)/s10*conv_),shift_hi=float(max(sh)/s10*conv_))
# significance of ridge vs lobe r31 difference
for mol in ('HCN','HCO+'):
    a=res[f'{mol}|OF-Red']; b=res[f'{mol}|South-OF-Blue']; d=np.log(a['r31']/b['r31']); e=np.hypot(a['er']/a['r31'],b['er']/b['r31'])
    P(f'{mol}: r31(ridge)/r31(lobe) = {a["r31"]/b["r31"]:.2f}; ln-difference {d:+.2f}+-{e:.2f} -> {abs(d)/e:.1f} sigma (1-0 errors from M18 table 3 included; positive-covariance of beam/aperture systematics not included)')
    res[f'{mol}|ratio_of_ratios']=dict(rr=float(a['r31']/b['r31']),nsig=float(abs(d)/e))
# ---- RADEX
d=np.load('/workspace/out/q3d/grid_lvg500.npz'); T0=np.log10(d['TK']); N0=d['LN']; C0=d['LNC']
def mk(key):
    a=d[key]; TR=a[...,0]; ok=np.isfinite(TR)&(TR>0)&(a[...,2]>0); return np.where(ok,TR,np.nan),np.where(ok,a[...,2],np.nan)
TRH,TAUH=mk('hcn'); TRO,TAUO=mk('hcop')       # [Tk,n,N,J]
rH=TRH[...,2]/TRH[...,0]; rO=TRO[...,2]/TRO[...,0]
Tk=np.broadcast_to(d['TK'][:,None,None],rH.shape); ln=np.broadcast_to(N0[None,:,None],rH.shape); lN=np.broadcast_to(C0[None,None,:],rH.shape)
dvG=150.
fH=(10**lN)/(TRH[...,2]*dvG*1.0645)   # N_HCN per unit 3-2 brightness (K km/s) [arbitrary common units]; mass at fixed flux ~ f
fO=(10**lN)/(TRO[...,2]*dvG*1.0645)
base=np.isfinite(rH)&np.isfinite(rO)&(Tk<=50)&(ln>=3.5)&(ln<=7)&(TAUH[...,2]<30)&(TAUO[...,2]<30)&(lN-np.log10(dvG)<=15.0)&(lN-np.log10(dvG)>=11.0)
def allowed(key):
    mH,eH=res[f'HCN|{key}']['r31'],res[f'HCN|{key}']['er']; mO,eO=res[f'HCO+|{key}']['r31'],res[f'HCO+|{key}']['er']
    chi=((rH-mH)/eH)**2+((rO-mO)/eO)**2; ok=base&(chi<4.0); return ok,chi
out={}
for key in ('OF-Red','South-OF-Blue'):
    ok,chi=allowed(key); 
    if ok.sum()==0: P(f'{key}: no RADEX solution with chi2<4 for both lines'); ok,chi=base&(chi<9),chi
    w=np.ones(ok.sum())
    nn=ln[ok]; tt=Tk[ok]; ff=fH[ok]; fo=fO[ok]
    P(f'RADEX {key}: allowed cells {ok.sum()}; log n(H2) 16/50/84% = {np.percentile(nn,16):.1f}/{np.percentile(nn,50):.1f}/{np.percentile(nn,84):.1f}; Tk 16/50/84 = {np.percentile(tt,16):.0f}/{np.percentile(tt,50):.0f}/{np.percentile(tt,84):.0f}; f_HCN(N per 3-2 K km/s) median {np.median(ff):.2e} [{np.percentile(ff,16):.1e},{np.percentile(ff,84):.1e}]')
    out[key]=dict(n=int(ok.sum()),logn=[float(np.percentile(nn,q)) for q in (16,50,84)],Tk=[float(np.percentile(tt,q)) for q in (16,50,84)],fH=[float(np.percentile(ff,q)) for q in (16,50,84)],fO=[float(np.percentile(fo,q)) for q in (16,50,84)])
res['radex']=out
for lab,k in (('HCN','fH'),('HCO+','fO')):
    rr=np.array(out['OF-Red'][k])/np.array(out['South-OF-Blue'][k][1]); P(f'  mass-per-unit-3-2-flux ({lab}): ridge/lobe median {out["OF-Red"][k][1]/out["South-OF-Blue"][k][1]:.2f} (16-84% ranges: ridge {out["OF-Red"][k][0]/out["South-OF-Blue"][k][1]:.2f}..{out["OF-Red"][k][2]/out["South-OF-Blue"][k][1]:.2f} of lobe median)')
# temperatures <15 K (Pereira-Santaella): what r31 does RADEX give for Tk 10-15 vs 40 at n 1e4..1e6, thin-moderate columns
for T_ in (10,14,40):
    i=np.argmin(abs(d['TK']-T_)); 
    for ln_ in (4.0,5.0,6.0):
        j=np.argmin(abs(N0-ln_)); k=np.argmin(abs(C0-(13+np.log10(dvG))))
        P(f'  RADEX Tk={d["TK"][i]:.0f} K, log n={ln_}, N/dv=1e13: r31(HCN)={rH[i,j,k]:.2f} r31(HCO+)={rO[i,j,k]:.2f}')
json.dump(res,open(OUT+'/s106d.json','w'),indent=1,default=float); open(OUT+'/s106d.txt','w').write('\n'.join(txt))
