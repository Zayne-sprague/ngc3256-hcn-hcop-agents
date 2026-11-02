# Q32 (round 2): ff / size of the CO(2-1) cores from (1) T_b(theta) at the core pixel for 3 velocity measures and 5 radial profile families,
# (2) direct 2D forward fits (core + extended + clumps + const, beam-convolved) to the native 0.18" map. Tk bounds per model.
import sys,json,numpy as np,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/scripts'); from q32_maps import *
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit,least_squares
H=6.62607e-27; KB=1.380649e-16; NU=230.538e9
Jf=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1); Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
tk=lambda T,tau,ff=1.0: float(Jinv(T/(ff*(1-np.exp(-tau)))+Jf(2.73)))
tt=json.load(open('/workspace/out/q30/tau_tk.json'))
O,(bmaj,bmin,pix)=load(); bg=np.sqrt(bmaj*bmin)
TH=np.array([bg,.26,.38,.57,.75,1.0,1.7])
def conv(m,th):
    if th==bg: return m
    return gaussian_filter(m,np.sqrt(th**2-bg**2)/pix/2.3548,mode='constant')
# ---- radial families: intrinsic circular profiles; value at centre of the map convolved with a Gaussian beam of FWHM th (total beam)
rr=np.linspace(0,12,6000)
def cen(prof,th):   # profile(r) -> central convolved value (unit-area Gaussian kernel, sigma=th/2.355)
    s=th/2.3548; k=np.exp(-rr**2/(2*s*s))*rr/s**2; return np.trapezoid(prof(rr)*k,rr)
def gauss_core(T0,s,th): return T0*s**2/(s**2+th**2)
core={'gauss':lambda r,T0,s: T0*np.exp(-4*np.log(2)*(r/s)**2),
      'exp':lambda r,T0,s: T0*np.exp(-np.log(2)*2*r/s),             # FWHM s
      'moffat2':lambda r,T0,s: T0*(1+(r/(s/2/np.sqrt(np.sqrt(2)-1)))**2)**-2}  # beta=2, FWHM s
def ffun(kind,T0,s,th): return gauss_core(T0,s,th) if kind=='gauss' else cen(lambda r:core[kind](r,T0,s),th)
MODELS={
 'M1 gauss+const':     (lambda th,T0,s,Te: np.array([ffun('gauss',T0,s,t) for t in th])+Te, [30,.4,4],([0,.03,0],[300,3,30]), 'gauss'),
 'M2 gauss+free-ext':  (lambda th,T0,s,B,s2: np.array([ffun('gauss',T0,s,t)+ffun('gauss',B,s2,t) for t in th]), [30,.4,5,3],([0,.03,0,1.2],[300,1.2,40,30]),'gauss'),
 'M3 single gauss':    (lambda th,T0,s: np.array([ffun('gauss',T0,s,t) for t in th]), [30,.5],([0,.03],[300,10]),'gauss'),
 'M4 exp+const':       (lambda th,T0,s,Te: np.array([ffun('exp',T0,s,t) for t in th])+Te, [30,.4,4],([0,.03,0],[300,3,30]),'exp'),
 'M5 moffat2+const':   (lambda th,T0,s,Te: np.array([ffun('moffat2',T0,s,t) for t in th])+Te, [30,.4,4],([0,.03,0],[300,3,30]),'moffat2'),
}
out={'th':TH.tolist(),'1d':{},'flux':{},'2d':{}}
RS=np.random.default_rng(1)
for nm in 'SN':
  for meas in ['peak','3ch','HCNwin']:
    m=O[nm]['maps'][meas]; HW=O[nm]['HW']
    T=np.array([conv(m,t)[HW,HW] for t in TH]); sig=0.05*T+0.3
    key=f'{nm}_{meas}'; out['1d'][key]={'T':T.tolist()}
    # explicit flux (peak flux density in the beam, Jy/beam proportional T*theta^2) 
    out['flux'][key]={'Tb':T.tolist(),'flux_ratio_0.26_over_1.0':float(T[1]*.26**2/(T[5]*1.0**2)),'flux_ratio_native_over_1.0':float(T[0]*bg**2/(T[5]*1.0)),'Tb_ratio_0.26/1.0':float(T[1]/T[5])}
    for mn,(f,p0,bd,kind) in MODELS.items():
        try:
            p,cv=curve_fit(f,TH,T,p0=p0,sigma=sig,bounds=bd,maxfev=20000)
        except Exception as e: out['1d'][key][mn]={'fail':str(e)}; print(key,mn,'FAIL',e); continue
        e=np.sqrt(np.diag(cv)); chi2=float(np.sum(((f(TH,*p)-T)/sig)**2))
        T0,s=p[0],p[1]; ff26=ffun(kind,1.0,s,.26)/1.0 if kind!='gauss' else s**2/(s**2+.26**2)
        ffn=s**2/(s**2+bg**2) if kind=='gauss' else ffun(kind,1.0,s,bg)
        rec={'p':p.tolist(),'e':e.tolist(),'chi2':chi2,'dof':len(TH)-len(p),'T0':float(T0),'s_fwhm':float(s),'ff_0.26':float(ff26),'ff_native':float(ffn)}
        # MC for ff(0.26) and s
        ss=[]
        for _ in range(150):
            Ti=T*(1+0.04*RS.standard_normal())*(1+0.04*RS.standard_normal(len(T)))+0.15*RS.standard_normal(len(T))
            try: q,_=curve_fit(f,TH,Ti,p0=p,sigma=sig,bounds=bd,maxfev=5000); ss.append([q[0],q[1],(q[1]**2/(q[1]**2+.26**2)) if kind=='gauss' else ffun(kind,1.,q[1],.26)])
            except Exception: pass
        ss=np.array(ss); rec['mc16_50_84']={'T0':np.percentile(ss[:,0],[16,50,84]).tolist(),'s':np.percentile(ss[:,1],[16,50,84]).tolist(),'ff26':np.percentile(ss[:,2],[16,50,84]).tolist()}
        if meas=='peak':
            tau={R:tt[f'{nm}_ratio_peak'][R]['tau'][1] for R in ('R40','R60','R90')}
            rec['Tk_core_T0_ff1']={R:tk(T0,t) for R,t in tau.items()}
        out['1d'][key][mn]=rec
        print(key,mn,'T0=%.1f s=%.2f ff26=%.2f ffnat=%.2f chi2=%.1f/%d'%(T0,s,ff26,ffn,chi2,rec['dof']),'s16-84=',np.round(rec['mc16_50_84']['s'],2),'ff26=',np.round(rec['mc16_50_84']['ff26'],2))
# ---- direct 2D forward fits to native peak and HCN-window maps
def blob(X,Y,A,x0,y0,s,b2):  # circular intrinsic Gaussian FWHM s convolved with circular beam b (Gaussian)
    S2=s**2+b2; return A*(s**2/S2)*np.exp(-4*np.log(2)*((X-x0)**2+(Y-y0)**2)/S2)
def fit2d(m,HW,R,nclump):
    Y,X=(np.mgrid[-HW:HW+1,-HW:HW+1])*pix; X=-X
    sel=np.hypot(X,Y)<R
    def model(p):
        r=blob(X,Y,p[0],p[1],p[2],p[3],bg**2)+blob(X,Y,p[4],p[5],p[6],p[7],bg**2)+p[8]
        for i in range(nclump): q=p[9+4*i:13+4*i]; r=r+blob(X,Y,q[0],q[1],q[2],q[3],bg**2)
        return r
    p=[m[HW,HW],0,0,.4, 4,0,0,2.5, 1.0]; L=[0,-.15,-.15,.05, 0,-1,-1,1.0, -3]; U=[300,.15,.15,2, 60,1,1,20, 20]
    for i in range(nclump):   # greedy: add a clump at the largest positive residual outside the core (r>0.3")
        res=m-model(np.array(p+[0]*0+[0,0,0,0]*(nclump-i))[:9+4*i] if False else np.array(p+[0,0,0,.3]*(nclump-i)))
        res=np.where(sel&(np.hypot(X,Y)>0.3),res,-9); jy,jx=np.unravel_index(np.argmax(res),res.shape)
        p+= [max(res[jy,jx],1.),X[jy,jx],Y[jy,jx],.3]; L+=[0,X[jy,jx]-.3,Y[jy,jx]-.3,.1]; U+=[60,X[jy,jx]+.3,Y[jy,jx]+.3,1.5]
        f=lambda q: (model(np.r_[q,[0,0,0,.3]*(nclump-i-1)])-m)[sel]
        r=least_squares(f,np.clip(p,np.array(L)+1e-9,np.array(U)-1e-9),bounds=(L,U),max_nfev=2000); p=list(r.x)
    f=lambda q: (model(q)-m)[sel]
    r=least_squares(f,np.clip(p,np.array(L)+1e-9,np.array(U)-1e-9),bounds=(L,U),max_nfev=4000)
    n=sel.sum(); k=len(p); rss=float(np.sum(r.fun**2)); bic=n*np.log(rss/n)+k*np.log(n/6.)  # n/6: ~pixels per beam correlation penalty
    return r.x,rss,bic,n,model
for nm in 'SN':
  for meas in ['peak','HCNwin']:
    m=O[nm]['maps'][meas]; HW=O[nm]['HW']
    for R in (1.0,1.5):
      for nc in (0,1,2,3):
        try: x,rss,bic,n,model=fit2d(m,HW,R,nc)
        except Exception as e: print('fail',nm,meas,R,nc,e); continue
        s=x[3]; rec={'T0':float(x[0]),'s_core':float(s),'ff26':float(s**2/(s**2+.26**2)),'core_xy':[float(x[1]),float(x[2])],'ext':{'A':float(x[4]),'s':float(x[7]),'xy':[float(x[5]),float(x[6])]},'const':float(x[8]),'rms_resid_K':float(np.sqrt(rss/n)),'BIC':float(bic),'nclump':nc,
             'clumps':[[float(x[9+4*i+j]) for j in range(4)] for i in range(nc)]}
        out['2d'][f'{nm}_{meas}_R{R}_nc{nc}']=rec
        print(nm,meas,R,nc,'T0=%.1f s=%.2f ff26=%.2f ext(A=%.1f,s=%.1f) c=%.1f rms=%.2f BIC=%.0f'%(x[0],s,rec['ff26'],x[4],x[7],x[8],rec['rms_resid_K'],bic))
json.dump(out,open('/workspace/out/q32/models.json','w'),indent=1)
