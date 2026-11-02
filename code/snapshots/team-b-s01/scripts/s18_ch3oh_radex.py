# CH3OH J_-1-J_0 E blend (165.05-165.37 GHz, Xb55 spw22): (1) unconstrained linear fit of J=1..5 with full covariance for several
# line-shape templates and fit windows (-> systematic spread of the blend integral); (2) non-LTE (RADEX, e-ch3oh.dat, LVG/slab) and LTE model grid
# fitted to the amplitude vector with the full covariance + blend systematic; (3) Tkin/n constraints as ranges.
import sys,subprocess,os,json,itertools; sys.path.insert(0,'/workspace/scripts')
from common import *
exec(open('/workspace/scripts/s12_line_table.py').read().split("rows=[]\nfor lab,f0,lst,note in LINES:")[0])
LJ={1:(165.050175,23.4,5.38627),2:(165.061130,28.0,8.94855),3:(165.099240,35.0,12.46907),4:(165.190475,44.3,15.92557),5:(165.369341,55.9,19.29095)}
def fitblend(reg,tmpl,win):
    T=make_template('Xb55',reg,'ap'); nu,y=spec('Xb55','spw22',reg,'ap')
    f0=LJ[2][0]; v=vel(nu,f0); o=np.argsort(v); v,y=v[o],y[o]; ok=np.isfinite(y); v,y=v[ok],y[ok]; dv=abs(np.median(np.diff(v)))
    vc=0.5*(T['hi']+T['lo']); half=0.5*(T['hi']-T['lo'])
    if tmpl=='hco+': tn=lambda vv: np.interp(vv,T['v'],T['t'],left=0,right=0)/T['t'].max()
    else:
        fw=float(tmpl); tn=lambda vv: np.exp(-4*np.log(2)*((vv-vc)/fw)**2)
    shift={J:C*(1-f*(1-vc/C)/f0) for J,(f,E,S_) in LJ.items()}
    free=abs(v-vc)<2500
    for J in LJ: free&=abs(v-shift[J])>half+100
    base=np.median(y[free]); yy=y-base; rms=np.std(yy[free])
    A=np.array([tn(v-shift[J]+vc) for J in LJ]).T; m=abs(v-vc)<win
    cov=np.linalg.inv(A[m].T@A[m])*rms**2*1.3**2     # x1.3: adjacent-channel correlation
    amp=cov@(A[m].T@yy[m])/rms**2/1.3**2
    norm=tn(np.arange(-1500,1500,dv)+vc).sum()*dv
    return amp*norm, cov*norm**2, rms
res={}
for reg in ('N','S'):
    fits=[]
    for tm in ('hco+','150','220','300'):
        for win in (600,900,1500):
            a,c,r=fitblend(reg,tm,win); fits.append((tm,win,a,c))
            print(reg,tm,win,'I_J=',np.round(a,2),'blend J1-3=%.2f+-%.2f'%(a[:3].sum(),np.sqrt(c[:3,:3].sum())),'J4=%.3f+-%.3f J5=%.3f+-%.3f'%(a[3],np.sqrt(c[3,3]),a[4],np.sqrt(c[4,4])))
    res[reg]=fits
# adopted: hco+ template, win 900 ; systematic of the blend sum from the spread
adopt={}
for reg in res:
    tm,win,a,c=[f for f in res[reg] if f[0]=='hco+' and f[1]==900][0]
    sums=np.array([f[2][:3].sum() for f in res[reg]]); 
    sysfrac=np.std(sums)/np.mean(sums)
    sig=np.sqrt(np.diag(c)); b=a[:3].sum(); be=np.sqrt(c[:3,:3].sum())
    ul4=max(a[3],0)+3*sig[3]; ul5=max(a[4],0)+3*sig[4]
    print(reg,'adopted blend %.2f+-%.2f stat, sys frac %.2f (range %.2f-%.2f); J4 3sig UL %.3f, J5 %.3f'%(b,be,sysfrac,sums.min(),sums.max(),ul4,ul5))
    adopt[reg]=dict(a=a.tolist(),cov=c.tolist(),blend=float(b),blend_stat=float(be),sysfrac=float(sysfrac),sums=sums.tolist(),ul4=float(ul4),ul5=float(ul5))
# ---- RADEX grid, fitted directly to the blend with ONE free scale (physical amplitude ratios), several templates
os.makedirs('/workspace/out/radex_ch3oh',exist_ok=True); os.chdir('/workspace/out/radex_ch3oh')
def radex(Tk,nH,N,dV=200.,geo='radex_lvg',tbg=2.73):
    inp=f"/opt/lamda/e-ch3oh.dat\nout.txt\n165.0 165.4\n{Tk}\n1\nH2\n{nH}\n{tbg}\n{N}\n{dV}\n0\n"
    subprocess.run([geo],input=inp,text=True,capture_output=True)
    out=np.zeros(5); tau=np.zeros(5)
    for l in open('out.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--':
            for J,(f0,E,S) in LJ.items():
                if abs(float(p[4])-f0)<1e-4: out[J-1]=float(p[8]); tau[J-1]=float(p[7])
    return out,tau
def lte(T):
    a=np.array([S*np.exp(-E/T) for J,(f,E,S) in LJ.items()]); return a/a.max()
def prep(reg,tmpl):
    T=make_template('Xb55',reg,'ap'); nu,y=spec('Xb55','spw22',reg,'ap')
    f0=LJ[2][0]; v=vel(nu,f0); o=np.argsort(v); v,y=v[o],y[o]; ok=np.isfinite(y); v,y=v[ok],y[ok]
    vc=0.5*(T['hi']+T['lo']); half=0.5*(T['hi']-T['lo'])
    tn=(lambda vv: np.interp(vv,T['v'],T['t'],left=0,right=0)/T['t'].max()) if tmpl=='hco+' else (lambda vv,fw=float(tmpl): np.exp(-4*np.log(2)*((vv-vc)/fw)**2))
    shift={J:C*(1-f*(1-vc/C)/f0) for J,(f,E,S_) in LJ.items()}
    free=abs(v-vc)<2500
    for J in LJ: free&=abs(v-shift[J])>half+100
    yy=y-np.median(y[free]); rms=np.std(yy[free])*1.3
    A=np.array([tn(v-shift[J]+vc) for J in LJ]).T; m=abs(v-vc)<900
    return A[m],yy[m],rms
def chi2_model(A,y,rms,amp):
    b=A@amp; s=(b@y)/(b@b) if b@b>0 else 0; s=max(s,0); return np.sum(((y-s*b)/rms)**2), s
Ts=[4,6,8,10,12,15,20,25,30,40,60,100,150]; nHs=[1e3,1e4,1e5,1e6,1e7]; Ns=[1e14,1e15,1e16]
out={}
for reg in ('N','S'):
    for tm in ('hco+','150','220'):
        A,y,rms=prep(reg,tm); chi0=np.sum((y/rms)**2)
        L=[]
        for T in Ts:
            c,s=chi2_model(A,y,rms,lte(T)); L.append(('LTE',T,None,None,c,s))
            for nH in nHs:
                for N in Ns:
                    amp,tau=radex(T,nH,N)
                    if amp.max()<=0: continue
                    c,s=chi2_model(A,y,rms,amp); L.append(('RADEX',T,nH,N,c,s))
        out[(reg,tm)]=(chi0,L)
        for mod in ('LTE','RADEX'):
            LL=[l for l in L if l[0]==mod]; best=min(LL,key=lambda l:l[4]); ok=[l for l in LL if l[4]<best[4]+4]   # ~2 sigma (dchi2<4)
            print(reg,tm,mod,'chi2_null=%.1f best chi2=%.1f (N_ch=%d) at T=%s n=%s N=%s ; T range (dchi2<4): %s-%s K'%(chi0,best[4],len(y),best[1],best[2],best[3],min(l[1] for l in ok),max(l[1] for l in ok)))
        # non-LTE only for sub-thermal n
        for nH in nHs:
            LL=[l for l in L if l[0]=='RADEX' and l[2]==nH]; ok=[l for l in LL if l[4]<best[4]+4] if False else None
        b0=min(l[4] for l in L)
        for nH in nHs:
            LL=[l for l in L if l[0]=='RADEX' and l[2]==nH]; ok=[l[1] for l in LL if l[4]<b0+4]
            print('    n=%g: allowed T(dchi2<4 wrt global best) = %s'%(nH,(min(ok),max(ok)) if ok else 'excluded'))
import pickle; pickle.dump({'adopt':adopt,'grid':{str(k):v for k,v in out.items()}},open('/workspace/out/s18_ch3oh_radex.pkl','wb'))
