# Q2 step e: nuclear line profiles are double-peaked (s75a), so decompose as 2 narrow + 1 broad Gaussian (broad sigma 60-300 km/s) and compare with 2 narrow only.
# Spectra: Xb4d 1.05" r=1" aperture (beam-averaged Jy) and Xb4b 0.5" peak pixel (Jy/beam). Noise per channel from line-free channels (-900..-450). Output out/q2/s75d.txt/json
import sys, json; sys.path.insert(0,'/workspace/scripts')
exec(open('/workspace/scripts/s75a_pv_spectra.py').read().split("def line_pos")[0].replace("SET={'Xb4d_1.05'","SET0={'Xb4d_1.05'").replace("SET0={'Xb4d_1.05':prep('Xb4d',1.05,24.0),'Xb4d_nat':prep('Xb4d',None,24.0),'Xb4b_0.5':prep('Xb4b',0.5,10.0)}","SET={'Xb4d_1.05':prep('Xb4d',1.05,24.0),'Xb4b_0.5':prep('Xb4b',0.5,10.0)}"))
from scipy.optimize import curve_fit
OUT='/workspace/out/q2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
def g1(v,a,c,s): return a*np.exp(-0.5*((v-c)/s)**2)
def m2(v,a1,c1,s1,a2,c2,s2): return g1(v,a1,c1,s1)+g1(v,a2,c2,s2)
def m3(v,a1,c1,s1,a2,c2,s2,a3,c3,s3): return m2(v,a1,c1,s1,a2,c2,s2)+g1(v,a3,c3,s3)
def peak_spec(R,pos,search=0.5):
    px,py=R['pxy'](*pos); n=R['n']; Y,X=np.mgrid[:n,:n]; m=np.hypot(X-px,Y-py)*R['ps']<=search
    M0=R['c'][np.abs(R['v'])<350].sum(0); i=np.argmax(np.where(m,M0,-1e9)); y,x=np.unravel_index(i,M0.shape); return R['c'][:,y,x]
res={}
for rn,pos in (('S',REGIONS['S'][:2]),('N',REGIONS['N'][:2])):
  for k,kind in (('Xb4d_1.05','ap'),('Xb4b_0.5','peak')):
    for mol in ('HCN','HCO+'):
        R=SET[k][mol]; v=R['v']; sp=aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],R['pxy'](*pos),1.0)) if kind=='ap' else peak_spec(R,pos)
        e=sp[(v>-900)&(v<-450)].std(); m=(v>-800)&(v<650); vv,y=v[m],sp[m]; a0=y.max()
        best2=None; best3=None
        for c1,c2 in ((-60,50),(-50,40),(-70,60)):
            try:
                p,cv=curve_fit(m2,vv,y,p0=[a0,c1,40,a0,c2,40],sigma=np.full(len(vv),e),absolute_sigma=True,bounds=([0,-150,15,0,-150,15],[3*a0,150,120,3*a0,150,120]),maxfev=20000)
                ch=np.sum(((y-m2(vv,*p))/e)**2)
                if best2 is None or ch<best2[0]: best2=(ch,p)
            except Exception: pass
        for cb in (-150,0,150):
          for sb in (80,150,250):
            try:
                p0=list(best2[1])+[0.1*a0,cb,sb]
                p,cv=curve_fit(m3,vv,y,p0=p0,sigma=np.full(len(vv),e),absolute_sigma=True,bounds=([0,-150,15,0,-150,15,0,-400,60],[3*a0,150,120,3*a0,150,120,a0,400,300]),maxfev=40000)
                ch=np.sum(((y-m3(vv,*p))/e)**2)
                if best3 is None or ch<best3[0]: best3=(ch,p,np.sqrt(np.diag(cv)))
            except Exception: pass
        ch2,p2=best2; ch3,p3,e3=best3
        Fb=np.sqrt(2*np.pi)*p3[6]*p3[8]; eFb=Fb*np.hypot(e3[6]/max(p3[6],1e-12),e3[8]/p3[8])
        Fnarrow=sum(np.sqrt(2*np.pi)*p3[i]*p3[i+2] for i in (0,3))
        tag=f'{rn}|{k}|{kind}|{mol}'
        res[tag]=dict(rms=float(e),chi2_2=float(ch2),chi2_3=float(ch3),narrow=p2.tolist(),full=p3.tolist(),err=e3.tolist(),Fbroad=float(Fb),Fnarrow=float(Fnarrow))
        P(f'{tag}: rms {e*1e3:.2f} mJy/ch; 2 narrow (c {p2[1]:.0f},{p2[4]:.0f}; FWHM {2.355*p2[2]:.0f},{2.355*p2[5]:.0f}) chi2={ch2:.0f}/{len(vv)-6}; +broad: A={p3[6]*1e3:.2f}+-{e3[6]*1e3:.2f} mJy c={p3[7]:.0f}+-{e3[7]:.0f} FWHM={2.355*p3[8]:.0f}+-{2.355*e3[8]:.0f} chi2={ch3:.0f}/{len(vv)-9}, dchi2={ch2-ch3:.0f}; broad flux {Fb:.3f}+-{eFb:.3f} (Jy or Jy/b)*km/s = {100*Fb/(Fb+Fnarrow):.1f}% of fitted line')
json.dump(res,open(OUT+'/s75d.json','w'),indent=1,default=float); open(OUT+'/s75d.txt','w').write('\n'.join(txt))
