# Q3 step 2: broad (|v-vc|>200 km/s) wing vs core at the nuclei: model-independent window fluxes and core+broad Gaussian fits, HCN and HCO+ 3-2.
# Spectra: Xb4d at matched 1.05" beam (and native peak pixel), Xb4b native 0.25" beam peak pixel. Positions: out/regions.csv (N, S).
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_lib import *
import pandas as pd
from scipy.optimize import curve_fit
g=lambda v,a,c,s: a*np.exp(-0.5*((v-c)/s)**2)
def fit(v,s,e,vmask):
    m=vmask; out={}
    f1=lambda v,a,c,sg: g(v,a,c,sg)
    f2=lambda v,a,c,sg,b,cb,sb: g(v,a,c,sg)+g(v,b,cb,sb)
    p1,_=curve_fit(f1,v[m],s[m],p0=[s[m].max(),VS-40,60],sigma=np.full(m.sum(),e),absolute_sigma=True)
    chi1=np.sum(((s[m]-f1(v[m],*p1))/e)**2)
    best=None
    for cb0 in (-100,0,100):
        try:
            p2,cv=curve_fit(f2,v[m],s[m],p0=[p1[0],p1[1],p1[2]*0.8,0.1*p1[0],p1[1]+cb0,150],sigma=np.full(m.sum(),e),absolute_sigma=True,
                bounds=([0,VS-300,10,0,VS-500,90],[np.inf,VS+300,200,np.inf,VS+500,400]))
            ch=np.sum(((s[m]-f2(v[m],*p2))/e)**2)
            if best is None or ch<best[0]: best=(ch,p2,cv)
        except Exception as ex: pass
    return p1,chi1,best,m.sum()
def run(ous,sp,fw,pos,tag):
    c=cube(sp,ous,crop=6. if ous=='Xb4b' else 14.,fw=fw); v=c['v']; rows=[]
    for nm,(x,y) in pos.items():
        # peak pixel of the core map within 0.5 beam of the nominal position
        m0,_=wsum(c,VS-200,VS+200); r=np.hypot(c['dx']-x,c['dy']-y); rr=0.6 if ous=='Xb4b' else 0.5
        j=np.argmax(np.where(r<rr,m0,-1e9)); iy,ix=np.unravel_index(j,m0.shape); s=c['a'][:,iy,ix]; px,py=c['dx'][iy,ix],c['dy'][iy,ix]
        fm=chan_noise(c); e=1.4826*np.median(np.abs(s[fm]-np.median(s[fm]))); base=np.median(s[fm]); s=s-base
        # windows relative to the line peak velocity (first moment of core +-150 around the peak channel)
        vp=v[np.argmax(np.convolve(s,np.ones(3)/3,'same'))]
        W={'blue':(vp-350,vp-200),'red':(vp+200,vp+350),'core':(vp-200,vp+200)}
        # window-sum noise per pixel from sliding line-free sums
        d={}
        for k,(lo,hi) in W.items():
            mm=(v>=lo)&(v<hi); n=mm.sum(); sg=wnoise(c,lo,hi)[0][iy,ix]; d[k]=(s[mm].sum()*c['dv'],sg,n)
        # vp-referenced windows must stay inside the frame
        p1,chi1,best,n=fit(v,s,e,(v>VS-900)&(v<VS+400)) if sp=='HCN' else fit(v,s,e,(v>VS-900)&(v<VS+330))
        row=dict(ous=ous,beam=tag,line=sp,region=nm,px=px,py=py,peak_mJy=s.max()*1e3,v_peak=vp,chan_rms_mJy=e*1e3,
          I_core=d['core'][0],e_core=d['core'][1],I_blue=d['blue'][0],e_blue=d['blue'][1],I_red=d['red'][0],e_red=d['red'][1],
          blue_over_core=d['blue'][0]/d['core'][0],e_blue_over_core=d['blue'][1]/d['core'][0],red_over_core=d['red'][0]/d['core'][0],e_red_over_core=d['red'][1]/d['core'][0],
          blue_snr=d['blue'][0]/d['blue'][1],red_snr=d['red'][0]/d['red'][1],
          core_sigma1=p1[2],chi2_1=chi1,ndat=n)
        if best:
            ch,p2,cv=best; er=np.sqrt(np.diag(cv)); row.update(chi2_2=ch,dchi2=chi1-ch,amp_core=p2[0],v_core=p2[1],sig_core=p2[2],amp_b=p2[3],e_amp_b=er[3],v_b=p2[4],sig_b=p2[5],e_sig_b=er[5],
              flux_b=p2[3]*p2[5]*np.sqrt(2*np.pi)*c['dv']/c['dv'],flux_core=p2[0]*p2[2]*np.sqrt(2*np.pi))
        rows.append(row)
    return rows
N,S=(0,2.55),(-0.15,-2.55)
rows=run('Xb4d','HCN',1.05,{'N':N,'S':S},'1.05"')+run('Xb4d','HCO+',1.05,{'N':N,'S':S},'1.05"')+run('Xb4d','HCN',None,{'N':N,'S':S},'native 0.95x0.75')+run('Xb4d','HCO+',None,{'N':N,'S':S},'native 0.95x0.75')
rows+=run('Xb4b','HCN',None,{'N':N,'S':S},'native 0.25"')+run('Xb4b','HCO+',None,{'N':N,'S':S},'native 0.25"')
T=pd.DataFrame(rows);T.to_csv('/workspace/out/q3/q3_nuclear_wings.csv',index=False)
pd.set_option('display.width',300);pd.set_option('display.max_columns',50)
print(T[['ous','beam','line','region','peak_mJy','chan_rms_mJy','v_peak','I_core','I_blue','e_blue','I_red','e_red','blue_over_core','red_over_core','blue_snr','red_snr']].round(3).to_string())
print(T[['ous','line','region','dchi2','chi2_1','chi2_2','ndat','v_core','sig_core','amp_b','v_b','sig_b','e_sig_b','flux_b','flux_core']].round(2).to_string())
