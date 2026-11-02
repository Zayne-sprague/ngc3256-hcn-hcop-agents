# Q14a: chemistry ratios to HCO+ 3-2 (and HCN 3-2) as a function of beam, all cubes convolved channel-wise to the SAME elliptical beam per column.
# Beams: Xb4d only: 1.05" circ, 1.7" circ (Harada beam); 4.45x2.95 PA87 (Band-7 beam+) with Xb4d + Xb47; 6.3x3.8 PA82 with Xb4d + Xb4f + Xb47.
# Line amplitude: Gaussian at the HCO+ (v0,FWHM) of the same beam/region + linear baseline (matched-template fit); 3 sigma limits where S/N<3.
# Noise: residual rms of the fit window (robust) x1.3 for channel correlation (assumption, as in team-b p0019). Output out/q14/beamseq.csv
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from q4_lib import spectrum,kb,vel,robust_sigma,C
INF=1.3
BEAMS={'1.05':(1.05,1.05,0.),'1.7':(1.7,1.7,0.),'4.45x2.95':(4.45,2.95,87.),'6.3x3.8':(6.3,3.8,82.)}
# line: (rest GHz, {OUS:(spw)}), Xb4d spw as in q1 LINES
L={'HCO+ 3-2':(267.557626,{'Xb4d':'spw25','Xb4f':'spw16'}),'HCN 3-2':(265.886434,{'Xb4d':'spw27','Xb4f':'spw18'}),
   'CH3OH 2(0)-1(-1)E':(254.015377,{'Xb4d':'spw29','Xb4f':'spw20'}),'HC3N 28-27+CH2NH 4_04-3_03 (blend)':(254.699500,{'Xb4d':'spw29','Xb4f':'spw20'}),
   'HOC+ 3-2':(268.451,{'Xb4d':'spw25','Xb4f':'spw16'}),'HCO+ 4-3':(356.734223,{'Xb47':'spw18'}),'HCN 4-3':(354.505478,{'Xb47':'spw16'}),
   'CS 7-6':(342.883000,{'Xb47':'spw22'}),'SO 8_8-7_7':(344.310612,{'Xb47':'spw20'})}
REGS=['N','S','Bridge','Wclump','Soutflow']
def gauss(v,A,v0,fw,a,b): return A/(fw/2.3548*2.5066)*np.exp(-0.5*((v-v0)/(fw/2.3548))**2)+a+b*(v-v0)/100.
def bright(v,y,v00=2780.,fw0=150.):
    m=np.abs(v-v00)<700
    p,cv=curve_fit(gauss,v[m],y[m],p0=[y[m].max()*fw0*2.5/2.3548*0.9,v00,fw0,0,0],maxfev=20000)
    r=y[m]-gauss(v[m],*p); far=np.abs(v[m]-p[1])>1.5*p[2]; sg=robust_sigma(r[far]) if far.sum()>8 else robust_sigma(r)
    dv=abs(np.median(np.diff(v))); return p,np.sqrt(abs(cv[0,0])),sg
def template(v,y,v0,fw,half=450.):
    m=np.abs(v-v0)<half; x=v[m]-v0; g=np.exp(-0.5*(x/(fw/2.3548))**2)/(fw/2.3548*2.5066)
    M=np.vstack([g,np.ones_like(x),x/100.]).T; use=np.ones(m.sum(),bool)
    for _ in range(4):
        p=np.linalg.lstsq(M[use],y[m][use],rcond=None)[0]; r=y[m]-M@p; far=np.abs(x)>fw
        sg=1.4826*np.median(np.abs(r[far&use]-np.median(r[far&use]))); use=np.abs(r)<4*sg
    cv=np.linalg.pinv(M[use].T@M[use])*sg**2; return p[0],np.sqrt(cv[0,0]),sg,x,r
rows=[]
for bn,tb in BEAMS.items():
    for reg in REGS:
        for ous in ('Xb4d','Xb4f','Xb47'):
            if bn in('1.05','1.7') and ous!='Xb4d': continue
            if bn=='4.45x2.95' and ous=='Xb4f': continue
            ref={}
            for ln,(f0,sp) in L.items():
                if ous not in sp: continue
                try: fr,s,beam,nb=spectrum(ous,sp[ous],reg,tb=tb)
                except ValueError as e: continue
                v=vel(fr,f0); y=s*kb(fr.mean(),*tb[:2]) # K
                if ln in('HCO+ 3-2','HCN 3-2','HCO+ 4-3','HCN 4-3'):
                    p,eA,sg=bright(v,y); I=p[0]; e=np.hypot(eA,0)*INF; rows.append(dict(beam=bn,ous=ous,region=reg,line=ln,I=I,eI=e,v0=p[1],fwhm=p[2],sig_ch=sg,kind='fit'))
                    ref[ln]=p
                else:
                    key='HCO+ 3-2' if '3-2' in ln or ln.startswith(('CH3','HC3N')) else 'HCO+ 4-3'
                    if key not in ref: continue
                    v0,fw=ref[key][1],ref[key][2]
                    A,eA,sg,x,r=template(v,y,v0,fw); rows.append(dict(beam=bn,ous=ous,region=reg,line=ln,I=A,eI=eA*INF,v0=v0,fwhm=fw,sig_ch=sg,kind='template'))
df=pd.DataFrame(rows); df['snr']=df.I/df.eI
# ratios to HCO+ 3-2 / HCN 3-2 / HCO+ 4-3 / HCN 4-3 (same beam, OUS)
out=[]
for (bn,ous,reg),g in df.groupby(['beam','ous','region']):
    g=g.set_index('line')
    for ln in g.index:
        for ref in ('HCO+ 3-2','HCN 3-2','HCO+ 4-3','HCN 4-3'):
            if ln==ref or ref not in g.index or ln in('HCO+ 3-2','HCN 3-2','HCO+ 4-3','HCN 4-3'): continue
            if ('4-3' in ref)!=(ous=='Xb47'): continue
            r=g.loc[ln,'I']/g.loc[ref,'I']; e=g.loc[ln,'eI']/g.loc[ref,'I']
            out.append(dict(beam=bn,ous=ous,region=reg,num=ln,den=ref,ratio=r,eratio=e,snr=r/e,lim3=max(3*e,r+3*e) if r/e<3 else np.nan))
R=pd.DataFrame(out); df.to_csv('/workspace/out/q14/beamseq_lines.csv',index=False); R.to_csv('/workspace/out/q14/beamseq.csv',index=False)
pd.set_option('display.width',220); pd.set_option('display.max_rows',500)
print(R[R.den.isin(['HCO+ 3-2','HCO+ 4-3'])].round(4).to_string())
