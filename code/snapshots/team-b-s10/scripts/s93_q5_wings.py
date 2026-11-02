# Q5: wing vs core line ratios as a shock test.  Windows relative to the HCO+ line centroid of each aperture (Xb4d: 0.94x0.74", Xb47: 4.3x2.8"):
# blue wing -400..-150, core |v|<100, red wing +150..+400 km/s; baseline from 500<|v|<1100 (linear); sigma from the same off-line channels (lag-1 inflated) and from the field.
# Xb4d: SO 6_5-5_4 and HCO+/HCN 3-2 at Sblue, S, N (the CH3OH A-ladder lines lie within +-110 km/s of SO, i.e. outside the wing windows except for their own wings).
# Xb47: CS 7-6 vs HCN 4-3 red wing (H2CS 10-9 and 29SiO 8-7 lie at -55 and -86 km/s from CS, so only the CS RED wing is blend-free).
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
def win(ous,spw,reg,f0,vc,lo,hi,kind='ap'):
    nu,y=spec(ous,spw,reg,kind); v=vel(nu,f0)-vc; o=np.argsort(v); v,y=v[o],y[o]; dv=abs(np.median(np.diff(v)))
    fl=(abs(v)>500)&(abs(v)<1100); pf=np.polyfit(v[fl],y[fl],1); yb=y-np.polyval(pf,v); res=yb[fl]; rms=res.std(); infl=lag1(res)
    m=(v>=lo)&(v<=hi); n=m.sum(); W=yb[m].sum()*dv*1e-3; s=rms*np.sqrt(n*infl)*dv*1e-3
    return W,s
rows=[]
for ous,reg_l,lines in (('Xb4d',['Sblue','S','N','Bridge'],[('SO 6_5-5_4','spw31',251.82577),('HCO+ 3-2','spw25',267.55763),('HCN 3-2','spw27',265.88618)]),
                        ('Xb47',['N','S','Bridge'],[('CS 7-6','spw22',342.88285),('HCN 4-3','spw16',354.50547),('HCO+ 4-3','spw18',356.73424)])):
  for reg in reg_l:
    shp,vc,fw,snr=line_shape(ous,reg,mode='gauss',tname='ref')
    for (nm,spw,f0) in lines:
      for wn,(lo,hi) in dict(blue=(-400,-150),core=(-100,100),red=(150,400)).items():
        W,s=win(ous,spw,reg,f0,vc,lo,hi); rows.append(dict(ous=ous,region=reg,line=nm,window=wn,W=W,sig=s,vc=vc))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/wings.csv',index=False)
pd.set_option('display.width',200); 
p=df.pivot_table(index=['ous','region','line'],columns='window',values=['W','sig']); print(p.round(3).to_string())
def ratio(ous,reg,a,b,w):
    A=df[(df.ous==ous)&(df.region==reg)&(df.line==a)&(df.window==w)].iloc[0]; B=df[(df.ous==ous)&(df.region==reg)&(df.line==b)&(df.window==w)].iloc[0]
    r=A.W/B.W; e=np.sqrt((A.sig/B.W)**2+(A.W*B.sig/B.W**2)**2+0.05**2*r**2); return r,e,A.W/A.sig
print('\nXb47 CS7-6/HCN4-3 (Jy ratio; x1.069 for K) by window')
for reg in ('N','S','Bridge'):
    print(reg,{w:tuple(round(x,3) for x in ratio('Xb47',reg,'CS 7-6','HCN 4-3',w)) for w in ('blue','core','red')})
print('Xb4d SO/HCO+ , SO/HCN windows')
for reg in ('Sblue','S','N','Bridge'):
    print(reg,'SO/HCO+',{w:tuple(round(x,3) for x in ratio('Xb4d',reg,'SO 6_5-5_4','HCO+ 3-2',w)) for w in ('blue','core','red')},'SO/HCN',{w:tuple(round(x,3) for x in ratio('Xb4d',reg,'SO 6_5-5_4','HCN 3-2',w)) for w in ('blue','core','red')})
