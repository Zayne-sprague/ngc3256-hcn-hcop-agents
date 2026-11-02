# Q14i (fix 3,4): thin-LTE columns and abundance proxies at the common 6.3x3.8" beam.
#  (1) rigid-rotor LTE line strengths (validated against LAMDA CS 7-6), CMB contrast included;
#  (2) N(SiO)/N(CS) limits from 29SiO 8-7 (28Si/29Si=20, range 15-30) and 30SiO 6-5 (28Si/30Si=30, range 25-35) vs Harada's SiO 2-1 / CS 2-1 (Table 13, 1.7") ;
#  (3) dust N(H2) at the matched beam from the Band-6 (258 GHz, Xb4f) and Band-7 (346 GHz, Xb47) continuum, free-free/synchrotron removed with the 258-346 GHz two-band split (MC), Td 20-50 K, kappa(346)=0.02 cm2/g(gas) (beta 1.5..2);
#  (4) N(CH3OH)/N(H2), N(CH3OH)/N(H13CN)(limit), N(CS)/N(H2); S/N contrast vs Tex.  Output out/q14e/abund_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
hP=6.62607e-27; cc=2.998e10; kB=1.380649e-16
def rr_line(B,mu,J):   # linear rotor J->J-1: nu GHz, Eu K, A s-1, gu
    nu=2*B*J; Eu=0.047992*B*J*(J+1); A=64*np.pi**4*(nu*1e9)**3*(mu*1e-18)**2/(3*hP*cc**3)*J/(2*J+1); return nu,Eu,A,2*J+1,B
def Wrr(B,mu,J,T,cmb=True):
    nu,Eu,A,gu,_=rr_line(B,mu,J); Q=kB*T/(hP*B*1e9)+1/3.; Nu=gu*np.exp(-Eu/T)/Q
    w=hP*cc**3*A*Nu/(8*np.pi*kB*(nu*1e9)**2)/1e5
    if cmb: x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; w*=1-(np.exp(x)-1)/(np.exp(xb)-1)
    return w
MOL={'CS76':(24.495,1.958,7),'CS21':(24.495,1.958,2),'SiO21':(21.7116,3.098,2),'29SiO87':(21.4363,3.098,8),'30SiO65':(21.1847,3.098,6),'H13CN43':(43.1699,2.985,4)}
W=lambda k,T:Wrr(*MOL[k],T)
# validation vs LAMDA CS 7-6 (q14c Mol class)
sys.argv=['x']
import importlib.util
def lamda_W(T):
    K1=1.4388; L=open('/opt/lamda/cs.dat').read().split('\n'); i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]; n=int(L[i+1])
    lv=np.array([[float(a.split()[1])*K1,float(a.split()[2])] for a in L[i+3:i+3+n]]); j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]
    tr=np.array([[float(x) for x in l.split()[:6]] for l in L[j+3:j+3+int(L[j+1])]]); r=tr[np.abs(tr[:,4]-342.883)<0.01][0]
    Q=np.sum(lv[:,1]*np.exp(-lv[:,0]/T)); up=int(r[1])-1; Nu=lv[up,1]*np.exp(-lv[up,0]/T)/Q; nu=r[4]
    w=hP*cc**3*r[3]*Nu/(8*np.pi*kB*(nu*1e9)**2)/1e5; x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; return w*(1-(np.exp(x)-1)/(np.exp(xb)-1))
print('validation CS 7-6 W(rigid)/W(LAMDA):',[round(W('CS76',T)/lamda_W(T),3) for T in (10,20,50)])

import q14f_lte as q
D=pd.read_csv('/workspace/out/q14e/deblend_all.csv'); RG=['N','S','Bridge','Wclump','Soutflow']
def dget(cl,ous,beam,reg,comp,model='all',tmpl='HCO+'):
    r=D[(D.cluster==cl)&(D.ous==ous)&(D.beam==beam)&(D.region==reg)&(D.comp==comp)&(D.model==model)&(D.template==tmpl)]
    return float(r.I.iloc[0]),float(r.eI.iloc[0])
def dspread(cl,ous,beam,reg,comp,model='all'):
    r=D[(D.cluster==cl)&(D.ous==ous)&(D.beam==beam)&(D.region==reg)&(D.comp==comp)&(D.model==model)]; return r.I.min(),r.I.max()
# ---- (A) CH3OH deblended intensities and (B) SiO limits
rows=[]
for reg in RG:
    for ous,beam in (('Xb4d','1.7'),('Xb4d','4.45x2.95'),('Xb4d','6.3x3.8'),('Xb4f','6.3x3.8')):
        I,e=dget('A',ous,beam,reg,'CH3OH'); lo,hi=dspread('A',ous,beam,reg,'CH3OH'); rows.append(dict(kind='CH3OH',ous=ous,beam=beam,region=reg,I=I,eI=e,tmpl_min=lo,tmpl_max=hi))
        I,e=dget('A',ous,beam,reg,'30SiO'); rows.append(dict(kind='30SiO 6-5',ous=ous,beam=beam,region=reg,I=I,eI=e,lim3=max(I,0)+3*e))
    for beam in ('4.45x2.95','6.3x3.8'):
        I,e=dget('B','Xb47',beam,reg,'CS'); rows.append(dict(kind='CS 7-6',ous='Xb47',beam=beam,region=reg,I=I,eI=e))
        I,e=dget('B','Xb47',beam,reg,'29SiO','all'); rows.append(dict(kind='29SiO 8-7 (joint CS+H2CS+29SiO)',ous='Xb47',beam=beam,region=reg,I=I,eI=e,lim3=max(I,0)+3*e))
        I,e=dget('B','Xb47',beam,reg,'29SiO','noH2CS'); rows.append(dict(kind='29SiO 8-7 (H2CS=0)',ous='Xb47',beam=beam,region=reg,I=I,eI=e,lim3=max(I,0)+3*e))
        I,e=dget('B','Xb47',beam,reg,'H2CS','noSiO'); rows.append(dict(kind='H2CS 10_0,10 (29SiO=0)',ous='Xb47',beam=beam,region=reg,I=I,eI=e,lim3=max(I,0)+3*e))
I_=pd.DataFrame(rows); I_.to_csv('/workspace/out/q14e/abund_intensities.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',500); print(I_[I_.beam=='6.3x3.8'].round(3).to_string())
# SiO/CS columns, Tex common; isotope ratios
out=[]
def get(kind,ous,beam,reg): r=I_[(I_.kind==kind)&(I_.ous==ous)&(I_.beam==beam)&(I_.region==reg)]; return r.iloc[0]
for reg in RG:
    for beam in ('4.45x2.95','6.3x3.8'):
        cs=get('CS 7-6','Xb47',beam,reg)
        for T in (10.,20.,30.,50.):
            NCS=cs.I/W('CS76',T)
            for kind,key,iso,isr in (('29SiO 8-7 (joint CS+H2CS+29SiO)','29SiO87',20.,(15,30)),('29SiO 8-7 (H2CS=0)','29SiO87',20.,(15,30))):
                l=get(kind,'Xb47',beam,reg).lim3; Nl=l/W(key,T)
                out.append(dict(region=reg,beam=beam,Tex=T,line=kind,N_CS=NCS,N_SiO_lim3=iso*Nl,N_SiO_lim3_iso15=isr[0]*Nl,N_SiO_lim3_iso30=isr[1]*Nl,SiO_over_CS_lim3=iso*Nl/NCS,SiO_over_CS_lim3_iso30=isr[1]*Nl/NCS))
O=pd.DataFrame(out); O.to_csv('/workspace/out/q14e/abund_sio_cs.csv',index=False); print(O[O.beam=='6.3x3.8'].round(4).to_string())
# 30SiO 6-5 (Xb4f 6.3x3.8, deblended) vs CS 7-6 (Xb47 6.3x3.8): different arrays/tunings, same beam
o30=[]
for reg in RG:
    cs=get('CS 7-6','Xb47','6.3x3.8',reg)
    for ous in ('Xb4f','Xb4d'):
        s30=get('30SiO 6-5',ous,'6.3x3.8',reg)
        for T in (10.,20.,30.,50.):
            for iso in (25.,30.,35.):
                o30.append(dict(region=reg,ous=ous,Tex=T,iso=iso,SiO_over_CS_lim3=iso*(s30.lim3/W('30SiO65',T))/(cs.I/W('CS76',T)),I30_lim3=s30.lim3,I30=s30.I,eI30=s30.eI))
O30=pd.DataFrame(o30); O30.to_csv('/workspace/out/q14e/abund_sio30_cs.csv',index=False); print(O30[O30.iso==30].round(4).to_string())
# Harada SiO/CS (2-1 lines, K km/s at 1.7", Table 13): N 6.41+-0.95 / 48.99+-0.64 ; S 2.75+-0.95 / 22.92+-0.64 ; TSW 3.95+-0.97/15.6+-0.66 ; OS 3.83+-0.97/10.34+-0.66
H=[]
for nm,si,esi,cs,ecs in (('N',6.41,0.95,48.99,0.64),('S',2.75,0.95,22.92,0.64),('C',5.26,0.95,29.89,0.66),('TSW',3.95,0.97,15.6,0.66),('OS',3.83,0.97,10.34,0.66),('SW',4.69,0.95,31.55,0.64)):
    for T in (10.,20.,30.):
        r=(si/W('SiO21',T))/(cs/W('CS21',T)); H.append(dict(pos=nm,Tex=T,SiO_over_CS=r,e=r*np.hypot(esi/si,ecs/cs),I_SiO=si,I_CS=cs))
H=pd.DataFrame(H); H.to_csv('/workspace/out/q14e/harada_sio_cs.csv',index=False); print(H.round(4).to_string())
