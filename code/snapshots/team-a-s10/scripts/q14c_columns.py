# Q14c: (1) Convert our 254 GHz CH3OH/HCO+ 3-2 into the Harada+2018 96.74 GHz (2k-1k, A+E sum) CH3OH/HCO+ 3-2 under thin LTE (LAMDA A, g, E);
#       (2) beam-averaged column ratios N(X)/N(HCO+) (thin LTE lower limit on HCO+ -> UPPER limit on X/HCO+ if HCO+ is thick; tau correction listed);
#       (3) significance of region contrasts from out/q14/beamseq.csv; (4) CS 7-6/HCO+ 4-3 double ratio N/S.
# Output out/q14/harada_convert.csv, columns.csv, contrasts.csv
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
K1=1.4388; h=6.62607e-27; c=2.998e10; kB=1.380649e-16
def lam(p):
    L=open(p).read().split('\n'); i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]; n=int(L[i+1])
    lv=np.array([[float(a.split()[1])*K1,float(a.split()[2])] for a in L[i+3:i+3+n]])
    j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]; m=int(L[j+1])
    tr=np.array([[float(x) for x in l.split()[:6]] for l in L[j+3:j+3+m]])  # id up low A f Eu
    return lv,tr
def bose(nu,T): x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; return 1-(np.exp(x)-1)/(np.exp(xb)-1)
class Mol:
    def __init__(s,f): s.lv,s.tr=lam(f)
    def Q(s,T): return float(np.sum(s.lv[:,1]*np.exp(-s.lv[:,0]/T)))
    def line(s,fGHz,tol=2e-3):
        r=s.tr[np.abs(s.tr[:,4]-fGHz)<tol][0]; up=int(r[1])-1; return dict(A=r[3],f=r[4],Eu=s.lv[up,0],gu=s.lv[up,1])
    def W(s,fGHz,T):   # thin K km/s per unit total column (RJ, Bose-corrected) at Tex=T
        l=s.line(fGHz); Nu=s.lv[int(s.tr[np.abs(s.tr[:,4]-fGHz)<2e-3][0][1])-1,1]*np.exp(-l['Eu']/T)/s.Q(T)
        return h*c**3*l['A']*Nu/(8*np.pi*kB*(l['f']*1e9)**2)*bose(l['f'],T)/1e5
E_=Mol('/opt/lamda/e-ch3oh.dat'); A_=Mol('/opt/lamda/ch3oh_a.dat'); HCO=Mol('/opt/lamda/hco+.dat'); CS=Mol('/opt/lamda/cs.dat')
# A+E with equal ortho/para (E:A=1:1 in total column; each ladder normalised on its own Q): W_AE(f)= 0.5*(W_E)+0.5*(W_A)
def W_AE(f,T,kind): return 0.5*(E_.W(f,T) if kind=='E' else A_.W(f,T))
rows=[]
for T in (8.,10.,15.,20.,30.,50.):
    w254=W_AE(254.015377,T,'E')
    w96=sum(W_AE(f,T,'E') for f in (96.739358,96.744545,96.755501))+0.5*A_.W(96.741375,T) if True else 0
    rows.append(dict(T=T,W254_per_Ntot=w254,W96p74_per_Ntot=w96,ratio_96p74_over_254=w96/w254))
H=pd.DataFrame(rows); H.to_csv('/workspace/out/q14/harada_convert.csv',index=False); print(H.to_string())
B=pd.read_csv('/workspace/out/q14/beamseq.csv')
def get(beam,ous,reg,num,den='HCO+ 3-2'):
    r=B[(B.beam==beam)&(B.ous==ous)&(B.region==reg)&(B.num==num)&(B.den==den)]; return (float(r.ratio.iloc[0]),float(r.eratio.iloc[0])) if len(r) else (np.nan,np.nan)
HAR={'N':(0.0337,0.0072),'S':(0.1424,0.0152)}   # Harada+2018 CH3OH 2k-1k/HCO+ 3-2 in K km/s ratio (q4c/harada_ratios.csv)
out=[]
for reg in ('N','S'):
    r,e=get('1.7','Xb4d',reg,'CH3OH 2(0)-1(-1)E')
    for T in (8.,10.,15.,20.,30.):
        f=float(H[H['T']==T].ratio_96p74_over_254.iloc[0])
        out.append(dict(region=reg,T=T,ours_254_over_HCOp32_1p7=r,e=e,pred_96p74_over_HCOp32=r*f,e_pred=e*f,harada=HAR[reg][0],e_harada=HAR[reg][1],pred_over_harada=r*f/HAR[reg][0],T_needed_note=''))
for reg in ('N','S'):    # same, with the Xb4f 6.3x3.8 beam (Xb4f recovers extended CH3OH; Xb4d resolves it out)
    r,e=get('6.3x3.8','Xb4f',reg,'CH3OH 2(0)-1(-1)E')
    for T in (8.,10.,15.,20.,30.):
        f=float(H[H['T']==T].ratio_96p74_over_254.iloc[0])
        out.append(dict(region=reg,T=T,ours_254_over_HCOp32_1p7=r,e=e,pred_96p74_over_HCOp32=r*f,e_pred=e*f,harada=HAR[reg][0],e_harada=HAR[reg][1],pred_over_harada=r*f/HAR[reg][0],T_needed_note='Xb4f 6.3x3.8 beam (not 1.7")'))
O=pd.DataFrame(out); O.to_csv('/workspace/out/q14/harada_compare.csv',index=False); print(O.round(4).to_string())
# contrasts
def ratio_err(a,b): return a[0]/b[0],np.hypot(a[1]/a[0],b[1]/b[0])*a[0]/b[0]
C=[]
for beam,ous in (('1.7','Xb4d'),('4.45x2.95','Xb4d'),('6.3x3.8','Xb4d'),('6.3x3.8','Xb4f')):
    n=get(beam,ous,'N','CH3OH 2(0)-1(-1)E'); s=get(beam,ous,'S','CH3OH 2(0)-1(-1)E'); b=get(beam,ous,'Bridge','CH3OH 2(0)-1(-1)E'); w=get(beam,ous,'Wclump','CH3OH 2(0)-1(-1)E')
    d=s[0]-n[0]; ed=np.hypot(s[1],n[1])
    C.append(dict(quantity='CH3OH254/HCO+3-2',beam=beam,ous=ous,N=n[0],eN=n[1],S=s[0],eS=s[1],Bridge=b[0],eB=b[1],Wclump=w[0],eW=w[1],S_minus_N=d,sig_S_gt_N=d/ed,S_over_N=s[0]/n[0] if n[0]>0 else np.nan))
for beam in ('4.45x2.95','6.3x3.8'):
    n=get(beam,'Xb47','N','CS 7-6','HCO+ 4-3'); s=get(beam,'Xb47','S','CS 7-6','HCO+ 4-3'); b=get(beam,'Xb47','Bridge','CS 7-6','HCO+ 4-3'); w=get(beam,'Xb47','Wclump','CS 7-6','HCO+ 4-3')
    rr=ratio_err(n,s); C.append(dict(quantity='CS7-6/HCO+4-3',beam=beam,ous='Xb47',N=n[0],eN=n[1],S=s[0],eS=s[1],Bridge=b[0],eB=b[1],Wclump=w[0],eW=w[1],N_over_S=rr[0],eN_over_S=rr[1],sig_N_gt_S=(n[0]-s[0])/np.hypot(n[1],s[1])))
    n=get(beam,'Xb47','N','CS 7-6','HCN 4-3'); s=get(beam,'Xb47','S','CS 7-6','HCN 4-3'); rr=ratio_err(n,s)
    C.append(dict(quantity='CS7-6/HCN4-3',beam=beam,ous='Xb47',N=n[0],eN=n[1],S=s[0],eS=s[1],N_over_S=rr[0],eN_over_S=rr[1],sig_N_gt_S=(n[0]-s[0])/np.hypot(n[1],s[1])))
C=pd.DataFrame(C); C.to_csv('/workspace/out/q14/contrasts.csv',index=False); pd.set_option('display.width',250); print(C.round(4).to_string())
# columns: thin-LTE beam-averaged N(X)/N(HCO+) at the common Xb4f 6.3x3.8 beam, Tex assumed (same for both)
Ld=pd.read_csv('/workspace/out/q14/beamseq_lines.csv')
def I(beam,ous,reg,line):
    r=Ld[(Ld.beam==beam)&(Ld.ous==ous)&(Ld.region==reg)&(Ld.line==line)]; return float(r.I.iloc[0]),float(r.eI.iloc[0])
rows=[]
for reg in ('N','S','Bridge','Wclump'):
    ih=I('6.3x3.8','Xb4f',reg,'HCO+ 3-2'); ic=I('6.3x3.8','Xb4f',reg,'CH3OH 2(0)-1(-1)E')
    i43=I('6.3x3.8','Xb47',reg,'HCO+ 4-3'); ics=I('6.3x3.8','Xb47',reg,'CS 7-6')
    for T in (10.,20.,30.,50.):
        NH=ih[0]/HCO.W(267.557626,T)           # per unit-column normalisation: N = I / (W per unit N)
        Nm=ic[0]/(0.5*E_.W(254.015377,T))      # E ladder carries half of total
        NH43=i43[0]/HCO.W(356.734223,T); Ncs=ics[0]/CS.W(342.883,T)
        rows.append(dict(region=reg,Tex=T,N_CH3OH_over_N_HCOp_thin=Nm/NH,e=ic[1]/ic[0]*Nm/NH if ic[0]>0 else np.nan,N_CS_over_N_HCOp_thin_from43=Ncs/NH43,N_HCOp_thin_32=NH,N_HCOp_thin_43=NH43,N_CH3OH=Nm,N_CS=Ncs,
            tau_corr_note='HCO+ is thick (tau32~6-30, Q2/q4c): true N(X)/N(HCO+) is the thin value times tau/(1-exp(-tau))^-1 -> divide by 6-30'))
K=pd.DataFrame(rows); K.to_csv('/workspace/out/q14/columns.csv',index=False); print(K.drop(columns='tau_corr_note').round(4).to_string())
