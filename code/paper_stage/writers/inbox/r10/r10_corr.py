"""r10: missing-flux bounds on 4-3 from Xb47 (4.5" beam, out/q60/imgplane.csv: Xb47 vs 00493/00223 same beam), applied to ladder.csv (r10_ladder.py).
(P) pedestal: add dT=I_Xb47-I_hr(4.5") (K km/s) to the 1"/1.5" integrals (uniform extended flux, same T_b at all beams);
(F) fractional: divide by ratio I_hr/I_Xb47 (all lost flux inside the core; maximal correction). Truth lies between. 3-2 uncorrected (Xb4b; Xb4d at 1.5" in brackets)."""
import pandas as pd, numpy as np
q=pd.read_csv('/workspace/out/q60/imgplane.csv'); q=q[q.beam==4.5]
g=lambda l,s,r: q[(q.line==l)&(q.src==s)&(q.reg==r)].I.iloc[0]
S={('HCN','N'):('Xb47','X1df1'),('HCN','S'):('Xb47','X1df1'),('HCO+','N'):('Xb47','X858'),('HCO+','S'):('Xb47','X858')}
L=pd.read_csv('ladder.csv'); L=L[L.ous32=='Xb4b']; out=[]
for _,r in L.iterrows():
    xb,hr=g(r.line,'Xb47',r.reg),g(r.line,S[(r.line,r.reg)][1],r.reg)
    out.append(dict(beam=r.beam,reg=r.reg,line=r.line,rec=hr/xb,ped=xb-hr,I43=r.I43,I43_P=r.I43+xb-hr,I43_F=r.I43*xb/hr,I32=r.I32))
D=pd.DataFrame(out)
for k in ('I43','I43_P','I43_F'): D['T'+k[3:]+'_43/32']=D[k]/D.I32
D.round(3).to_csv('ladder_corr.csv',index=False); pd.set_option('display.width',200); print(D.round(3))
print('\nHCN/HCO+ and N/S contrast')
for bm in (1.0,1.5):
  d=D[D.beam==bm].set_index(['reg','line']); row={}
  for k in ('I43','I43_P','I43_F','I32'):
    rN=d.loc[('N','HCN'),k]/d.loc[('N','HCO+'),k]; rS=d.loc[('S','HCN'),k]/d.loc[('S','HCO+'),k]
    row[k]=(round(rN,3),round(rS,3),round(rN/rS,2))
  print(bm,row, 'contrast ratio 4-3/3-2: raw %.2f P %.2f F %.2f'%tuple(row[k][2]/row['I32'][2] for k in ('I43','I43_P','I43_F')))
