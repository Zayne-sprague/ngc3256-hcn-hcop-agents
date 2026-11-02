# Q4d: RADEX (LVG, dv=150 km/s FWHM, Tbg=2.73 K; LAMDA hcn.dat, hco+.dat, h13cn@xpol.dat) consequences of the two readings of the N-core R_K(3-2)~1 (Xb4b 0.25", s43) vs the outer R_K~0.4:
# (a) abundance reading: optically thin-ish gas, N(HCN)/N(HCO+)=A>1;  (b) opacity reading: both lines thick and at similar Tex -> R_K->1 independent of A.
# For each (Tkin, n, N(HCN), A) tabulate R_K(3-2), R_K(4-3), tau(HCN 3-2), tau(HCO+ 3-2), and the predicted I(HCN 4-3)/I(H13CN 4-3) for 12C/13C = 30, 50, 100 (equal beam filling; H13CN column = N(HCN)/X).
# Data constraints they are compared with: R_K = R_flux * (nu_HCO+/nu_HCN)^2 = 1.0127 R_flux; N-core R(0.25")=0.984+-0.015 (s43); H13CN limits from s41; vib limit from s42.
import subprocess,os,json,itertools,numpy as np,pandas as pd
os.makedirs('/workspace/out/radex_q4',exist_ok=True); os.chdir('/workspace/out/radex_q4')
LO={'hcn':{2:177.2611,3:265.8864,4:354.5055},'hco+':{2:178.3751,3:267.5576,4:356.7342},'h13cn@xpol':{3:259.0118,4:345.3398}}
def radex(mol,Tk,nH,N,dV=150.):
    inp=f"/opt/lamda/{mol}.dat\nout_{mol}.txt\n170 370\n{Tk}\n1\nH2\n{nH}\n2.73\n{N}\n{dV}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    I={};tau={};tex={}
    for l in open(f'out_{mol}.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--':
            for J,f in LO[mol].items():
                if abs(float(p[4])-f)<2e-3: I[J]=float(p[11]); tau[J]=float(p[7]); tex[J]=float(p[6])
    return I,tau,tex
rows=[]
Ts=[30,60,120]; ns=[1e5,1e6,1e7]; NH=10**np.arange(13,17.51,0.5); As=[1,3,10]; Xs=[30,50,100]
for T in Ts:
    for n in ns:
        H={N:radex('hcn',T,n,N) for N in NH}
        for N in NH:
            for A in As:
                Io,to,_=radex('hco+',T,n,N/A)
                for X in Xs:
                    I13,t13,_=radex('h13cn@xpol',T,n,N/X)
                    if len(H[N][0])<3 or len(Io)<3 or 4 not in I13 or min(Io[3],Io[4],Io[2],I13[4])<=0: continue
                    rows.append(dict(Tk=T,n=n,logN_HCN=np.log10(N),A=A,X=X,R32=H[N][0][3]/Io[3],R43=H[N][0][4]/Io[4],R21=H[N][0][2]/Io[2],tau32_HCN=H[N][1][3],tau32_HCO=to[3],tau43_HCN=H[N][1][4],
                                     HCN_H13CN_43=H[N][0][4]/I13[4],Tpk32_HCN=H[N][0][3]/150./1.0645))   # K (approx peak=flux/(1.0645*FWHM))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s44_radex_opacity.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',400)
print('R_K(3-2) vs N(HCN), A=1,3,10 (X=50) at n=1e6, T=60:'); print(df[(df.Tk==60)&(df.n==1e6)&(df.X==50)].pivot(index='logN_HCN',columns='A',values='R32').round(2).to_string())
print('tau32_HCN at the same cells:'); print(df[(df.Tk==60)&(df.n==1e6)&(df.X==50)&(df.A==1)][['logN_HCN','tau32_HCN','tau32_HCO','HCN_H13CN_43']].round(2).to_string())
core=df[(df.R32>0.90)&(df.R32<1.07)&(df.X==50)]
print('\nCells with R_K(3-2)=0.90-1.07 (N-core, within 1 sigma+systematic): count',len(core),'of',len(df[df.X==50]))
print(core.groupby('A')[['tau32_HCN','tau32_HCO']].agg(['min','median','max']).round(2).to_string())
out=df[(df.R32>0.33)&(df.R32<0.47)&(df.X==50)]
print('\nCells with R_K(3-2)=0.33-0.47 (outer N annulus / S):',len(out)); print(out.groupby('A')[['tau32_HCN','tau32_HCO']].agg(['min','median','max']).round(2).to_string())
print('\nPredicted HCN/H13CN(4-3) in the R_K~1 cells, per X and per A:'); print(df[(df.R32>0.90)&(df.R32<1.07)].groupby(['X','A']).HCN_H13CN_43.agg(['min','median','max']).round(1).to_string())
print('\nPredicted HCN/H13CN(4-3) in the R_K~0.4 cells:'); print(df[(df.R32>0.33)&(df.R32<0.47)].groupby(['X','A']).HCN_H13CN_43.agg(['min','median','max']).round(1).to_string())
