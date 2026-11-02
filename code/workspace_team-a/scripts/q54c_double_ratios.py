"""Q54: N/S double ratios D = (X/HCN)_N/(X/HCN)_S (flux-scale independent: same two cubes at both nuclei) with window/baseline variants.
Variants: velocity half-window HW 250/350/450 (HCN, HNC) ; CN frequency window pad -0.04/0/+0.06 GHz (3-2) or -0.03/0/+0.04 (1-0); baseline (450-950) or (550-1100).
Output out/q54/double_ratios.csv (median and range over variants; stat error from the fiducial)."""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, itertools, copy
import q54_chem as c, q54b_cn_ratios as b
W0=copy.deepcopy(b.W)
PAIRS=[('CN1-0_X4d','HCN1-0_00993',2.5),('CN1-0_X4d','HCN1-0_00412',2.5),('CN1-0_main','HCN1-0_00993',2.5),
       ('HNC1-0_00993','HCN1-0_00993',2.5),('HNC1-0_00412','HCN1-0_00412',1.4),('HNC1-0_00412','HCN1-0_00412',2.5),
       ('CN3-2_X1df1','HCN3-2_Xb4d',1.4),('CN3-2_X1df1','HCN3-2_Xb4d',2.5),('CN3-2_X1df1','HCN4-3_X1df1',1.4),('CN3-2_X1df1','HCN4-3_X1df1',2.5),
       ('HNC3-2_00412','HCN3-2_Xb4d',1.4),('HNC3-2_00412','HCN3-2_Xb4d',2.5)]
def val(key,bm,reg,hw,pad):
    if key in b.W:
        w=b.W[key]; w['win']=(W0[key]['win'][0]-pad[0],W0[key]['win'][1]+pad[1]); return b.cn(key,bm,reg)
    return c.meas(key,bm,reg,hw)[:2]
rows=[]
for a,d,bm in PAIRS:
    D=[];fid=None
    for hw,bl,k in itertools.product([250,350,450],[(450,950),(550,1100)],[-1,0,1]):
        c.BL=bl
        for kk in b.W: b.W[kk]['bl']=W0[kk]['bl']
        if bl[0]==550:
            for kk in b.W: b.W[kk]['bl']=[(lo-0.03,hi) if lo<hi else (lo,hi) for lo,hi in W0[kk]['bl']]  # modestly shifted baseline windows
        p=(0.04,0.06) if '3-2' in a else (0.03,0.04); pad=(-0.5*p[0]*(k==-1)+p[0]*(k==1)*0,p[1]*(k==1)-p[0]*(k==-1)) if False else ((-p[0] if k==-1 else 0),(p[1] if k==1 else (-p[0] if k==-1 else 0)))
        try:
            r={}
            for reg in 'NS':
                A=val(a,bm,reg,hw,pad); B=val(d,bm,reg,hw,pad); r[reg]=(A[0]/B[0],abs(A[0]/B[0])*np.hypot(A[1]/A[0],B[1]/B[0]))
        except Exception as ex: print('skip',a,d,ex); break
        Dv=r['N'][0]/r['S'][0]; D.append((Dv,r['N'][0],r['S'][0]))
        if hw==350 and bl==(450,950) and k==0: fid=(Dv,Dv*np.hypot(r['N'][1]/r['N'][0],r['S'][1]/r['S'][0]),r['N'],r['S'])
    D=np.array(D)
    if fid is None: continue
    rows.append(dict(num=a,den=d,beam=bm,RN=fid[2][0],eRN=fid[2][1],RS=fid[3][0],eRS=fid[3][1],D=fid[0],eD_stat=fid[1],Dmin=D[:,0].min(),Dmax=D[:,0].max(),RNmin=D[:,1].min(),RNmax=D[:,1].max(),RSmin=D[:,2].min(),RSmax=D[:,2].max()))
    print(rows[-1],flush=True)
pd.DataFrame(rows).round(3).to_csv('/workspace/out/q54/double_ratios.csv',index=False)
