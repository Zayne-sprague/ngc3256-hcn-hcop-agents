import pandas as pd,glob,numpy as np
for f in sorted(glob.glob('/workspace/out/q12b/runs/base_*_r8_tbg2.73_s1.csv')):
    r=pd.read_csv(f).iloc[0]
    p=[r[f'pull{i}'] for i in range(8)]
    print(f.split('/')[-1],'chi2prof %.2f sum pull2 %.2f dof %d'%(r.chi2_prof,sum(x*x for x in p),r.dof),'g %.2f h %.2f'%(r.g_best,r.h_best),np.round(p,1), 'Tc',round(r.Tc_50),round(r.tau32H_c_50) if 'tau32H_c_50' in r else '')
d=np.load(glob.glob('/workspace/out/q12b/runs/base_1z_N_r8_tbg2.73_s1.npz')[0])['dobs']; print(np.round(d,3))
# algebraic redundancy check: R32 vs R10*(H32/10)/(O32/10)
print('R32',d[4,0],'R10*ratio',d[5,0]*d[6,0]/d[7,0])
