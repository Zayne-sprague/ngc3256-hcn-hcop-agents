# Q14o (fix 1): CH3OH 254.015 and HCO+ 3-2 aperture flux (Jy km/s, circular apertures of radius R around N, S and the N-S midpoint) at NATIVE resolution from Xb4d (12-m long baselines) and Xb4f (compact), Harada window (2775+-150 km/s),
# local per-pixel linear baseline; noise from the same aperture sums in 6 shifted line-free 300 km/s windows (offsets +-[550..850] km/s; the CH3OH window offsets avoid SiS/30SiO at +104/+237 km/s).  Output out/q14e/aperture_recovery.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
CENS={'N':(REG.loc['N','RA_deg'],REG.loc['N','Dec_deg']),'S':(REG.loc['S','RA_deg'],REG.loc['S','Dec_deg']),'mid':CEN}
SP={'HCO+':({'Xb4d':'spw25','Xb4f':'spw16'},267.557626),'CH3OH':({'Xb4d':'spw29','Xb4f':'spw20'},254.015377)}
OFFS=[-850,-700,-550,550,700,850]; RAD=[1.5,3,4.5,6,9,12]; rows=[]
for ln,(spw,f0) in SP.items():
    for ous in ('Xb4d','Xb4f'):
        h=fits.open(f(ous,spw[ous],'cube.I.pbcor'),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9; v=vel(fr,f0)
        nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600; barea=OM(nb)/pix**2
        px,py=[float(z) for z in w.world_to_pixel_values(*CEN)]; n=int(24/pix); x0=max(0,int(px)-n); y0=max(0,int(py)-n)
        d=h[0].data; d=d[0] if d.ndim==4 else d; a=np.nan_to_num(np.array(d[:,y0:int(py)+n+1,x0:int(px)+n+1],np.float32))
        yy,xx=np.mgrid[0:a.shape[1],0:a.shape[2]]
        maps={0:m0map(v,a,2775.,300.,excl=(2625,2925))}
        for o in OFFS:
            c=2775.+o
            if c-150<v.min()+100 or c+150>v.max()-100: continue
            maps[o]=m0map(v,a,c,300.,excl=(2625-50,2925+50))
        for cn,(ra,de) in CENS.items():
            cx,cy=[float(z) for z in w.world_to_pixel_values(ra,de)]; cx-=x0; cy-=y0; rr=np.hypot(xx-cx,yy-cy)*pix
            for R in RAD:
                m=rr<R; fl={o:float(mp[m].sum()/barea) for o,mp in maps.items()}
                nul=np.array([fl[o] for o in fl if o!=0])
                rows.append(dict(line=ln,ous=ous,centre=cn,R=R,flux_Jykms=fl[0],null_rms=np.std(nul) if len(nul)>2 else np.nan,null_mean=np.mean(nul) if len(nul)>2 else np.nan,n_null=len(nul)))
A=pd.DataFrame(rows); P=A.pivot_table(index=['line','centre','R'],columns='ous',values=['flux_Jykms','null_rms']); P.columns=['_'.join(c) for c in P.columns]; P=P.reset_index()
P['rec_4d_over_4f']=P.flux_Jykms_Xb4d/P.flux_Jykms_Xb4f; P['e_rec']=np.abs(P.rec_4d_over_4f)*np.hypot(P.null_rms_Xb4d/P.flux_Jykms_Xb4d,P.null_rms_Xb4f/P.flux_Jykms_Xb4f)
P['sig_4f']=P.flux_Jykms_Xb4f/P.null_rms_Xb4f; P['sig_4d']=P.flux_Jykms_Xb4d/P.null_rms_Xb4d
P.to_csv('/workspace/out/q14e/aperture_recovery.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_rows',200); print(P.round(3).to_string())
