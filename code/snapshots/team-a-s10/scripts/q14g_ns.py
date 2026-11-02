# Q14g (fix 2): S-N contrast of CH3OH 254 / HCO+ 3-2 with a paired empirical null: the SAME estimator (local baseline, HCO+ (v0,FWHM) window per region) applied at the same N and S pixels to shifted line-free windows of the CH3OH cube
# (null_dR = [x_S(off)/H_S - x_N(off)/H_N]) and to 60 random off-source position pairs (9-16" ring, same window); significance = dR_obs / rms(null). Also the S/N and Bridge-vs-N contrasts.
# Output out/q14e/ns_contrast.csv. python3 -I q14g_ns.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
F0={'HCO+':267.557626,'CH3OH':254.015377}; SPW={'Xb4d':{'HCO+':'spw25','CH3OH':'spw29'},'Xb4f':{'HCO+':'spw16','CH3OH':'spw20'}}
BEAMS=[('Xb4d','1.7',(1.7,1.7,0.)),('Xb4d','3.0',(3.0,3.0,0.)),('Xb4d','4.45x2.95',(4.45,2.95,87.)),('Xb4d','6.3x3.8',(6.3,3.8,82.)),('Xb4f','6.3x3.8',(6.3,3.8,82.))]
rng=np.random.default_rng(7); out=[]
for ous,bn,tb in BEAMS:
    D={}
    for ln in F0:
        fr,a,cw,pix,beam,nb=cut(ous,SPW[ous][ln],tb,half=20.); D[ln]=(vel(fr,F0[ln]),a,cw)
    vh,ah,cwh=D['HCO+']; vc,ac,cwc=D['CH3OH']
    for wname in ('FWHM','Harada'):
        win={}
        for reg in ('N','S'):
            ra,de=REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']; p=hco_fit(vh,sample(ah,cwh,ra,de))
            win[reg]=(p[1],p[2]) if wname=='FWHM' else (2775.,300.)
        def est(reg,off=0.,pos=None):
            ra,de=(REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']) if pos is None else pos
            c,W=win[reg]; sh=sample(ah,cwh,REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']); sc=sample(ac,cwc,ra,de)
            H,_=m0(vh,sh,c,W,excl=(c-W/2,c+W/2))
            x,_=m0(vc,sc,c+off,W,excl=(c-W,c+W)) if off else m0(vc,sc,c,W,excl=(c-W/2,c+W/2))
            return x,H
        xN,HN=est('N'); xS,HS=est('S'); dR=xS/HS-xN/HN
        nd=[];nx={'N':[],'S':[]}
        for off in np.arange(-1500,1501,150):
            if abs(off)<max(win['N'][1],win['S'][1])*1.3: continue
            a1,_=est('N',off); a2,_=est('S',off)
            if np.isfinite(a1) and np.isfinite(a2): nd.append(a2/HS-a1/HN); nx['N'].append(a1); nx['S'].append(a2)
        nd=np.array(nd); nb_=[]
        def P(r,t): return (CEN[0]+r*np.sin(t)/3600/np.cos(np.deg2rad(CEN[1])),CEN[1]+r*np.cos(t)/3600)
        for k in range(60):
            r1,r2=rng.uniform(9,16,2); t1,t2=rng.uniform(0,2*np.pi,2)
            a1,_=est('N',pos=P(r1,t1)); a2,_=est('S',pos=P(r2,t2)); nb_.append(a2/HS-a1/HN)
        nb_=np.array(nb_); sa=np.std(nd); sb=np.std(nb_); s=max(sa,sb)
        out.append(dict(ous=ous,beam=bn,window=wname,I_N=xN,H_N=HN,I_S=xS,H_S=HS,R_N=xN/HN,R_S=xS/HS,dR=dR,n_null_win=len(nd),null_win_rms=sa,null_win_mad=robust_sigma(nd),null_pos_rms=sb,
                        sig_win=dR/sa,sig_pos=dR/sb,sig_conservative=dR/s,frac_null_win_ge=np.mean(nd>=dR),frac_null_pos_ge=np.mean(nb_>=dR),
                        S_sig_vs_null=xS/np.std(nx['S']),N_sig_vs_null=xN/np.std(nx['N']),R_N_err=np.std(nx['N'])/HN,R_S_err=np.std(nx['S'])/HS))
    print(ous,bn,flush=True)
O=pd.DataFrame(out); O.to_csv('/workspace/out/q14e/ns_contrast.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_columns',40); print(O.round(4).to_string())
