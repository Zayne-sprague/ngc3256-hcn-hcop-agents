# Q14f (fix 1): Harada+2018 matched test. (i) LTE conversion conventions (CMB contrast on/off), 241.79/96.74 Trot from Harada's own Table 13;
# (ii) Harada-window (+-150 km/s about 2775) moment-0 maps: Xb4d at 1.7", Xb4d & Xb4f at 6.3x3.8; "feathered" 1.7" estimate = Xb4d(1.7") + [Xb4f - Xb4d](6.3x3.8) (adds the low-spatial-frequency deficit of the 12-m data; HCO+ used as the control);
# (iii) positional test: max of feathered/Xb4d maps within 0.5/1" of Harada's N/S (our N,S positions) ; (iv) aperture recovery Xb4d/Xb4f vs radius from native-resolution moment maps.
# Output out/q14f/*.csv. python3 -I q14f_harada.py
import sys,os; sys.path.insert(0,'/workspace/scripts')
os.makedirs('/workspace/out/q14f',exist_ok=True)
from q14e_lib import *
sys.argv=['x']
import importlib.util
K1=1.4388; hP=6.62607e-27; cc=2.998e10; kB=1.380649e-16
def lam(p):
    L=open(p).read().split('\n'); i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]; n=int(L[i+1])
    lv=np.array([[float(a.split()[1])*K1,float(a.split()[2])] for a in L[i+3:i+3+n]])
    j=[k for k,l in enumerate(L) if 'NUMBER OF RADIATIVE' in l][0]; m=int(L[j+1])
    tr=np.array([[float(x) for x in l.split()[:6]] for l in L[j+3:j+3+m]]); return lv,tr
def Wtab(path,sel,T,cmb=True):
    lv,tr=lam(path); Q=np.sum(lv[:,1]*np.exp(-lv[:,0]/T)); s=0.
    for r in tr:
        if not sel(r[4]): continue
        up=int(r[1])-1; Nu=lv[up,1]*np.exp(-lv[up,0]/T)/Q; nu=r[4]
        w=hP*cc**3*r[3]*Nu/(8*np.pi*kB*(nu*1e9)**2)/1e5
        if cmb: x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; w*=1-(np.exp(x)-1)/(np.exp(xb)-1)
        s+=w
    return s,lv,tr
E_='/opt/lamda/e-ch3oh.dat'; A_='/opt/lamda/ch3oh_a.dat'
def grp(sel,T,cmb=True):   # thin-LTE K km/s per unit TOTAL column (N_A=N_E=N/2, each ladder on its own Q)
    return 0.5*Wtab(E_,sel,T,cmb)[0]+0.5*Wtab(A_,sel,T,cmb)[0]
S254=lambda nu:abs(nu-254.015)<0.01; S96=lambda nu:96.73<nu<96.76; S241=lambda nu:241.69<nu<241.93
rows=[]
for T in (3,4,5,6,8,10,15,20,30,50,100):
    for cmb in (True,False):
        w254=0.5*Wtab(E_,S254,T,cmb)[0]; w96=grp(S96,T,cmb); w241=grp(S241,T,cmb)
        rows.append(dict(T=T,cmb_contrast=cmb,R96_254=w96/w254,R241_96=w241/w96,R241_254=w241/w254))
H=pd.DataFrame(rows); H.to_csv('/workspace/out/q14f/lte_ratios.csv',index=False)
lv,tr=lam(E_); print('lines in 96.73-96.76 (E file):',[(r[4],r[3]) for r in tr if S96(r[4])]); lv2,tr2=lam(A_); print('A file:',[(r[4],r[3]) for r in tr2 if S96(r[4])])
print('241 group E:',[round(r[4],3) for r in tr if S241(r[4])],'A:',[round(r[4],3) for r in tr2 if S241(r[4])])
print(H[H.cmb_contrast].round(3).to_string())
# Harada's own 241.79/96.74 Trot (K km/s; Table 13): N 5.84+-0.57 / 3.09+-0.66, S 4.35+-0.57(0.6)/6.22+-0.66
rng=np.random.default_rng(3)
Ts=np.linspace(2,60,300); R241=np.interp(Ts,H[H.cmb_contrast].T,H[H.cmb_contrast].R241_96) if False else None
Tg=np.array(sorted(set(H['T']))); Hc=H[H.cmb_contrast].set_index('T')
from scipy.interpolate import PchipInterpolator
f241=PchipInterpolator(np.log(Tg),Hc.loc[Tg,'R241_96'].values); f96=PchipInterpolator(np.log(Tg),Hc.loc[Tg,'R96_254'].values)
tg=np.exp(np.linspace(np.log(3),np.log(100),600))
out=[]
for reg,(i241,e241,i96,e96) in {'N':(5.84,0.57,3.09,0.66),'S':(4.35,0.57,6.22,0.66)}.items():
    r=rng.normal(i241,e241,40000)/np.clip(rng.normal(i96,e96,40000),0.3,None)
    T=np.exp(np.interp(r,f241(np.log(tg)),np.log(tg),left=np.nan,right=np.nan))
    ok=np.isfinite(T); pc=np.percentile(T[ok],[16,50,84,97.5]); Rp=f96(np.log(T[ok]))
    out.append(dict(region=reg,R241_96=i241/i96,Tmed=pc[1],T16=pc[0],T84=pc[2],T975=pc[3],frac_in_range=ok.mean(),R96_254_med=np.median(Rp),R96_254_16=np.percentile(Rp,16),R96_254_84=np.percentile(Rp,84),
                    pred_254_Kkms_med=i96/np.median(Rp),pred_254_lo=np.percentile(i96/Rp,16),pred_254_hi=np.percentile(i96/Rp,84)))
HT=pd.DataFrame(out); HT.to_csv('/workspace/out/q14f/harada_trot.csv',index=False); print(HT.round(3).to_string())
# Harada Table 13 other positions with both 254 and 96.74: TNE 254 0.46+-0.72 / 4.33+-0.68; TSW 0.57+-0.69 / 15.4+-0.68  (K km/s, 1.7")
for nm,i254,e254,i96,e96 in (('TNE',0.46,0.72,4.33,0.68),('TSW',0.57,0.69,15.4,0.68)):
    print(nm,'254/96.74 = %.3f +- %.3f'%(i254/i96,e254/i96),' LTE R96_254>=: 3sigma-ish limit on 254/96.74 <',(i254+3*e254)/i96, '-> R96/254 > ',i96/(i254+3*e254))
# --- maps
HW=(2775.,300.)
def maps(ous,tb,line):
    spw={'HCO+':{'Xb4d':'spw25','Xb4f':'spw16'},'CH3OH':{'Xb4d':'spw29','Xb4f':'spw20'}}[line][ous]
    f0={'HCO+':267.557626,'CH3OH':254.015377}[line]
    fr,a,cw,pix,beam,nb=cut(ous,spw,tb,half=22.); v=vel(fr,f0)
    return m0map(v,a,HW[0],HW[1],excl=(HW[0]-HW[1]/2,HW[0]+HW[1]/2)),cw,pix,beam,nb,v,a
M={}
for ln in ('HCO+','CH3OH'):
    M[(ln,'4d17')]=maps('Xb4d',(1.7,1.7,0.),ln); M[(ln,'4d63')]=maps('Xb4d',(6.3,3.8,82.),ln); M[(ln,'4f63')]=maps('Xb4f',(6.3,3.8,82.),ln)
def regval(key,reg): m,cw=M[key][:2]; return val(m,cw,REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg'])
res=[]
for reg in ('N','S','Bridge','Wclump'):
    d={}
    for ln in ('HCO+','CH3OH'):
        a=regval((ln,'4d17'),reg); b=regval((ln,'4d63'),reg); c=regval((ln,'4f63'),reg)
        d[ln]=dict(x4d17=a,x4d63=b,x4f63=c,deficit=c-b,feather=a+(c-b),rec4d_4f=b/c)
    res.append(dict(region=reg,**{f'{ln}_{k}':v for ln in d for k,v in d[ln].items()}))
RS=pd.DataFrame(res); RS.to_csv('/workspace/out/q14f/feather_regions.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_columns',40); print(RS.round(3).T.to_string())
# positional test: max within r of the region for Xb4d 1.7 (feather offset constant)
pos=[]
for reg in ('N','S'):
    for r in (0.3,0.5,1.0):
        out_={}
        for ln in ('HCO+','CH3OH'):
            m,cw,pix=M[(ln,'4d17')][:3]; dfc=regval((ln,'4f63'),reg)-regval((ln,'4d63'),reg)
            ny,nx=m.shape; yy,xx=np.mgrid[0:ny,0:nx]; wc=cw.pixel_to_world_values(xx,yy)
            ra,de=REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']; dr=np.hypot((wc[0]-ra)*np.cos(np.deg2rad(de))*3600,(wc[1]-de)*3600)
            out_[ln]=np.max(np.where(dr<r,m,-1e9))+dfc; out_[ln+'_min']=np.min(np.where(dr<r,m,1e9))+dfc
        pos.append(dict(region=reg,r=r,HCO_max=out_['HCO+'],CH3OH_max=out_['CH3OH'],CH3OH_min=out_['CH3OH_min'],ratio_max=out_['CH3OH']/out_['HCO+']))
PS=pd.DataFrame(pos); PS.to_csv('/workspace/out/q14f/position_test.csv',index=False); print(PS.round(3).to_string())
# --- aperture recovery from native-resolution moment maps (Jy km/s within radius r of the N-S midpoint)
ap=[]
for ln in ('HCO+','CH3OH'):
    for ous in ('Xb4d','Xb4f'):
        spw={'HCO+':{'Xb4d':'spw25','Xb4f':'spw16'},'CH3OH':{'Xb4d':'spw29','Xb4f':'spw20'}}[ln][ous]
        f0={'HCO+':267.557626,'CH3OH':254.015377}[ln]
        h=fits.open(f(ous,spw,'cube.I.pbcor'),memmap=True); hd=h[0].header; fr=freqaxis(hd)/1e9
        nb=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']); w=WCS(hd).celestial; pix=abs(hd['CDELT2'])*3600
        px,py=[float(z) for z in w.world_to_pixel_values(*CEN)]; n=int(26/pix); x0=int(px)-n; y0=int(py)-n
        d=h[0].data; d=d[0] if d.ndim==4 else d; a=np.nan_to_num(np.array(d[:,max(0,y0):y0+2*n+1,max(0,x0):x0+2*n+1],np.float32)); v=vel(fr,f0)
        mp=m0map(v,a,HW[0],HW[1],excl=(HW[0]-HW[1]/2,HW[0]+HW[1]/2)); yy,xx=np.mgrid[0:mp.shape[0],0:mp.shape[1]]
        cx,cy=px-max(0,x0),py-max(0,y0); rr=np.hypot(xx-cx,yy-cy)*pix; barea=OM(nb)/pix**2
        for R in (1.5,3,4.5,6,9,12,16,20):
            ap.append(dict(line=ln,ous=ous,R=R,flux_Jykms=float(mp[rr<R].sum()/barea)))
AP=pd.DataFrame(ap).pivot_table(index=['line','R'],columns='ous',values='flux_Jykms').reset_index(); AP['rec_4d_over_4f']=AP.Xb4d/AP.Xb4f
AP.to_csv('/workspace/out/q14f/aperture_recovery.csv',index=False); print(AP.round(3).to_string())
