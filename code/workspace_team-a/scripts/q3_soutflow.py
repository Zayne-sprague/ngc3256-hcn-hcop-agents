# Q3 step 3: Soutflow (blue lobe S of S nucleus) and OF-Red lobe: velocity-resolved HCN/HCO+ 3-2 ratio (Xb4d, 1.05"), spectral shape vs scaled S-nucleus spectrum,
# 3-2 vs Michiyama+2018 HCN/HCO+ 1-0 at a matched 2.1" beam, dense-gas luminosity.
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_pos import *
import pandas as pd
nu={'HCN':265.886434/(1+VS/C),'HCO+':267.557626/(1+VS/C)}
c1={sp:cube(sp,fw=1.05) for sp in nu}
cT=Tb_factor(nu['HCN'],1.05)/Tb_factor(nu['HCO+'],1.05)
rows=[]
for pn in ('Soutflow','M18 OF-Red','S','N'):
    x,y=POS[pn]; S={sp:spec(c1[sp],x,y) for sp in nu}
    for lo in range(2175,3225,100):
        r=dict(pos=pn,v_lo=lo,v_hi=lo+100)
        for sp in nu:
            c=c1[sp]; m=(c['v']>=lo)&(c['v']<lo+100); n=m.sum(); fm=chan_noise(c)
            e=1.4826*np.median(np.abs(S[sp][fm]-np.median(S[sp][fm])))   # per-channel (Jy/beam)
            r['I_'+sp]=S[sp][m].sum()*c['dv']; r['e_'+sp]=e*c['dv']*np.sqrt(n)
        a,b,ea,eb=r['I_HCN'],r['I_HCO+'],r['e_HCN'],r['e_HCO+']
        if a/ea>3 and b/eb>3: r['R']=a/b*cT; r['eR']=r['R']*np.hypot(ea/a,eb/b); r['kind']='value'
        elif a/ea>3: r['R']=a/(3*eb)*cT; r['kind']='lower limit'
        elif b/eb>3: r['R']=3*ea/b*cT; r['kind']='upper limit'
        rows.append(r)
V=pd.DataFrame(rows); V.round(3).to_csv('/workspace/out/q3/q3_velocity_ratio.csv',index=False)
pd.set_option('display.width',250); print(V[V.pos.isin(['Soutflow','M18 OF-Red'])].round(3).to_string())
# --- blend limit: CH2NH 4(1,3)-3(1,2) (263.8008 GHz sky at Vsys = HCN-frame v=%.0f) and CH3OH 3(2,1)-3(0,3) E etc. measured at S and N at the same sky frequency
f0=265.886434
for nm,fs in (('CH2NH 4_13-3_12',263.8008),('CH3NH2 263.641',263.6413),('CH3OH 3(-2)-3(0)E',263.72785)):
    vh=(1-fs/f0)*C
    for pn in ('S','N','Soutflow'):
        c=c1['HCN']; s=spec(c,*POS[pn]); fm=chan_noise(c); e=1.4826*np.median(np.abs(s[fm]-np.median(s[fm])))
        # line at the position's own systemic velocity: shift by (v_local-VS)
        vloc=2722. if pn=='S' else (2800. if pn=='N' else 2650.)
        m=(c['v']>=vh+(vloc-VS)-60)&(c['v']<=vh+(vloc-VS)+60)
        print(f'{nm} @{pn}: expected HCN-frame v={vh+vloc-VS:.0f}; mean in +-60 km/s = {s[m].mean()*1e3:.2f} mJy/beam (chan rms {e*1e3:.2f}, mean err {e*1e3/np.sqrt(m.sum()):.2f}); HCN peak {s.max()*1e3:.1f}')
# --- 3-2 vs Michiyama 1-0 at 2.1" beam: surface brightness in Jy/beam km/s in the M18 blue range
c2={sp:cube(sp,fw=2.1,crop=14.) for sp in nu}
R32={}
for vs in (2740,2775,2810):
    for pn in ('M18 South-OF-Blue','M18 OF-Red'):
        x,y=POS[pn]; rng=(vs-360,vs-180) if 'Blue' in pn else (vs+240,vs+360); r=dict(vsys=vs,pos=pn,range=rng)
        for sp in nu:
            c=c2[sp]; m0,_=wsum(c,*rng); sg,_=wnoise(c,*rng); r['I32_'+sp]=samp(c,m0,x,y); r['e32_'+sp]=samp(c,sg,x,y)
        R32[(vs,pn)]=r
M18={'M18 South-OF-Blue':(0.07,0.02,0.10,0.02),'M18 OF-Red':(0.08,0.02,0.08,0.02)}   # Michiyama+2018 Table 3: HCN 1-0, HCO+ 1-0 (Jy km/s) with errors
F10={'HCN':88.6316,'HCO+':89.1885}
out=[]
for (vs,pn),r in R32.items():
    h,eh,o,eo=M18[pn]
    for sp,(a,ea) in (('HCN',(h,eh)),('HCO+',(o,eo))):
        I32,e32=r['I32_'+sp],r['e32_'+sp]
        q=(I32/a)*(F10[sp]/nu[sp])**2   # Tb(3-2)/Tb(1-0) for equal beams (2.1")
        eq=q*np.hypot(e32/I32,ea/a) if I32>0 else np.nan
        out.append(dict(vsys=vs,pos=pn,line=sp,range=str(r['range']),I32_Jybeamkms=I32,e32=e32,I10_M18=a,e10=ea,r31=q,e_r31=eq))
O=pd.DataFrame(out); O.round(3).to_csv('/workspace/out/q3/q3_r31_vs_M18.csv',index=False); print(O.round(3).to_string())
