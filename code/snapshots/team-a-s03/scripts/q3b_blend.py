# Q3 round 2, fix 1: is the HCN blue-side excess at Soutflow a blend?  (a) HCN = k*HCO+ residual; (b) independent LTE forward-model of candidate blenders (CH3NH2, CH3OH, CH2NH, CH3CHO, CH3SH, HDO) fitted to ALL other Band-6 channels at the same position; (c) the same on Xb4b smoothed to 1.05" and on the tapered Xb4d cubes.
# Outputs: out/q3b/blend_*.csv, blend.log, blend_fig.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
from ident import cat
import pandas as pd
POSN={'Soutflow':POS['Soutflow'],'M18 South-OF-Blue':POS['M18 South-OF-Blue']}
RAD=0.5     # aperture radius (arcsec) on the 1.05" smoothed cube (~1" box mean as in the review check)
HW=(2375,2575)   # HCN-frame bins in question
# ---------- (a) k*HCO+ model ------------
def prof(sp,x,y,ous='Xb4d',fw=1.05):
    c=vcube(sp,ous,fw,crop=9.); s=aper(c,x,y,RAD)*1e3; v=c['v']
    s,e=baseline(v,s,FREE[sp],1); return v,s,e,c
def binned_noise(v,s,free,width,dv):
    """empirical sigma of a boxcar mean of `width` km/s on line-free channels (includes channel correlation)"""
    n=max(1,int(round(width/dv))); m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    k=np.ones(n)/n; sm=np.convolve(s,k,'valid'); mv=np.convolve(m.astype(float),k,'valid')>0.999
    return robust_sigma(sm[mv])
def bins_on(v,s,edges): return np.array([s[(v>=a)&(v<b)].mean() for a,b in zip(edges[:-1],edges[1:])])
def part_a(pn,ous='Xb4d'):
    x,y=POSN[pn]; vh,sh,eh,ch=prof('HCN',x,y,ous); vo,so,eo,co=prof('HCO+',x,y,ous)
    edges=np.arange(2175,2925,50.); cen=0.5*(edges[1:]+edges[:-1])
    H=bins_on(vh,sh,edges); O=bins_on(vo,so,edges)
    nh=binned_noise(vh,sh,FREE['HCN'],50,ch['dv']); no=binned_noise(vo,so,FREE['HCO+'],50,co['dv'])
    m=(cen>2575)&(cen<2875)
    # weighted LS for k on the blend-free range (errors on both axes: ODR-like via effective variance, iterated)
    k=np.sum(H[m]*O[m])/np.sum(O[m]**2)
    for _ in range(20):
        w=1/(nh**2+(k*no)**2); k=np.sum(w*H[m]*O[m])/np.sum(w*O[m]**2)
    ek=1/np.sqrt(np.sum(O[m]**2/(nh**2+(k*no)**2)))
    res=H-k*O; eres=np.full(len(H),np.sqrt(nh**2+(k*no)**2))
    chi2=np.sum((res[m]/eres[m])**2)/(m.sum()-1)
    return dict(pn=pn,ous=ous,cen=cen,H=H,O=O,res=res,eres=eres,k=k,ek=ek,chi2red_fit=chi2,nh=nh,no=no)
# ---------- (b) LTE forward model ------------
SPECIES={'CH3NH2':['CH3NH2'],'CH3OH':['CH3OH'],'CH2NH':['CH2NH'],'CH3CHO':['CH3CHO'],'CH3SH':['CH3SH']}   # HDO (single in-band line, Eu=157 K) has no other line to calibrate it: not testable, discussed in the text
SPWS=[('Xb4d','spw31'),('Xb4d','spw29'),('Xb4d','spw27'),('Xb4d','spw25')]
def load_spectra(x,y,ous='Xb4d',fw=1.05,rad=RAD):
    out={}
    for o,sw in [(ous,s) for _,s in SPWS]:
        c=fcube(o,sw,fw,9.); s=aper(c,x,y,rad)*1e3; out[sw]=(c['fr'],s,c)
    return out
def template(v,s,width=30.):
    """normalised gas-velocity profile from the HCO+ spectrum at the same position (smoothed, positive part, zero outside 2150..2950)"""
    n=int(round(width/4.4)); k=np.ones(n)/n; sm=np.convolve(s,k,'same'); sm=np.clip(sm,0,None); sm[(v<2150)|(v>2950)]=0
    sm/=sm.max(); o=np.argsort(v); return lambda vg: np.interp(vg,v[o],sm[o],left=0,right=0)
def box_template(lo=2375,hi=2875,edge=40.):
    return lambda vg: np.clip(np.minimum((vg-lo)/edge,(hi-vg)/edge)+0.5,0,1)
MASK={'spw27':[(2000,3400)],'spw25':[(1650,1900),(2000,3400)],'spw29':[(1900,3400)],'spw31':[(1900,3400)]}
def fit_species(sp,T,P,spec,cats):
    L=cats[cats.sp.str.split().str[0].str.strip(',').isin(SPECIES[sp])&(cats.Eu<300)].copy()
    L['w']=L.sijmu2*L.f0**3*np.exp(-L.Eu/T); L=L[L.w>0]
    mods={};datas={};wts={}
    for sw,(fr,s,c) in spec.items():
        v=np.zeros(len(fr)); M=np.zeros(len(fr))
        sel=L[(L.f0*(1-3300/C)<fr.max()+0.01)&(L.f0*(1-1500/C)>fr.min()-0.01)]
        for r in sel.itertuples(): M+=r.w*P((1-fr/r.f0)*C)
        mods[sw]=M
    wmax=L.w.max(); 
    num=0;den=0;resid={}
    # baseline: linear fit per spw to unmasked channels, simultaneously with A (iterate: fit baseline to s - A M, then A)
    A=0.; E={}
    for it in range(3):
        num=0;den=0
        for sw,(fr,s,c) in spec.items():
            v=vel(fr,F0['HCN']) if sw!='spw25' else vel(fr,F0['HCO+'])
            # use own velocity convention just for masks: spw29/31 mask by SO-ish window around the band's line centres is approximated by whole-window cut 1900..3400 of HCN-like frame (conservative)
            if sw in('spw29','spw31'): vv=vel(fr,F0['HCN']*0+251.825770 if sw=='spw29' else 249.9)
            else: vv=v
            mk=np.ones(len(fr),bool)
            for a,b in MASK[sw]: mk&=~((vv>=a)&(vv<=b))
            y=s-A*mods[sw]; X=np.vander(fr-fr.mean(),3); coef=np.linalg.lstsq(X[mk],y[mk],rcond=None)[0]; bl=X@coef
            r=s-bl; sig=robust_sigma((r-A*mods[sw])[mk]); E[sw]=(mk,sig,bl)
            num+=np.sum(r[mk]*mods[sw][mk])/sig**2; den+=np.sum(mods[sw][mk]**2)/sig**2
        A=num/den if den>0 else 0
    eA=1/np.sqrt(den) if den>0 else np.nan
    return A,eA,mods,L,E
def autocorr_factor(spec):
    """variance inflation for correlated channel noise: 1+2*sum rho_k (k<=6), from masked residuals of the HCN-spw baseline"""
    fr,s,c=spec['spw27']; v=vel(fr,F0['HCN']); m=(v>1400)&(v<2050); r=s[m]-np.polyval(np.polyfit(v[m],s[m],1),v[m]); r-=r.mean()
    rho=[np.sum(r[:-k]*r[k:])/np.sum(r*r) for k in range(1,7)]; return max(1.0,1+2*np.sum(rho)),rho
def part_b(pn,tag,temps=(15,30,60,120),ous='Xb4d',fw=1.05,use_ident=True):
    x,y=POSN[pn] if pn in POSN else POS[pn]
    spec=load_spectra(x,y,ous,fw)
    cats=cat(300); infl,rho=autocorr_factor(spec)
    fr,s,c=spec['spw25']; vo=vel(fr,F0['HCO+']); so,_=baseline(vo,s,FREE['HCO+'],1)
    P1=template(vo,so); P2=box_template()
    fh,sh,_=spec['spw27']; vh=vel(fh,F0['HCN'])
    rows=[]
    for spn in SPECIES:
        for T in temps:
            for pname,P in (('HCO+prof',P1),('box',P2)):
                A,eA,mods,L,E=fit_species(spn,T,P,spec,cats); eA*=np.sqrt(infl)
                Aul=max(A,0)+2*eA     # ~97.7% one-sided upper bound
                Mh=mods['spw27']
                row=dict(pos=pn,species=spn,T=T,templ=pname,A=A,eA=eA,A_over_eA=A/eA if eA>0 else np.nan,nlines=len(L))
                for lo in (2275,2375,2475,2575):
                    k=(vh>=lo)&(vh<lo+100); row[f'blend_ul_{lo}']=Aul*Mh[k].mean(); row[f'blend_best_{lo}']=A*Mh[k].mean()
                rows.append(row)
    return pd.DataFrame(rows),infl,rho,spec
def part_c(pn,tag,loader,edges,free_h,free_o,order=0):
    """k*HCO+ residual in an independent data product; loader(sp)->(v,spec mJy at position,chan_dv)"""
    x,y=POSN[pn]; out={}
    for sp,fr_ in (('HCN',free_h),('HCO+',free_o)):
        v,s0=loader(sp,x,y); s1,e=baseline(v,s0,fr_,order); out[sp]=(v,s1,e)
    cen=0.5*(edges[1:]+edges[:-1]); H=bins_on(out['HCN'][0],out['HCN'][1],edges); O=bins_on(out['HCO+'][0],out['HCO+'][1],edges)
    nh=out['HCN'][2]/np.sqrt(len(edges)); no=out['HCO+'][2]/np.sqrt(len(edges))
    n_ch=np.median([((out['HCN'][0]>=a)&(out['HCN'][0]<b)).sum() for a,b in zip(edges[:-1],edges[1:])])
    nh=out['HCN'][2]/np.sqrt(n_ch); no=out['HCO+'][2]/np.sqrt(n_ch)
    m=(cen>2575)&(cen<2875); k=np.sum(H[m]*O[m])/np.sum(O[m]**2)
    for _ in range(20):
        w=1/(nh**2+(k*no)**2); k=np.sum(w*H[m]*O[m])/np.sum(w*O[m]**2)
    res=H-k*O; e=np.sqrt(nh**2+(k*no)**2)
    return pd.DataFrame(dict(v=cen,HCN=H,HCOp=O,resid=res,e=e,sig=res/e)),k
def loader_native(ous):
    def f_(sp,x,y):
        c=vcube(sp,ous,1.05,crop=7. if ous=='Xb4b' else 9.); return c['v'],aper(c,x,y,RAD)*1e3
    return f_
def loader_taper(sp,x,y):
    c=fcube('Xb4d',SPW[sp],None,9.,'cube.I.tapered.image.pbcor','_'); v=vel(c['fr'],F0[sp])
    return v,aper(c,x,y,RAD)*1e3
if __name__=='__main__':
    log=open(OUT+'blend.log','w')
    def pr(*a):
        s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
    from scipy.optimize import curve_fit
    g=lambda v,a,mu,sg: a*np.exp(-0.5*((v-mu)/sg)**2)
    summ=[]
    for pn in POSN:
        r=part_a(pn,'Xb4d')
        pr(f'== (a) {pn} Xb4d 1.05": HCN = k*HCO+ fitted on 2575-2875: k={r["k"]:.3f}+-{r["ek"]:.3f}, chi2_red(fit range)={r["chi2red_fit"]:.2f} (>1: the HCN/HCO+ ratio is not constant even there); 50 km/s-bin noise HCN {r["nh"]:.2f}, HCO+ {r["no"]:.2f} mJy/beam (empirical, includes channel correlation)')
        d=pd.DataFrame(dict(v=r['cen'],HCN=r['H'],HCOp=r['O'],resid=r['res'],e=r['eres'],sig=r['res']/r['eres'])); pr(d.round(2).to_string(index=False))
        d.round(3).to_csv(OUT+f'blend_a_{pn.replace(" ","_")}_Xb4d.csv',index=False)
        cen=r['cen']; m=(cen>=HW[0])&(cen<HW[1]); ex=r['res'][m].sum()*50/1e3; eex=np.sqrt(m.sum())*r['eres'][m].mean()*50/1e3
        hcn=r['H'][m].sum()*50/1e3; pr(f'   HCN flux in {HW}: {hcn:.3f} Jy/beam km/s; residual excess over k*HCO+: {ex:.3f}+-{eex:.3f}  ({ex/hcn:.2f} of HCN)')
        mm=(cen>2150)&(cen<2900)
        p,cv=curve_fit(g,cen[mm],r['res'][mm],p0=[1,2500,80],sigma=r['eres'][mm],absolute_sigma=True,maxfev=5000)
        chi_g=np.sum(((r['res'][mm]-g(cen[mm],*p))/r['eres'][mm])**2); chi_0=np.sum((r['res'][mm]/r['eres'][mm])**2)
        pr(f'   residual Gaussian: amp {p[0]:.2f} mJy/beam, centre {p[1]:.0f}+-{np.sqrt(cv[1,1]):.0f}, sigma {abs(p[2]):.0f}+-{np.sqrt(cv[2,2]):.0f} km/s (FWHM {2.355*abs(p[2]):.0f}); chi2 {chi_g:.1f} vs {chi_0:.1f} (no excess), {mm.sum()} bins. HCO+ line FWHM at the position is ~250-300 km/s, so a single line-width blend and a continuous excess are indistinguishable by shape alone -> need (b)')
        summ.append(dict(pn=pn,k=r['k'],excess_Jykms=ex,e_excess=eex,hcn_win=hcn,res_centre=p[1],res_fwhm=2.355*abs(p[2])))
    pd.DataFrame(summ).round(4).to_csv(OUT+'blend_a_summary.csv',index=False)
    # (c) independent data products
    for tag,ld,edges,fh,fo,order in (('Xb4b smoothed to 1.05"',loader_native('Xb4b'),np.arange(2175,2925,50.),FREE['HCN'],FREE['HCO+'],1),
                                    ('Xb4d tapered (26 MHz channels)',loader_taper,np.arange(2355,3000,60.),[(2950,3213)],[(2239,2300),(2950,3103)],0)):
        for pn in POSN:
            d,k=part_c(pn,tag,ld,edges,fh,fo,order); pr(f'== (c) {pn}, {tag}: k={k:.3f} (fit 2575-2875)'); pr(d.round(2).to_string(index=False))
            d.round(3).to_csv(OUT+f'blend_c_{pn.replace(" ","_")}_{tag.split()[0]}.csv',index=False)
            m=(d.v>=2375)&(d.v<2575); step=d.v.iloc[1]-d.v.iloc[0]
            pr(f'   excess in 2375-2575: {d.resid[m].sum()*step/1e3:.3f}+-{np.sqrt(m.sum())*d.e[m].mean()*step/1e3:.3f} Jy/beam km/s; HCN there {d.HCN[m].sum()*step/1e3:.3f}')
    # (b)
    allb=[]
    for pn in ('Soutflow','M18 South-OF-Blue','S','N'):
        dfb,infl,rho,spec=part_b(pn,pn); allb.append(dfb)
        pr(f'== (b) {pn}: noise autocorr rho(1..6)={np.round(rho,2)} -> variance inflation {infl:.2f}')
    B=pd.concat(allb); B.round(5).to_csv(OUT+'blend_b_stack.csv',index=False)
    pd.set_option('display.width',250)
    for pn in ('Soutflow','M18 South-OF-Blue'):
        b=B[B.pos==pn]; pr(f'-- (b) {pn}: A/eA (best-fit amplitude over its error) per species and temperature; predicted mean blend (2-sigma upper bound) in mJy/beam per HCN-frame 100 km/s bin')
        pr(b[(b.templ=='HCO+prof')][['species','T','nlines','A_over_eA','blend_ul_2275','blend_ul_2375','blend_ul_2475','blend_ul_2575']].round(3).to_string(index=False))
        # envelope: per species max over T and template, summed over species
        env=b.groupby('species')[[f'blend_ul_{lo}' for lo in (2275,2375,2475,2575)]].max().clip(lower=0).sum()
        x,y=POSN[pn]; c=vcube('HCN','Xb4d',1.05,9.); s,_=baseline(c['v'],aper(c,x,y,RAD)*1e3,FREE['HCN'],1)
        for lo in (2275,2375,2475,2575):
            k=(c['v']>=lo)&(c['v']<lo+100); hb=s[k].mean(); pr(f'   HCN-frame {lo}-{lo+100}: HCN mean {hb:.2f} mJy/beam; summed envelope of blend upper bounds {env[f"blend_ul_{lo}"]:.3f} -> blend fraction <= {env[f"blend_ul_{lo}"]/hb:.3f}')
