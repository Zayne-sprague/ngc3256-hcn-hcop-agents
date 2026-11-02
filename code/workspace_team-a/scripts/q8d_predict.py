# Q8d: predictions of the uv forward model (posterior chains from q8c_uvmodel.py): recovery of Xb4d/Xb4b/Xb47 relative to Xb4f at N,S,Bridge,
# posterior-predictive check of the calibration data, Xb47-vs-Xb4f recovery of the 4-3 / 3-2 peak and the correction to the N/S T43/T32 contrast.
# r_c (compact-component T43/T32) ~ U(0.4,1) prior; r_e (extended) solved so the model reproduces the observed ratio at the matched beam. 4-3 sizes = q x 3-2 sizes, q in {1,0.7}.
import sys,importlib,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
OUT='/workspace/out/q8c/'
VARS=sys.argv[1:] or ['taper_phys','taper_raw_phys','sharp_phys','taper','taper_raw','sharp']
OBS={'HCN':{'N':0.425,'S':0.319,'Bridge':0.364},'HCO+':{'N':0.547,'S':0.483,'Bridge':0.4}}   # T43/T32 (Tb) at the Xb4f beam, Xb4f 3-2 (out/q8/T43_T32_contrast.csv)
BRIDGE_HCOP=None
rng=np.random.default_rng(7); rows=[]; pp=[]; ctab=[]; FAC={}; CMP={}
for V in VARS:
    sys.argv=['x',V]; import q8c_uvmodel as M; importlib.reload(M)
    Z=np.load(OUT+f'chain_{V}.npz',allow_pickle=True); flat=Z['flat']; lp=Z['lp']
    idx=rng.choice(len(flat),4000,replace=False); th=flat[idx]; ix=M.ix; peak=M.peak; unf=M.unf
    BM=(6.2,3.7)
    def LL(a,cont=False):
        if cont and a=='47': return M.arrL(th,'f',True)*(346.5/257.7*0+ (1/346.5)/(1/257.7))/np.exp(th[:,ix['lnrB47']])
        return M.arrL(th,a,cont)
    # ---- posterior predictive: recoveries of the 3-2 peak at the Xb4f beam (filtered/unfiltered) per array
    rec={}
    for reg in M.REG:
        for l in M.LIN:
            s,se,Ac,Ae=M.comp(th,reg,l); U=unf(Ac,s,BM)+unf(Ae,se,BM)
            for a in ('f','d','b'):
                rec[(reg,l,a)]=(peak(Ac,s,BM,LL(a),M.NATB[a])+peak(Ae,se,BM,LL(a),M.NATB[a]))/U
            for q in (1.0,0.7):
                U4=unf(Ac,s*q,BM)+unf(Ae,se,BM)
                rec[(reg,l,'47',q)]=(peak(Ac,s*q,BM,LL('47'),M.NATB['47'])+peak(Ae,se,BM,LL('47'),M.NATB['47']))/U4
    for reg in M.REG:
        for l in M.LIN:
            f=rec[(reg,l,'f')]
            for a in ('d','b'):
                r=rec[(reg,l,a)]/f*M.gscale(th,a); pp.append(dict(variant=V,region=reg,line=l,quantity=f'pred {a}/f at 6.2x3.7',p16=np.percentile(r,16),p50=np.median(r),p84=np.percentile(r,84)))
            for q in (1.0,0.7):
                pp.append(dict(variant=V,region=reg,line=l,quantity=f'rec Xb47(4-3,q={q})/rec Xb4f(3-2)',p16=np.percentile(rec[(reg,l,'47',q)]/f,16),p50=np.median(rec[(reg,l,'47',q)]/f),p84=np.percentile(rec[(reg,l,'47',q)]/f,84)))
            pp.append(dict(variant=V,region=reg,line=l,quantity='rec Xb4f (3-2 peak, abs.)',p16=np.percentile(f,16),p50=np.median(f),p84=np.percentile(f,84)))
    cont=pd.read_csv('/workspace/out/q8/continuum_matched.csv'); cp=cont[cont.beam=='Xb4f(6.2x3.7,82)'].pivot_table(index='array',columns='region',values='mJy_beam')
    ALPHA={}
    for reg in M.REG:
        s_,se_,rE=[np.exp(th[:,ix[f'cont_{reg}_{k}']]) for k in ('lns','lnse','lnrE')]
        U=unf(1.,s_,BM)+unf(rE,se_,BM); rf=(peak(1.,s_,BM,LL('f',True),M.NATB['f'])+peak(rE,se_,BM,LL('f',True),M.NATB['f']))/U
        r47=(peak(1.,s_,BM,LL('47',True),M.NATB['47'])+peak(rE,se_,BM,LL('47',True),M.NATB['47']))/U
        ALPHA[reg]=np.log(cp.loc['Xb47',reg]/cp.loc['Xb4f',reg]/(r47/rf))/np.log(346.5/257.7)
    W=np.ones(len(th))
    for reg in M.REG: W*=((ALPHA[reg]>1.5)&(ALPHA[reg]<4.5))
    ess=W.sum(); rows.append(dict(variant=V,param='alpha-bound weight: fraction of draws with 1.5<alpha<4.5 at N,S,Bridge',p16=ess/len(th),p50=np.nan,p84=np.nan))
    # ---- correction to N/S: equal-r (r_c=r_e) and conditional r_c scan
    for l in M.LIN:
        for q in (1.0,0.7):
            c={reg:rec[(reg,l,'f')]/rec[(reg,l,'47',q)] for reg in ('N','S','Bridge')}   # multiply observed T43/T32 by c to remove the loss (r_c=r_e)
            CNS=c['N']/c['S']; obsNS=OBS[l]['N']/OBS[l]['S']
            ctab.append(dict(variant=V,line=l,q43=q,mode='r_c=r_e',C_NS_p16=np.percentile(CNS,16),C_NS_p50=np.median(CNS),C_NS_p84=np.percentile(CNS,84),C_NS_p2_5=np.percentile(CNS,2.5),C_NS_p97_5=np.percentile(CNS,97.5),
                 c_N=np.median(c['N']),c_S=np.median(c['S']),c_B=np.median(c['Bridge']),obs_NS=obsNS,corr_NS_p16=obsNS*np.percentile(CNS,16),corr_NS_p50=obsNS*np.median(CNS),corr_NS_p84=obsNS*np.percentile(CNS,84)))
            # compact/extended split: r_c prior, r_e solved
            rc=rng.uniform(0.4,1.0,len(th)); tr={}; ok=np.ones(len(th),bool); re_={}; piece={}
            for reg in ('N','S'):
                s,se,Ac,Ae=M.comp(th,reg,l); Lf=LL('f'); L47=LL('47')
                den=peak(Ac,s,BM,Lf,M.NATB['f'])+peak(Ae,se,BM,Lf,M.NATB['f'])
                numc=peak(Ac,s*q,BM,L47,M.NATB['47']); nume=peak(Ae,se,BM,L47,M.NATB['47'])
                piece[reg]=(den,numc,nume,unf(Ac,s*q,BM),unf(Ae,se,BM),unf(Ac,s,BM)+unf(Ae,se,BM)); re=(OBS[l][reg]*den-rc*numc)/nume; ok&=(re>0.03)&(re<1.3); re_[reg]=re
                true=(rc*unf(Ac,s*q,BM)+re*unf(Ae,se,BM))/(unf(Ac,s,BM)+unf(Ae,se,BM)); tr[reg]=true
            for re0 in (0.30,0.36,0.45):   # extended gas with the SAME T43/T32 at N and S: what the compact component must do
                rcN=(OBS[l]['N']*piece['N'][0]-re0*piece['N'][2])/piece['N'][1]; rcS=(OBS[l]['S']*piece['S'][0]-re0*piece['S'][2])/piece['S'][1]
                good=(rcN>0)&(rcS>0)&(W>0)
                if good.sum()>20: CMP[(V,l,q,re0)]=(rcN[good],rcS[good],np.mean((piece['N'][3]/(piece['N'][3]+piece['N'][4]))[good]),np.mean((piece['S'][3]/(piece['S'][3]+piece['S'][4]))[good]))
            if ok.sum()>50:
                ratio=(tr['N']/tr['S'])[ok]
                FAC[(V,l,q)]=(((tr['N']/tr['S'])/(OBS[l]['N']/OBS[l]['S']))[ok],W[ok],rc[ok],re_['N'][ok],re_['S'][ok])
                ctab.append(dict(variant=V,line=l,q43=q,mode='r_c~U(0.4,1), r_e solved',frac_valid=ok.mean(),alpha_ok_frac=W[ok].mean(),corr_NS_p16=np.percentile(ratio,16),corr_NS_p50=np.median(ratio),corr_NS_p84=np.percentile(ratio,84),corr_NS_p2_5=np.percentile(ratio,2.5),corr_NS_p97_5=np.percentile(ratio,97.5),obs_NS=obsNS,
                    reN_p16=np.percentile(re_['N'][ok],16),reN_p50=np.median(re_['N'][ok]),reN_p84=np.percentile(re_['N'][ok],84),reS_p16=np.percentile(re_['S'][ok],16),reS_p50=np.median(re_['S'][ok]),reS_p84=np.percentile(re_['S'][ok],84),
                    trueN_p50=np.median(tr['N'][ok]),trueS_p50=np.median(tr['S'][ok]),rc_p50=np.median(rc[ok])))
    # ---- array parameters
    for n in ('lnLf','lnLd','lnLb','lngd','lngb','lnrB47','lnrB55'):
        v=np.exp(flat[:,ix[n]]); rows.append(dict(variant=V,param=n[2:],p16=np.percentile(v,16),p50=np.median(v),p84=np.percentile(v,84)))
    Lf=np.exp(flat[:,ix['lnLf']]); rows.append(dict(variant=V,param='L47 (derived)',p16=np.percentile(Lf*0.75/np.exp(flat[:,ix['lnrB47']]),16),p50=np.median(Lf*0.75/np.exp(flat[:,ix['lnrB47']])),p84=np.percentile(Lf*0.75/np.exp(flat[:,ix['lnrB47']]),84)))
    rows.append(dict(variant=V,param='chi2_min / ndata',p16=-2*lp.max(),p50=len(M.NAT)+len(M.RAT)+len(M.CRAT),p84=np.nan))
    # ---- continuum: implied spectral index from the Xb47/Xb4f continuum ratio after the model recovery ratio
    cont=pd.read_csv('/workspace/out/q8/continuum_matched.csv'); cp=cont[cont.beam=='Xb4f(6.2x3.7,82)'].pivot_table(index='array',columns='region',values='mJy_beam')
    for reg in M.REG:
        s,se,rE=[np.exp(th[:,ix[f'cont_{reg}_{k}']]) for k in ('lns','lnse','lnrE')]
        U=unf(1.,s,BM)+unf(rE,se,BM); rf=(peak(1.,s,BM,LL('f',True),M.NATB['f'])+peak(rE,se,BM,LL('f',True),M.NATB['f']))/U
        r47=(peak(1.,s,BM,LL('47',True),M.NATB['47'])+peak(rE,se,BM,LL('47',True),M.NATB['47']))/U
        obs=cp.loc['Xb47',reg]/cp.loc['Xb4f',reg]; al=np.log(obs/(r47/rf))/np.log(346.5/257.7)
        pp.append(dict(variant=V,region=reg,line='continuum',quantity='rec Xb47/rec Xb4f (cont.)',p16=np.percentile(r47/rf,16),p50=np.median(r47/rf),p84=np.percentile(r47/rf,84)))
        pp.append(dict(variant=V,region=reg,line='continuum',quantity='implied spectral index alpha (Xb47/Xb4f, 258-347 GHz)',p16=np.percentile(al,16),p50=np.median(al),p84=np.percentile(al,84)))
    print(V,'done',flush=True)
import pickle; pickle.dump((FAC,CMP),open(OUT+'fac_samples.pkl','wb'))
A=pd.DataFrame(rows); P=pd.DataFrame(pp); Cc=pd.DataFrame(ctab)
A.to_csv(OUT+'array_params.csv',index=False); P.to_csv(OUT+'predictions.csv',index=False); Cc.to_csv(OUT+'NS_correction.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500,'display.max_columns',40)
print(A.round(3).to_string()); print(P[P.quantity.str.contains('rec|alpha')].round(3).to_string()); print(Cc.round(3).to_string())
