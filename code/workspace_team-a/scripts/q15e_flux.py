# Q15e: ridge vs blue-lobe wing fluxes, Xb4d vs Xb4b (both convolved to 1.05" circular), HCN and HCO+ 3-2, with a red-edge-anchored linear baseline (Q15d)
# Windows (radio km/s own frame; v0=2730): red +230..+470 (2960-3200), blue mirror -470..-230 (2260-2500); team-b-like windows 150..400 also reported.
# Baseline: linear fit to blue line-free (Q3 FREE) + red edge >3290 km/s. Systematic = spread over {free-blue order0, order1, red-edge anchored order1/order2}.
# Noise: scatter of the same window in random r=1" apertures over the field (>=3.5" from N,S,ridge,lobe,bridge,W clump; core(2575-2875) aperture mean < 0.25 Jy/b km/s).
# Output: out/q15/q15e.log, q15e_table.csv, q15e_maps.npz, q15e_fig.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
from q3_pos import POS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/workspace/out/q15/'; log=open(O+'q15e.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
Sx,Sy=POS['S']; N=POS['N']; LOBE=POS['Soutflow']; RIDGE=(0.3,Sy+3.3)
WRED=(2960,3200); WBLU=(2260,2500); WCORE=(2575,2875)
VAR={'A free-blue order1':(0,1),'B free-blue order0':(0,0),'C blue+edge order1':(1,1),'D blue+edge order2':(1,2)}
def reduce(sp,ous,var):
    c=vcube(sp,ous,1.05,9.); v=c['v']; a=c['a']; edge,od=VAR[var]
    free=list(FREE[sp])+([(3290,v.max())] if edge else [])
    m=np.zeros(len(v),bool)
    for p,q in free: m|=(v>=p)&(v<=q)
    X=np.vander(v-v[m].mean(),od+1); A=a[m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
    r=a-(X@coef).reshape(len(v),*a.shape[1:]); return c,r,v
def apm(M,dx,dy,x,y,rad=1.0): return M[np.hypot(dx-x,dy-y)<=rad].mean()
rows=[]; rng=np.random.default_rng(15); store={}
for ous in ('Xb4d','Xb4b'):
    for sp in ('HCN','HCO+'):
        for var in VAR:
            c,r,v=reduce(sp,ous,var); dx,dy,dv=c['dx'],c['dy'],c['dv']
            Mred=r[(v>=WRED[0])&(v<WRED[1])].sum(0)*dv; Mblu=r[(v>=WBLU[0])&(v<WBLU[1])].sum(0)*dv; Mcore=r[(v>=WCORE[0])&(v<WCORE[1])].sum(0)*dv
            if var.startswith('C'): store[(ous,sp)]=dict(Mred=Mred,Mblu=Mblu,dx=dx,dy=dy)
            ex=[N,POS['S'],RIDGE,LOBE,POS['Bridge'],(-1.8,-2.4),(5.2,-2.4)]
            cands=[]
            for _ in range(6000):
                x,y=rng.uniform(-7.5,7.5,2)
                if np.hypot(x,y)>7.5 or min(np.hypot(x-a_,y-b_) for a_,b_ in ex)<3.5: continue
                if apm(Mcore,dx,dy,x,y)>0.25: continue
                cands.append((x,y))
            cands=cands[:400]
            nr=np.array([apm(Mred,dx,dy,x,y) for x,y in cands]); nb=np.array([apm(Mblu,dx,dy,x,y) for x,y in cands])
            row=dict(ous=ous,line=sp,baseline=var,n_null=len(cands),null_sd_red=nr.std(),null_mean_red=nr.mean(),null_sd_blue=nb.std())
            for nm,(x,y) in (('ridge',RIDGE),('lobe',LOBE),('N',N),('S',POS['S']),('E-clump',(5.2,-2.4))):
                row[nm+'_red']=apm(Mred,dx,dy,x,y); row[nm+'_blue']=apm(Mblu,dx,dy,x,y)
            # mean velocity in 150..500 over ridge/lobe aperture (positive channels weighted), with blue/red sign
            for nm,(x,y),sgn in (('ridge',RIDGE,1),('lobe',LOBE,-1)):
                ap=np.hypot(dx-x,dy-y)<=1.0; s=r[:,ap].mean(1); w=(v>=2730+sgn*150)&(v<=2730+sgn*500) if sgn>0 else (v>=2730-500)&(v<=2730-150)
                row[nm+'_meanv']=(s[w]*(v[w]-2730)).sum()/s[w].sum() if s[w].sum()>0 else np.nan; row[nm+'_F150_500']=s[w].sum()*dv
            rows.append(row)
T=pd.DataFrame(rows); T.round(4).to_csv(O+'q15e_table.csv',index=False)
pd.set_option('display.width',300); pd.set_option('display.max_columns',50)
pr(T[['ous','line','baseline','null_sd_red','ridge_red','ridge_blue','lobe_red','lobe_blue','N_red','N_blue','E-clump_red','ridge_meanv','lobe_meanv','ridge_F150_500','lobe_F150_500']].round(3).to_string(index=False))
# summary with systematics
for ous in ('Xb4d','Xb4b'):
    for sp in ('HCN','HCO+'):
        t=T[(T.ous==ous)&(T.line==sp)&(T.baseline.str[0].isin(list('ABC')))]; c=T[(T.ous==ous)&(T.line==sp)&(T.baseline.str[0]=='C')].iloc[0]
        sd=c.null_sd_red
        pr(f'{ous} {sp}: ridge red(+230..+470) = {c.ridge_red:.3f} +- {sd:.3f} (stat, random apertures) +- {t.ridge_red.max()-t.ridge_red.min():.3f} (baseline range A-C: {t.ridge_red.min():.3f}..{t.ridge_red.max():.3f}); z_stat={c.ridge_red/sd:.1f}, z_incl_sys={c.ridge_red/np.hypot(sd,(t.ridge_red.max()-t.ridge_red.min())/2):.1f}; ridge blue-mirror {c.ridge_blue:.3f}; lobe blue(-470..-230) {c.lobe_blue:.3f}+-{c.null_sd_blue:.3f} (range {t.lobe_blue.min():.3f}..{t.lobe_blue.max():.3f}); ridge/lobe = {c.ridge_red/c.lobe_blue:.2f}; lobe red {c.lobe_red:.3f}; N red {c.N_red:.3f}')
        pr(f'   <v> wings: ridge {c.ridge_meanv:+.0f}, lobe {c.lobe_meanv:+.0f} km/s rel 2730')
# figure: red wing and blue wing HCN, HCO+ maps (Xb4d)
fig,ax=plt.subplots(2,2,figsize=(9,8),sharex=True,sharey=True)
for i,sp in enumerate(('HCN','HCO+')):
    s=store[('Xb4d',sp)]
    for k,(key,tt) in enumerate((('Mblu','blue -470..-230'),('Mred','red +230..+470'))):
        a=ax[i,k]; im=a.imshow(s[key],origin='lower',extent=[s['dx'].min(),s['dx'].max(),s['dy'].min(),s['dy'].max()] if False else [s['dx'][0,0],s['dx'][0,-1],s['dy'][0,0],s['dy'][-1,0]],vmin=-0.05,vmax=0.3,cmap='viridis')
        a.invert_xaxis(); a.set_title(f'{sp} 3-2 {tt} km/s (Xb4d, 1.05")'); a.plot(*N,'w^'); a.plot(*POS['S'],'wv'); a.plot(*RIDGE,'rx'); a.plot(*LOBE,'rx'); a.plot(*POS['M18 OF-Red'],'mo',mfc='none'); a.set_xlim(4,-4); a.set_ylim(-8,6)
        plt.colorbar(im,ax=a,label='Jy/beam km/s')
fig.supxlabel('dRA (arcsec)'); fig.supylabel('dDec (arcsec)'); plt.tight_layout(); plt.savefig(O+'q15e_fig.png',dpi=110)
