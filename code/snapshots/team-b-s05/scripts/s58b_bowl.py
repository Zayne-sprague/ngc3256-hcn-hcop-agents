# Q12 round 2, step b: empirical filtering of the three ACA sets from their own maps. For each array/line (and the line-free continuum) at a FIXED 10" circular beam:
#  - radial profile and curve of growth of the PB-UNcorrected map (pbcor map x fitted Gaussian PB) about the N-S midpoint out to the full map (circle r<=25" common to all
#    three arrays; larger circles and the full square for Xb4f/Xb55 separately);
#  - bowl depth (annulus mean of the pbcor map, minimum beyond R0) and depth/peak, bowl radius, zero-crossing radius R0, F_in=cumulative flux to R0 (the positive source),
#    F_edge = cumulative flux to the edge, F_neg = F_edge - F_in (integrated bowl).  Zero-spacing identity: a TP-less dirty image integrates to 0 => F_edge/F_in -> 0.
#    A cleaned image does not obey it exactly (clean components carry flux), so F_neg/F_in is an empirical *index* of the filtering, compared between the arrays.
# Errors: 16 sign-flipped line-free noise realisations pushed through the same pipeline.  Fluxes: Jy km/s = sum(Jy/beam km/s)*pixarea/beam_area(10").
# usage: python3 -I s58b_bowl.py -> out/q12b/bowl.json, bowl.txt, profiles.csv
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from common import *
import s56_lib as L
CEN=L.CEN; B10=(10.,10.,0.)
RPB={'Xb55':56.1,'Xb4f':37.5,'Xb47':28.2}      # mean fitted PB FWHM (s58a) of HCO+/HCN; 7-m expectation 56.9/37.6/28.4
RMAX={'Xb55':50.,'Xb4f':33.,'Xb47':25.}
log=[]
def P(s=''): print(s,flush=True); log.append(s)
def setup(ous,ps,wc,n):
    Y,X=np.mgrid[:n,:n]; cx,cy=wc.world_to_pixel_values(*CEN); rr=np.hypot(X-cx,Y-cy)*ps
    ra,dec=wc.pixel_to_world_values(X,Y); dx=(ra-156.9635)*np.cos(np.deg2rad(-43.9046))*3600; dy=(dec+43.9046)*3600
    pb=np.exp(-4*np.log(2)*(dx**2+dy**2)/RPB[ous]**2); return rr,pb
def metrics(m,rr,pb,ps,b0,tgt,rmax):
    """m: Jy/beam km/s map at the target beam (pbcor).  Returns dict"""
    area=1.1331*tgt[0]*tgt[1]/ps**2; u=m*pb; edges=np.arange(0,rmax+1e-6,max(1.0,ps)); rc=0.5*(edges[1:]+edges[:-1])
    prof=np.array([m[(rr>=a)&(rr<b)].mean() for a,b in zip(edges[:-1],edges[1:])]); uprof=np.array([u[(rr>=a)&(rr<b)].mean() for a,b in zip(edges[:-1],edges[1:])])
    cog=np.array([u[rr<=r].sum()/area for r in edges[1:]]); pk=m[rr<4].max()
    # zero crossing of the PB-uncorrected annular mean profile (first annulus beyond the peak with mean<0)
    neg=np.where(uprof<0)[0]; i0=neg[0] if len(neg) else len(uprof)-1; R0=edges[i0]
    F_in=cog[max(i0-1,0)]; F_edge=cog[-1]; F_max=cog.max(); Rmax=edges[1:][cog.argmax()]
    # bowl: minimum of the pbcor annular mean beyond R0
    j=i0+np.argmin(prof[i0:]) if i0<len(prof) else len(prof)-1
    return dict(peak=float(pk),R0=float(R0),F_in=float(F_in),F_edge=float(F_edge),F_max=float(F_max),R_Fmax=float(Rmax),F_neg=float(F_edge-F_max),bowl=float(prof[j]),R_bowl=float(rc[j]),bowl_over_peak=float(prof[j]/pk),
                bowl_unc=float(uprof[j]),prof=prof.tolist(),uprof=uprof.tolist(),cog=cog.tolist(),r=rc.tolist())
OUT={}
for line in ('HCO+','HCN'):
    for ous in ('Xb55','Xb4f','Xb47'):
        d=np.load(f'/workspace/out/q12b/full_{ous}_{line}.npz'); ps=float(d['ps']); b0=tuple(d['beam']); J={'Xb55':2,'Xb4f':3,'Xb47':4}[ous]
        spw=dict(Xb55={'HCN':'spw16','HCO+':'spw18'},Xb4f={'HCO+':'spw16','HCN':'spw18'},Xb47={'HCN':'spw16','HCO+':'spw18'})[ous][line]
        wc=WCS(fits.getheader(files(f'*{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0])).celestial
        n=d['m0'].shape[0]; rr,pb=setup(ous,ps,wc,n)
        M=L.conv(np.concatenate([d['m0'][None],d['noise']]),ps,b0,B10)
        for tag,rmax in (('r25',25.),('rfull',RMAX[ous])):
            r=metrics(M[0],rr,pb,ps,b0,B10,rmax); nz=[metrics(q,rr,pb,ps,b0,B10,rmax) for q in M[1:]]
            area=1.1331*100/ps**2; sq=float((M[0]*pb).sum()/area); sqe=float(np.std([(q*pb).sum()/area for q in M[1:]]))
            r['sq_edge']=sq; r['err_sq']=sqe
            for k in ('F_in','F_edge','F_max','F_neg','bowl','peak'): r['e_'+k]=float(np.std([x[k] for x in nz]))
            OUT[f'{line}|{ous}|{J}|{tag}']=r
            P(f'{line:5s} {ous} J={J} 10" {tag:5s}: peak {r["peak"]:.1f} Jy/b km/s; R0={r["R0"]:.0f}" F_in={r["F_in"]:.1f} F_max={r["F_max"]:.1f}@{r["R_Fmax"]:.0f}" F_edge={r["F_edge"]:.1f}+-{r["e_F_edge"]:.1f}  F_neg(edge-max)={r["F_neg"]:.1f}+-{r["e_F_neg"]:.1f} '
              f'[F_neg/F_max={r["F_neg"]/r["F_max"]:.2f}, F_edge/F_max={r["F_edge"]/r["F_max"]:.2f}]  bowl={r["bowl"]:.2f}+-{r["e_bowl"]:.2f} at {r["R_bowl"]:.0f}" (bowl/peak {r["bowl_over_peak"]:.3f}); square-field integral {sq:.1f}+-{sqe:.1f}')
# continuum at 10" (Xb55,Xb4f,Xb47 and the 12-m Xb4d for reference): same metrics, line-free continuum maps in Jy/beam
P('\n## continuum (line-free spw16_18_20_22 cont.I.pbcor) at 10"')
for ous in ('Xb55','Xb4f','Xb47'):
    d=np.load(f'/workspace/out/q12b/cont_{ous}.npz'); ps=float(d['ps']); b0=tuple(d['beam']); n=d['img'].shape[0]
    wc=WCS(fits.getheader(files(f'*{ous}.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits')[0])).celestial
    rr,pb=setup(ous,ps,wc,n); M=L.conv(d['img'],ps,b0,B10)[0]
    for tag,rmax in (('r25',25.),):
        r=metrics(M,rr,pb,ps,b0,B10,rmax); area=1.1331*100/ps**2
        OUT[f'cont|{ous}|{tag}']=r
        P(f'cont {ous}: peak {r["peak"]*1e3:.1f} mJy/b; R0={r["R0"]:.0f}" F_in={r["F_in"]*1e3:.1f} F_max={r["F_max"]*1e3:.1f}@{r["R_Fmax"]:.0f}" F_edge={r["F_edge"]*1e3:.1f} mJy; F_neg/F_max={r["F_neg"]/r["F_max"]:.2f}; bowl/peak {r["bowl_over_peak"]:.3f} at {r["R_bowl"]:.0f}"')
json.dump(OUT,open('/workspace/out/q12b/bowl.json','w')); open('/workspace/out/q12b/bowl.txt','w').write('\n'.join(log))
