# Q9b library: smoothed (cached) HCN/HCO+ 3-2 cubes, aperture spectra, window-integrated fluxes with null-based errors
# (random positions x sliding line-free windows) + baseline-fit variance, and ratio with censored-limit handling.
import sys,os; sys.path.insert(0,'/workspace/scripts')
from q7_lib import *
CD='/workspace/out/q9b/cache/'
TBF=(265.886434/267.557626)**-2      # Tb_HCN/Tb_HCO+ = (S_HCN/S_HCO+)*(nu_O/nu_H)^2 = 1.0126*Jy ratio ... (nu_O/nu_H)^2
REG={'N':(0.0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
_cubes={}
def cube(ous,ln):
    if (ous,ln) not in _cubes: _cubes[(ous,ln)]=Cube(ous,ln)
    return _cubes[(ous,ln)]
def smcube(ous,ln,tb):
    """baseline-subtracted cube convolved to circular beam tb (Jy/beam(tb)); cached .npy; tb=None -> native"""
    c=cube(ous,ln)
    if tb is None: return c.cube
    fn=CD+f'{ous}_{ln}_{tb[0]:.2f}.npy'
    if os.path.exists(fn): return np.load(fn)
    out=np.empty_like(c.cube)
    for i in range(c.cube.shape[0]): out[i]=smooth(c.cube[i],c,tb)
    np.save(fn,out); return out
def centre(ous,tb,reg,line='HCO+'):
    c=cube(ous,line); s=smcube(ous,line,tb); img=s[(c.v>=-350)&(c.v<=350)].sum(0)*c.dv
    return peak_near(c,img,*REG[reg])
def apmask(c,x0,y0,r1,r2):
    r=np.hypot(c.X-x0,c.Y-y0); return (r>=r1)&(r<r2)
def apspec(c,s,m,tb):
    """Jy per channel summed in aperture mask (flux density)"""
    om=OM(tb) if tb is not None else OM(c.beam)
    return s[:,m].sum(1)/(om/c.pix**2)
def wweights(v,dv,lo,hi):
    """fractional channel overlap weights with window [lo,hi)"""
    a=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,None); return a/dv
def winflux(c,spec,lo,hi):
    return (spec*wweights(c.v,c.dv,lo,hi)).sum()*c.dv     # Jy km/s
def nullstats(c,s,tb,r1,r2,widths,rng,ndraw=160,box=5.0):
    """per-channel and window-sum noise from random aperture positions; windows slide over line-free channels (|v|>550, contiguous).
    returns dict width(km/s)->sigma(Jy km/s), and sigma_chan (Jy)"""
    free=np.where(np.abs(c.v)>550)[0]; seg=np.split(free,np.where(np.diff(free)>1)[0]+1)
    segs=[g for g in seg if len(g)>3]
    sp=[]
    cx,cy=0.0,0.0
    while len(sp)<ndraw:
        x,y=rng.uniform(-box,box,2); m=apmask(c,x,y,r1,r2)
        if m.sum()<1: continue
        sp.append(apspec(c,s,m,tb))
    sp=np.array(sp)
    out={}
    for w in widths:
        n=max(1,int(round(w/c.dv))); vals=[]
        for g in segs:
            if len(g)<n: continue
            cs=np.cumsum(np.pad(sp[:,g],((0,0),(1,0))),axis=1); vals.append((cs[:,n:]-cs[:,:-n]).ravel())
        vals=np.concatenate(vals)*c.dv*(w/(n*c.dv))   # rescale to the nominal width
        out[w]=1.4826*np.median(np.abs(vals-np.median(vals)))
    schan=1.4826*np.median(np.abs(sp[:,free]-np.median(sp[:,free])))
    return out,schan
def basevar(c,lo,hi):
    """variance factor (in units of channel noise^2 * dv^2) of the linear-baseline-fit contribution to a window sum over [lo,hi): a^T (X^T X)^-1 a, with a = sum of weighted design rows"""
    X=np.vstack([np.ones(c.bm.sum()),c.v[c.bm]/500.]).T; G=np.linalg.inv(X.T@X)
    w=wweights(c.v,c.dv,lo,hi); A=np.vstack([np.ones(len(c.v)),c.v/500.]).T; a=(A*w[:,None]).sum(0); return a@G@a
def window_table(ous,tb,reg,r1,r2,wins,rng,ndraw=160):
    """HCN & HCO+ window fluxes (Jy km/s) with errors (null + baseline fit) in an aperture about the HCO+ peak. wins: list of (lo,hi)"""
    x0,y0=centre(ous,tb,reg); rows=[]; res={}
    for ln in ('HCN','HCO+'):
        c=cube(ous,ln); s=smcube(ous,ln,tb); m=apmask(c,x0,y0,r1,r2); sp=apspec(c,s,m,tb)
        ns,sch=nullstats(c,s,tb,r1,r2,sorted({round(b-a) for a,b in wins}),rng,ndraw)
        for lo,hi in wins:
            F=winflux(c,sp,lo,hi); e0=ns[round(hi-lo)]; eb=sch*c.dv*np.sqrt(basevar(c,lo,hi))
            res[(ln,lo,hi)]=(F,np.hypot(e0,eb),eb)
    return (x0,y0),res
def ratio(A,eA,B,eB,nsig=3.,lim=2.):
    """HCN/HCO+ Jy-flux ratio with censored handling. returns (R,eR,kind,limit) kind: 'det','lo' (HCN-only),'hi' (HCO+-only),'none'"""
    dA,dB=A>nsig*eA,B>nsig*eB
    R=A/B if B!=0 else np.nan; eR=abs(R)*np.hypot(eA/A,eB/B) if (A!=0 and B!=0) else np.nan
    if dA and dB: return R,eR,'det',np.nan
    if dA and not dB: return R,eR,'lo',(A-lim*eA)/(max(B,0)+lim*eB)
    if dB and not dA: return R,eR,'hi',(max(A,0)+lim*eA)/(B-lim*eB)
    return R,eR,'none',np.nan
