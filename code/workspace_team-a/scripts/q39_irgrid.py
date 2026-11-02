# Q39/Q24: HCN excitation with IR pumping through the nu2=1 bending mode (14.0 um, 712 cm^-1).  Extends the q36_rt escape-probability solver with
# a second vibrational ladder (nu2=1, J'=1..JV; l-type doublets lumped, g=6(2J'+1) incl. nuclear-spin 3).  Rotational A and collision rates of nu2=1 are copied from nu=0
# (standard approximation); vibrational A: A0=1.7 s^-1 band total (Lahuis/Boonman: A(14um) ~1.7 s^-1), split P/Q/R with Honl-London factors S_R=(J+2)/2, S_Q=(2J+1)/2, S_P=(J-1)/2
# (J=lower), A(J'->J)=A0*S(J->J')/(2J'+1)*(nu/nu0)^3.  The 14um radiation field is J_nu = W*B_nu(T_IR) (photon occupation W/(exp(hnu/kT_IR)-1)); T_IR=0 -> pumping off (reproduces q36_rt).
# The rotational background in the cores is the q36 J_bg(nu) of the chosen scenario.  usage: q39_irgrid.py <TIR K> <W> <tag>  -> out/q39/grid_hcn_{N,S}_{scen}_<tag>.npz ; q39_irgrid.py test
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, q36_rt as R, q36_bg as B
from q36_rt import H,KB,C
from concurrent.futures import ProcessPoolExecutor
TK=np.array([10,15,20,25,30,40,50,65,80,100,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25)
NU0=712.0*C; A0=1.7; JV=14
def extend(M):
    nl=M['nl']; Ev=712.0+M['E'][1:JV+1]   # nu2=1 J'=1..JV (index J'): rotational ladder shifted (B' ~ B)
    E=np.concatenate([M['E'],Ev]); g=np.concatenate([M['g'],6*(2*np.arange(1,JV+1)+1.)])
    up=list(M['up']); lo=list(M['lo']); A=list(M['A']); nu=list(M['nu'])
    for J in range(1,JV):    # rotational J+1->J in nu2=1 (same A, nu as v=0)
        k=J   # LAMDA transition k: level k+1 -> k (J=k)  ; v=0 transition J+1->J is index J
        up.append(nl+J); lo.append(nl+J-1); A.append(M['A'][J]); nu.append(M['nu'][J])
    for Jp in range(1,JV+1):
        for J,S in ((Jp-1,(Jp+1)/2),(Jp,(2*Jp+1)/2),(Jp+1,Jp/2)):
            if S<=0 or J>25: continue
            f=(E[nl+Jp-1]-E[J])*C
            up.append(nl+Jp-1); lo.append(J); A.append(A0*S/(2*Jp+1)*(f/NU0)**3); nu.append(f)
    # collisions: copy v=0 rates to v2=1 (both levels in ladder)
    cu=list(M['cu']); cl=list(M['cl']); K=[r for r in M['K']]
    for u,l,r in zip(M['cu'],M['cl'],M['K']):
        if 1<=u<=JV and 1<=l<=JV: cu.append(nl+u-1); cl.append(nl+l-1); K.append(r)   # v=0 level index u ~ J=u (levels 1..), J=u>=1 maps to nl+u-1
    return dict(nl=nl+JV,E=E,g=g,up=np.array(up),lo=np.array(lo),A=np.array(A),nu=np.array(nu),cT=M['cT'],cu=np.array(cu),cl=np.array(cl),K=np.array(K))
def solve(M,T,n,N,dv,jbg_rot,TIR,Wd,niter=2000):
    M2=M
    def jb(f):
        f=np.atleast_1d(f); out=np.where(f<2e12,jbg_rot(f),0.)
        if TIR>0: out=np.where(f>=2e12,Wd*R.Jp(np.maximum(f,1e12),TIR),out)
        return out
    return R.solve(M2,T,n,N,dv,jb,niter)
def job(a):
    reg,scen,TIR,Wd=a; M=extend(R.parse('/opt/lamda/hcn.dat')); f=(lambda nu: B.Jbg(reg,scen,nu/1e9))
    out=np.full((len(TK),len(LN),len(LNC),4,4),np.nan)
    for i,T in enumerate(TK):
        for j,ln in enumerate(LN):
            for k,lN in enumerate(LNC):
                o=solve(M,T,10**ln,10**lN,150.,f,TIR,Wd); out[i,j,k,:,0]=o['TR'][:4]; out[i,j,k,:,1]=o['tau'][:4]; out[i,j,k,:,2]=o['Tex'][:4]; out[i,j,k,:,3]=o['x'][26:].sum()/o['x'][:26].sum()
    return a,out
if __name__=='__main__':
    if sys.argv[1]=='test':
        M=extend(R.parse('/opt/lamda/hcn.dat')); M0=R.parse('/opt/lamda/hcn.dat'); f=lambda nu: B.Jbg('N','cmb',nu/1e9)
        print('vib A (J1->J0,J1->J1,J1->J2):',[(round(a,3),round(n/C,1)) for u,l,a,n in zip(M['up'],M['lo'],M['A'],M['nu']) if u==M['nl']-14]);
        import time
        for (T,n,N) in ((30,1e5,1e15),(30,1e6,1e17),(100,1e7,1e17),(100,1e5,1e16)):
            t=time.time(); a=R.solve(M0,T,n,N,150.,f); b=solve(M,T,n,N,150.,f,0,1); c=solve(M,T,n,N,150.,f,150.,1.); d=solve(M,T,n,N,150.,f,300.,1.)
            print(T,n,N,'TR[1-0..4-3] v0 only',np.round(a['TR'][:4],4),'| ext, no IR',np.round(b['TR'][:4],4),'| TIR150',np.round(c['TR'][:4],4),'| TIR300',np.round(d['TR'][:4],4),'%.1fs'%(time.time()-t))
        sys.exit()
    TIR=float(sys.argv[1]); Wd=float(sys.argv[2]); tag=sys.argv[3]; scens=sys.argv[4].split(',') if len(sys.argv)>4 else ['cmb','nominal','dusty','high']
    jobs=[(r,s,TIR,Wd) for r in 'NS' for s in scens]
    with ProcessPoolExecutor(8) as ex:
        for a,out in ex.map(job,jobs):
            np.savez(f'/workspace/out/q39/grid_hcn_{a[0]}_{a[1]}_{tag}.npz',TR=out[...,0],tau=out[...,1],Tex=out[...,2],xvib=out[...,0,3],TK=TK,LN=LN,LNC=LNC); print('done',a,flush=True)
