# Q12 step c: Fourier-domain transfer-function model of the nested Band 6 set, fitted to whole maps (not single positions).
# Model: an array a sees the sky with a smooth uv hole, T_a(q)=1-exp(-(q/q_a)^p), flux scale g_a (g_Xb4f=1 by definition: all recoveries are RELATIVE to Xb4f).
# At a common circular beam B (6.5", 10"): I_d = g_d F^-1[ T_d/T_f F I_f ],  I_b = g_b F^-1[ T_b/T_f F I_f ];  and for the finer beams (2,3,4.5,6.5,10"): I_b = (g_b/g_d) F^-1[ T_b/T_d F I_d ].
# The effective hole radius is quoted as a baseline L_a=q_a*lambda (m) for HCO+/HCN 3-2 (lambda=1.128 mm). Xb4f is a 7-m (ACA) data set (see s55a PB: noise-map FWHM 37"), hole fixed by scenario L_f (not
# measurable from the nested set); p in {1,2,4} (shape of the uv taper). Pixels: r<15" of the N-S midpoint on the 0.7" Xb4f grid; sigma^2 = noise_a^2 + (r*noise_ref)^2 + (4% I)^2 + (0.3 Jy/b km/s)^2.
# usage: python3 -I s56c_transfer_fit.py -> out/q12/transfer_fit.json/.txt, out/q12/transfer_pred.json
import sys,json; sys.path.insert(0,'/workspace/scripts')
import numpy as np, s56_lib as L
from scipy.optimize import least_squares
PAD=128
def filt(I,R,ps):
    n=I.shape[0]; a=np.zeros((PAD,PAD)); a[:n,:n]=I; fy=np.fft.fftfreq(PAD,ps); q=np.hypot(*np.meshgrid(fy,fy)); return np.fft.ifft2(np.fft.fft2(a)*R(q)).real[:n,:n],q
def T(q,qa,p): return -np.expm1(-(np.maximum(q,1e-9)/qa)**p)
def build(line,beams_f=((6.5,6.5,0),(10.,10.,0)),beams_d=((2.,2.,0),(3.,3.,0),(4.5,4.5,0),(6.5,6.5,0),(10.,10.,0))):
    ra,dec,ps=L.common_grid(line); n=ra.shape[0]
    A={o:L.Arr(o,line) for o in ('Xb4f','Xb4d','Xb4b')}
    cx,cy=A['Xb4f'].sky2idx(*L.CEN); Y,X=np.mgrid[:n,:n]; rr=np.hypot(X-cx,Y-cy)*ps; mask=rr<15
    G={}
    for t in set(beams_f)|set(beams_d):
        for o in A:
            if o=='Xb4f' and t[0]<5.8: continue
            G[(o,t[0])]=L.maps_on_grid(o,line,t,A[o])[:2]
    return G,mask,ps,beams_f,beams_d
def make_resid(G,mask,ps,beams_f,beams_d,Lf,p,line,detail=False):
    qf=L.L2q(Lf)
    def sig(o,B,pred_from=None):
        M,N=G[(o,B)]; s=N.std(0)
        return M,s
    def resid(th):
        qd,qb=np.exp(th[0:2]); gd,gb=th[2:4]; out=[]; det={}
        for t in beams_f:
            B=t[0]; If,sf=G[('Xb4f',B)][0],G[('Xb4f',B)][1].std(0)
            for o,qa,ga in (('Xb4d',qd,gd),('Xb4b',qb,gb)):
                pred,_=filt(If,lambda q:ga*T(q,qa,p)/T(q,qf,p),ps); M,N=G[(o,B)]; s=np.sqrt(N.std(0)**2+(pred.sum()*0+ (ga*0))**2+(abs(pred)*0.04)**2+0.3**2+ (ga*sf)**2)
                out.append(((M-pred)/s)[mask]); det[(o+'<-Xb4f',B)]=(pred,M)
        for t in beams_d:
            B=t[0]; Id=G[('Xb4d',B)][0]; sd=G[('Xb4d',B)][1].std(0)
            pred,_=filt(Id,lambda q:(gb/gd)*T(q,qb,p)/T(q,qd,p),ps); M,N=G[('Xb4b',B)]; s=np.sqrt(N.std(0)**2+(abs(pred)*0.04)**2+0.3**2+sd**2)
            out.append(((M-pred)/s)[mask]); det[('Xb4b<-Xb4d',B)]=(pred,M)
        return (np.concatenate(out),det) if detail else np.concatenate(out)
    return resid
if __name__=='__main__':
    res={}; txt=[]; preds={}
    for line in ('HCO+','HCN'):
        D=build(line)
        for p in (1,2,4):
            for Lf in (6.,8.7,12.,16.):
                rs=make_resid(*D,Lf,p,line)
                best=None
                for q0 in ((20.,60.),(30.,100.),(15.,40.)):
                    th0=[np.log(L.L2q(q0[0])),np.log(L.L2q(q0[1])),0.9,0.9]
                    r=least_squares(rs,th0,bounds=([np.log(L.L2q(3)),np.log(L.L2q(10)),0.3,0.3],[np.log(L.L2q(200)),np.log(L.L2q(400)),1.6,1.6]))
                    if best is None or r.cost<best.cost: best=r
                r=best; n=len(r.fun); chi2=2*r.cost; red=chi2/(n-4)
                J=r.jac; cov=np.linalg.pinv(J.T@J)*max(red,1); er=np.sqrt(np.diag(cov))
                Ld,Lb=L.q2L(np.exp(r.x[:2])); 
                key=f'{line}|p={p}|Lf={Lf:g}'
                res[key]=dict(L_d=float(Ld),L_b=float(Lb),g_d=float(r.x[2]),g_b=float(r.x[3]),eL_d=float(Ld*er[0]),eL_b=float(Lb*er[1]),eg_d=float(er[2]),eg_b=float(er[3]),chi2_red=float(red),npix=n)
                txt.append(f'{key:22s} L_d={Ld:6.1f}+-{Ld*er[0]:.1f} m  L_b={Lb:6.1f}+-{Lb*er[1]:.1f} m  g_d={r.x[2]:.3f}+-{er[2]:.3f} g_b={r.x[3]:.3f}+-{er[3]:.3f}  chi2_red={red:.2f} (N={n})'); print(txt[-1],flush=True)
                if p==2 and Lf==8.7:
                    rs2=make_resid(*D,Lf,p,line,detail=True); _,det=rs2(r.x); preds[line]={f'{k[0]}|{k[1]}':[v[0].tolist(),v[1].tolist()] for k,v in det.items()}
    json.dump(res,open('/workspace/out/q12/transfer_fit.json','w'),indent=1); open('/workspace/out/q12/transfer_fit.txt','w').write('\n'.join(txt))
    json.dump(preds,open('/workspace/out/q12/transfer_pred.json','w'))
