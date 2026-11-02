"""Q29 s188: opacity and abundance quantities from the weak lines (s185) at the 2.2" circular beam.
(1) Direct tau(HCO+ 1-0) and tau(HCN 1-0) from the isotopologue integral ratios R=I(12)/I(13) for 12C/13C = 40, 70, 100:  R = (1-e^-tau)/(1-e^-tau/X)  (same Tex, filling, width).
(2) Column densities from optically thin RADEX (LVG, dv=150 km/s) conversions K(Tk,n) = I_model/N for H13CO+ (HCO+ file, scaled A: <10%), H13CN (hcn.dat), SiO 2-1 (sio-h2.dat), HC3N 10-9 (hc3n-h2.dat),
over Tk 25-100 K, n(H2) 1e4-1e6; N(HCO+)=X*N(H13CO+), N(HCN)=X*N(H13CN) (or limit). Abundance ratios relative to HCO+ (common filling, same gas):  SiO/HCO+, HC3N/HCO+, HCN/HCO+.
usage: python3 -I s188_q29_abundances.py -> out/q29/s188.txt/json"""
import json, subprocess, os, tempfile, numpy as np
from scipy.optimize import brentq
O='/workspace/out/q29/'; J=json.load(open(O+'s185.json')); T=J['table']; B=json.load(open(O+'s185_box.json')); log=[]
def P(s=''): print(s,flush=True); log.append(s)
def Rm(tau,X): return (1-np.exp(-tau))/(1-np.exp(-tau/X))
def tau_of(R,X):
    if R>=Rm(1e-4,X): return 0.0
    if R<=Rm(60,X): return 60.0
    return brentq(lambda t:Rm(t,X)-R,1e-4,60)
def tauband(R,eR,X):
    return tau_of(R,X),tau_of(R*(1+eR),X),tau_of(R/(1+eR),X)      # central, high-R(low tau), low-R(high tau) -> returns (c, lo_tau, hi_tau)
res={}
P('== (1) direct opacity from isotopologue ratios, 2.2" beam ==')
rows=[]
for reg in ('N','S','Bridge'):
    a=T[f'{reg}|H13COp|comb']; I12=J['HCOp_int'][reg][0]; R=I12/a['amp']; eR=a['err']/a['amp']; eR=np.hypot(eR,0.08)
    rows.append((reg,'HCO+/H13CO+ (template integral)',R,eR))
    b=B[f'pm150|{reg}|H13COp']; rows.append((reg,'HCO+/H13CO+ (+-150 aperture)',b['HCOp_int']/b['amp'],np.hypot(b['err']/b['amp'],0.08)))
    b=B[f'pm150pk|{reg}|H13COp']; rows.append((reg,'HCO+/H13CO+ (+-150 peak pixel of 2.2" map)',b['HCOp_int']/b['amp'],np.hypot(b['err']/b['amp'],0.08)))
    b=B[f'pm150|{reg}|H13CN']; rows.append((reg,'HCN/H13CN (+-150 aperture)',b['HCN_int']/b['amp'] if b['amp']>0 else np.inf,np.hypot(b['err']/max(b['amp'],1e-9),0.08)))
    b=B[f'pm150pk|{reg}|H13CN']; rows.append((reg,'HCN/H13CN (+-150 peak pixel)',b['HCN_int']/b['amp'] if b['amp']>0 else np.inf,np.hypot(b['err']/max(b['amp'],1e-9),0.08)))
for reg,lab,R,eR in rows:
    s=f'{reg:6s} {lab:44s} R={R:6.1f} (x/ {1+eR:.2f}) '
    out={}
    for X in (40,70,100):
        c,lo,hi=tauband(R,eR,X); out[X]=(c,lo,hi)
        s+=f'| X={X}: tau={c:.2f} [{lo:.2f},{hi:.2f}] '
    P(s); res[f'tau|{reg}|{lab}']=dict(R=R,eR=eR,tau=out)
# (2) thin RADEX conversions
def radex(mol,Tk,n,N,fr,dv=150.,lo=50,hi=400):
    d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); cwd=os.getcwd(); os.chdir(d)
    try:
        inp=f"/opt/lamda/{mol}.dat\no.txt\n{lo} {hi}\n{Tk}\n1\nH2\n{n}\n2.73\n{N}\n{dv}\n0\n"
        subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
        for l in open('o.txt'):
            p=l.split()
            if len(p)>=12 and p[1]=='--' and abs(float(p[4])-fr)<0.02: return dict(Tex=float(p[6]),tau=float(p[7]),TR=float(p[8]),flux=float(p[11]))
    finally: os.chdir(cwd)
LINES=dict(H13COp=('hco+',86.7543),H13CN=('h13cn@xpol',86.3402),SiO=('sio-h2',86.8470),HC3N=('hc3n-h2',90.9790),HCOp=('hco+',89.1885),HCN=('hcn',88.6318))
def K(line,Tk,n):
    mol,fr=LINES[line]
    if line=='H13COp': fr=89.1885      # no H13CO+ LAMDA file: use the HCO+ file in the thin limit (same collisional rates), flux scaled x0.90 for (A(H13CO+)/A(HCO+))*(nu ratio)
    r=radex(mol,Tk,n,1e12,fr)
    return r['flux']/1e12*(0.90 if line=='H13COp' else 1.0), r['tau']
P('\n== (2) thin-line RADEX flux per unit column K(Tk,n) [K km/s per cm^-2, x1e-12 -> K km/s for N=1e12], tau at 1e12 ==')
Kt={}
for Tk in (25,50,100):
    for n in (1e4,1e5,1e6):
        row=f'Tk={Tk:3d} n=1e{int(np.log10(n))}: '
        for l in ('H13COp','H13CN','SiO','HC3N'):
            k,t=K(l,Tk,n); Kt[(l,Tk,n)]=k; row+=f'{l} {k*1e12:.3g} (tau {t:.1e})  '
        P(row)
P('\n== column densities and abundance ratios relative to HCO+ (N(HCO+)=X N(H13CO+), X=40-100), spans over Tk 25-100 K, n 1e4-1e6 ==')
def amp(reg,l,box=None): a=T[f'{reg}|{l}|'+('comb' if l=='H13COp' else {'SiO':'X95d_spw1','HC3N':'X961_spw1','H13CN':'X95d_spw0'}[l])]; return a['amp'],a['err'],a['snr'],a['cls']
for reg in ('N','S','Bridge'):
    h,eh,sh,ch=amp(reg,'H13COp'); 
    for l in ('SiO','HC3N','H13CN'):
        a,e,sn,cl=amp(reg,l); vals=[]
        for (ll,Tk,n),k in Kt.items():
            if ll!='H13COp': continue
            kl=Kt[(l,Tk,n)]
            # N(l)/N(H13CO+) = (a/kl)/(h/kH) ; then /X13 for relative abundance to HCO+ (SiO/HC3N) or direct ratio  N(H13CN)/N(H13CO+) = N(HCN)/N(HCO+) (same X)
            vals.append(((a/kl)/(h/k)))
        vals=np.array(vals)
        if l=='H13CN':
            lab='N(HCN)/N(HCO+) = N(H13CN)/N(H13CO+)  [X13 cancels if the same]'; v=vals; sc=1.0
        else: lab=f'N({l})/N(HCO+) x1e-3 [X13=40..100 -> /{40}..{100}]'; v=vals; sc=1.0
        if l=='H13CN':
            P(f'{reg:6s} {lab}: {np.min(v):.2f}-{np.max(v):.2f} median {np.median(v):.2f} (line {sn:.1f} sig {cl}; 3sig lim {3*e:.2f} K km/s -> ratio lim {np.median(3*e/a*v) if a>0 else np.median(v*0+3*e/np.maximum(a,3*e)*0):.2f})')
            res[f'abund|{reg}|HCN/HCO+']=dict(min=float(v.min()),max=float(v.max()),med=float(np.median(v)),snr=sn,cls=cl)
        else:
            lo=np.min(v)/100*1e3; hi=np.max(v)/40*1e3; mid=np.median(v)/70*1e3
            P(f'{reg:6s} X({l})/X(HCO+) = {lo:.2f}-{hi:.2f} x1e-3 (median-model & X=70: {mid:.2f}) ; line {a:.2f}+-{e:.2f} ({sn:.1f} sig {cl}); thin-model spread of N(l)/N(H13CO+): {v.min():.2f}-{v.max():.2f}')
            res[f'abund|{reg}|{l}/HCO+']=dict(lo=float(lo),hi=float(hi),mid=float(mid),snr=sn,cls=cl,spread=[float(v.min()),float(v.max())])
# column density of HCO+ itself (thin conversion of H13CO+)
for reg in ('N','S','Bridge'):
    h,eh,sh,ch=amp(reg,'H13COp'); Ns=[]
    for Tk in (25,50,100):
        for n in (1e4,1e5,1e6): Ns.append(h/Kt[('H13COp',Tk,n)])
    Ns=np.array(Ns); P(f'{reg:6s} N(H13CO+) = {Ns.min():.2e}-{Ns.max():.2e} cm-2 (beam-avg., 2.2" aperture, +-RJ); N(HCO+) = {Ns.min()*40:.1e}-{Ns.max()*100:.1e}')
    res[f'N|{reg}|HCO+']=dict(H13=[float(Ns.min()),float(Ns.max())],HCO=[float(Ns.min()*40),float(Ns.max()*100)],H13amp=h,snr=sh)
open(O+'s188.txt','w').write('\n'.join(log)); json.dump(res,open(O+'s188.json','w'),indent=1,default=float)
