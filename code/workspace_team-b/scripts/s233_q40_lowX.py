# Q40: extend s230 to lower 12C/13C (X=25, 30) for variants H and A at N, S, C (Harada R(1-0)=22 at N suggests low X or thick lines). usage: python3 -I s233_q40_lowX.py
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s230_q40_isotopologue_fit as s,time
for v in ('H','A'):
  for X in (25,30):
    for pos in 'NSC':
        o=s.run_one(pos,np.log10(X),v); tag=f'{pos}_{v}_X{X}'; np.savez(s.O+f'fit_{tag}.npz',**o)
        cm=o['chimin']; Tk=s.m.TK; print(tag,f"chi2min {o['chibest']:.2f} cold-hot dchi2 {cm[Tk<=30].min()-cm[Tk>=100].min():.2f}",flush=True)
