import subprocess,sys,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S
F=S.Fine()
def radex(mol,fr,T,n,N,dv=150.):
    f={'hcn':'hcn.dat','hcop':'hco+.dat'}[mol]
    inp=f"{f}\nout.txt\n{fr[0]} {fr[1]}\n{T}\n1\nH2\n{n}\n2.73\n{N}\n{dv}\n0\n"
    inp=f"/opt/lamda/{f}\nout.txt\n{fr[0]} {fr[1]}\n{T}\n1\nH2\n{n}\n2.73\n{N}\n{dv}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    L=[l.split() for l in open('out.txt') if l.strip() and not l.startswith('*') and l.split()[0].isdigit()]
    return L
# grid at off-node points
for T,n,N in [(150,10**4.62,10**13.5),(150,10**4.62,10**14.3),(60,10**4.3,10**14),(35,10**4.8,10**13.8)]:
    for mol,fr in (('hcn',(80,400)),):
        L=radex(mol,fr,T,n,N)
        print(T,n,N,[(r[0],r[1],r[7],r[8],r[11]) for r in L[:5]][:4])
