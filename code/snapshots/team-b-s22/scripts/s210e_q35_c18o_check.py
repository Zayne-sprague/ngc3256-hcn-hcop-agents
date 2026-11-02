# Q35: independent check of the N-vs-S excitation contrast with C18O / 13CO 2-1/1-0 (Harada 1.7"): RADEX (LAMDA c18o/co) at the (Tk,n) posteriors of s210c ('ladder+4.3': N 199 K / 10^4.5; S 100 K / 10^4.6; and the HCN+HCO+-only
# (hh) solutions N 22 K / 10^4.0, S 50 K / 10^3.9). Predicted K-km/s ratio I(2-1)/I(1-0) for thin C18O (N=1e15, dv 150) vs observed N 1.64, S 1.11.  usage: python3 -I s210e_q35_c18o_check.py
import subprocess,os,tempfile,numpy as np
def run(mol,Tk,ln,N,fs,dv=150.):
    d=tempfile.mkdtemp(prefix='rdx_',dir='/tmp'); os.chdir(d)
    inp=f"/opt/lamda/{mol}.dat\no.txt\n80 370\n{Tk}\n1\nH2\n{10**ln:.6g}\n2.73\n{N:.6g}\n{dv}\n0\n"; subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    out={}
    for l in open('o.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--':
            for k,f in fs.items():
                if abs(float(p[4])-f)<2e-3: out[k]=1.0645*dv*float(p[8])
    return out
T=[]
fs={'1-0':109.782,'2-1':219.560}; fs13={'1-0':110.201,'2-1':220.399}
for lab,Tk,ln in (('N ladder+4.3',199,4.5),('S ladder+4.3',100,4.6),('N hh',22,4.0),('S hh',50,3.9),('N at Tk=100, same n',100,4.5)):
    r=run('c18o',Tk,ln,1e15,fs); r2=run('c18o',Tk,ln,3e16,fs)
    T.append(f'C18O {lab}: Tk {Tk} n 10^{ln}: I21/I10 = {r["2-1"]/r["1-0"]:.2f} (N=1e15), {r2["2-1"]/r2["1-0"]:.2f} (N=3e16)')
print('\n'.join(T)); open('/workspace/out/q35/s210e.txt','w').write('\n'.join(T)+'\nobserved (Harada 1.7"): C18O 2-1/1-0 N 1.64 S 1.11; 13CO 2-1/1-0 N 1.67 S 1.25\n')
