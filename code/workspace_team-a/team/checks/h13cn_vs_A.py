# reviewer check: what H13CN 4-3/HCN 4-3 does scenario A (HCN thick, sub-thermal) predict at S? vs observed limit alpha<~0.12 (3 sigma, team-b/ours)
import subprocess,re,os
def run(mol,T,n,N,dv=150.):
    wd='/tmp/revchk'; os.makedirs(wd,exist_ok=True)
    subprocess.run(['radex'],input=f'/opt/lamda/{mol}.dat\n{wd}/o.txt\n150 400\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n{dv}\n0\n',text=True,capture_output=True)
    r={}
    for l in open(f'{wd}/o.txt'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit(): r[int(p[1])]=(float(re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',p[7])),float(re.sub(r'(\d)([+-]\d{3})$',r'\1E\2',p[6])))
    return r
for (T,ln,lN) in [(40,4.25,16.75),(40,4.25,16.25),(40,4.5,16.5),(30,4.5,17.0),(40,3.0,17.0)]:
    h=run('hcn',T,10**ln,10**lN)
    for A in (30,50,90):
        x=run('h13cn@xpol',T,10**ln,10**lN/A)
        print(T,ln,lN,'A=',A,'tau43 HCN %.1f'%h[3][1],'tau43 H13CN %.2f'%x[3][1],'H13CN/HCN TR(4-3)=%.3f'%(x[3][0]/h[3][0]))
