import subprocess,os,sys
def run(f,Tk,n,N,dv=150):
    inp=f"/opt/lamda/{f}\no.txt\n80 370\n{Tk}\n1\nH2\n{n}\n2.73\n{N}\n{dv}\n0\n"
    subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True)
    r={}
    for l in open('o.txt'):
        p=l.split()
        if len(p)>=12 and p[1]=='--': r[round(float(p[4]),1)]=(float(p[6]),float(p[7]),float(p[8]))
    return r
for Tk,ln,lN,X in [(29,5.4,16.15,32),(300,4.6,15.5,32),(29,5.4,16.3,45)]:
    a=run('hcn@xpol.dat',Tk,10**ln,10**lN); b=run('h13cn@xpol.dat',Tk,10**ln,10**lN/X)
    print(Tk,ln,lN,X)
    for k in a:
        if abs(k-354.5)<1 or abs(k-88.6)<1: print(' HCN',k,a[k])
    for k in b:
        if abs(k-345.3)<1 or abs(k-86.3)<1: print(' H13CN',k,b[k])
