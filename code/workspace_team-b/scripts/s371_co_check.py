import numpy as np
L=open('/opt/lamda/co.dat').read().split('\n')
nl=41; E=[float(L[i].split()[1]) for i in range(7,7+nl)]; g=[float(L[i].split()[2]) for i in range(7,7+nl)]
A={}
for l in L[51:91]:
    p=l.split(); A[(int(p[1])-1,int(p[2])-1)]=float(p[3])
idx=[i for i,l in enumerate(L) if l.startswith('!TRANS+')]
temps=[float(x) for x in L[idx[0]-1].split()]
print(len(temps),temps[:3],temps[-3:],idx)
def block(i0,n):
    d={}
    for l in L[i0+1:i0+1+n]:
        p=l.split(); d[(int(p[1])-1,int(p[2])-1)]=np.array(p[3:],float)
    return d
bl=[block(idx[0],820),block(idx[1],820)]
def rates(T,x):
    it=np.argmin(abs(np.array(temps)-T)); C=np.zeros((nl,nl))
    for (u,l),v in bl[0].items():
        r=x[0]*v[it]+x[1]*bl[1][(u,l)][it]
        C[u,l]=r; C[l,u]=r*g[u]/g[l]*np.exp(-(E[u]-E[l])*1.4388/T)
    return C
for n,T in((1e3,100),(1e4,100),(1e4,20),(1e6,100)):
    C=rates(T,(0.25*n,0.75*n)); M=C.T.copy()
    for (u,l),a in A.items(): M[l,u]+=a
    Mm=M-np.diag(M.sum(0)); Mm[-1,:]=1; b=np.zeros(nl); b[-1]=1
    x=np.linalg.solve(Mm,b)
    tex=lambda u,l: -(E[u]-E[l])*1.4388/np.log(x[u]*g[l]/(x[l]*g[u]))
    print(n,T,[round(tex(i+1,i),1) for i in range(3)])
