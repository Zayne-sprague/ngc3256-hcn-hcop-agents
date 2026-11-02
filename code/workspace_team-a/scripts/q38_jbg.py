# Q38 step 3: envelope J_bg upper bound from global FIR-SED dust disk (Td, beta, F) vs the q36 'envhi' scenario (70 mK @258 GHz, beta 1.8)
import numpy as np,sys
sys.path.insert(0,'/snapshots/team-a-s22/scripts')
import q36_bg as B
H=6.62607e-27;KB=1.380649e-16;C=2.998e10;arc=np.pi/180/3600
Bn=lambda nu,T: 2*H*nu**3/C**2/np.expm1(H*nu/(KB*T)); nu160=C/160e-4
nus=np.array([88.63,177.26,265.89,354.51])
def tb(nu,F,Td,beta,fdiff=1.0):
    S=fdiff*93.48*(nu/nu160)**beta*Bn(nu,Td)/Bn(nu160,Td); om=np.pi*F**2/(4*np.log(2))*arc**2
    return S*1e-23/om*C**2/(2*KB*nu**2)
L=['Envelope J_bg = J_CMB + eta*Tb_RJ,dust(F,Td,beta) [eta=1 gas behind the whole disk; 0.5 embedded at mid-plane]. Columns: HCN/HCO+ 1-0,2-1,3-2,4-3. Tbg_eq in K.']
L.append('q36 envhi: Tbg_eq='+' '.join(f'{B.Tbg_eq(n,float(B.Jbg("E","envhi",n))):.2f}' for n in nus)+'   (Tb258=70 mK)')
for Td,beta in ((30,1.8),(30,1.5),(34,2.0),(25,1.5)):
  for F in (10,15,20,30):
    for eta in (1.0,0.5):
      J=np.array([float(B.Jbg('E','envmid',n))-0.02*(n/258.)**1.8+eta*tb(n*1e9,F,Td,beta) for n in nus])  # CMB + dust only
      L.append(f'Td={Td} beta={beta} F={F}" eta={eta}: Tb258={eta*tb(257.7e9,F,Td,beta)*1e3:.0f} mK  Tbg_eq='+' '.join(f'{B.Tbg_eq(n,j):.2f}' for n,j in zip(nus,J)))
print('\n'.join(L)); open('/workspace/out/q38/q38_jbg.log','w').write('\n'.join(L))
