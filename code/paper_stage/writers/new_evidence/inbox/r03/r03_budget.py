"""r03: adopted N/S budget (quadrature). Inputs are from team outputs cited in answer.md; no new data."""
import math
q=lambda *a: math.sqrt(sum(x*x for x in a))
# component: (up, down)
B={"1.0\"":dict(c=2.1,stat=.06,base=.01,array=.04,resid=(0,.08),uvform=(.10,.10)),
   "0.3\"":dict(c=2.2,stat=.09,base=.07,array=.05,resid=(0,.08),uvform=(.10,0))}
for k,b in B.items():
    up=q(b['stat'],b['base'],b['array'],b['resid'][0],b['uvform'][0]); dn=q(b['stat'],b['base'],b['array'],b['resid'][1],b['uvform'][1])
    print(k,'central',b['c'],'+%.2f/-%.2f'%(up,dn),'sym',round(max(up,dn),2))
