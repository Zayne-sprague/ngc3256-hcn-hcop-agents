"""r11: N/S line-ratio table from team-b's single pipeline (scripts/s320_q52_chem_ratios.py -> out/q52/s320.json; s321_q52_ratios.py -> s321.json).
usage: python3 -I r11_table.py -> r11_table.md, r11_table.csv"""
import json
d=json.load(open('/workspace/out/q52/s321.json'))
rows=[('HNC/HCO+'),('CN/HCO+'),('HCN/HCO+'),('HNC/HCN'),('CN/HCN')]
cols=[('10',2.5,'1-0 2.5"'),('32',1.4,'3-2 1.4"'),('32',2.5,'3-2 2.5"')]
L=['| ratio (T_b) | '+' | '.join(f'{c[2]} N | S | N/S' for c in cols)+' |','|'+'---|'*(1+3*len(cols))]; C=['ratio,'+','.join(f'{c[2]} N,S,N/S' for c in cols)]
for r in rows:
    a=[];b=[]
    for J,t,_ in cols:
        x=d[f'{J}|{t}|{r}']; a+= [f"{x['N']:.3f}±{x['eN']:.3f}",f"{x['S']:.3f}±{x['eS']:.3f}",f"**{x['NS']:.2f}±{x['eNS']:.2f}**"]
        b+= [f"{x['N']:.4f}",f"{x['S']:.4f}",f"{x['NS']:.3f}"]
    L.append(f'| {r} | '+' | '.join(a)+' |'); C.append(r+','+','.join(b))
open('r11_table.md','w').write('\n'.join(L)+'\n'); open('r11_table.csv','w').write('\n'.join(C)+'\n'); print('\n'.join(L))
