"""Q60 step 6: feed uv-plane 4-3 HCN/HCO+ (Tb ratio, out/q60/summary.csv, uvset all, model 2, hw 350) into q55e_joint (same likelihood/grid).
Replacement R43 sets: (a) 00493+00223, (b) Xb47, (c) pooled per-EB with inter-EB+inter-project scatter as error; for the 2.5" (A) and core (B_core3) base sets.
Run: python3 -I q60h_joint.py > ../out/q60/joint43.txt"""
import sys,json; sys.path.insert(0,'/workspace/scripts')
import q55e_joint as Q
A=Q.SETS['A_full']; B=Q.SETS['B_core3']
vals={'a_00493+00223':(0.71,0.10,0.52,0.05),'b_Xb47':(0.64,0.07,0.39,0.08),'c_pooled':(0.675,0.09,0.455,0.09)}
S={}
for base,bn in ((A,'A'),(B,'B')):
    for k,(n,en,s,es) in vals.items():
        d=dict(base); d['R43']=('N',n,en,'S',s,es); d['sig']={**base['sig'],'R43':0.10}; S[f'{bn}_43uv_{k}']=d
S['A_old_43']=A; S['B_old_43']=B; S['A_43Xb47_img']=Q.SETS['A_43Xb47']
ALL={}
for nm,s in S.items():
    for sh in (True,False):
        rows=Q.summarise(nm,s,sh); ALL[f'{nm}|{sh}']=rows
        print(f'\n== {nm} X13 {"shared" if sh else "free"}')
        for d in rows: print(f"  {d['hyp']:38s} lnZ {d['lnZ']:8.2f}  lnBF vs best {d['lnBF_vs_best']:7.2f}")
json.dump(ALL,open('/workspace/out/q60/joint43.json','w'),indent=1,default=float)
