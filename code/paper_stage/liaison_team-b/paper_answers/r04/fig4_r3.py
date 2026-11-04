"""R04 / Figure 4: observed N/S HCN/HCO+ contrast at 2.2" vs excitation-only / column-only / abundance-only / full hypotheses.
Inputs: out/q57/s356_s*_m*.txt (profile likelihood of 1-0 model contrast, BF), out/q57/s354_1.6-2.0.txt (best-fit contrasts, 3-2),
out/q29/s183.json, out/q34/s192.json (observed 2.2" data, +-150 km/s), out/q59/s350.txt (13C ratios), p0062 for A=N(HCN)/N(HCO+).
usage: python3 -I fig4_r3.py"""
import re, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
Q='/workspace/out/q57/'
def blocks(fn):
    B=[]; cur=None
    for l in open(fn):
        if l.startswith('## '): cur=dict(head=l[3:].strip(),prof={}); B.append(cur)
        elif cur is not None:
            m=re.search(r'BF\(Habun/Hexc\)=([\d.]+)',l)
            if m: cur['bf']=float(m.group(1))
            m=re.search(r'profile 1-0 model contrast (\w+) .*min ([\d.]+) at c=([\d.]+)\): (.*)',l)
            if m:
                cur['prof'][m.group(1)]=(float(m.group(2)),sorted([(float(a),float(b)) for a,b in re.findall(r'([\d.]+):([\d.]+)',m.group(4))]+[(float(m.group(3)),float(m.group(2)))]))
            m=re.search(r'delta\(-2lnL\) vs Hfull: (.*)',l)
            if m: cur['dchi']={k:float(v) for k,v in re.findall(r'(\w+) ([\d.]+)',m.group(1))}
    return B
EM={'halved':'s0.025_m0.05','default':'s0.05_m0.1','doubled':'s0.1_m0.2'}
BL={k:blocks(Q+f's356_{v}.txt') for k,v in EM.items()}
def pick(B,xr,phys):
    for b in B:
        if f'X {xr}' in b['head'] and (('L<=100pc (N core), x<=1e-8' in b['head'] and 'log n>=3.8' in b['head']) if phys=='phys' else b['head'].startswith('unconstrained')): return b
# observed
r183=json.load(open('/workspace/out/q29/s183.json')); r192=json.load(open('/workspace/out/q34/s192.json'))['res']
def obs(pos):
    a=r183['2015.1.00993 X961|2.2|'+pos]; c=r192['2p2|comb|'+pos]
    h1,eh1=a['HCN']['h150_ap']; o1,eo1=a['HCOp']['h150_ap']; h3,eh3=c['HCN32|ap|h150']; o3,eo3=c['HCOp32|ap|h150']
    f=lambda h,eh,o,eo:(h/o,h/o*np.hypot(eh/h,eo/o))
    return f(h1,eh1,o1,eo1),f(h3,eh3,o3,eo3)
(N1,N3),(S1,S3)=obs('N'),obs('S')
def con(n,s): c=n[0]/s[0]; return c,c*np.hypot(n[1]/n[0],s[1]/s[0])
c10=con(N1,S1); c32=con(N3,S3); print('obs contrast 1-0',c10,'3-2',c32)
# plot
plt.rcParams.update({'font.size':10,'axes.labelsize':11})
fig=plt.figure(figsize=(12,4.8)); ax=fig.add_axes([0.07,0.13,0.45,0.8]); ax2=fig.add_axes([0.60,0.13,0.37,0.8])
col={'Hexc':'#0072B2','Hcol':'#E69F00'}; nm={'Hexc':'excitation only (Hexc)','Hcol':'column only (Hcol)'}
for xr,ls in (('40-100','-'),('20-200','--')):
    b=pick(BL['default'],xr,'unc')
    for h in ('Hexc','Hcol'):
        mn,pts=b['prof'][h]; ax.plot([p[0] for p in pts],[p[1] for p in pts],ls,marker='o',ms=4,color=col[h],label=f'{nm[h]}, X={xr}' )
b=pick(BL['default'],'40-100','unc')
ax.plot([1.33],[b['dchi']['Habun']],marker='s',ms=9,color='#009E73',ls='',label='abundance only (Habun): best fit, contrast 1.33')
ax.plot([1.29],[0],marker='D',ms=8,clip_on=False,color='k',ls='',label='full (Hfull): best fit, contrast 1.29')
ax.axvspan(c10[0]-c10[1],c10[0]+c10[1],color='#CC79A7',alpha=0.35,lw=0)
ax.axvline(c10[0],color='#CC79A7'); ax.text(c10[0]+0.01,18.8,'observed 1-0\n%.2f±%.2f (stat)'%c10,color='#8a3a6e',fontsize=9,va='top')
ax.set_xlabel('model HCN/HCO$^+$ N/S contrast at 1-0 (equal abundance ratio at N and S)'); ax.set_ylabel(r'$\Delta\chi^2$ relative to Hfull (profile)')
ax.set_ylim(-1.2,19); ax.set_xlim(0.95,1.65); ax.legend(fontsize=7.5,loc='upper left',frameon=False,bbox_to_anchor=(0,0.93))
ax.set_title('2.2" beam, ±150 km/s. Default error model (5% scale + 10% model)',fontsize=9)
t=("Best-fit 3-2 contrast (profile not computed):\n obs %.2f±%.2f; Hexc 1.23, Hcol 1.35, Habun 1.57, Hfull 1.53\n"
   "Best-fit 1-0: Hexc 1.09, Hcol 1.13, Habun 1.33, Hfull 1.29")%c32
pass
# right: ln BF
cases=[('X 40-100','40-100','unc'),('X 40-100\nphys.','40-100','phys'),('X 20-200','20-200','unc'),('X 20-200\nphys.','20-200','phys')]
w=0.26; cc={'halved':'#56B4E9','default':'#0072B2','doubled':'#001f4d'}
for i,(lab,xr,ph) in enumerate(cases):
    for j,(k,c) in enumerate(cc.items()):
        v=np.log(pick(BL[k],xr,ph)['bf']); ax2.bar(i+(j-1)*w,v,w,color=c,label=f'errors {k}' if i==0 else None)
        ax2.text(i+(j-1)*w,v+0.05,'%.1f'%v,ha='center',fontsize=7)
ax2.axhline(0,color='k',lw=.8); ax2.axhline(np.log(3),color='0.5',ls=':',lw=1); ax2.text(-0.45,np.log(3)+0.05,'BF=3',fontsize=7,ha='left',color='0.4')
ax2.set_xticks(range(4)); ax2.set_xticklabels([c[0] for c in cases],fontsize=8); ax2.set_ylabel(r'ln BF (Habun / Hexc)'); ax2.legend(fontsize=8,frameon=False,loc='upper center')
ax2.set_ylim(0,8); ax2.set_title('Evidence for abundance-only over excitation-only',fontsize=9)
fig.savefig('/workspace/paper_answers/r04/fig4_contrast_hypotheses.png',dpi=180)
# table text
L=['ln BF(Habun/Hexc) (from s356 files):']
for lab,xr,ph in cases: L.append('  %-14s '%lab.replace('\n',' ')+'  '.join('%s %.2f (BF %.1f)'%(k,np.log(pick(BL[k],xr,ph)['bf']),pick(BL[k],xr,ph)['bf']) for k in EM))
open('/workspace/paper_answers/r04/lnBF_table.txt','w').write('\n'.join(L)+'\n'); print('\n'.join(L))
