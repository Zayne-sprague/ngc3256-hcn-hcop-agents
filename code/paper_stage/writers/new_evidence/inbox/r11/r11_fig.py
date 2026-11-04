"""r11 figure: N/S HCN/HCO+ 3-2 contrast vs estimator. Reads r11_apertures.json, r11_xb4b_imgtaper.json. Usage: python3 -I r11_fig.py"""
import json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
a=json.load(open('r11_apertures.json'))['tab']; x=json.load(open('r11_xb4b_imgtaper.json'))
labs=['peak','r0.25','r0.5','r0.75','rbeam']; xt=['peak\npixel','r=0.25"','r=0.5"','r=0.75"','r=beam\nFWHM\n(1.5" / 2.2")']
fig,ax=plt.subplots(1,2,figsize=(10,4.2),sharey=True)
for i,(bn,B) in enumerate((('1p5',1.5),('2p2',2.2))):
    for vn,c,m,lbl in (('comb','#0072B2','o','12m+7m uv-taper'),('c12','#D55E00','s','12m-only uv-taper')):
        y=[a[f'{bn}|{vn}|{l}']['NS'][0] for l in labs]; e=[a[f'{bn}|{vn}|{l}']['NS'][1] for l in labs]
        ax[i].errorbar(range(5),y,e,color=c,marker=m,capsize=3,label=lbl,ms=6)
    if B==1.5:
        y=[x[f'1.5|{l}']['NS'][0] for l in labs]; e=[x[f'1.5|{l}']['NS'][1] for l in labs]
        ax[i].errorbar(range(5),y,e,color='#009E73',marker='^',ls='--',capsize=3,label='Xb4b image convolved (0.24"→1.5")',ms=6)
    ax[i].set_xticks(range(5)); ax[i].set_xticklabels(xt,fontsize=9); ax[i].set_title(f'{B}" circular beam',fontsize=11); ax[i].grid(alpha=.3)
    ax[i].axvspan(3.5,4.5,color='0.92',zorder=0); ax[i].legend(fontsize=8,loc='lower left'); ax[i].set_xlabel('estimator',fontsize=10)
ax[0].set_ylabel('N/S contrast of HCN/HCO$^+$ 3–2 (ratio of T$_b$ ratios)',fontsize=10); ax[0].set_ylim(1.3,2.2)
plt.tight_layout(); plt.savefig('estimator_spread.png',dpi=170)
