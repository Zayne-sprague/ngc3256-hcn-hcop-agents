# Q35: N/S double ratios vs J-pair (HCN/CS, HCO+/CS, HCN/HCO+): Harada 1.7" (1-0/2-1, 3-2/5-4) and Xb47 4.3" deblend (4-3/7-6, MC errors s210b).  usage: python3 -I s211_q35_plot.py -> out/q35/s211_ns_vs_J.png
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
x=[0,1,2]; lab=['1-0/2-1\n(1.7")','3-2/5-4\n(1.7")','4-3/7-6\n(4.3" deblend)']
D={'HCN/CS':([1.63,1.20,1.10],[.17,.21,.47]),'HCO+/CS':([1.09,.63,.68],[.11,.11,.26]),'HCN/HCO+':([1.49,1.92,1.62],[.11,.16,.46])}
col={'HCN/CS':'#1f77b4','HCO+/CS':'#d62728','HCN/HCO+':'#444444'}
fig,ax=plt.subplots(figsize=(6,4))
for i,(k,(v,e)) in enumerate(D.items()): ax.errorbar(np.array(x)+0.06*(i-1),v,e,fmt='o-',color=col[k],label=k,capsize=3)
ax.axhline(1,color='k',lw=.5,ls=':'); ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel('N/S double ratio'); ax.legend(); ax.set_title('NGC 3256 N/S: normalised line ratios vs J')
fig.tight_layout(); fig.savefig('/workspace/out/q35/s211_ns_vs_J.png',dpi=130)
