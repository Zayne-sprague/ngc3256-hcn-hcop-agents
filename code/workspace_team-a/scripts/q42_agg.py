# Q42: across-seed aggregation of out/q42/q42_summary_runs2.csv (run q42_post2.py first) -> out/q42/q42_agg.csv
import pandas as pd,numpy as np
d=pd.read_csv('/workspace/out/q42/q42_summary_runs2.csv'); d=d[d.env=='sh']
c=['logZ','TcN_16','TcN_50','TcN_84','TcS_50','Te_50','P50','P_TcN_gt_TcS','dlogT_50','P_TcN_le30','P_TcN_lt20']
g=d.groupby(['mode','scen'])[c]; m=g.mean(); s=g.std(); n=g.size()
o=m.round(2).join(s.round(2),rsuffix='_sd').join(n.rename('nseed')); o.to_csv('/workspace/out/q42/q42_agg.csv')
pd.set_option('display.width',300); pd.set_option('display.max_columns',50); print(o.T.to_string())
