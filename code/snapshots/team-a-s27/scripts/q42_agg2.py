# Q42: pooled across all runs (runs2+runs3, seeds), per (mode,scen): mean and sd of per-run statistics -> out/q42/q42_agg_all.csv
import pandas as pd
d=pd.concat([pd.read_csv(f'/workspace/out/q42/q42_summary_{k}.csv') for k in ('runs2','runs3')]); d=d[d.env=='sh']
c=['logZ','TcN_16','TcN_50','TcN_84','TcS_50','Te_50','P50','P_TcN_gt_TcS','dlogT_50','P_TcN_le30','P_TcN_lt20']
g=d.groupby(['mode','scen'])[c]; o=g.mean().round(2).join(g.std().round(2),rsuffix='_sd').join(g.size().rename('n')); o.to_csv('/workspace/out/q42/q42_agg_all.csv')
pd.set_option('display.width',300); print(o.T.to_string())
