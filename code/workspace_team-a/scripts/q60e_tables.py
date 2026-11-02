"""Q60 step 4: tables from out/q60/uvfit.csv, cog.csv"""
import pandas as pd, numpy as np
pd.set_option('display.width',250)
d=pd.read_csv('/workspace/out/q60/uvfit.csv',keep_default_na=False,na_values=['nan'])
def tab(us,mod,hw):
    s=d[(d.hw==hw)&(d.uvset==us)&(d['model']==mod)]
    p=s.pivot_table(index=['proj','eb','line'],columns='comp',values='I').reset_index()
    e=s.pivot_table(index=['proj','eb','line'],columns='comp',values='err').reset_index()
    p['eN']=e['N'];p['eS']=e['S']; p['N/S']=p['N']/p['S']; return p
if __name__=='__main__':
    for us,mod,hw in [('all',2,350),('all',3,350),('B',2,350),('C',2,350)]:
        print('\n==',us,'model',mod,'hw',hw); print(tab(us,mod,hw).round(2).to_string())
    print(d[(d.hw==350)&(d.uvset=='all')&(d['model']==2)&(d.comp=='N')][['proj','eb','line','vmin','vmax','nbase','r50']].to_string())
