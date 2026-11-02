import re
t=open('lit/dl/b21_ar5iv.html').read()
i=t.find('Mass-weighted percentiles'); seg=t[i:i+80000]
def dd(s):
    s=s.strip()
    n=len(s)
    return s[:n//2] if n%2==0 and s[:n//2]==s[n//2:] else s
for r in re.findall(r'<tr.*?</tr>',seg,re.S)[3:20]:
    c=[dd(re.sub(r'<[^>]+>','',x)) for x in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>',r,re.S)]
    print(c)
