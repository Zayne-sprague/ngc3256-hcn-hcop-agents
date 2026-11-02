# compare Table 2 of Brunetti+2021 (arXiv:2011.01250): numbers parsed from the PDF text layer vs from the ar5iv HTML rendering
import re,pdfplumber,json
pdf=pdfplumber.open('/workspace/lit/dl/b21.pdf').pages[6].extract_text().splitlines()
P=[]
for l in pdf:
    t=l.split()
    if re.match(r'^(Non-nuclear|Nuclei|Combined|Northern|Nucleus|Southern)?\s*',l) and len(t)>=16:
        nums=[x for x in t if re.fullmatch(r'\d+(\.\d+)?',x)]
        if len(nums)==16 and nums[0] in('55','80','120'): P.append(nums)
        elif len(nums)==17 and nums[0] in ('55','80','120'): P.append(nums)
print(len(P),'PDF rows')
H=[]
for l in open('/workspace/out/q23_table2_html.txt'):
    c=eval(l); n=[x for x in c if re.fullmatch(r'\d+(\.\d+)?',x)]
    if n[0]=='5': n[0]='55'   # repeated-digit artefact of the HTML de-duplication ('5555'->'5')
    H.append(n)
print(len(H),'HTML rows')
bad=0
for i,(a,b) in enumerate(zip(P,H)):
    if a[:16]!=b[:16]: bad+=1; print('MISMATCH row',i,a,b)
print('rows compared',min(len(P),len(H)),'mismatches',bad)
