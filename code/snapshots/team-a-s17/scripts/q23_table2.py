import sys,pdfplumber
p=pdfplumber.open(sys.argv[1])
for i,pg in enumerate(p.pages):
    t=pg.extract_text() or ''
    if 'Mass-weighted percentiles' in t:
        print('page',i+1)
        for tb in pg.extract_tables():
            for r in tb: print(r)
        print(t[:4000])
