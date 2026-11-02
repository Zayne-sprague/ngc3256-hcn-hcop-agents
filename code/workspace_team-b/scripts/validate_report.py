import json,sys
from pydantic import BaseModel,ValidationError,TypeAdapter
import jsonschema
s=json.load(open('/task/output_schema.json'))
d=json.load(open(sys.argv[1] if len(sys.argv)>1 else '/output/report.json'))
v=jsonschema.Draft202012Validator(s)
errs=list(v.iter_errors(d))
print(len(errs))
for e in errs: print(list(e.absolute_path), e.message[:150])
