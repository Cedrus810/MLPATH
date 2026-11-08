import json,sys
from pathlib import Path
from types import SimpleNamespace
sys.path[:0]=['/home/ruigengji/MLPATH/explore/src','/home/ruigengji/MLPATH/src']
from explore.d2 import _classify,_points,_hf_bracket,_intent
x=json.loads(Path('/tmp/mlpath_hf_completed_audit.json').read_text());out={}
for seed in ['sA','sB']:
 r=x['b05_hf_v11_'+seed];p=Path(r['run']);ids={e['task'][3:] for e in r['errors']};need={prefix+'|'+iid for prefix in ['L1','L2'] for iid in ids};needles=['"task": "'+t+'"' for t in need];settled={}
 with (p/'journal.jsonl').open() as f:
  for line in f:
   head=line[:2200]
   if '"event": "settle"' not in head or not any(n in head for n in needles):continue
   e=json.loads(line)
   if e['event']=='settle' and e['task'] in need:settled[e['task']]=e
 rows=[json.loads(line) for line in (p/'population.jsonl').read_text().splitlines()]
 run=SimpleNamespace(settled=settled,settings=SimpleNamespace(noise_eV=r['manifest']['settings']['noise_eV']))
 result=[]
 for row in rows:
  if row['intent_id'] not in ids:continue
  cls=_classify(run,row);points=_points(settled['L1|'+row['intent_id']]['payload']);a,b=_hf_bracket(cls,points,_intent(row));result.append({'intent_id':row['intent_id'],'class':cls['class'],'q_targets':[a.q_target,b.q_target],'width':abs(a.q_target-b.q_target)})
 out[seed]=result;print(seed,[(z['q_targets'],z['width']) for z in result],flush=True)
Path('/tmp/mlpath_hf_bracket_audit.json').write_text(json.dumps(out,indent=2))
