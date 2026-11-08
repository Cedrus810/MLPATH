import json, collections, hashlib
from pathlib import Path
base=Path('/home/ruigengji/MLPATH/benchmarks/explore_d2/results')
out={}
for date,prefix in [('2026-10-26','b05_hf_'),('2026-10-28','b05_hf_v11_')]:
 for seed in ['sA','sB']:
  p=base/date/(prefix+seed); m=json.loads((p/'manifest.json').read_text()); report=json.loads((p/'report.json').read_text())
  first=last=None; hs={};hc={}; errors=[]; settled=collections.Counter(); skip=collections.Counter();cost=collections.Counter(); dispatch={}; duplicates=[]; rows=[]
  with (p/'journal.jsonl').open() as f:
   for line in f:
    ev=json.loads(line); first=ev['t'] if first is None else first;last=ev['t'];kind=ev.get('event');task=ev.get('task','');stage=task.split('|')[0]
    if kind=='dispatch': dispatch[task]=ev['t']
    if kind=='skip':skip[ev['reason']]+=1
    if kind!='settle':continue
    settled[stage]+=1; cost[stage]+=ev['spent']
    if stage not in ('HS','HC'):continue
    z={k:ev.get(k) for k in ['task','status','spent','cumulative_spent','t']};z['elapsed_s']=ev['t']-dispatch[task];pl=ev['payload'];z['cluster']=pl.get('start_cluster');z['class']=pl.get('class');z['payload']=pl if stage=='HC' or ev['status'].startswith('error:') else None
    (hs if stage=='HS' else hc)[task]=z
    if ev['status'].startswith('error:'):errors.append(z)
  ranks=collections.Counter('first' if z['cluster']['rank']==0 else 'second' for z in hs.values() if z['cluster'])
  stage_round={}
  for r in ['first','second']:
   group=[z for z in hs.values() if z['cluster'] and (z['cluster']['rank']==0)==(r=='first')]
   stage_round[r]={'items':len(group),'forwards':sum(z['spent'] for z in group),'statuses':dict(collections.Counter(z['status'] for z in group)),'zero_cost':sum(z['spent']==0 for z in group)}
  ledger=collections.Counter();open_seg={}
  with (p/'ledger.jsonl').open() as f:
   for line in f:
    ev=json.loads(line)
    if ev['event']=='close':ledger[ev['task'].split('|')[0]]+=ev['forwards']
  compact={'run':str(p),'host':m['host'],'wall_journal_s':last-first,'report':{k:v for k,v in report.items() if k not in ('c1','c2')},'settled':dict(settled),'spent_by_stage':dict(cost),'ledger_by_stage':dict(ledger),'errors':errors,'rounds':stage_round,'covered_clusters':len({z['cluster']['leader'] for z in hs.values() if z['cluster']}),'hs':hs,'hc':hc,'manifest':m}
  out[p.name]=compact
  print('\nRUN',p.name,'wall_h',round((last-first)/3600,3),'cost',dict(cost),'ledger',dict(ledger))
  print('rounds',stage_round,'covered',compact['covered_clusters'])
  for z in errors[:1]:print('ERROR_EXAMPLE',json.dumps(z))
  print('ERROR_COUNTS',dict(collections.Counter(z['payload']['error'] for z in errors)))
  print('CONNECTED',[(z['cluster'],z['spent']) for z in hc.values() if z['status']=='connected'])
for s in ['sA','sB']:
 a=out['b05_hf_'+s];b=out['b05_hf_v11_'+s];am=a['manifest'];bm=b['manifest']
 print('COMPARE',s,'settings_equal',am['settings']==bm['settings'],'hf_equal',am['hf']==bm['hf'],'config_equal',am['config']==bm['config'],'l0_equal',am['identity']['l0']==bm['identity']['l0'])
 x=a['hs'];y=b['hs'];common=x.keys()&y.keys();new=y.keys()-x.keys();lost=x.keys()-y.keys()
 print('HS_OVERLAP',len(common),'new',len(new),'old_not_run',len(lost),'common_changes',dict(collections.Counter((x[t]['status'],y[t]['status']) for t in common if x[t]['status']!=y[t]['status'])))
 print('NEW_STATUSES',dict(collections.Counter(y[t]['status'] for t in new)))
Path('/tmp/mlpath_hf_completed_audit.json').write_text(json.dumps(out))
