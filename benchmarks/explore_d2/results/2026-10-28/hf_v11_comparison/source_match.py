import sys,json
from pathlib import Path
import numpy as np
sys.path[:0]=['/home/ruigengji/MLPATH/explore/src','/home/ruigengji/MLPATH/src']
from ase import Atoms
from prrs.config import SearchConfig
from prrs.chemistry import automorphisms,symmetric_rmsd
x=json.loads(Path('/tmp/mlpath_hf_completed_audit.json').read_text());out={}
for seed in ['sA','sB']:
 old=x['b05_hf_'+seed];new=x['b05_hf_v11_'+seed];p=Path(new['run']);pose=json.loads((p/'poses.jsonl').read_text().splitlines()[0]);template=Atoms(numbers=pose['numbers'],positions=pose['positions']);c=SearchConfig.from_dict(new['manifest']['config']);group,overflow=automorphisms(template,c.bond_scale,c.active_atoms,c.automorphism_limit)
 oe=[v for v in old['hc'].values() if v['status']=='connected'];ne=[v for v in new['hc'].values() if v['status']=='connected'];mat=[]
 for e in oe:
  a=template.copy();a.positions=e['payload']['saddle_positions'];row=[]
  for f in ne:
   b=template.copy();b.positions=f['payload']['saddle_positions'];row.append(symmetric_rmsd(a,b,group,c.active_atoms))
  mat.append(row)
 matches=[]
 for i,e in enumerate(oe):
  j=int(np.argmin(mat[i]));f=ne[j]
  matches.append({'old_task':e['task'],'new_task':f['task'],'rmsd_A':mat[i][j],'energy_delta_eV':abs(e['payload']['saddle_energy_eV']-f['payload']['saddle_energy_eV']),'chemical_pair_equal':e['payload']['chemical_pair']==f['payload']['chemical_pair']})
 print(seed,'matrix',np.round(mat,4).tolist(),'MATCH',[(round(z['rmsd_A'],5),z['energy_delta_eV'],z['chemical_pair_equal']) for z in matches])
 first=[]
 for e in ne:
  cl=e['cluster'];leader=new['hs']['HS|'+cl['leader']];hc=new['hc'].get('HC|'+cl['leader']);z={'event_task':e['task'],'rank':cl['rank'],'leader_hs_status':leader['status'],'leader_hc_status':hc['status'] if hc else None};first.append(z)
 print('event leader status',first)
 genuine_new=[t for t in new['hs'] if t not in old['hs'] and new['hs'][t]['spent']>0]
 print('NEW_NONZERO_COST_HS',len(genuine_new),'NEW_TOTAL_HS',len(new['hs'])-len(new['errors']))
 out[seed]={'group_size':len(group),'overflow':overflow,'old_to_new_rmsd_A':mat,'matches':matches,'new_events_cluster_leaders':first,'newly_run_nonzero_cost_hs':len(genuine_new)}
Path('/tmp/mlpath_hf_event_match.json').write_text(json.dumps(out,indent=2))
