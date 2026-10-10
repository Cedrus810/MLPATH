import collections
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path('/home/ruigengji/MLPATH')
sys.path[:0] = [str(ROOT / 'src'), str(ROOT)]
from prrs.config import SearchConfig
from prrs.perturbations import propose
from prrs.search import _refuse_unaffordable_enumeration
from prrs.state import aligned_rmsd
from prrs.io import structure_hash
import runpy
implementation_sha256 = runpy.run_path(str(ROOT / "scripts/implementation_hash.py"))["implementation_sha256"]

run = ROOT / 'runs/b09_probe1'
manifest = json.loads((run / 'config.json').read_text())
network = json.loads((run / 'network.json').read_text())
records = [json.loads(line) for line in (run / 'attempts.jsonl').read_text().splitlines()]
cfg = SearchConfig.from_dict(manifest['config'])
proposal = network['amplitude_scans'][0]['proposal']['budget']
state = network['chemical_nodes'][0]['microstates'][0]
atoms = read(run / state['structure'])
heavy = np.flatnonzero(atoms.numbers > 1).tolist()
source_e = state['energy_eV']
selected = set(proposal['selected_ids'])
executed = set(network['direction_ledger']['executed_ids'])
uncovered = set(network['direction_ledger']['uncovered_ids'])
assert len(records) == network['termination_detail']['trials_spent'] == 32
assert selected == executed | set(network['direction_ledger']['not_executed_ids'])
assert not selected & uncovered
assert implementation_sha256() == manifest['implementation_sha256']
assert structure_hash(read(run / 'input.extxyz')) == manifest['input_structure_sha256']
assert (run / 'input.extxyz').read_bytes() == (ROOT / 'benchmarks/carbene_ambimodal/probe_b09/input.extxyz').read_bytes()

exhaustive = replace(cfg, direction_mode='exhaustive')
full_report = {}
full = propose(atoms, exhaustive, cfg.seed, report=full_report)
assert set(full_report['candidate_ids']) == set(proposal['candidate_ids'])
lengths = [len(batch) for batch in full]
coarse = sum(lengths)
# Every adjacent coarse pair could conservatively get a refinement loop. In fact,
# returned -> departed brackets cannot be adjacent; deliberately use a loose bound.
refine_bound = cfg.refinement_steps * sum(max(0, count - 1) for count in lengths)
family_lengths = {'stretch': len(cfg.geometry_amplitudes_A), 'compress': len(cfg.geometry_amplitudes_A),
                  'kick': len(cfg.kick_energies_eV), 'bend': len(cfg.bend_amplitudes_rad)}
control_lengths = [family_lengths[cid.split(':')[0]] for cid in proposal['selected_ids']]
control_coarse = sum(control_lengths)
control_bound = control_coarse + cfg.refinement_steps * sum(count - 1 for count in control_lengths)
configs = {}
for name, updates in [
    ('config_control24.json', dict(max_trials=control_bound, chemical_trials_per_node=control_bound,
                                   conformer_trials_per_microstate=control_bound)),
    ('config_exhaustive.json', dict(direction_mode='exhaustive', max_trials=coarse + refine_bound,
                                    chemical_trials_per_node=coarse + refine_bound,
                                    conformer_trials_per_microstate=coarse + refine_bound)),
]:
    new = replace(cfg, **updates)
    path = ROOT / 'benchmarks/carbene_ambimodal/probe_b09' / name
    assert not path.exists(), path
    path.write_text(json.dumps(new.to_dict(), ensure_ascii=False, indent=2) + '\n')
    configs[name] = {key: {'before': cfg.to_dict()[key], 'after': value} for key, value in updates.items()}
    if new.direction_mode == 'exhaustive':
        assert _refuse_unaffordable_enumeration(read(run / 'input.extxyz'), new) == coarse

saddles = []
for ts in network['ts_candidates']:
    geom = read(run / ts['structure'])
    saddles.append(dict(id=ts['id'], status=ts['status'], imaginary_wavenumbers_icm=ts['imaginary_wavenumbers_icm'],
                        rise_meV=(ts['energy_eV']-source_e)*1000,
                        heavy_rmsd_A=aligned_rmsd(atoms,geom,active=heavy),
                        r_Ccarb_Ca_A=geom.get_distance(2,13), r_Ccarb_Cb_A=geom.get_distance(2,15),
                        r_Ccarb_Oald_A=geom.get_distance(2,1), connects=ts['connects']))
targets = {}
for pair in ['2-13', '2-15']:
    for family, sign in [('compress','-'), ('kick','-'), ('kick','+')]:
        cid=f'{family}:{pair}:{sign}'
        targets[cid] = 'executed' if cid in executed else 'selected_not_executed' if cid in selected else 'candidate_not_selected' if cid in uncovered else 'basis_uncovered'
        assert targets[cid] == 'candidate_not_selected'
report = dict(date='2026-10-09', run=str(run), runtime_host_recorded=manifest['host'],
              recorded_implementation_sha256=manifest['implementation_sha256'], current_implementation_sha256=implementation_sha256(),
              source_files_sha256={name: hashlib.sha256((run/name).read_bytes()).hexdigest() for name in ['config.json','network.json','attempts.jsonl','input.extxyz']},
              status=network['status'], termination=network['termination'], termination_detail=network['termination_detail'],
              reaction_count=len(network['reactions']), chemical_node_count=len(network['chemical_nodes']),
              outcomes=dict(collections.Counter(r['status'] for r in records)), stages=dict(collections.Counter(r['stage'] for r in records)),
              ranking=proposal['ranked_by'], ranking_scores='Not saved in the historical run; no Hessian recomputation performed.',
              proposal_families=proposal['families'], selected_ids=proposal['selected_ids'],
              executed_ids=sorted(executed), target_direction_audit=targets,
              source_distances_A={'Ccarb_Ca':atoms.get_distance(2,13),'Ccarb_Cb':atoms.get_distance(2,15),'Ccarb_Oald':atoms.get_distance(2,1)},
              largest_injected_kinetic_eV=max(r.get('injected_kinetic_eV') or 0 for r in records),
              rejections=[{'attempt':r['attempt_id'],'failure':r['failure']} for r in records if r['status']=='rejected'],
              ts_candidates=saddles, enumeration={'directions':len(full),'coarse_trials':coarse,
                 'refinement_loose_upper_bound':refine_bound,'source_trial_budget_upper_bound':coarse+refine_bound,
                 'control24_coarse_trials':control_coarse,'control24_trial_budget_upper_bound':control_bound,
                 'scope':'One source proposal with the current primitive basis/grid only; gates, deeper states and scientific acceptance remain separate.'},
              prepared_config_differences=configs)
out = ROOT / 'benchmarks/carbene_ambimodal/results/2026-10-09'
out.mkdir(exist_ok=True)
p=out/'b09_probe1_audit.json'
assert not p.exists(), p
p.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:report[k] for k in ['outcomes','stages','source_distances_A','enumeration','prepared_config_differences']},ensure_ascii=False,indent=2))
print('TS rise meV',min(t['rise_meV'] for t in saddles),max(t['rise_meV'] for t in saddles))
print('TS heavy RMSD',min(t['heavy_rmsd_A'] for t in saddles),max(t['heavy_rmsd_A'] for t in saddles))
print('TS imaginary',min(v for t in saddles for v in t['imaginary_wavenumbers_icm']),max(v for t in saddles for v in t['imaginary_wavenumbers_icm']))
