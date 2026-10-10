import json,shutil,subprocess,tempfile
from pathlib import Path
root=Path('/home/ruigengji/MLPATH'); python='/home/ruigengji/miniforge3/envs/openmm_dev/bin/python'; results=[]
for mutation,test in [('interleaved','test_coarse_first.py'),('normalized_scraps','test_restricted_modes.py'),('unconstrained_ring','test_constrained_delivery.py')]:
 case=Path(tempfile.mkdtemp(prefix='mlpath-mechanism-mutation-')); shutil.copytree(root/'src/prrs',case/'src/prrs',ignore=shutil.ignore_patterns('__pycache__')); (case/'tests').mkdir();shutil.copy2(root/'tests'/test,case/'tests'/test);shutil.copy2(root/'pyproject.toml',case/'pyproject.toml')
 if mutation=='interleaved':
  p=case/'src/prrs/search.py';s=p.read_text();assert 'pending_scans.append((probes, coarse, direction))' in s;s=s.replace('pending_scans.append((probes, coarse, direction))','finish_scan(node, state, probes, coarse, proposal, schedule, direction)');p.write_text(s);extra=['-k','exact_coarse_budget']
 elif mutation=='normalized_scraps':
  p=case/'src/prrs/runner.py';s=p.read_text();a=s.index('def _restricted_soft_modes(');b=s.index('\n\ndef soft_internal_modes(',a);s=s[:a]+'''def _restricted_soft_modes(eigenvalues, vectors, excluded, limit):
    q = np.linalg.qr(excluded.T)[0]
    result = []
    for value, mode in zip(eigenvalues, vectors.T):
        residual = mode - q @ (q.T @ mode)
        norm = np.linalg.norm(residual)
        if norm >= 1e-8:
            result.append((float(value), residual / norm))
        if len(result) >= limit:
            break
    return tuple(result)
'''+s[b:];p.write_text(s);extra=[]
 else:
  p=case/'src/prrs/perturbations.py';s=p.read_text();old='internal.constrained_displace(\n            atoms, probe.kind, probe.indices, requested, graph.edges\n        )';assert old in s;s=s.replace(old,'internal.displace(atoms, probe.kind, probe.indices, requested)');p.write_text(s);extra=['-k','ring_delivery']
 out=subprocess.run([python,'-m','pytest','-q','tests/'+test,*extra],cwd=case,capture_output=True,text=True)
 (case/'result.log').write_text(out.stdout+out.stderr)
 row={'mutation':mutation,'exit_code':out.returncode,'summary':out.stdout.strip().splitlines()[-1] if out.stdout else out.stderr[-500:],'log':str(case/'result.log')};results.append(row);print(json.dumps(row),flush=True)
 assert out.returncode==1, row
Path('/tmp/mlpath-mechanism-mutations.json').write_text(json.dumps(results,indent=2)+'\n')
