#!/usr/bin/env python3
"""Contract/syntax checks plus a real execution of the production aggregation functions.

Analytical unit inputs are software test values, not simulated study results.
No network request or third-party dataset download is made.
"""
import ast, json, math, subprocess, sys
from pathlib import Path
import numpy as np
import yaml

root=Path(__file__).resolve().parents[1]
config=yaml.safe_load((root/'config/config.example.yaml').read_text())
assert config['frozen_rank']=={'TXNL1':1,'MAPK3':2,'FGFR3':3}
stages=yaml.safe_load((root/'workflow/stages.yaml').read_text())
assert len(stages)==11
count=0
for p in (root/'scripts').rglob('*.py'):
    ast.parse(p.read_text(encoding='utf-8'),filename=str(p.relative_to(root))); count+=1
for spec in stages.values():
    assert (root/'envs'/f"{spec['environment']}.yaml").exists()
    for command in spec['commands']:
        path=command.split()[1]
        assert (root/path).is_file(),path
for p in (root/'scripts').rglob('*'):
    assert not p.name.endswith(('.pyc','.xlsx','.gz','.rds'))

source=root/'scripts/07_gsmap/parse_phase3_gsmap.py'
tree=ast.parse(source.read_text())
funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'bh','acat'}]
namespace={'np':np}
exec(compile(ast.Module(body=funcs,type_ignores=[]),'production_aggregation_functions','exec'),namespace)
assert math.isclose(namespace['acat']([.5,.5]),.5,abs_tol=1e-14)
assert np.allclose(namespace['bh']([.01,.04,.03]),[.03,.04,.04])
result={'status':'PASS','python_syntax_files':count,'configuration':'PASS','stage_script_resolution':'PASS','production_ACAT_identity':'PASS','production_BH_known_order':'PASS','full_data_workflow_rerun':False}
out=root/'results/smoke';out.mkdir(parents=True,exist_ok=True)
(out/'smoke_output.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
