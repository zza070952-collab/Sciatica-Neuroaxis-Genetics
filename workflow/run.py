#!/usr/bin/env python3
"""Sequential entry point with explicit environment and prerequisite checks."""
import argparse, json, os, subprocess
from pathlib import Path
import yaml

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(); ap.add_argument('--stage',default='all'); ap.add_argument('--dry-run',action='store_true'); ap.add_argument('--conda',action='store_true')
    a=ap.parse_args(); os.chdir(root)
    stages=yaml.safe_load((root/'workflow/stages.yaml').read_text())
    names=list(stages) if a.stage=='all' else a.stage.split(',')
    if any(n not in stages for n in names): ap.error('Unknown stage')
    report=[]
    for name in names:
        spec=stages[name]; missing=[p for p in spec['prerequisites'] if not (root/p).exists()]
        report.append({'stage':name,'environment':spec['environment'],'missing_inputs':missing,'commands':spec['commands']})
        if a.dry_run: continue
        if missing: raise SystemExit(f'{name}: missing external or upstream inputs: '+', '.join(missing))
        (root/'logs').mkdir(exist_ok=True)
        for cmd in spec['commands']:
            args=['bash','-lc',cmd]
            if a.conda: args=['conda','run','--no-capture-output','-n','sciatica-'+spec['environment']]+args
            env=dict(os.environ,SCIATICA_PROJECT_ROOT=str(root))
            with (root/'logs'/f'release_{name}.log').open('a') as log:
                log.write(cmd+'\n'); log.flush()
                subprocess.run(args,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
