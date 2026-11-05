"""PYTHONPATH=../DynaComp/src python examples/dynacomp_return.py [--replay file]."""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.verification.dynacomp_bridge import sampled_step, periodic_certificate, replay_periodic


def examples():
    cases = {'scalar': [([['-1']], '1')],
             'nonnormal': [([['-1','3'],['0','-2']], '1/2')],
             'noncommuting': [([['-1','2'],['0','-2']], '1/2'),
                              ([['-2','0'],['1','-1']], '1/3')]}
    return {name: periodic_certificate([sampled_step(a,h)[1] for a,h in steps])
            for name,steps in cases.items()}

if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--replay':
        data = json.loads(Path(sys.argv[2]).read_text())
        if not data or not all(replay_periodic(r) for r in data.values()):
            raise SystemExit('replay failed')
        print('all certificates replayed')
    else:
        print(json.dumps(examples(), indent=2, sort_keys=True))
