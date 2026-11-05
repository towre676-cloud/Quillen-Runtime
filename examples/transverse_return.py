"""Reproduce or replay source-bound transverse certificates."""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.verification.dynacomp_bridge import sampled_step
from src.verification.transverse import sampled_return_certificate, replay


def examples():
    _, scalar = sampled_step([[-1, 0], [0, 0]], 1)
    a = [[-1, 1, 0], [0, -1, 0], [0, 0, 0]]
    b = [[-1, 0, 0], [1, -1, 0], [0, 0, 0]]
    steps = [sampled_step(x, 1)[1] for x in (a, b)]
    return dict(
        diagonal=sampled_return_certificate([scalar], [[1], [0]], [[1, 0]], neutral=[[0], [1]]),
        oblique=sampled_return_certificate([sampled_step([[-1,0],[-1,0]],1)[1]],
                                          [[2],[2]], [[1,'-1/2']], neutral=[[0],[1]]),
        noncommuting=sampled_return_certificate(steps, [[1,0],[0,1],[0,0]],
                                               [[1,0,0],[0,1,0]], neutral=[[0],[0],[1]]))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--replay':
        data = json.loads(Path(sys.argv[2]).read_text())
        if set(data) != {'diagonal','oblique','noncommuting'} or not all(replay(c) for c in data.values()):
            raise SystemExit('replay rejected')
        print('all three transverse certificates replayed')
    elif len(sys.argv) == 1:
        print(json.dumps(examples(), indent=2, sort_keys=True))
    else:
        raise SystemExit('usage: transverse_return.py [--replay RECEIPT]')
