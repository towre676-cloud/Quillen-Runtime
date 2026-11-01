"""Run: python examples/validated_return_map.py [--replay receipt.json]."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.verification.matrix import MatrixEnclosure, certify, replay, periodic_return_enclosure


def example():
    # Coupled nonnormal monodromy: eigenvalues 1/2, singular values differ.
    # These step boxes stand for externally validated integration outputs.
    step = MatrixEnclosure([['499/1000', '1999/1000'], [0, '499/1000']],
                           [['501/1000', '2001/1000'], [0, '501/1000']])
    return certify(periodic_return_enclosure([step]), candidate=[1, 0], rhs=['1/2', 0])


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--replay':
        valid = replay(json.loads(Path(sys.argv[2]).read_text()))
        print('VALID' if valid else 'INVALID')
        raise SystemExit(0 if valid else 1)
    if len(sys.argv) != 1:
        raise SystemExit(__doc__)
    print(json.dumps(example(), indent=2))
