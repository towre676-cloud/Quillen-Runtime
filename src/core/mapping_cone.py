"""Linear initial-value and periodic defects; affine forcing is separate."""

def defect_matrix(steps, dimension, periodic=False):
    """steps=[S1,...,Sk]. IV has k+1 blocks; periodic has k blocks.

    Entries may be Fraction objects. Fixed u0 is not part of this linear map.
    """
    if dimension < 1 or any(len(s) != dimension or any(len(r) != dimension for r in s) for s in steps):
        raise ValueError('invalid block dimensions')
    if periodic and not steps:
        raise ValueError('periodic closure requires a step')
    blocks = len(steps) if periodic else len(steps) + 1
    a = [[int(i == j) for j in range(blocks*dimension)] for i in range(blocks*dimension)]
    for k, s in enumerate(steps):
        row = (k+1) % blocks if periodic else k+1
        for i in range(dimension):
            for j in range(dimension):
                a[row*dimension+i][k*dimension+j] -= s[i][j]
    return a
