"""Recover affine symmetries and an asymmetric unit from complete sites.

The example uses metric diag(1, 2, 3) and the eight diagonal sign matrices.
Two full-occupancy X sites at x=0 and x=1/2 admit both zero and half-cell
translations for every sign matrix: 16 operations in total.  The lex-smallest
point keeps its one source record as the compressed asymmetric unit.  The
published document contains process-compatible operations/sites and a complete
per-source trace.

Run from the repository root:

    PYTHONPATH=src python3 examples/run_recover.py
"""
import json
from fractions import Fraction

from site_orbits import process, recover_asymmetric_unit
from site_orbits.io import publish


def fraction_pair(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


METRIC = [
    [fraction_pair(1), [0, 1], [0, 1]],
    [[0, 1], fraction_pair(2), [0, 1]],
    [[0, 1], [0, 1], fraction_pair(3)],
]

SIGN_GROUP = [
    [[a, 0, 0], [0, b, 0], [0, 0, c]]
    for a in (-1, 1)
    for b in (-1, 1)
    for c in (-1, 1)
]

COMPLETE_SITES = [
    {
        'id': 'origin-x',
        'species': 'X',
        'point': [[0, 1], [0, 1], [0, 1]],
        'occupancy': [1, 1],
    },
    {
        'id': 'half-x',
        'species': 'X',
        'point': [[1, 2], [0, 1], [0, 1]],
        'occupancy': [1, 1],
    },
]


if __name__ == '__main__':
    recovered = recover_asymmetric_unit(COMPLETE_SITES, METRIC, SIGN_GROUP)
    expanded = process({
        'operations': recovered['operations'],
        'sites': recovered['sites'],
    })
    publish('run/recover-asymmetric-unit.json', recovered)
    print(json.dumps({
        'format': recovered['format'],
        'operations': len(recovered['operations']),
        'asymmetric_sites': len(recovered['sites']),
        'source_map_records': len(recovered['source_map']),
        'expanded_positions': len(expanded['positions']),
        'cell_composition': expanded['cell_composition'],
        'output': 'run/recover-asymmetric-unit.json',
    }, sort_keys=True))
