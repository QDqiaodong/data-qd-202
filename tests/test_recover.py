import copy
import json
import unittest
from fractions import Fraction

from site_orbits import process, recover_asymmetric_unit
from site_orbits.affine import compose
from site_orbits.model import parse_operations


def pair(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


def metric_diagonal(*values):
    return [
        [pair(values[i]) if i == j else [0, 1] for j in range(3)]
        for i in range(3)
    ]


def sign_matrices():
    return [
        [[a, 0, 0], [0, b, 0], [0, 0, c]]
        for a in (-1, 1)
        for b in (-1, 1)
        for c in (-1, 1)
    ]


def site(site_id, species, point, occupancy=(1, 1)):
    if not isinstance(point[0], list):
        point = [pair(value) for value in point]
    return {
        'id': site_id,
        'species': species,
        'point': point,
        'occupancy': list(occupancy),
    }


IDENTITY_OPERATION = {
    'id': 'identity',
    'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    'translation': [[0, 1], [0, 1], [0, 1]],
}


def expanded_species(data):
    result = process(data)
    points = {}
    for position in result['positions']:
        coordinate = tuple(Fraction(value) for value in position['point'])
        totals = {}
        for contribution in position['contributions']:
            species = contribution['species']
            totals[species] = (
                totals.get(species, Fraction(0))
                + Fraction(contribution['occupancy'])
            )
        points[coordinate] = totals
    return points, result['cell_composition']


class RecoverAsymmetricUnitTests(unittest.TestCase):
    def setUp(self):
        self.metric = metric_diagonal(1, 2, 3)
        self.matrices = sign_matrices()
        self.sites = [
            site('a', 'X', [0, 0, 0]),
            site('b', 'X', [Fraction(1, 2), 0, 0]),
        ]

    def test_finds_all_sixteen_translations_and_compresses_sources(self):
        original = copy.deepcopy(self.sites)
        result = recover_asymmetric_unit(self.sites, self.metric, self.matrices)

        self.assertEqual(result['format'], 'asymmetric-unit/1')
        self.assertEqual(len(result['operations']), 16)
        self.assertEqual(len(result['sites']), 1)
        self.assertEqual(result['sites'][0]['id'], 'a')
        self.assertEqual(result['sites'][0]['point'],
                         [[0, 1], [0, 1], [0, 1]])

        identity = [op for op in result['operations'] if op['id'] == 'identity']
        self.assertEqual(len(identity), 1)
        self.assertEqual(identity[0]['matrix'], IDENTITY_OPERATION['matrix'])
        self.assertEqual(identity[0]['translation'], IDENTITY_OPERATION['translation'])

        # Operations are sorted by flattened matrix and exact translation.
        sort_keys = []
        for op in result['operations']:
            matrix_key = tuple(value for row in op['matrix'] for value in row)
            translation_key = tuple(Fraction(n, d) for n, d in op['translation'])
            sort_keys.append((matrix_key, translation_key))
        self.assertEqual(sort_keys, sorted(sort_keys))

        # The emitted operation set passes the downstream closure contract.
        parsed = parse_operations({'operations': result['operations']})
        operation_keys = {key for _, key in parsed}
        for _, left in parsed:
            for _, right in parsed:
                self.assertIn(compose(left, right), operation_keys)

        compressed_input = {
            'operations': result['operations'],
            'sites': result['sites'],
        }
        compressed = process(compressed_input)
        self.assertEqual(compressed['cell_composition'], {'X': '2'})
        self.assertEqual(len(compressed['positions']), 2)
        self.assertEqual(
            [(p['point'], p['occupancy']) for p in compressed['positions']],
            [(['0', '0', '0'], '1'), (['1/2', '0', '0'], '1')],
        )

        original_points, original_composition = expanded_species({
            'operations': [IDENTITY_OPERATION],
            'sites': original,
        })
        compressed_points, compressed_composition = expanded_species(compressed_input)
        self.assertEqual(compressed_points, original_points)
        self.assertEqual(compressed_composition, original_composition)

        original_process = process({
            'operations': [IDENTITY_OPERATION],
            'sites': original,
        })
        compressed_process = process(compressed_input)
        self.assertEqual(
            [(p['point'], p['occupancy'], p['vacancy'])
             for p in compressed_process['positions']],
            [(p['point'], p['occupancy'], p['vacancy'])
             for p in original_process['positions']],
        )

        records = {row['id']: row for row in result['source_map']}
        self.assertEqual(set(records), {'a', 'b'})
        trace = records['b']
        self.assertEqual(trace['point'], [[1, 2], [0, 1], [0, 1]])
        self.assertEqual(trace['representative_point'], [[0, 1], [0, 1], [0, 1]])

        # The recorded ID is exactly the first sorted operation mapping the
        # representative to the original coordinate.
        matching_ids = []
        for op in result['operations']:
            matrix = tuple(tuple(row) for row in op['matrix'])
            translation = tuple(Fraction(n, d) for n, d in op['translation'])
            mapped = tuple(
                (sum(matrix[i][j] * Fraction(trace['representative_point'][j][0],
                                             trace['representative_point'][j][1])
                     for j in range(3)) + translation[i]) % 1
                for i in range(3)
            )
            if mapped == (Fraction(1, 2), Fraction(0), Fraction(0)):
                matching_ids.append(op['id'])
        self.assertTrue(matching_ids)
        self.assertEqual(trace['operation_id'], matching_ids[0])
        self.assertEqual(trace['operation_id'], 'op-2')

        # Reordering either input list does not change the deterministic value.
        reordered = recover_asymmetric_unit(
            list(reversed(self.sites)),
            self.metric,
            list(reversed(self.matrices)),
        )
        self.assertEqual(reordered, result)

        # JSON round trip preserves the complete contract and process input.
        round_tripped = json.loads(json.dumps(result, sort_keys=True))
        self.assertEqual(round_tripped, result)
        process({
            'operations': round_tripped['operations'],
            'sites': round_tripped['sites'],
        })

        # The caller's data is not modified.
        self.assertEqual(self.sites, original)

    def test_species_difference_prevents_orbit_merge(self):
        sites = [
            site('a', 'X', [0, 0, 0]),
            site('b', 'Y', [Fraction(1, 2), 0, 0]),
        ]
        result = recover_asymmetric_unit(sites, self.metric, self.matrices)

        self.assertEqual(len(result['operations']), 8)
        self.assertEqual(len(result['sites']), 2)
        self.assertEqual([item['id'] for item in result['sites']], ['a', 'b'])
        self.assertTrue(all(
            all(component == [0, 1] for component in op['translation'])
            for op in result['operations']
        ))

        output = process({
            'operations': result['operations'],
            'sites': result['sites'],
        })
        self.assertEqual(output['cell_composition'], {'X': '1', 'Y': '1'})
        by_point = {p['point'][0]: p for p in output['positions']}
        self.assertEqual(by_point['0']['contributions'][0]['species'], 'X')
        self.assertEqual(by_point['1/2']['contributions'][0]['species'], 'Y')

    def test_equal_aggregate_vectors_with_different_source_splits_merge(self):
        sites = [
            site('a1', 'X', [0, 0, 0], (1, 3)),
            site('a2', 'X', [0, 0, 0], (2, 3)),
            site('b1', 'X', [Fraction(1, 2), 0, 0], (1, 2)),
            site('b2', 'X', [Fraction(1, 2), 0, 0], (1, 2)),
        ]
        result = recover_asymmetric_unit(sites, self.metric, self.matrices)

        self.assertEqual(len(result['operations']), 16)
        self.assertEqual([item['id'] for item in result['sites']], ['a1', 'a2'])
        self.assertEqual(len(result['source_map']), 4)
        records = {row['id']: row for row in result['source_map']}
        self.assertEqual(records['b1']['occupancy'], [1, 2])
        self.assertEqual(records['b2']['occupancy'], [1, 2])
        self.assertEqual(records['b1']['operation_id'], records['b2']['operation_id'])

        compressed_input = {
            'operations': result['operations'],
            'sites': result['sites'],
        }
        original_points, original_composition = expanded_species({
            'operations': [IDENTITY_OPERATION],
            'sites': sites,
        })
        compressed_points, compressed_composition = expanded_species(compressed_input)
        self.assertEqual(compressed_points, original_points)
        self.assertEqual(compressed_composition, original_composition)

    def test_closed_group_not_preserving_metric_is_rejected(self):
        swap_xy = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, self.metric,
                                    [IDENTITY_OPERATION['matrix'], swap_xy])

    def test_anchor_only_candidate_translation_is_not_listed(self):
        sites = [
            site('a', 'X', [0, 0, 0]),
            site('b', 'X', [Fraction(1, 2), 0, 0]),
            site('c', 'X', [Fraction(1, 4), 0, 0]),
        ]
        result = recover_asymmetric_unit(sites, self.metric, self.matrices)
        x_reflection = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
        translations = [
            tuple(Fraction(n, d) for n, d in op['translation'])
            for op in result['operations']
            if op['matrix'] == x_reflection
        ]
        self.assertIn((Fraction(1, 2), Fraction(0), Fraction(0)), translations)
        self.assertNotIn((Fraction(1, 4), Fraction(0), Fraction(0)), translations)

        compressed_input = {
            'operations': result['operations'],
            'sites': result['sites'],
        }
        original_points, _ = expanded_species({
            'operations': [IDENTITY_OPERATION], 'sites': sites,
        })
        compressed_points, _ = expanded_species(compressed_input)
        self.assertEqual(compressed_points, original_points)

    def test_invalid_metric_and_candidate_groups_are_rejected(self):
        asymmetric = metric_diagonal(1, 2, 3)
        asymmetric[0][1] = [1, 1]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, asymmetric, self.matrices)

        non_positive = metric_diagonal(0, 2, 3)
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, non_positive, self.matrices)

        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, [[[1, 1]] * 3] * 2, self.matrices)

        duplicate_matrices = self.matrices + [self.matrices[0]]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, self.metric, duplicate_matrices)

        nonclosed = [IDENTITY_OPERATION['matrix'],
                     [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
                     [[0, 0, 1], [0, 1, 0], [1, 0, 0]]]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, self.metric, nonclosed)

        bad_determinant = [[[2, 0, 0], [0, 1, 0], [0, 0, 1]]]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(self.sites, self.metric, bad_determinant)

    def test_overoccupied_complete_decoration_is_rejected(self):
        sites = [
            site('a', 'X', [0, 0, 0], (3, 4)),
            site('b', 'Y', [0, 0, 0], (3, 4)),
        ]
        with self.assertRaises(ValueError):
            recover_asymmetric_unit(sites, self.metric, [IDENTITY_OPERATION['matrix']])

    def test_operation_limit_is_enforced(self):
        metric = metric_diagonal(1, 1, 1)
        matrices = [IDENTITY_OPERATION['matrix']]

        acceptable = [
            site(f's{i:02d}', 'X', [Fraction(i, 24), 0, 0])
            for i in range(24)
        ]
        result = recover_asymmetric_unit(acceptable, metric, matrices)
        self.assertEqual(len(result['operations']), 24)
        process({'operations': result['operations'], 'sites': result['sites']})

        too_many = [
            site(f's{i:02d}', 'X', [Fraction(i, 25), 0, 0])
            for i in range(25)
        ]
        with self.assertRaisesRegex(ValueError, '24'):
            recover_asymmetric_unit(too_many, metric, matrices)


if __name__ == '__main__':
    unittest.main()
