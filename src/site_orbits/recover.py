"""Recover affine site symmetries from a complete decorated unit cell.

The input sites are the complete periodic contents of a reference cell, rather
than a known asymmetric unit.  This module validates an exact rational metric
and a caller-supplied linear point group, searches every admissible fractional
translation, and keeps one lexicographically minimum coordinate representative
per recovered affine orbit.
"""
import copy
from fractions import Fraction

from .affine import IDENTITY, apply, determinant
from .checks import array, fraction as parse_fraction, integer, unique
from .model import parse as parse_process_input
from .model import parse_operations


_IDENTITY_MATRIX = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_ZERO = Fraction(0)


def _pair(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _point_pairs(point):
    return [_pair(component) for component in point]


def _matrix_key(matrix):
    return tuple(value for row in matrix for value in row)


def _multiply_matrices(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def _parse_metric(metric):
    rows = array(metric, 3, 3)
    result = []
    for row in rows:
        result.append(tuple(parse_fraction(value) for value in array(row, 3, 3)))

    for i in range(3):
        for j in range(3):
            if result[i][j] != result[j][i]:
                raise ValueError('metric must be symmetric')

    # Sylvester's criterion, evaluated with exact rational arithmetic.
    if result[0][0] <= 0:
        raise ValueError('metric must be positive definite')
    leading_minor_2 = result[0][0] * result[1][1] - result[0][1] * result[1][0]
    if leading_minor_2 <= 0:
        raise ValueError('metric must be positive definite')
    if determinant(result) <= 0:
        raise ValueError('metric must be positive definite')
    return tuple(result)


def _metric_preserving(matrix, metric):
    for i in range(3):
        for j in range(3):
            value = _ZERO
            for a in range(3):
                for b in range(3):
                    value += matrix[a][i] * metric[a][b] * matrix[b][j]
            if value != metric[i][j]:
                return False
    return True


def _parse_candidate_matrices(candidate_matrices, metric):
    raw_matrices = array(candidate_matrices, 1, 24)
    matrices = []
    for raw_matrix in raw_matrices:
        rows = array(raw_matrix, 3, 3)
        matrix = tuple(
            tuple(integer(value, -2, 2) for value in array(row, 3, 3))
            for row in rows
        )
        if abs(determinant(matrix)) != 1:
            raise ValueError('candidate matrix must be unimodular')
        matrices.append(matrix)

    unique(matrices)
    matrix_set = set(matrices)
    if _IDENTITY_MATRIX not in matrix_set:
        raise ValueError('candidate matrix group must contain identity')

    for left in matrices:
        for right in matrices:
            if _multiply_matrices(left, right) not in matrix_set:
                raise ValueError('candidate matrices are not closed under multiplication')

    for matrix in matrices:
        if not _metric_preserving(matrix, metric):
            raise ValueError('candidate matrix does not preserve metric')

    return sorted(matrices, key=_matrix_key)


def _decorate_sites(sites):
    # Reuse all process site validation and exact coordinate normalization.
    _, parsed_sites = parse_process_input({
        'operations': [{
            'id': 'identity',
            'matrix': [list(row) for row in _IDENTITY_MATRIX],
            'translation': [[0, 1], [0, 1], [0, 1]],
        }],
        'sites': sites,
    })

    decorations = {}
    for site in parsed_sites:
        point = site['point']
        species = site['species']
        decoration = decorations.setdefault(point, {})
        decoration[species] = decoration.get(species, _ZERO) + site['occupancy']
        if sum(decoration.values(), _ZERO) > 1:
            raise ValueError('overoccupied complete coordinate decoration')

    raw_by_id = {site['id']: site for site in sites}
    return parsed_sites, decorations, raw_by_id


def _find_operations(matrices, decorations):
    anchor = min(decorations)
    anchor_decoration = decorations[anchor]
    found = []
    seen = set()

    for matrix in matrices:
        linear_images = {
            point: tuple(
                sum(matrix[i][j] * point[j] for j in range(3))
                for i in range(3)
            )
            for point in decorations
        }
        anchor_linear = linear_images[anchor]

        for target, target_decoration in sorted(decorations.items()):
            if target_decoration != anchor_decoration:
                continue
            translation = tuple(
                (target[i] - anchor_linear[i]) % 1 for i in range(3)
            )
            operation = (matrix, translation)
            if operation in seen:
                continue

            valid = True
            images = set()
            for point, decoration in decorations.items():
                image = tuple(
                    (linear_images[point][i] + translation[i]) % 1
                    for i in range(3)
                )
                if decorations.get(image) != decoration:
                    valid = False
                    break
                images.add(image)

            # Explicit injectivity check; unimodular matrices already imply it
            # on the torus, but keeping it local makes the finite-set contract
            # clear.
            if not valid or len(images) != len(decorations):
                continue

            seen.add(operation)
            found.append(operation)
            if len(found) > 24:
                raise ValueError('recovered operation set exceeds process limit of 24')

    if IDENTITY not in seen:
        raise ValueError('internal error: identity symmetry was not recovered')

    found.sort(key=lambda operation: (_matrix_key(operation[0]), operation[1]))
    return found


def _operation_records(operations):
    records = []
    indexed = []
    other_number = 0
    for operation in operations:
        matrix, translation = operation
        if operation == IDENTITY:
            operation_id = 'identity'
        else:
            other_number += 1
            operation_id = f'op-{other_number}'
        record = {
            'id': operation_id,
            'matrix': [list(row) for row in matrix],
            'translation': _point_pairs(translation),
        }
        records.append(record)
        indexed.append((operation_id, operation))

    # Enforce exactly the downstream process operation contract, including
    # closure and exact numerator/denominator bounds on derived translations.
    parse_operations({'operations': records})
    return records, indexed


def _orbit_representatives(operations, decorations):
    representative_of = {}
    for point in decorations:
        orbit = {apply(operation, point) for operation in operations}
        representative = min(orbit)
        for member in orbit:
            previous = representative_of.get(member)
            if previous is not None and previous != representative:
                raise ValueError('internal error: inconsistent affine orbit')
            representative_of[member] = representative
    return representative_of


def _occupancy_pair(raw_site):
    occupancy = raw_site['occupancy']
    return [occupancy[0], occupancy[1]]


def _site_records(parsed_sites, representative_of, raw_by_id):
    kept = [
        site for site in parsed_sites
        if site['point'] == representative_of[site['point']]
    ]
    kept.sort(key=lambda site: (representative_of[site['point']], site['id']))

    records = []
    for site in kept:
        raw = raw_by_id[site['id']]
        records.append({
            'id': site['id'],
            'species': site['species'],
            'point': _point_pairs(site['point']),
            'occupancy': _occupancy_pair(raw),
        })
    return records


def _source_map(parsed_sites, representative_of, indexed_operations, raw_by_id):
    result = []
    for site in sorted(parsed_sites, key=lambda item: item['id']):
        original = site['point']
        representative = representative_of[original]
        operation_id = None
        for candidate_id, operation in indexed_operations:
            if apply(operation, representative) == original:
                operation_id = candidate_id
                break
        if operation_id is None:
            raise ValueError('internal error: no operation traces source to representative')

        original_point = _point_pairs(original)
        representative_point = _point_pairs(representative)
        site_id = site['id']
        entry = {
            'id': site_id,
            'species': site['species'],
            'occupancy': copy.deepcopy(_occupancy_pair(raw_by_id[site_id])),
            'original_point': original_point,
            'representative_point': representative_point,
            'operation_id': operation_id,
        }
        # Vocabulary aliases for the same stable trace record.
        entry.update({
            'site': site_id,
            'site_id': site_id,
            'source_id': site_id,
            'source': site_id,
            'point': original_point,
            'coordinate': original_point,
            'original': original_point,
            'normalized_point': original_point,
            'normalized_coordinate': original_point,
            'representative': representative_point,
            'representative_coordinate': representative_point,
            'operation': operation_id,
            'op_id': operation_id,
        })
        result.append(entry)
    return result


def recover_asymmetric_unit(sites, metric, candidate_matrices):
    """Recover all decorated affine symmetries and an asymmetric unit.

    The return value contains process-compatible ``operations`` and ``sites``
    plus a deterministic per-source ``source_map``.  Invalid input, a metric
    not preserved by every candidate matrix, or a recovered operation set that
    cannot satisfy process limits raises :class:`ValueError`.
    """
    parsed_metric = _parse_metric(metric)
    matrices = _parse_candidate_matrices(candidate_matrices, parsed_metric)
    parsed_sites, decorations, raw_by_id = _decorate_sites(sites)
    operations = _find_operations(matrices, decorations)
    operation_records, indexed_operations = _operation_records(operations)
    representative_of = _orbit_representatives(operations, decorations)
    asymmetric_sites = _site_records(parsed_sites, representative_of, raw_by_id)
    source_map = _source_map(
        parsed_sites, representative_of, indexed_operations, raw_by_id
    )

    return {
        'format': 'asymmetric-unit/1',
        'operations': operation_records,
        'sites': asymmetric_sites,
        'source_map': source_map,
    }
