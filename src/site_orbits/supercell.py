"""超胞分片生成与离线拼接。

在 process 的原胞展开语义之上，按新基底把超胞切成胞片（旧原胞平移模新基底
整数张成的陪集）分别交付；接收方不持有原胞输入，仅凭分片包离线合并，跨胞片
的有向连接以出端/入端记录按实例配对。全部算术为精确有理数，不随未选中胞片
数或基本平行六面体的外接整数盒增长。
"""
import hashlib
import json
import re
from fractions import Fraction

from .affine import determinant
from .checks import array, fields, fraction, ident, integer, unique
from .expand import expand
from .model import parse

FORMAT = 'supercell-shard/1'
BOUND = 10 ** 9
CELL_BOUND = 3 * 10 ** 9
MAX_CELLS = 64
MAX_LINKS = 1000
MAX_POSITIONS = 4800
MAX_RECORDS = MAX_LINKS * MAX_CELLS
_VERSION = re.compile(r'sha256:[0-9a-f]{64}')
_RATIONAL = re.compile(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?')


def _adjugate(matrix):
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return ((e * i - f * h, c * h - b * i, b * f - c * e),
            (f * g - d * i, a * i - c * g, c * d - a * f),
            (d * h - e * g, b * g - a * h, a * e - b * d))


def _image(adj, det, vec):
    """B⁻¹·vec 的精确有理数三元组。"""
    return tuple(Fraction(sum(adj[i][j] * vec[j] for j in range(3)), det) for i in range(3))


def _floor(value):
    return value.numerator // value.denominator


def _canon(basis, adj, det, translation):
    """旧胞平移的规范胞片代表：basis⁻¹c 各分量落在 [0,1) 的唯一整数向量。"""
    f = tuple(_floor(u) for u in _image(adj, det, translation))
    return tuple(translation[i] - sum(basis[i][j] * f[j] for j in range(3)) for i in range(3))


def _check_basis(value):
    basis = tuple(tuple(integer(v, -BOUND, BOUND) for v in array(row, 3, 3)) for row in array(value, 3, 3))
    det = determinant(basis)
    if not 1 <= det <= BOUND:
        raise ValueError('basis determinant must be in 1..10^9')
    return basis, det


def _cell(value, bound):
    return tuple(integer(v, -bound, bound) for v in array(value, 3, 3))


def _rational(value):
    if not isinstance(value, str) or not _RATIONAL.fullmatch(value):
        raise ValueError('canonical rational string required')
    result = Fraction(value)
    if str(result) != value:
        raise ValueError('canonical rational string required')
    return result


def _unit(value):
    result = _rational(value)
    if not 0 <= result < 1:
        raise ValueError('rational in [0,1) required')
    return result


def _positive(value):
    result = _rational(value)
    if not 0 < result <= 1:
        raise ValueError('rational in (0,1] required')
    return result


def _integer(value):
    if type(value) is not int:
        raise ValueError('integer required')
    return value


def _check_links(value, expanded):
    links = []
    for raw in array(value, 0, MAX_LINKS):
        fields(raw, ('id', 'from', 'to', 'shift'))
        lid = ident(raw['id'])
        ends = []
        for key in ('from', 'to'):
            point = tuple(fraction(v) % 1 for v in array(raw[key], 3, 3))
            if point not in expanded:
                raise ValueError('link endpoint is not an expanded position')
            ends.append(point)
        shift = tuple(integer(v, -BOUND, BOUND) for v in array(raw['shift'], 3, 3))
        links.append({'id': lid, 'from': ends[0], 'to': ends[1], 'shift': shift})
    unique([link['id'] for link in links])
    return links


def _version(basis, positions, links):
    declared = [{'id': link['id'],
                 'from': [str(v) for v in link['from']],
                 'to': [str(v) for v in link['to']],
                 'shift': list(link['shift'])} for link in sorted(links, key=lambda link: link['id'])]
    payload = {'basis': [list(row) for row in basis], 'links': declared, 'positions': positions}
    blob = json.dumps(payload, sort_keys=True, separators=(',', ':')).encode('utf-8')
    return 'sha256:' + hashlib.sha256(blob).hexdigest()


def _winding(adj, det, start, shift, source, target):
    head = _image(adj, det, tuple(start[i] + shift[i] + target[i] for i in range(3)))
    tail = _image(adj, det, tuple(start[i] + source[i] for i in range(3)))
    return tuple(_floor(head[i]) - _floor(tail[i]) for i in range(3))


def _record(lid, start, source, end, target, winding):
    return {'link': lid,
            'cell': list(start),
            'from': {'cell': list(start), 'point': [str(v) for v in source]},
            'to': {'cell': list(end), 'point': [str(v) for v in target]},
            'winding': list(winding)}


def build_supercell_shard(data, basis, links, cells):
    """生成超胞分片包；输入校验与语义见 docs/contract.md。失败抛 ValueError。"""
    operations, sites = parse(data)
    positions, _, _ = expand(operations, sites)
    basis, det = _check_basis(basis)
    adj = _adjugate(basis)
    requested = [_cell(raw, BOUND) for raw in array(cells, 1, MAX_CELLS)]
    coverage = sorted({_canon(basis, adj, det, c) for c in requested})
    if len(coverage) != len(requested):
        raise ValueError('duplicate cell shard in request')
    expanded = {tuple(Fraction(v) for v in rec['point']) for rec in positions}
    declarations = _check_links(links, expanded)
    version = _version(basis, positions, declarations)
    points = []
    for cell in coverage:
        for rec in positions:
            source = tuple(Fraction(v) for v in rec['point'])
            position = tuple(v % 1 for v in _image(adj, det, tuple(cell[i] + source[i] for i in range(3))))
            points.append({'cell': list(cell),
                           'point': list(rec['point']),
                           'position': [str(v) for v in position],
                           'occupancy': rec['occupancy'],
                           'vacancy': rec['vacancy'],
                           'contributions': [{'site': c['site'], 'species': c['species'],
                                              'occupancy': c['occupancy'],
                                              'operations': list(c['operations'])}
                                             for c in rec['contributions']]})
    selected = set(coverage)
    complete, outward, inward = [], [], []
    for declaration in declarations:
        lid = declaration['id']
        source, target, shift = declaration['from'], declaration['to'], declaration['shift']
        for cell in coverage:
            end = _canon(basis, adj, det, tuple(cell[i] + shift[i] for i in range(3)))
            record = _record(lid, cell, source, end, target, _winding(adj, det, cell, shift, source, target))
            (complete if end in selected else outward).append(record)
        for cell in coverage:
            start = _canon(basis, adj, det, tuple(cell[i] - shift[i] for i in range(3)))
            if start not in selected:
                inward.append(_record(lid, start, source, cell, target, _winding(adj, det, start, shift, source, target)))
    key = lambda record: (record['link'], record['cell'])
    return {'format': FORMAT,
            'version': version,
            'basis': [list(row) for row in basis],
            'cell_count': det,
            'cells': [list(cell) for cell in coverage],
            'points': points,
            'links': sorted(complete, key=key),
            'out_links': sorted(outward, key=key),
            'in_links': sorted(inward, key=key)}


def _read_endpoint(raw, basis, adj, det):
    fields(raw, ('cell', 'point'))
    cell = _cell(raw['cell'], CELL_BOUND)
    if _canon(basis, adj, det, cell) != cell:
        raise ValueError('endpoint cell is not the canonical representative')
    point = tuple(_unit(v) for v in array(raw['point'], 3, 3))
    return cell, point


def _read_point(raw, adj, det, coverage):
    fields(raw, ('cell', 'point', 'position', 'occupancy', 'vacancy', 'contributions'))
    cell = _cell(raw['cell'], CELL_BOUND)
    if cell not in coverage:
        raise ValueError('point cell outside shard coverage')
    point = tuple(_unit(v) for v in array(raw['point'], 3, 3))
    position = tuple(_unit(v) for v in array(raw['position'], 3, 3))
    expect = tuple(v % 1 for v in _image(adj, det, tuple(cell[i] + point[i] for i in range(3))))
    if position != expect:
        raise ValueError('position inconsistent with cell and point')
    occupancy = _positive(raw['occupancy'])
    vacancy = _unit(raw['vacancy'])
    if occupancy + vacancy != 1:
        raise ValueError('occupancy and vacancy do not sum to one')
    contributions = []
    for raw_contribution in array(raw['contributions'], 1, 200):
        fields(raw_contribution, ('site', 'species', 'occupancy', 'operations'))
        operations = [ident(op) for op in array(raw_contribution['operations'], 1, 24)]
        unique(operations)
        contributions.append({'site': ident(raw_contribution['site']),
                              'species': ident(raw_contribution['species']),
                              'occupancy': _positive(raw_contribution['occupancy']),
                              'operations': operations})
    unique([c['site'] for c in contributions])
    if sum((c['occupancy'] for c in contributions), Fraction(0)) != occupancy:
        raise ValueError('contributions do not sum to point occupancy')
    record = {'cell': list(cell),
              'point': [str(v) for v in point],
              'position': [str(v) for v in position],
              'occupancy': str(occupancy),
              'vacancy': str(vacancy),
              'contributions': [{'site': c['site'], 'species': c['species'],
                                 'occupancy': str(c['occupancy']),
                                 'operations': sorted(c['operations'])}
                                for c in sorted(contributions, key=lambda c: c['site'])]}
    return (cell, point), record


def _read_links(raw_list, kind, basis, adj, det, coverage, points):
    records = {}
    direction = {'complete': (True, True), 'out': (True, False), 'in': (False, True)}[kind]
    for raw in array(raw_list, 0, MAX_RECORDS):
        fields(raw, ('link', 'cell', 'from', 'to', 'winding'))
        lid = ident(raw['link'])
        cell = _cell(raw['cell'], CELL_BOUND)
        if _canon(basis, adj, det, cell) != cell:
            raise ValueError('start cell is not the canonical representative')
        from_end = _read_endpoint(raw['from'], basis, adj, det)
        to_end = _read_endpoint(raw['to'], basis, adj, det)
        if from_end[0] != cell:
            raise ValueError('instance cell must be the start cell')
        winding = tuple(_integer(v) for v in array(raw['winding'], 3, 3))
        if (from_end[0] in coverage, to_end[0] in coverage) != direction:
            raise ValueError('link record has wrong direction for its class')
        for end, local in ((from_end, direction[0]), (to_end, direction[1])):
            if local and end not in points:
                raise ValueError('link endpoint is not a listed point')
        iid = (lid, cell)
        if iid in records:
            raise ValueError('duplicate link instance')
        records[iid] = {'link': lid,
                        'cell': list(cell),
                        'from': {'cell': list(from_end[0]), 'point': [str(v) for v in from_end[1]]},
                        'to': {'cell': list(to_end[0]), 'point': [str(v) for v in to_end[1]]},
                        'winding': list(winding)}
    return records


def _read_package(pkg):
    fields(pkg, ('format', 'version', 'basis', 'cell_count', 'cells', 'points', 'links', 'out_links', 'in_links'))
    if pkg['format'] != FORMAT:
        raise ValueError('unsupported shard format')
    if not isinstance(pkg['version'], str) or not _VERSION.fullmatch(pkg['version']):
        raise ValueError('invalid structure version tag')
    basis, det = _check_basis(pkg['basis'])
    if type(pkg['cell_count']) is not int or pkg['cell_count'] != det:
        raise ValueError('cell count does not match basis determinant')
    adj = _adjugate(basis)
    cells = []
    for raw in array(pkg['cells'], 1, MAX_CELLS):
        cell = _cell(raw, CELL_BOUND)
        if _canon(basis, adj, det, cell) != cell:
            raise ValueError('cell is not the canonical representative')
        cells.append(cell)
    unique(cells)
    coverage = set(cells)
    points = {}
    for raw in array(pkg['points'], 1, MAX_CELLS * MAX_POSITIONS):
        key, record = _read_point(raw, adj, det, coverage)
        if key in points:
            raise ValueError('duplicate point')
        points[key] = record
    if {key[0] for key in points} != coverage:
        raise ValueError('covered cell carries no points')
    complete = _read_links(pkg['links'], 'complete', basis, adj, det, coverage, points)
    outward = _read_links(pkg['out_links'], 'out', basis, adj, det, coverage, points)
    inward = _read_links(pkg['in_links'], 'in', basis, adj, det, coverage, points)
    if set(complete) & set(outward) or set(complete) & set(inward) or set(outward) & set(inward):
        raise ValueError('link instance recorded in conflicting states')
    return {'version': pkg['version'], 'basis': basis, 'det': det,
            'cells': cells, 'points': points, 'complete': complete, 'out': outward, 'in': inward}


def assemble_supercell_shards(shards, *, finalize=False):
    """离线合并分片包；finalize=True 仅在覆盖齐全且无未配边界记录时返回。"""
    if type(finalize) is not bool:
        raise ValueError('finalize must be a bool')
    if not isinstance(shards, (list, tuple)) or not shards:
        raise ValueError('at least one shard package required')
    parsed = [_read_package(shard) for shard in shards]
    head = parsed[0]
    for other in parsed[1:]:
        if other['version'] != head['version'] or other['basis'] != head['basis']:
            raise ValueError('shards must share structure version and basis')
    cells, seen = [], set()
    points, complete, outward, inward = {}, {}, {}, {}
    for part in parsed:
        for cell in part['cells']:
            if cell in seen:
                raise ValueError('overlapping shard coverage')
            seen.add(cell)
            cells.append(cell)
        if len(cells) > MAX_CELLS:
            raise ValueError('merged coverage exceeds 64 cells')
        for key, record in part['points'].items():
            points[key] = record
        for target, source in ((complete, part['complete']), (outward, part['out']), (inward, part['in'])):
            for iid, record in source.items():
                if iid in target:
                    raise ValueError('duplicate link instance record')
                target[iid] = record
    if set(complete) & (set(outward) | set(inward)):
        raise ValueError('link instance recorded in conflicting states')
    links = dict(complete)
    out_rest, in_rest = {}, {}
    for iid, record in outward.items():
        mate = inward.get(iid)
        if mate is None:
            out_rest[iid] = record
        elif mate != record:
            raise ValueError('inconsistent link instance halves')
        else:
            links[iid] = record
    for iid, record in inward.items():
        if iid not in outward:
            in_rest[iid] = record
    for record in out_rest.values():
        if tuple(record['to']['cell']) in seen:
            raise ValueError('missing inward half inside merged coverage')
    for record in in_rest.values():
        if tuple(record['from']['cell']) in seen:
            raise ValueError('missing outward half inside merged coverage')
    result = {'format': FORMAT,
              'version': head['version'],
              'basis': [list(row) for row in head['basis']],
              'cell_count': head['det'],
              'cells': [list(cell) for cell in sorted(cells)],
              'points': [points[key] for key in sorted(points)],
              'links': [links[iid] for iid in sorted(links)],
              'out_links': [out_rest[iid] for iid in sorted(out_rest)],
              'in_links': [in_rest[iid] for iid in sorted(in_rest)]}
    if finalize and (len(cells) != head['det'] or result['out_links'] or result['in_links']):
        raise ValueError('finalize requires full coverage and no unpaired boundary records')
    return result
