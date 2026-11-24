"""占位区间草稿求解与确值结构交付。

草稿沿用 process 的显式仿射操作与源位点（ID、物种、分数坐标），把每条来源的
标量占位换成整数单位区间。调用方给出整数刻度 K（1..10⁹）：来源 s 分到 u_s 个
单位时，其轨道的每个点占位都是 u_s/(mK)，m 为轨道的互异展开点数。约束为：每条
来源单位落在各自非空闭区间内（端点 1..mK）；每个轨道的来源单位之和不超过 mK
（容量），并可加闭区间观测；每个物种的晶胞单位总数（该物种各来源单位之和）可
加闭区间观测（端点 0..K·M，M 为全草稿互异展开点数）。所有来源保留，不把缺测
来源默认丢掉。

轨道按完整展开点集合等价划分（与源坐标写在轨道哪个点无关，同一点由多个操作
产生只贡献一次），物种按标签划分；每条来源恰好是（轨道, 物种）二部图上的一条
边，上述约束即该图上的整数可行循环流（轨道节点、物种节点各带总量区间，边带来
源区间）。边单位数与组内单位和的实际最小/最大值由残余网络的最小割精确给出：
对可行流 x，边 e=(u,v) 的最大值为 x_e+min(h_e-x_e, F(v→u)-(x_e-l_e))，最小值
对称；同轨道同物种的平行边成组时，组和的最值把组边聚合成一条边后同样成立，
且所需残余最大流与单边相同，按割的修正是精确的。全部算术为整数，不枚举可行
分配，时间与空间不随 K 的数值铺开。

交付请求把每个选中轨道的全部来源划分成互不重叠的非空组（同组只能同物种），
并指定组内一个来源作为输出代表。每组单位和固定、且换成每点占位约为最简分数
后落在 process 有理数格式（分子分母 1..100000）内时，交付可直接传给 process
的输入（原操作 + 每组代表来源的定点占位）以及代表到成员的映射；否则返回各阻
塞组的原因、精确范围与达到两个端点的完整可行分配见证，不给部分结构。输入错
误抛 ValueError；正常响应可经 JSON 往返；仅重排操作、来源、观测、组或成员的
同一输入得到相同响应；不修改调用方数据。
"""
from collections import deque
from fractions import Fraction

from .affine import apply
from .checks import array, fields, fraction, ident, integer, unique
from .model import parse_operations

FORMAT = 'occupancy-draft/1'
MAX_SCALE = 10 ** 9
MAX_FRACTION = 100000
_MAX_OPERATIONS = 24
_MAX_SOURCES = 200


def _maxflow(node_count, arcs, source, sink):
    """Dinic 最大流；返回 (流量值, 每条输入弧的流量)。邻接按 arcs 顺序，结果确定。"""
    graph = [[] for _ in range(node_count)]
    positions = []
    for u, v, capacity in arcs:
        positions.append((u, len(graph[u])))
        graph[u].append([v, capacity, len(graph[v])])
        graph[v].append([u, 0, len(graph[u]) - 1])
    total = 0
    infinity = sum(capacity for _, _, capacity in arcs) or 1
    while True:
        level = [-1] * node_count
        level[source] = 0
        queue = deque([source])
        while queue:
            x = queue.popleft()
            for edge in graph[x]:
                if edge[1] > 0 and level[edge[0]] < 0:
                    level[edge[0]] = level[x] + 1
                    queue.append(edge[0])
        if level[sink] < 0:
            break
        it = [0] * node_count

        def push(x, flow):
            if x == sink:
                return flow
            while it[x] < len(graph[x]):
                edge = graph[x][it[x]]
                if edge[1] > 0 and level[edge[0]] == level[x] + 1:
                    sent = push(edge[0], min(flow, edge[1]))
                    if sent:
                        edge[1] -= sent
                        graph[edge[0]][edge[2]][1] += sent
                        return sent
                it[x] += 1
            return 0

        while True:
            pushed = push(source, infinity)
            if not pushed:
                break
            total += pushed
    flows = [capacity - graph[u][index][1]
             for (u, index), (_, _, capacity) in zip(positions, arcs)]
    return total, flows


def _feasible_circulation(node_count, edges):
    """下界/上界循环流可行性；返回与 edges 对齐的整数流量列表，不可行返回 None。"""
    demand = [0] * node_count
    arcs = []
    for u, v, lo, hi, _ in edges:
        demand[u] -= lo
        demand[v] += lo
        arcs.append((u, v, hi - lo))
    super_source, super_sink = node_count, node_count + 1
    required = 0
    for node in range(node_count):
        if demand[node] > 0:
            arcs.append((super_source, node, demand[node]))
            required += demand[node]
        elif demand[node] < 0:
            arcs.append((node, super_sink, -demand[node]))
    value, flows = _maxflow(node_count + 2, arcs, super_source, super_sink)
    if value != required:
        return None
    return [lo + f for (_, _, lo, _, _), f in zip(edges, flows[:len(edges)])]


def _residual_maxflow(node_count, edges, flows, source, sink):
    """当前可行循环流的残余网络上 source→sink 的最大流（用于精确范围端点）。"""
    arcs = []
    for (u, v, lo, hi, _), x in zip(edges, flows):
        if x < hi:
            arcs.append((u, v, hi - x))
        if lo < x:
            arcs.append((v, u, x - lo))
    return _maxflow(node_count, arcs, source, sink)[0]


def _interval(value, low, high, what):
    pair = array(value, 2, 2)
    lo = integer(pair[0], low, high)
    hi = integer(pair[1], low, high)
    if lo > hi:
        raise ValueError(f'{what} interval is empty')
    return lo, hi


def _parse_draft(data):
    fields(data, ('scale', 'operations', 'sites'), ('orbit_observations', 'species_observations'))
    scale = integer(data['scale'], 1, MAX_SCALE)
    operations = parse_operations(data)
    bound = _MAX_OPERATIONS * MAX_SCALE
    sites = []
    for raw in array(data['sites'], 1, _MAX_SOURCES):
        fields(raw, ('id', 'species', 'point', 'units'))
        lo, hi = _interval(raw['units'], 1, bound, 'source units')
        sites.append({'id': ident(raw['id']), 'species': ident(raw['species']),
                      'point': tuple(fraction(v) % 1 for v in array(raw['point'], 3, 3)),
                      'lo': lo, 'hi': hi})
    unique([site['id'] for site in sites])
    sites.sort(key=lambda site: site['id'])
    orbit_observations = []
    for raw in array(data.get('orbit_observations', []), 0, _MAX_SOURCES):
        fields(raw, ('site', 'units'))
        orbit_observations.append({'site': ident(raw['site']),
                                   'units': _interval(raw['units'], 0, bound, 'orbit units')})
    species_bound = MAX_SCALE * _MAX_OPERATIONS * _MAX_SOURCES
    species_observations = []
    for raw in array(data.get('species_observations', []), 0, _MAX_SOURCES):
        fields(raw, ('species', 'units'))
        species_observations.append({'species': ident(raw['species']),
                                     'units': _interval(raw['units'], 0, species_bound, 'species units')})
    unique([observation['species'] for observation in species_observations])
    return operations, sites, orbit_observations, species_observations, scale


def _orbits(operations, sites):
    """按完整展开点集合划分轨道；轨道按最小成员 ID 排序，成员按 ID 排序。"""
    expansions = [frozenset(apply(op, site['point']) for _, op in operations) for site in sites]
    groups = {}
    for index, points in enumerate(expansions):
        groups.setdefault(points, []).append(index)
    orbits = []
    for members in groups.values():
        members = sorted(members, key=lambda i: sites[i]['id'])
        points = expansions[members[0]]
        if len(operations) % len(points):
            raise ValueError('orbit size invariant')
        orbits.append({'members': members, 'points': points, 'm': len(points)})
    orbits.sort(key=lambda orbit: sites[orbit['members'][0]]['id'])
    for orbit_index, orbit in enumerate(orbits):
        for index in orbit['members']:
            sites[index]['orbit'] = orbit_index
    return orbits


def _prepare(data):
    operations, sites, orbit_observations, species_observations, scale = _parse_draft(data)
    orbits = _orbits(operations, sites)
    by_id = {site['id']: index for index, site in enumerate(sites)}
    capacities = [orbit['m'] * scale for orbit in orbits]
    for site in sites:
        if site['hi'] > capacities[site['orbit']]:
            raise ValueError('source units interval exceeds orbit capacity mK')
    orbit_bounds = {}
    for observation in orbit_observations:
        ref = observation['site']
        if ref not in by_id:
            raise ValueError(f'orbit observation references unknown source: {ref}')
        orbit_index = sites[by_id[ref]]['orbit']
        if orbit_index in orbit_bounds:
            raise ValueError('duplicate orbit observation')
        lo, hi = observation['units']
        if hi > capacities[orbit_index]:
            raise ValueError('orbit observation interval exceeds mK')
        orbit_bounds[orbit_index] = (lo, hi)
    cell_points = sum(orbit['m'] for orbit in orbits)
    species_bound = scale * cell_points
    known_species = sorted({site['species'] for site in sites})
    species_bounds = {}
    for observation in species_observations:
        species = observation['species']
        if species not in known_species:
            raise ValueError(f'species observation references unknown species: {species}')
        lo, hi = observation['units']
        if hi > species_bound:
            raise ValueError('species observation interval exceeds cell bound')
        species_bounds[species] = (lo, hi)
    return {'operations': operations, 'sites': sites, 'orbits': orbits, 'scale': scale,
            'cell_points': cell_points, 'species_list': known_species,
            'orbit_bounds': orbit_bounds, 'species_bounds': species_bounds}


def _network(prepared):
    """构建 (轨道, 物种) 二部循环流网络；边为 (u, v, lo, hi, tag)，顺序确定。"""
    sites, orbits, species_list = prepared['sites'], prepared['orbits'], prepared['species_list']
    scale = prepared['scale']
    node_count = 2 + len(orbits) + len(species_list)
    onode = [2 + index for index in range(len(orbits))]
    pnode = {species: 2 + len(orbits) + index for index, species in enumerate(species_list)}
    infinity = (1 + sum(site['hi'] for site in sites)
                + sum(orbit['m'] * scale for orbit in orbits)
                + sum(hi for _, hi in prepared['species_bounds'].values()))
    edges = []
    for index, orbit in enumerate(orbits):
        lo, hi = prepared['orbit_bounds'].get(index, (0, orbit['m'] * scale))
        edges.append((0, onode[index], lo, hi, ('orbit', index)))
    for index, site in enumerate(sites):
        edges.append((onode[site['orbit']], pnode[site['species']],
                      site['lo'], site['hi'], ('source', index)))
    for species in species_list:
        lo, hi = prepared['species_bounds'].get(species, (0, infinity))
        edges.append((pnode[species], 1, lo, hi, ('species', species)))
    edges.append((1, 0, 0, infinity, ('total', None)))
    return node_count, edges, onode, pnode


def _solve(prepared):
    """求一个可行循环流并给出每条来源边的精确单位范围；不可行返回 None。"""
    node_count, edges, onode, pnode = _network(prepared)
    flows = _feasible_circulation(node_count, edges)
    if flows is None:
        return None
    sites = prepared['sites']
    source_flow = {}
    for (_, _, _, _, tag), x in zip(edges, flows):
        if tag[0] == 'source':
            source_flow[tag[1]] = x
    pairs = {}
    for index, site in enumerate(sites):
        pairs.setdefault((site['orbit'], site['species']), []).append(index)
    pair_flow = {}
    ranges = {}
    for (orbit_index, species), members in sorted(pairs.items()):
        f_back = _residual_maxflow(node_count, edges, flows, pnode[species], onode[orbit_index])
        f_fwd = _residual_maxflow(node_count, edges, flows, onode[orbit_index], pnode[species])
        pair_flow[(orbit_index, species)] = (f_back, f_fwd)
        for index in members:
            x, lo, hi = source_flow[index], sites[index]['lo'], sites[index]['hi']
            ranges[index] = (x - min(x - lo, f_fwd - (hi - x)),
                             x + min(hi - x, f_back - (x - lo)))
    return {'node_count': node_count, 'edges': edges, 'onode': onode, 'pnode': pnode,
            'flows': flows, 'source_flow': source_flow, 'pair_flow': pair_flow, 'ranges': ranges}


def _parse_request(request, sites, orbits):
    fields(request, ('groups',))
    by_id = {site['id']: index for index, site in enumerate(sites)}
    seen = set()
    covered = {}
    groups = []
    for raw in array(request['groups'], 1, _MAX_SOURCES):
        fields(raw, ('members', 'representative'))
        member_ids = [ident(member) for member in array(raw['members'], 1, _MAX_SOURCES)]
        unique(member_ids)
        representative = ident(raw['representative'])
        if representative not in member_ids:
            raise ValueError('representative must be a group member')
        members = []
        for member_id in member_ids:
            if member_id not in by_id:
                raise ValueError(f'group references unknown source: {member_id}')
            if member_id in seen:
                raise ValueError('source assigned to more than one group')
            seen.add(member_id)
            members.append(by_id[member_id])
        orbit_ids = {sites[index]['orbit'] for index in members}
        if len(orbit_ids) != 1:
            raise ValueError('group members span multiple orbits')
        species_ids = {sites[index]['species'] for index in members}
        if len(species_ids) != 1:
            raise ValueError('group members span multiple species')
        groups.append({'members': sorted(members, key=lambda i: sites[i]['id']),
                       'representative': by_id[representative],
                       'orbit': orbit_ids.pop(), 'species': species_ids.pop()})
        covered.setdefault(groups[-1]['orbit'], set()).update(member_ids)
    for orbit_index, ids in covered.items():
        if ids != {sites[index]['id'] for index in orbits[orbit_index]['members']}:
            raise ValueError('grouping omits sources of a selected orbit')
    groups.sort(key=lambda group: sites[group['representative']]['id'])
    return groups


def _group_range(group, sites, context):
    """组内单位和在全部完整可行分配上的实际最小、最大值（不按成员界分别相加）。"""
    members = group['members']
    source_flow = context['source_flow']
    total = sum(source_flow[index] for index in members)
    slack = sum(sites[index]['hi'] - source_flow[index] for index in members)
    drop = sum(source_flow[index] - sites[index]['lo'] for index in members)
    f_back, f_fwd = context['pair_flow'][(group['orbit'], group['species'])]
    return total - min(drop, f_fwd - slack), total + min(slack, f_back - drop)


def _witness(prepared, context, group, target):
    """达到组内单位和 target 的完整可行分配（含未选中轨道），按 ID 排序返回。"""
    sites = prepared['sites']
    member = set(group['members'])
    edges = [edge for edge in context['edges']
             if not (edge[4][0] == 'source' and edge[4][1] in member)]
    edges.append((context['onode'][group['orbit']], context['pnode'][group['species']],
                  target, target, ('group', None)))
    flows = _feasible_circulation(context['node_count'], edges)
    assert flows is not None, 'witness target outside computed group range'
    assignment = {}
    for (_, _, _, _, tag), x in zip(edges, flows):
        if tag[0] == 'source':
            assignment[sites[tag[1]]['id']] = x
    remaining = target
    members = group['members']
    for position, index in enumerate(members):
        rest = sum(sites[other]['lo'] for other in members[position + 1:])
        lo, hi = sites[index]['lo'], sites[index]['hi']
        value = max(lo, min(hi, remaining - rest))
        assignment[sites[index]['id']] = value
        remaining -= value
    return assignment


def _structure(prepared, fixed_groups):
    operations = [{'id': name,
                   'matrix': [list(row) for row in matrix],
                   'translation': [[t.numerator, t.denominator] for t in translation]}
                  for name, (matrix, translation) in sorted(prepared['operations'])]
    sites = []
    for group in fixed_groups:
        representative = prepared['sites'][group['representative_index']]
        sites.append({'id': representative['id'], 'species': representative['species'],
                      'point': [[f.numerator, f.denominator] for f in representative['point']],
                      'occupancy': list(group['occupancy'])})
    return {'operations': operations, 'sites': sites}


def _deliver(prepared, context, groups):
    sites, orbits, scale = prepared['sites'], prepared['orbits'], prepared['scale']
    blocking = []
    fixed = []
    for group in groups:
        lo, hi = _group_range(group, sites, context)
        base = {'representative': sites[group['representative']]['id'],
                'members': [sites[index]['id'] for index in group['members']],
                'orbit': group['orbit'], 'species': group['species']}
        if lo != hi:
            blocking.append({**base, 'units': [lo, hi], 'reason': 'not_fixed',
                             'witnesses': {'min': _witness(prepared, context, group, lo),
                                           'max': _witness(prepared, context, group, hi)}})
            continue
        occupancy = Fraction(lo, orbits[group['orbit']]['m'] * scale)
        numerator, denominator = occupancy.numerator, occupancy.denominator
        if numerator > MAX_FRACTION or denominator > MAX_FRACTION:
            witness = _witness(prepared, context, group, lo)
            blocking.append({**base, 'units': [lo, hi], 'occupancy': [numerator, denominator],
                             'reason': 'occupancy_unrepresentable',
                             'witnesses': {'min': witness, 'max': witness}})
            continue
        fixed.append({'representative_index': group['representative'],
                      'occupancy': [numerator, denominator],
                      'entry': {**base, 'units': lo, 'occupancy': [numerator, denominator]}})
    if blocking:
        return {'status': 'blocked', 'blocking': blocking}
    return {'status': 'delivered', 'structure': _structure(prepared, fixed),
            'groups': [group['entry'] for group in fixed]}


def solve_draft(data, delivery=None):
    """求解占位区间草稿；给出 delivery 请求时同时判定能否交付确值结构。

    无 delivery 时只做草稿求解。格式非法抛 ValueError；约束无解返回
    {'format': 'occupancy-draft/1', 'status': 'infeasible'}；可行时返回每条
    来源在全晶胞约束下实际可达的单位范围。delivery 为 {'groups': [...]}，
    每组 {'members': [...], 'representative': id}；结果附加 'delivery' 字段，
    或为可直接传给 process 的 'delivered' 结构，或为带原因、精确范围与端点
    见证的 'blocked' 报告。不修改传入的数据。
    """
    prepared = _prepare(data)
    groups = _parse_request(delivery, prepared['sites'], prepared['orbits']) \
        if delivery is not None else None
    context = _solve(prepared)
    if context is None:
        return {'format': FORMAT, 'status': 'infeasible'}
    sites, orbits = prepared['sites'], prepared['orbits']
    result = {'format': FORMAT,
              'status': 'feasible',
              'scale': prepared['scale'],
              'cell_points': prepared['cell_points'],
              'orbits': [{'index': index, 'multiplicity': orbit['m'],
                          'sites': [sites[member]['id'] for member in orbit['members']]}
                         for index, orbit in enumerate(orbits)],
              'sources': [{'id': site['id'], 'species': site['species'], 'orbit': site['orbit'],
                           'units': list(context['ranges'][index])}
                          for index, site in enumerate(sites)]}
    if groups is not None:
        result['delivery'] = _deliver(prepared, context, groups)
    return result
