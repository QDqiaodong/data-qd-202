"""周期结构整组抽取与离线作业交付。

在 finalize 验收的超胞包（supercell-shard/1）之上，把包沿基底三个方向重复
(nx, ny, nz) 次：包内点（胞片代表 cell、原胞展开点 point）加各分量落在
[0, n) 的重复像编号 replica 标识本次作业的一个实际位点。包中每条完整有向
连接从任一重复像 r 出发给出一条指向 (r+winding) mod n 的有向边；连通组按
无向可达定义。调用方按优先顺序给出请求点、可给保留标记并给作业额度：同组
请求归并为一项，按请求首次出现顺序逐组判定——组内任一实际位点带保留标记
则整组保留（不输出、不扣额度，即使原本超额）；否则仅当整组互异实际位点都
能放进剩余额度时才交付（按位点数扣额度），不够则记额度不足并继续后面的组；
不截部分组。

重复像上的图是包内点图（顶点为 (cell, point)，边为完整连接记录、带整向量
winding）在 Z_n = Z_n1×Z_n2×Z_n3 上的提升。固定包内连通组 C 中任一树根，
记 σ(key) 为根到 key 的树路径 winding 和；闭合路径的 winding 全体构成子群
H ≤ Z_n，用格 L = ⟨各记录的 σ(from)+w−σ(to), n1e1, n2e2, n3e3⟩（Z³ 子格，
H ≅ L/(n1Z×n2Z×n3Z)）精确表示。两个实际位点同组当且仅当包内点同组且
r−σ(key) 模 L 同余；组大小为 |C|·|H|，其中 |H| = n1n2n3/det L。保留判定
与额度判定只做同余类算术：不展开重复网格，不扫描无请求的连通组，不靠遍历
超额组判断同组；只有确定交付的组才枚举其 |C|·|H| 个位点（不超过额度上限）
与连接实例。全部算术为精确整数/有理数。

交付组自包含：原结构版本、包基底、重复次数、全部实际位点（身份、来源贡献、
占位、空位、重复后基底中的规范分数坐标 (p+r)/n）与全部连接实例（包内实例
身份加起点重复像、两个实际端点、重复后基底中的整向量绕行量，满足
终点规范坐标 + 绕行量 − 起点规范坐标 = 包内有向位移逐分量除以重复次数，
不由最短距离猜测；取模后两端相同的边不合并）。输入错误抛 ValueError，不
修改调用方数据；结果可经 JSON 往返；包内记录与保留标记重排不改变结果，
请求顺序具有业务含义。
"""
import re
from collections import deque
from fractions import Fraction

from .checks import array, fields, ident, integer, unique
from .supercell import CELL_BOUND, assemble_supercell_shards

FORMAT = 'supercell-extract/1'
MAX_REPEATS = 10 ** 9
MAX_QUOTA = 10000
MAX_REQUESTS = 10000
MAX_RESERVES = 10000
_RATIONAL = re.compile(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?')


def _unit(value):
    if not isinstance(value, str) or not _RATIONAL.fullmatch(value):
        raise ValueError('canonical rational string required')
    result = Fraction(value)
    if str(result) != value or not 0 <= result < 1:
        raise ValueError('canonical rational in [0,1) required')
    return result


def _xgcd(a, b):
    """返回 (g, s, t)：g = gcd(a, b) >= 0 且 s·a + t·b = g。"""
    r0, r1, s0, s1, t0, t1 = a, b, 1, 0, 0, 1
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if r0 < 0:
        r0, s0, t0 = -r0, -s0, -t0
    return r0, s0, t0


class _Lattice:
    """含 repeats 对角的 Z³ 子格；列式下三角基、正对角，非对角按所在行对角约化。"""

    def __init__(self, repeats):
        self.cols = [[repeats[0], 0, 0], [0, repeats[1], 0], [0, 0, repeats[2]]]

    def insert(self, vector):
        v = list(vector)
        for j in range(3):
            b = v[j]
            if b == 0:
                continue
            col = self.cols[j]
            g, s, t = _xgcd(col[j], b)
            self.cols[j] = [s * col[k] + t * v[k] for k in range(3)]
            f, q = col[j] // g, b // g
            v = [f * v[k] - q * col[k] for k in range(3)]

    def canonicalize(self):
        for j in range(2):
            for i in range(j + 1, 3):
                q = self.cols[j][i] // self.cols[i][i]
                if q:
                    self.cols[j] = [self.cols[j][k] - q * self.cols[i][k] for k in range(3)]

    def determinant(self):
        return self.cols[0][0] * self.cols[1][1] * self.cols[2][2]

    def reduce(self, vector):
        """vector 模本格的规范代表：第 i 分量落在 [0, h_ii)。"""
        v = list(vector)
        for j in range(3):
            q = v[j] // self.cols[j][j]
            if q:
                v = [v[k] - q * self.cols[j][k] for k in range(3)]
        return tuple(v)

    def coset(self, c, repeats):
        """枚举 c + L 落在 [0, n) 盒中的全部元素（恰 |H| 个）。"""
        h = (self.cols[0][0], self.cols[1][1], self.cols[2][2])
        for a0 in range(repeats[0] // h[0]):
            r0 = c[0] + a0 * h[0]
            shift1 = c[1] + a0 * self.cols[0][1]
            base1 = shift1 % h[1]
            q1 = (shift1 - base1) // h[1]
            for k1 in range(repeats[1] // h[1]):
                shift2 = c[2] + a0 * self.cols[0][2] + (k1 - q1) * self.cols[1][2]
                base2 = shift2 % h[2]
                for k2 in range(repeats[2] // h[2]):
                    yield (r0, base1 + k1 * h[1], base2 + k2 * h[2])


def _key(cell, point):
    return tuple(cell), tuple(Fraction(v) for v in point)


def _parse_entry(raw, repeats, points):
    fields(raw, ('id', 'cell', 'point', 'replica'))
    entry_id = ident(raw['id'])
    cell = tuple(integer(v, -CELL_BOUND, CELL_BOUND) for v in array(raw['cell'], 3, 3))
    point = tuple(_unit(v) for v in array(raw['point'], 3, 3))
    replica = tuple(integer(v, 0, repeats[i] - 1) for i, v in enumerate(array(raw['replica'], 3, 3)))
    if (cell, point) not in points:
        raise ValueError('entry references an unknown package point')
    return entry_id, (cell, point), replica


def _explore(start, adjacency, repeats):
    """包内点图一个连通组：树势 σ、闭合 winding 格 L（已约化）、组内连接记录。"""
    sigma = {start: (0, 0, 0)}
    lattice = _Lattice(repeats)
    records = []
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v, w, record, forward in adjacency.get(u, ()):
            if forward:
                records.append(record)
            if v not in sigma:
                sigma[v] = tuple(sigma[u][i] + w[i] for i in range(3))
                queue.append(v)
            else:
                lattice.insert(tuple(sigma[u][i] + w[i] - sigma[v][i] for i in range(3)))
    lattice.canonicalize()
    return {'sigma': sigma, 'lattice': lattice, 'records': records}


def _site(record, cell, point, replica, repeats):
    position = tuple(Fraction(v) for v in record['position'])
    return {'cell': list(cell),
            'point': [str(v) for v in point],
            'replica': list(replica),
            'position': [str((position[i] + replica[i]) / repeats[i]) for i in range(3)],
            'occupancy': record['occupancy'],
            'vacancy': record['vacancy'],
            'contributions': [{'site': c['site'], 'species': c['species'],
                               'occupancy': c['occupancy'], 'operations': list(c['operations'])}
                              for c in record['contributions']]}


def _deliver(comp, inv, points, repeats, final, request_ids):
    lattice, sigma = comp['lattice'], comp['sigma']
    sites = []
    for key in sigma:
        cell, point = key
        c = lattice.reduce(tuple(inv[i] + sigma[key][i] for i in range(3)))
        for replica in lattice.coset(c, repeats):
            sites.append((cell, point, replica))
    sites.sort()
    out_sites = [_site(points[key[:2]], key[0], key[1], key[2], repeats) for key in sites]
    links = []
    for record in comp['records']:
        source = _key(record['from']['cell'], record['from']['point'])
        target = _key(record['to']['cell'], record['to']['point'])
        w = tuple(record['winding'])
        c = lattice.reduce(tuple(inv[i] + sigma[source][i] for i in range(3)))
        for replica in lattice.coset(c, repeats):
            landed = tuple((replica[i] + w[i]) % repeats[i] for i in range(3))
            winding = tuple((replica[i] + w[i] - landed[i]) // repeats[i] for i in range(3))
            links.append({'link': record['link'],
                          'cell': list(source[0]),
                          'replica': list(replica),
                          'from': {'cell': list(source[0]), 'point': [str(v) for v in source[1]],
                                   'replica': list(replica)},
                          'to': {'cell': list(target[0]), 'point': [str(v) for v in target[1]],
                                 'replica': list(landed)},
                          'winding': list(winding)})
    links.sort(key=lambda link: (link['link'], link['cell'], link['replica']))
    return {'requests': list(request_ids),
            'status': 'delivered',
            'version': final['version'],
            'basis': [list(row) for row in final['basis']],
            'repeats': list(repeats),
            'sites': out_sites,
            'links': links}


def extract_supercell_groups(package, repeats, requests, reserves=None, quota=None):
    """从 finalize 验收的超胞包整组抽取本次作业交付；语义见 docs/contract.md。"""
    final = assemble_supercell_shards([package], finalize=True)
    repeats = tuple(integer(v, 1, MAX_REPEATS) for v in array(repeats, 3, 3))
    quota = integer(quota, 1, MAX_QUOTA)
    points = {}
    for record in final['points']:
        points[_key(record['cell'], record['point'])] = record
    adjacency = {}
    for record in final['links']:
        source = _key(record['from']['cell'], record['from']['point'])
        target = _key(record['to']['cell'], record['to']['point'])
        w = tuple(record['winding'])
        adjacency.setdefault(source, []).append((target, w, record, True))
        adjacency.setdefault(target, []).append((source, tuple(-v for v in w), record, False))
    request_entries = [_parse_entry(raw, repeats, points) for raw in array(requests, 1, MAX_REQUESTS)]
    unique([entry[0] for entry in request_entries])
    if reserves is None:
        reserves = []
    reserve_entries = [_parse_entry(raw, repeats, points) for raw in array(reserves, 0, MAX_RESERVES)]
    unique([entry[0] for entry in reserve_entries])
    components = []
    comp_of = {}
    groups = []
    group_index = {}
    for entry_id, key, replica in request_entries:
        comp_id = comp_of.get(key)
        if comp_id is None:
            comp_id = len(components)
            comp = _explore(key, adjacency, repeats)
            components.append(comp)
            for member in comp['sigma']:
                comp_of[member] = comp_id
        comp = components[comp_id]
        inv = comp['lattice'].reduce(tuple(replica[i] - comp['sigma'][key][i] for i in range(3)))
        gkey = (comp_id, inv)
        if gkey in group_index:
            groups[group_index[gkey]]['requests'].append(entry_id)
        else:
            group_index[gkey] = len(groups)
            groups.append({'requests': [entry_id], 'comp': comp_id, 'inv': inv, 'reserves': []})
    for entry_id, key, replica in reserve_entries:
        comp_id = comp_of.get(key)
        if comp_id is None:
            continue
        comp = components[comp_id]
        inv = comp['lattice'].reduce(tuple(replica[i] - comp['sigma'][key][i] for i in range(3)))
        owner = group_index.get((comp_id, inv))
        if owner is not None:
            groups[owner]['reserves'].append(entry_id)
    remaining = quota
    volume = repeats[0] * repeats[1] * repeats[2]
    out_groups = []
    for group in groups:
        reserved = sorted(group['reserves'])
        if reserved:
            out_groups.append({'requests': group['requests'], 'status': 'reserved', 'reserves': reserved})
            continue
        comp = components[group['comp']]
        size = len(comp['sigma']) * (volume // comp['lattice'].determinant())
        if size > remaining:
            out_groups.append({'requests': group['requests'], 'status': 'quota_exceeded'})
            continue
        remaining -= size
        out_groups.append(_deliver(comp, group['inv'], points, repeats, final, group['requests']))
    return {'format': FORMAT,
            'version': final['version'],
            'basis': [list(row) for row in final['basis']],
            'repeats': list(repeats),
            'groups': out_groups}
