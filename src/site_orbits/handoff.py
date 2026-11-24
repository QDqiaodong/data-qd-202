"""周期迁移轨迹的定位校准与时间窗交接。

在 finalize 验收的超胞包（supercell-shard/1）之上解释一份压缩轨迹：轨迹由步
（单条包内连接实例）、串接（先左后右）与重复（同一子程序连续执行至少两次）
组成有向无环图，根不指定绝对起始像编号。调用方一次给出若干时间窗，每窗带
两端独立实测像编号与三轴闭区间作业区域。处理顺序固定为：

1. 以左端实测像为绝对起点，按轨迹压缩算术推得右端像；与右端实测不符即
   「测量冲突」，不查出界、不交接，绝不丢弃任一实测值。
2. 两端相符后，对窗内全部步后状态做精确区域判定；越界时返回最早越界步的
   根步号、包内位点身份、实际像编号与精确实际坐标，不交付越界前缀。
3. 全窗通过才输出可机读、自包含的交接路线：共享的压缩定义表只含获准路线
   实际执行的定义，接收方无需原包、原轨迹或进程缓存即可重建有向路线。

全部算术为精确整数/有理数。单窗的定位、冲突与越界判定工作量只受图深与
整数位长约束：不展开轨迹、不按重复次数逐轮推进、不扫描到首个越界步，也
不为每个窗口重扫整张结点表。重复结点内同一子程序状态随迭代轮号线性变化，
整轮可行的轮号构成闭整数区间，越界轮号边界由子程序状态的相对坐标最值
（自底向上一次性汇总）配合整除直接算出。输入错误抛 ValueError，不修改
调用方数据；业务结果（冲突、越界）互不影响；结果可经 JSON 往返。
"""
import hashlib
import json
import re
from fractions import Fraction

from .checks import array, fields, ident, integer, unique
from .supercell import CELL_BOUND, assemble_supercell_shards

FORMAT = 'orbit-handoff/1'
MAX_NODES = 5000
MAX_DEPTH = 48
MAX_STEPS = 10 ** 18
MAX_WINDOWS = 2000
IMAGE_BOUND = 10 ** 60
_RATIONAL = re.compile(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?')
_STEP = 'step'
_SEQUENCE = 'sequence'
_REPEAT = 'repeat'


def _rational(value):
    if not isinstance(value, str) or not _RATIONAL.fullmatch(value):
        raise ValueError('canonical rational string required')
    result = Fraction(value)
    if str(result) != value:
        raise ValueError('canonical rational string required')
    return result


def _plain_int(value, low=None, high=None):
    if type(value) is not int:
        raise ValueError('integer required')
    if low is not None and value < low:
        raise ValueError('integer out of range')
    if high is not None and value > high:
        raise ValueError('integer out of range')
    return value


def _image_triple(value):
    triple = tuple(_plain_int(v) for v in array(value, 3, 3))
    if any(abs(v) > IMAGE_BOUND for v in triple):
        raise ValueError('image component exceeds 10^60 in absolute value')
    return triple


def _cell_triple(value):
    return tuple(integer(v, -CELL_BOUND, CELL_BOUND) for v in array(value, 3, 3))


def _floor(value):
    return value.numerator // value.denominator


def _ceil(value):
    return -((-value.numerator) // value.denominator)


class _Node:
    """轨迹结点的内部表示；summary 为局部状态 p(site)+累计 winding 的逐轴最值。"""

    def __init__(self, nid, kind, steps, start, end, winding, rel_min, rel_max, height):
        self.id = nid
        self.kind = kind
        self.steps = steps
        self.start = start
        self.end = end
        self.winding = winding
        self.rel_min = rel_min
        self.rel_max = rel_max
        self.height = height


class _StepNode(_Node):
    def __init__(self, nid, instance, record, positions):
        self.instance = instance
        self.record = record
        start = (tuple(record['from']['cell']), tuple(Fraction(v) for v in record['from']['point']))
        end = (tuple(record['to']['cell']), tuple(Fraction(v) for v in record['to']['point']))
        winding = tuple(record['winding'])
        p_from = positions[start]
        p_to = positions[end]
        rel_min = tuple(min(p_from[i], p_to[i] + winding[i]) for i in range(3))
        rel_max = tuple(max(p_from[i], p_to[i] + winding[i]) for i in range(3))
        super().__init__(nid, _STEP, 1, start, end, winding, rel_min, rel_max, 1)


class _SequenceNode(_Node):
    def __init__(self, nid, left, right):
        if left.end != right.start:
            raise ValueError('sequence children do not close on an in-package site')
        steps = left.steps + right.steps
        winding = tuple(left.winding[i] + right.winding[i] for i in range(3))
        rel_min = tuple(min(left.rel_min[i], right.rel_min[i] + left.winding[i]) for i in range(3))
        rel_max = tuple(max(left.rel_max[i], right.rel_max[i] + left.winding[i]) for i in range(3))
        super().__init__(nid, _SEQUENCE, steps, left.start, right.end, winding,
                         rel_min, rel_max, 1 + max(left.height, right.height))
        self.left = left
        self.right = right


class _RepeatNode(_Node):
    def __init__(self, nid, child, count):
        if child.start != child.end:
            raise ValueError('repeat child does not close on an in-package site')
        self.child = child
        steps = count * child.steps
        winding = tuple(count * w for w in child.winding)
        span = (count - 1)
        rel_min = tuple(child.rel_min[i] + min(0, span * child.winding[i]) for i in range(3))
        rel_max = tuple(child.rel_max[i] + max(0, span * child.winding[i]) for i in range(3))
        super().__init__(nid, _REPEAT, steps, child.start, child.end, winding,
                         rel_min, rel_max, 1 + child.height)


def _parse_trajectory(raw, links_by_instance, positions):
    fields(raw, ('root', 'nodes'))
    root_id = ident(raw['root'])
    entries = array(raw['nodes'], 1, MAX_NODES)
    node_ids = []
    parsed_raw = {}
    for entry in entries:
        fields(entry, ('id', 'kind'), ('link', 'cell', 'left', 'right', 'child', 'count'))
        nid = ident(entry['id'])
        kind = entry['kind']
        if kind not in (_STEP, _SEQUENCE, _REPEAT):
            raise ValueError('unknown trajectory node kind')
        if kind == _STEP:
            fields(entry, ('id', 'kind', 'link', 'cell'))
            lid = ident(entry['link'])
            cell = _cell_triple(entry['cell'])
            instance = (lid, cell)
            if instance not in links_by_instance:
                raise ValueError('step node references an unknown link instance')
            parsed_raw[nid] = ('step', instance)
        elif kind == _SEQUENCE:
            fields(entry, ('id', 'kind', 'left', 'right'))
            parsed_raw[nid] = ('sequence', ident(entry['left']), ident(entry['right']))
        else:
            fields(entry, ('id', 'kind', 'child', 'count'))
            parsed_raw[nid] = ('repeat', ident(entry['child']),
                               _plain_int(entry['count'], 2, MAX_STEPS))
        node_ids.append(nid)
    unique(node_ids)
    if root_id not in parsed_raw:
        raise ValueError('root is not among node definitions')

    built = {}
    visiting = set()

    def build(nid):
        if nid in built:
            return built[nid]
        if nid not in parsed_raw:
            raise ValueError(f'node references unknown definition: {nid}')
        if nid in visiting:
            raise ValueError('trajectory node graph contains a cycle')
        visiting.add(nid)
        spec = parsed_raw[nid]
        if spec[0] == 'step':
            node = _StepNode(nid, spec[1], links_by_instance[spec[1]], positions)
        elif spec[0] == 'sequence':
            node = _SequenceNode(nid, build(spec[1]), build(spec[2]))
        else:
            node = _RepeatNode(nid, build(spec[1]), spec[2])
        visiting.discard(nid)
        if node.steps > MAX_STEPS:
            raise ValueError('trajectory exceeds 10^18 root steps')
        built[nid] = node
        return node

    root = build(root_id)
    if len(built) != len(parsed_raw):
        raise ValueError('not all nodes are reachable from the root')
    if root.height > MAX_DEPTH:
        raise ValueError('root-to-leaf path exceeds 48 nodes')
    return root


def _prefix_winding(node, t):
    """结点局部第 t 步（0..steps）相对结点起点的累计绕行。"""
    if node.kind == _STEP:
        return tuple(node.winding) if t == 1 else (0, 0, 0)
    if node.kind == _SEQUENCE:
        cut = node.left.steps
        if t <= cut:
            return _prefix_winding(node.left, t)
        tail = _prefix_winding(node.right, t - cut)
        return tuple(node.left.winding[i] + tail[i] for i in range(3))
    q, u = divmod(t, node.child.steps)
    tail = _prefix_winding(node.child, u)
    return tuple(q * node.child.winding[i] + tail[i] for i in range(3))


def _state_at(node, t):
    """结点局部第 t 步的包内位点身份。"""
    if node.kind == _STEP:
        return node.end if t == 1 else node.start
    if node.kind == _SEQUENCE:
        cut = node.left.steps
        if t <= cut:
            return _state_at(node.left, t)
        return _state_at(node.right, t - cut)
    _, u = divmod(t, node.child.steps)
    return _state_at(node.child, u)


def _segment_winding(node, a, b):
    """局部步区间 (a, b] 的净绕行。"""
    wa = _prefix_winding(node, a)
    wb = _prefix_winding(node, b)
    return tuple(wb[i] - wa[i] for i in range(3))


def _inside(coordinate, region):
    return all(region[i][0] <= coordinate[i] <= region[i][1] for i in range(3))


def _violation(node, t, base, positions):
    site = _state_at(node, t)
    image = tuple(base[i] + _prefix_winding(node, t)[i] for i in range(3))
    position = positions[site]
    coordinate = tuple(position[i] + image[i] for i in range(3))
    return site, image, coordinate


def _earliest_violation(node, a, b, base, region, positions):
    """局部步闭区间 [a,b] 内最早越界状态；无越界返回 None。

    base 为结点第 0 步的绝对像编号。重复结点的整轮区间用线性最值与整除
    直接定界，不逐轮推进、不逐状态扫描。
    """
    if node.kind == _STEP:
        for t in range(a, b + 1):
            site, image, coordinate = _violation(node, t, base, positions)
            if not _inside(coordinate, region):
                return t, site, image, coordinate
        return None
    if node.kind == _SEQUENCE:
        cut = node.left.steps
        if b <= cut:
            return _earliest_violation(node.left, a, b, base, region, positions)
        if a >= cut:
            right_base = tuple(base[i] + node.left.winding[i] for i in range(3))
            found = _earliest_violation(node.right, a - cut, b - cut, right_base, region, positions)
            if found is None:
                return None
            return found[0] + cut, found[1], found[2], found[3]
        found = _earliest_violation(node.left, a, cut, base, region, positions)
        if found is not None:
            return found
        right_base = tuple(base[i] + node.left.winding[i] for i in range(3))
        found = _earliest_violation(node.right, 0, b - cut, right_base, region, positions)
        if found is None:
            return None
        return found[0] + cut, found[1], found[2], found[3]

    child = node.child
    m = child.steps
    w = node.child.winding
    k0, u0 = divmod(a, m)
    k1, u1 = divmod(b, m)
    if k0 == k1:
        child_base = tuple(base[i] + k0 * w[i] for i in range(3))
        found = _earliest_violation(child, u0, u1, child_base, region, positions)
        if found is None:
            return None
        return found[0] + k0 * m, found[1], found[2], found[3]

    # 首轮：只可能是局部切片。
    if u0 > 0:
        child_base = tuple(base[i] + k0 * w[i] for i in range(3))
        found = _earliest_violation(child, u0, m, child_base, region, positions)
        if found is not None:
            return found[0] + k0 * m, found[1], found[2], found[3]
    full_lo = k0 + (1 if u0 > 0 else 0)
    full_hi = k1 - 1
    if full_lo <= full_hi:
        anchor_base = tuple(base[i] + full_lo * w[i] for i in range(3))
        # 先确认区间内第一轮本身全程可行（O(图深)，与轮数无关）。
        found = _earliest_violation(child, 0, m, anchor_base, region, positions)
        if found is not None:
            return found[0] + full_lo * m, found[1], found[2], found[3]
        k_min, k_max = None, None
        for i in range(3):
            lo, hi = region[i]
            if w[i] == 0:
                continue
            base_lo = base[i] + child.rel_min[i]
            base_hi = base[i] + child.rel_max[i]
            if w[i] > 0:
                lower = _ceil(Fraction(lo - base_lo, w[i]))
                upper = _floor(Fraction(hi - base_hi, w[i]))
            else:
                lower = _ceil(Fraction(hi - base_hi, w[i]))
                upper = _floor(Fraction(lo - base_lo, w[i]))
            k_min = lower if k_min is None else max(k_min, lower)
            k_max = upper if k_max is None else min(k_max, upper)
        # full_lo 已验证可行，故 k_min <= full_lo <= k_max。
        if k_max is not None and k_max < full_hi:
            bad_k = k_max + 1
            bad_base = tuple(base[i] + bad_k * w[i] for i in range(3))
            found = _earliest_violation(child, 0, m, bad_base, region, positions)
            assert found is not None, 'feasibility boundary must contain a violation'
            return found[0] + bad_k * m, found[1], found[2], found[3]
    # 末轮：局部切片（u1 可为 0，即只含该轮起点状态）。
    child_base = tuple(base[i] + k1 * w[i] for i in range(3))
    found = _earliest_violation(child, 0, u1, child_base, region, positions)
    if found is None:
        return None
    return found[0] + k1 * m, found[1], found[2], found[3]


class _HandoffBuilder:
    """把获准窗口的步区间编译成内容寻址、跨窗共享的压缩定义表。"""

    def __init__(self, records, positions):
        self._records = records
        self._positions = positions
        self._segments = {}
        self._terms = {}
        self._used_instances = set()

    def _intern(self, term):
        payload = json.dumps(term, sort_keys=True, separators=(',', ':')).encode('utf-8')
        digest = 'h_' + hashlib.sha256(payload).hexdigest()[:16]
        if digest in self._terms and self._terms[digest] != term:
            raise ValueError('handoff definition hash collision')
        self._terms[digest] = term
        return digest

    def emit(self, node, a, b):
        key = (node.id, a, b)
        cached = self._segments.get(key)
        if cached is not None:
            return cached
        if node.kind == _STEP:
            if (a, b) != (0, 1):
                raise ValueError('step segment must cover the single hop')
            self._used_instances.add(node.instance)
            result = self._intern(('hop', node.instance[0], list(node.instance[1])))
        elif node.kind == _SEQUENCE:
            cut = node.left.steps
            parts = []
            if a < cut:
                parts.append(self.emit(node.left, a, min(b, cut)))
            if b > cut:
                parts.append(self.emit(node.right, max(0, a - cut), b - cut))
            if len(parts) == 1:
                result = parts[0]
            else:
                result = self._intern(('sequence', parts[0], parts[1]))
        else:
            m = node.child.steps
            k0, u0 = divmod(a, m)
            k1, u1 = divmod(b, m)
            if k0 == k1:
                # 单轮内部的切片（含整轮）直接由子程序定义承担。
                result = self.emit(node.child, u0, u1)
                self._segments[key] = result
                return result
            parts = []
            if u0 > 0:
                parts.append(self.emit(node.child, u0, m))
            full_lo = k0 + (1 if u0 > 0 else 0)
            full_hi = k1 - 1
            if full_lo <= full_hi:
                whole = self._whole(node.child)
                count = full_hi - full_lo + 1
                parts.append(whole if count == 1 else self._intern(('repeat', whole, count)))
            if u1 > 0:
                parts.append(self.emit(node.child, 0, u1))
            if not parts:
                raise ValueError('empty repeat segment')
            result = parts[0]
            for extra in parts[1:]:
                result = self._intern(('sequence', result, extra))
        self._segments[key] = result
        return result

    def _whole(self, node):
        return self.emit(node, 0, node.steps)

    def definitions(self):
        out = []
        for digest, term in sorted(self._terms.items()):
            if term[0] == 'hop':
                out.append({'id': digest, 'kind': 'hop', 'link': term[1], 'cell': term[2]})
            elif term[0] == 'sequence':
                out.append({'id': digest, 'kind': 'sequence', 'left': term[1], 'right': term[2]})
            else:
                out.append({'id': digest, 'kind': 'repeat', 'child': term[1], 'count': term[2]})
        return out

    def link_table(self):
        out = []
        for lid, cell in sorted(self._used_instances):
            record = self._records[(lid, cell)]
            from_site = (tuple(record['from']['cell']), tuple(Fraction(v) for v in record['from']['point']))
            to_site = (tuple(record['to']['cell']), tuple(Fraction(v) for v in record['to']['point']))
            out.append({'link': lid, 'cell': list(cell),
                        'from': {'cell': record['from']['cell'], 'point': record['from']['point'],
                                 'position': [str(v) for v in self._positions[from_site]]},
                        'to': {'cell': record['to']['cell'], 'point': record['to']['point'],
                               'position': [str(v) for v in self._positions[to_site]]},
                        'winding': list(record['winding'])})
        return out


def _site_ref(site, positions):
    cell, point = site
    return {'cell': list(cell), 'point': [str(v) for v in point],
            'position': [str(v) for v in positions[site]]}


def _parse_region(raw):
    pairs = array(raw, 3, 3)
    region = []
    for pair in pairs:
        bounds = array(pair, 2, 2)
        lo, hi = _rational(bounds[0]), _rational(bounds[1])
        if lo > hi:
            raise ValueError('region interval is empty')
        region.append((lo, hi))
    return tuple(region)


def handoff_windows(package, trajectory, windows):
    """核对时间窗测量并判定交接；输入契约与输出语义见 docs/contract.md。

    失败输入抛 ValueError；测量冲突与越界是逐窗业务结果，不影响其他窗口。
    不修改调用方传入的任何对象，返回值仅含 JSON 可序列化类型。
    """
    final = assemble_supercell_shards([package], finalize=True)
    positions = {}
    for record in final['points']:
        site = (tuple(record['cell']), tuple(Fraction(v) for v in record['point']))
        positions[site] = tuple(Fraction(v) for v in record['position'])
    records = {}
    for record in final['links']:
        records[(record['link'], tuple(record['cell']))] = record
    root = _parse_trajectory(trajectory, records, positions)
    total = root.steps

    raw_windows = array(windows, 1, MAX_WINDOWS)
    parsed_windows = []
    for entry in raw_windows:
        fields(entry, ('id', 'l', 'r', 'left_image', 'right_image', 'region'))
        wid = ident(entry['id'])
        left_step = _plain_int(entry['l'], 0, total)
        right_step = _plain_int(entry['r'], 0, total)
        if not left_step < right_step:
            raise ValueError('window requires 0 <= l < r <= total steps')
        left_image = _image_triple(entry['left_image'])
        right_image = _image_triple(entry['right_image'])
        region = _parse_region(entry['region'])
        parsed_windows.append({'id': wid, 'l': left_step, 'r': right_step,
                               'left_image': left_image, 'right_image': right_image,
                               'region': region})
    unique([w['id'] for w in parsed_windows])

    results = []
    builder = _HandoffBuilder(records, positions)
    routes = []
    for window in parsed_windows:
        l, r = window['l'], window['r']
        net = _segment_winding(root, l, r)
        derived = tuple(window['left_image'][i] + net[i] for i in range(3))
        if derived != window['right_image']:
            results.append({'id': window['id'], 'status': 'measurement_conflict',
                            'l': l, 'r': r,
                            'measured_left': list(window['left_image']),
                            'measured_right': list(window['right_image']),
                            'derived_right': list(derived)})
            continue
        root_base = tuple(window['left_image'][i] - _prefix_winding(root, l)[i] for i in range(3))
        found = _earliest_violation(root, l, r, root_base, window['region'], positions)
        if found is not None:
            step, site, image, coordinate = found
            results.append({'id': window['id'], 'status': 'out_of_region', 'l': l, 'r': r,
                            'violation': {'step': step,
                                          'site': {'cell': list(site[0]),
                                                   'point': [str(v) for v in site[1]],
                                                   'position': [str(v) for v in positions[site]]},
                                          'image': list(image),
                                          'coordinate': [str(v) for v in coordinate]}})
            continue
        program = builder.emit(root, l, r)
        start_site = _state_at(root, l)
        end_site = _state_at(root, r)
        end_image = tuple(window['left_image'][i] + net[i] for i in range(3))
        route = {'id': window['id'], 'l': l, 'r': r, 'hops': r - l,
                 'program': program,
                 'start': _site_ref(start_site, positions),
                 'image': list(window['left_image']),
                 'end_site': _site_ref(end_site, positions),
                 'end_image': list(end_image)}
        routes.append(route)
        results.append({'id': window['id'], 'status': 'handoff', 'l': l, 'r': r, 'route': route})

    output = {'format': FORMAT, 'version': final['version'],
              'basis': final['basis'], 'total_steps': total,
              'windows': sorted(results, key=lambda item: item['id'])}
    if routes:
        output['handoff'] = {'definitions': builder.definitions(),
                             'links': builder.link_table(),
                             'routes': sorted(routes, key=lambda item: item['id'])}
    return output


class RebuiltHandoff:
    """仅凭交接文档重建获准路线；不持有原包、原轨迹或任何进程缓存。"""

    def __init__(self, document):
        if not isinstance(document, dict) or document.get('format') != FORMAT:
            raise ValueError('not an orbit-handoff/1 document')
        handoff = document.get('handoff')
        if not isinstance(handoff, dict):
            raise ValueError('document carries no handoff data')
        fields(handoff, ('definitions', 'links', 'routes'))
        self._defs = {}
        for raw in array(handoff['definitions'], 1, 10 ** 9):
            kind = raw.get('kind')
            if kind == 'hop':
                fields(raw, ('id', 'kind', 'link', 'cell'))
                term = ('hop', ident(raw['link']), tuple(_cell_triple(raw['cell'])))
            elif kind == 'sequence':
                fields(raw, ('id', 'kind', 'left', 'right'))
                term = ('sequence', ident(raw['left']), ident(raw['right']))
            elif kind == 'repeat':
                fields(raw, ('id', 'kind', 'child', 'count'))
                term = ('repeat', ident(raw['child']), _plain_int(raw['count'], 1, MAX_STEPS))
            else:
                raise ValueError('unknown handoff definition kind')
            did = ident(raw['id'])
            if did in self._defs:
                raise ValueError('duplicate handoff definition')
            self._defs[did] = term
        self._links = {}
        for raw in array(handoff['links'], 1, 10 ** 9):
            fields(raw, ('link', 'cell', 'from', 'to', 'winding'))
            cell = tuple(_cell_triple(raw['cell']))
            lid = ident(raw['link'])
            for side in ('from', 'to'):
                fields(raw[side], ('cell', 'point', 'position'))
            winding = tuple(_plain_int(v) for v in array(raw['winding'], 3, 3))
            iid = (lid, cell)
            if iid in self._links:
                raise ValueError('duplicate handoff link record')
            self._links[iid] = raw
        self._routes = {}
        for raw in array(handoff['routes'], 1, MAX_WINDOWS):
            fields(raw, ('id', 'l', 'r', 'hops', 'program', 'start', 'image', 'end_site', 'end_image'))
            rid = ident(raw['id'])
            if rid in self._routes:
                raise ValueError('duplicate handoff route')
            if raw['program'] not in self._defs:
                raise ValueError('route references an unknown definition')
            l_step = _plain_int(raw['l'], 0, MAX_STEPS)
            r_step = _plain_int(raw['r'], 0, MAX_STEPS)
            if l_step >= r_step:
                raise ValueError('route step interval is empty')
            self._routes[rid] = raw
        self._validate_refs()
        for raw in self._routes.values():
            if _plain_int(raw['hops'], 1, MAX_STEPS) != self._steps_of(raw['program']):
                raise ValueError('route hop count does not match its program')
            for key in ('image', 'end_image'):
                if len(_image_triple(raw[key])) != 3:
                    raise ValueError('route image must be a triple')

    def _validate_refs(self):
        reachable = set()
        stack = [route['program'] for route in self._routes.values()]
        while stack:
            did = stack.pop()
            if did in reachable:
                continue
            if did not in self._defs:
                raise ValueError('definition references an unknown id')
            reachable.add(did)
            term = self._defs[did]
            if term[0] == 'sequence':
                stack.extend([term[1], term[2]])
            elif term[0] == 'repeat':
                stack.append(term[1])
            else:
                if (term[1], term[2]) not in self._links:
                    raise ValueError('hop references an unknown link record')
        if reachable != set(self._defs):
            raise ValueError('handoff table contains an unreachable definition')

    def route_ids(self):
        return sorted(self._routes)

    def _steps_of(self, did):
        term = self._defs[did]
        if term[0] == 'hop':
            return 1
        if term[0] == 'sequence':
            return self._steps_of(term[1]) + self._steps_of(term[2])
        return term[2] * self._steps_of(term[1])

    def _state_at(self, did, t):
        """路线局部第 t 步的 (位点引用, 相对像编号)。"""
        term = self._defs[did]
        if term[0] == 'hop':
            record = self._links[(term[1], term[2])]
            if t == 0:
                return record['from'], (0, 0, 0)
            return record['to'], tuple(record['winding'])
        if term[0] == 'sequence':
            cut = self._steps_of(term[1])
            if t <= cut:
                return self._state_at(term[1], t)
            site, rel = self._state_at(term[2], t - cut)
            left_w = self._net_winding(term[1])
            return site, tuple(left_w[i] + rel[i] for i in range(3))
        child_steps = self._steps_of(term[1])
        q, u = divmod(t, child_steps)
        site, rel = self._state_at(term[1], u)
        child_w = self._net_winding(term[1])
        return site, tuple(q * child_w[i] + rel[i] for i in range(3))

    def _net_winding(self, did):
        term = self._defs[did]
        if term[0] == 'hop':
            return tuple(self._links[(term[1], term[2])]['winding'])
        if term[0] == 'sequence':
            wl = self._net_winding(term[1])
            wr = self._net_winding(term[2])
            return tuple(wl[i] + wr[i] for i in range(3))
        child_w = self._net_winding(term[1])
        return tuple(term[2] * w for w in child_w)

    def state_at(self, window_id, step):
        """路线第 step 步（0..hops）的位点与绝对像编号，O(定义树深)。"""
        route = self._routes.get(window_id)
        if route is None:
            raise ValueError(f'unknown route: {window_id}')
        _plain_int(step, 0, route['hops'])
        site, rel = self._state_at(route['program'], step)
        base = tuple(route['image'])
        return {'site': {'cell': list(site['cell']), 'point': list(site['point']),
                         'position': list(site['position'])},
                'image': [base[i] + rel[i] for i in range(3)]}

    def end_state(self, window_id):
        route = self._routes[window_id]
        return {'site': route['end_site'], 'image': list(route['end_image'])}

    def iter_hops(self, window_id):
        """按次序惰性产出路线的每一跳动（重复以计数器展开，不预物化）。"""
        route = self._routes.get(window_id)
        if route is None:
            raise ValueError(f'unknown route: {window_id}')
        image = list(route['image'])
        index = 0
        stack = [(route['program'], 0)]
        while stack:
            did, cursor = stack.pop()
            term = self._defs[did]
            if term[0] == 'hop':
                record = self._links[(term[1], term[2])]
                index += 1
                landed = tuple(image[i] + record['winding'][i] for i in range(3))
                yield {'index': index, 'link': term[1], 'cell': list(term[2]),
                       'winding': list(record['winding']),
                       'from': {'cell': record['from']['cell'], 'point': record['from']['point'],
                                'position': record['from']['position']},
                       'to': {'cell': record['to']['cell'], 'point': record['to']['point'],
                              'position': record['to']['position']},
                       'from_image': list(image), 'to_image': list(landed)}
                image = list(landed)
            elif term[0] == 'sequence':
                stack.append((term[2], 0))
                stack.append((term[1], 0))
            else:
                remaining = term[2] - cursor
                if remaining <= 0:
                    continue
                if remaining > 1:
                    stack.append((did, cursor + 1))
                stack.append((term[1], 0))
        if index != route['hops']:
            raise ValueError('rebuilt route hop count does not match handoff record')
        if image != route['end_image']:
            raise ValueError('rebuilt end image does not match handoff record')


def rebuild_handoff(document):
    """校验交接文档并返回可离线重建路线的 RebuiltHandoff。"""
    return RebuiltHandoff(document)
