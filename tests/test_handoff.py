import copy
import json
import unittest
from fractions import Fraction

from site_orbits import (assemble_supercell_shards, build_supercell_shard,
                         handoff_windows, rebuild_handoff)

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'translation': [[0, 1], [0, 1], [0, 1]]}
UNIT = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
A_POINT = ['1/4', '0', '0']
B_POINT = ['3/4', '0', '0']


def ab_data():
    return {'operations': [IDENTITY],
            'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
                      {'id': 'b', 'species': 'Y', 'point': [[3, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]}]}


def ab_links():
    return [
        {'id': 'E', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]},
        {'id': 'H', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [0, 0, 0]},
        {'id': 'F', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [-1, 0, 0]},
    ]


def package(data=ab_data(), basis=UNIT, links=ab_links(), cells=((0, 0, 0),)):
    shards = [build_supercell_shard(data, basis, links, [list(cell)]) for cell in cells]
    return assemble_supercell_shards(shards, finalize=True)


def step(nid, lid, cell=(0, 0, 0)):
    return {'id': nid, 'kind': 'step', 'link': lid, 'cell': list(cell)}


def sequence(nid, left, right):
    return {'id': nid, 'kind': 'sequence', 'left': left, 'right': right}


def repeat(nid, child, count):
    return {'id': nid, 'kind': 'repeat', 'child': child, 'count': count}


def cycle_m_trajectory(count):
    """M = E 后接 F（净绕行 (+1,0,0)），根为 M 的 count 次重复。"""
    nodes = [step('e', 'E'), step('f', 'F'), sequence('M', 'e', 'f')]
    if count == 1:
        return {'root': 'M', 'nodes': nodes}
    nodes.append(repeat('R', 'M', count))
    return {'root': 'R', 'nodes': nodes}


def window(wid, l, r, left_image, right_image, region):
    return {'id': wid, 'l': l, 'r': r,
            'left_image': list(left_image), 'right_image': list(right_image),
            'region': [list(pair) for pair in region]}


REGION_NARROW = [['-1', '2'], ['0', '0'], ['0', '0']]


def naive_states(pkg, raw_trajectory, start_image):
    """从绝对起始像出发的朴素逐步状态序列（位点身份, 像编号），长度为步数+1。"""
    defs = {n['id']: n for n in raw_trajectory['nodes']}
    link_records = {(r['link'], tuple(r['cell'])): r for r in pkg['links']}
    hops = []

    def site_of(endpoint):
        return tuple(endpoint['cell']), tuple(Fraction(v) for v in endpoint['point'])

    def walk(nid):
        node = defs[nid]
        if node['kind'] == 'step':
            rec = link_records[(node['link'], tuple(node['cell']))]
            hops.append((node['link'], tuple(node['cell']), tuple(rec['winding']),
                         site_of(rec['from']), site_of(rec['to'])))
        elif node['kind'] == 'sequence':
            walk(node['left'])
            walk(node['right'])
        else:
            for _ in range(node['count']):
                walk(node['child'])

    walk(raw_trajectory['root'])
    image = tuple(start_image)
    site = hops[0][3]
    states = [(site, image)]
    for lid, cell, winding, from_site, to_site in hops:
        assert site == from_site
        image = tuple(image[i] + winding[i] for i in range(3))
        site = to_site
        states.append((site, image))
    return hops, states


def position_of(pkg, site):
    cell, point = site
    for record in pkg['points']:
        if tuple(record['cell']) == cell and tuple(Fraction(v) for v in record['point']) == point:
            return tuple(Fraction(v) for v in record['position'])
    raise KeyError(site)


def naive_violation(pkg, states, l, r, region):
    for t in range(l, r + 1):
        site, image = states[t]
        position = position_of(pkg, site)
        coordinate = tuple(position[i] + image[i] for i in range(3))
        if not all(region[i][0] <= coordinate[i] <= region[i][1] for i in range(3)):
            return t, site, image, coordinate
    return None


class AcceptanceOneTests(unittest.TestCase):
    """N=2 的循环段 M=E·F，三个同批窗口：越界、可交接、测量冲突。"""

    def setUp(self):
        self.pkg = package()
        self.trajectory = cycle_m_trajectory(2)

    def test_batch_outcome_for_three_windows(self):
        windows = [
            window('w-out', 0, 2, (0, 0, 0), (1, 0, 0), REGION_NARROW),
            window('w-ok', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW),
            window('w-bad', 0, 2, (-1, 0, 0), (1, 0, 0), REGION_NARROW),
        ]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        self.assertEqual(result['format'], 'orbit-handoff/1')
        self.assertEqual(result['version'], self.pkg['version'])
        self.assertEqual(result['basis'], UNIT)
        self.assertEqual(result['total_steps'], 4)
        by_id = {w['id']: w for w in result['windows']}
        self.assertEqual([w['id'] for w in result['windows']], ['w-bad', 'w-ok', 'w-out'])

        out = by_id['w-out']
        self.assertEqual(out['status'], 'out_of_region')
        self.assertEqual(out['violation']['step'], 1)
        self.assertEqual(out['violation']['site']['cell'], [0, 0, 0])
        self.assertEqual(out['violation']['site']['point'], B_POINT)
        self.assertEqual(out['violation']['image'], [2, 0, 0])
        self.assertEqual(out['violation']['coordinate'], ['11/4', '0', '0'])
        self.assertNotIn('route', out)

        bad = by_id['w-bad']
        self.assertEqual(bad['status'], 'measurement_conflict')
        self.assertEqual(bad['measured_left'], [-1, 0, 0])
        self.assertEqual(bad['measured_right'], [1, 0, 0])
        self.assertEqual(bad['derived_right'], [0, 0, 0])
        self.assertNotIn('violation', bad)
        self.assertNotIn('route', bad)

        ok = by_id['w-ok']
        self.assertEqual(ok['status'], 'handoff')
        route = ok['route']
        self.assertEqual(route['l'], 0)
        self.assertEqual(route['r'], 2)
        self.assertEqual(route['hops'], 2)
        self.assertEqual(route['image'], [-1, 0, 0])
        self.assertEqual(route['start']['point'], A_POINT)
        self.assertEqual(route['end_site']['point'], A_POINT)
        self.assertEqual(route['end_image'], [0, 0, 0])

    def test_handoff_table_excludes_unhandled_and_second_round(self):
        windows = [
            window('w-out', 0, 2, (0, 0, 0), (1, 0, 0), REGION_NARROW),
            window('w-ok', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW),
            window('w-bad', 0, 2, (-1, 0, 0), (1, 0, 0), REGION_NARROW),
        ]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        handoff = result['handoff']
        # 只用 E、F 两个实例；同端点平行边 H 不得代用，失败窗口不占用定义
        self.assertEqual(sorted((l['link'], tuple(l['cell'])) for l in handoff['links']),
                         [('E', (0, 0, 0)), ('F', (0, 0, 0))])
        kinds = sorted(d['kind'] for d in handoff['definitions'])
        self.assertEqual(kinds, ['hop', 'hop', 'sequence'])
        hop_links = sorted(d['link'] for d in handoff['definitions'] if d['kind'] == 'hop')
        self.assertEqual(hop_links, ['E', 'F'])
        # 路线只含一轮 M，不附第二轮让接收方跳过
        self.assertNotIn('repeat', kinds)
        # 链接表携带两端包内规范坐标
        e_record = next(l for l in handoff['links'] if l['link'] == 'E')
        self.assertEqual(e_record['from']['position'], A_POINT)
        self.assertEqual(e_record['to']['position'], B_POINT)
        self.assertEqual(e_record['winding'], [2, 0, 0])
        # 唯一起作用的路线
        self.assertEqual([r['id'] for r in handoff['routes']], ['w-ok'])

    def test_rebuild_offline_matches_calibrated_route(self):
        windows = [window('w-ok', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW)]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        rebuilt = rebuild_handoff(result)
        self.assertEqual(rebuilt.route_ids(), ['w-ok'])
        hops = list(rebuilt.iter_hops('w-ok'))
        self.assertEqual([(h['link'], h['from_image'], h['to_image']) for h in hops],
                         [('E', [-1, 0, 0], [1, 0, 0]), ('F', [1, 0, 0], [0, 0, 0])])
        self.assertEqual(rebuilt.state_at('w-ok', 0)['image'], [-1, 0, 0])
        self.assertEqual(rebuilt.state_at('w-ok', 1)['site']['point'], B_POINT)
        self.assertEqual(rebuilt.state_at('w-ok', 2)['image'], [0, 0, 0])
        self.assertEqual(rebuilt.end_state('w-ok'),
                         {'site': result['handoff']['routes'][0]['end_site'], 'image': [0, 0, 0]})

    def test_conflict_takes_precedence_over_region(self):
        windows = [window('w', 0, 2, (0, 0, 0), (9, 0, 0), [['0', '0'], ['0', '0'], ['0', '0']])]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        self.assertEqual(result['windows'][0]['status'], 'measurement_conflict')
        self.assertNotIn('handoff', result)

    def test_no_handoff_section_when_all_fail(self):
        windows = [window('w-bad', 0, 2, (-1, 0, 0), (1, 0, 0), REGION_NARROW)]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        self.assertNotIn('handoff', result)

    def test_parallel_edge_never_substituted(self):
        # 只含 E 的轨迹（H 同端点但零 winding），交接路线必须引用 E
        trajectory = {'root': 'e', 'nodes': [step('e', 'E')]}
        region = [['-2', '3'], ['0', '0'], ['0', '0']]
        result = handoff_windows(self.pkg, trajectory,
                                 [window('w', 0, 1, (0, 0, 0), (2, 0, 0), region)])
        self.assertEqual(result['handoff']['links'][0]['link'], 'E')
        hops = list(rebuild_handoff(result).iter_hops('w'))
        self.assertEqual([h['link'] for h in hops], ['E'])


class AcceptanceTwoTests(unittest.TestCase):
    """N=10¹²：越界定位不逐轮扫描；中段交接的定义表不随轮数增长。"""

    def setUp(self):
        self.n = 10 ** 12
        self.k = 500_000_000_000
        self.pkg = package()
        self.trajectory = cycle_m_trajectory(self.n)

    def test_earliest_violation_at_exact_step(self):
        win = window('big', 0, 2 * self.n, (0, 0, 0), (self.n, 0, 0),
                     [['0', str(self.k + 2)], ['0', '0'], ['0', '0']])
        result = handoff_windows(self.pkg, self.trajectory, [win])
        status = result['windows'][0]
        self.assertEqual(status['status'], 'out_of_region')
        violation = status['violation']
        self.assertEqual(violation['step'], 2 * self.k + 1)
        self.assertEqual(violation['site']['point'], B_POINT)
        self.assertEqual(violation['image'], [self.k + 2, 0, 0])
        self.assertEqual(Fraction(violation['coordinate'][0]), self.k + Fraction(11, 4))
        self.assertEqual(violation['coordinate'][0], f'{4 * self.k + 11}/4')

    def test_middle_window_handoff_structure(self):
        win = window('mid', 1, 2 * self.n - 1, (2, 0, 0), (self.n + 1, 0, 0),
                     [['0', str(self.n + 2)], ['0', '0'], ['0', '0']])
        result = handoff_windows(self.pkg, self.trajectory, [win])
        route = result['windows'][0]['route']
        self.assertEqual(route['hops'], 2 * self.n - 2)
        self.assertEqual(route['start']['point'], B_POINT)
        self.assertEqual(route['image'], [2, 0, 0])
        self.assertEqual(route['end_site']['point'], B_POINT)
        self.assertEqual(route['end_image'], [self.n + 1, 0, 0])
        definitions = result['handoff']['definitions']
        # F、E 两个 hop，M 串接，一次 (N-2) 轮重复，外加两片串接：共 6 条定义
        self.assertEqual(len(definitions), 6)
        repeat_defs = [d for d in definitions if d['kind'] == 'repeat']
        self.assertEqual(len(repeat_defs), 1)
        self.assertEqual(repeat_defs[0]['count'], self.n - 2)
        hops = sorted(d['link'] for d in definitions if d['kind'] == 'hop')
        self.assertEqual(hops, ['E', 'F'])
        # 定义表不得包含原根（N 轮重复）或被裁掉的跳转以外的轨迹残余
        self.assertTrue(all(d.get('count', 0) != self.n for d in definitions))
        rebuilt = rebuild_handoff(result)
        first = next(rebuilt.iter_hops('mid'))
        self.assertEqual((first['link'], first['from_image'], first['to_image']),
                         ('F', [2, 0, 0], [1, 0, 0]))
        self.assertEqual(rebuilt.state_at('mid', 2 * self.n - 2)['image'], [self.n + 1, 0, 0])
        self.assertEqual(rebuilt.state_at('mid', 2 * self.n - 2)['site']['point'], B_POINT)
        # 末跳为 E
        last_state = rebuilt.state_at('mid', 2 * self.n - 3)
        self.assertEqual(last_state['site']['point'], A_POINT)

    def test_many_identical_windows_share_definitions(self):
        region = [['0', str(self.n + 2)], ['0', '0'], ['0', '0']]
        windows = [window(f'w{i:04d}', 1, 2 * self.n - 1, (2, 0, 0), (self.n + 1, 0, 0), region)
                   for i in range(200)]
        result = handoff_windows(self.pkg, self.trajectory, windows)
        self.assertEqual(len(result['handoff']['definitions']), 6)
        self.assertEqual(len(result['handoff']['links']), 2)
        programs = {r['program'] for r in result['handoff']['routes']}
        self.assertEqual(len(programs), 1)
        self.assertEqual(len(result['handoff']['routes']), 200)
        # 定义表与链接表规模与窗口份数、10¹² 轮数均无关
        blob = json.dumps(result['handoff']['definitions'] + result['handoff']['links'])
        self.assertNotIn(str(self.n - 1), blob)
        self.assertIn(str(self.n - 2), blob)

    def test_output_independent_of_input_ordering(self):
        win = window('mid', 1, 2 * self.n - 1, (2, 0, 0), (self.n + 1, 0, 0),
                     [['0', str(self.n + 2)], ['0', '0'], ['0', '0']])
        first = handoff_windows(self.pkg, self.trajectory, [win])
        shuffled_pkg = copy.deepcopy(self.pkg)
        shuffled_pkg['points'] = shuffled_pkg['points'][::-1]
        shuffled_pkg['links'] = shuffled_pkg['links'][::-1]
        shuffled_traj = copy.deepcopy(self.trajectory)
        shuffled_traj['nodes'] = shuffled_traj['nodes'][::-1]
        second = handoff_windows(shuffled_pkg, shuffled_traj, [win])
        self.assertEqual(json.loads(json.dumps(first)), json.loads(json.dumps(second)))


class SmallNCrossCheckTests(unittest.TestCase):
    """小 N 同形输入：朴素逐步序列与交接结果逐项比对。"""

    def _region(self, n):
        return [['0', str(n + 2)], ['0', '0'], ['0', '0']]

    def test_full_and_middle_windows_match_naive_simulation(self):
        pkg = package()
        for n in (2, 3, 5, 8):
            with self.subTest(n=n):
                trajectory = cycle_m_trajectory(n)
                hops_flat, states = naive_states(pkg, trajectory, (0, 0, 0))
                self.assertEqual(len(hops_flat), 2 * n)
                # 整窗交接
                win = window('full', 0, 2 * n, (0, 0, 0), (n, 0, 0), self._region(n))
                result = handoff_windows(pkg, trajectory, [win])
                self.assertEqual(result['windows'][0]['status'], 'handoff')
                rebuilt = list(rebuild_handoff(result).iter_hops('full'))
                naive_slice = hops_flat[:]
                self.assertEqual([(h['link'], h['cell'], h['winding']) for h in rebuilt],
                                 [(h[0], list(h[1]), list(h[2])) for h in naive_slice])
                for t, (site, image) in enumerate(states):
                    got = rebuild_handoff(result).state_at('full', t)
                    self.assertEqual(tuple(got['image']), image)
                    self.assertEqual(tuple(Fraction(v) for v in got['site']['point']), site[1])
                # 中段窗口 [1, 2n-1]，以 B 起止
                win_mid = window('mid', 1, 2 * n - 1, (2, 0, 0), (n + 1, 0, 0), self._region(n))
                result_mid = handoff_windows(pkg, trajectory, [win_mid])
                self.assertEqual(result_mid['windows'][0]['status'], 'handoff')
                rebuilt_mid = list(rebuild_handoff(result_mid).iter_hops('mid'))
                self.assertEqual([h['link'] for h in rebuilt_mid],
                                 [h[0] for h in hops_flat[1:2 * n - 1]])
                self.assertEqual(rebuilt_mid[0]['link'], 'F')
                self.assertEqual(rebuilt_mid[-1]['link'], 'E')

    def test_violation_step_matches_naive_simulation(self):
        pkg = package()
        for n in (2, 3, 7):
            with self.subTest(n=n):
                trajectory = cycle_m_trajectory(n)
                _, states = naive_states(pkg, trajectory, (0, 0, 0))
                cap = 3
                region = [['0', str(cap)], ['0', '0'], ['0', '0']]
                win = window('w', 0, 2 * n, (0, 0, 0), (n, 0, 0), region)
                result = handoff_windows(pkg, trajectory, [win])
                expected = naive_violation(pkg, states, 0, 2 * n,
                                           tuple((Fraction(a), Fraction(b)) for a, b in region))
                self.assertIsNotNone(expected)
                violation = result['windows'][0]['violation']
                self.assertEqual(violation['step'], expected[0])
                self.assertEqual(tuple(violation['image']), expected[2])
                self.assertEqual(tuple(Fraction(v) for v in violation['coordinate']), expected[3])

    def test_violation_under_negative_calibration_matches_naive(self):
        pkg = package()
        n = 4
        trajectory = cycle_m_trajectory(n)
        _, states = naive_states(pkg, trajectory, (-5, -1, 2))
        region_str = [['-6', '0'], ['-1', '0'], ['0', '3']]
        region = tuple((Fraction(a), Fraction(b)) for a, b in region_str)
        win = window('w', 0, 2 * n, (-5, -1, 2), (n - 5, -1, 2), region_str)
        result = handoff_windows(pkg, trajectory, [win])
        expected = naive_violation(pkg, states, 0, 2 * n,
                                   tuple((Fraction(a), Fraction(b)) for a, b in region))
        if expected is None:
            self.assertEqual(result['windows'][0]['status'], 'handoff')
        else:
            self.assertEqual(result['windows'][0]['violation']['step'], expected[0])
            self.assertEqual(tuple(result['windows'][0]['violation']['image']), expected[2])


class SliceAndBoundaryTests(unittest.TestCase):
    """单个重复段中间切片、跨共享结点边界、左端已越界。"""

    def setUp(self):
        self.pkg = package()
        self.trajectory = cycle_m_trajectory(3)
        self.wide = [['-10', '10'], ['0', '0'], ['0', '0']]

    def test_slice_inside_one_repeat_round(self):
        # [1,2]：B →(F) A，完全落在第 0 轮 M 内
        win = window('slice', 1, 2, (1, 0, 0), (0, 0, 0), self.wide)
        result = handoff_windows(self.pkg, self.trajectory, [win])
        self.assertEqual(result['windows'][0]['status'], 'handoff')
        defs = result['handoff']['definitions']
        self.assertEqual(len(defs), 1)
        self.assertEqual(defs[0]['kind'], 'hop')
        self.assertEqual(defs[0]['link'], 'F')
        hops = list(rebuild_handoff(result).iter_hops('slice'))
        self.assertEqual([(h['link'], h['from_image'], h['to_image']) for h in hops],
                         [('F', [1, 0, 0], [0, 0, 0])])

    def test_window_crossing_shared_node_boundary(self):
        # [1,3]：B →(F) A →(E) B，跨过第 0/1 轮边界（步 2）
        win = window('cross', 1, 3, (1, 0, 0), (2, 0, 0), self.wide)
        result = handoff_windows(self.pkg, self.trajectory, [win])
        self.assertEqual(result['windows'][0]['status'], 'handoff')
        hops = list(rebuild_handoff(result).iter_hops('cross'))
        self.assertEqual([h['link'] for h in hops], ['F', 'E'])
        self.assertEqual([h['from_image'] for h in hops], [[1, 0, 0], [0, 0, 0]])
        # 与朴素展开一致
        flat, _ = naive_states(self.pkg, self.trajectory, (0, 0, 0))
        self.assertEqual([h[0] for h in flat[1:3]], ['F', 'E'])

    def test_left_state_already_out_of_region(self):
        # 左端 B 像 (1,0,0)，x=7/4，区域 x∈[0,1]：第 l 步即越界，不交前缀
        region = [['0', '1'], ['0', '0'], ['0', '0']]
        win = window('leftbad', 1, 2, (1, 0, 0), (0, 0, 0), region)
        result = handoff_windows(self.pkg, self.trajectory, [win])
        status = result['windows'][0]
        self.assertEqual(status['status'], 'out_of_region')
        self.assertEqual(status['violation']['step'], 1)
        self.assertEqual(status['violation']['site']['point'], B_POINT)
        self.assertEqual(status['violation']['image'], [1, 0, 0])
        self.assertNotIn('handoff', result)

    def test_repeat_node_with_nonzero_net_winding_is_legal(self):
        # M 净绕行 (+1,0,0)，重复它本身合法（只要求包内位点身份闭合）
        result = handoff_windows(self.pkg, self.trajectory,
                                 [window('w', 0, 6, (0, 0, 0), (3, 0, 0), self.wide)])
        self.assertEqual(result['windows'][0]['status'], 'handoff')

    def test_self_loop_step_repeat_nonzero_winding(self):
        data = {'operations': [IDENTITY],
                'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                           'occupancy': [1, 1]}]}
        links = [{'id': 'L', 'from': [[1, 4], [0, 1], [0, 1]],
                  'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 1, 0]}]
        pkg = package(data, UNIT, links)
        trajectory = {'root': 'R', 'nodes': [step('s', 'L'), repeat('R', 's', 5)]}
        result = handoff_windows(pkg, trajectory,
                                 [window('w', 0, 5, (0, 0, 0), (10, 5, 0),
                                         [['0', '12'], ['0', '6'], ['0', '0']])])
        self.assertEqual(result['windows'][0]['status'], 'handoff')
        hops = list(rebuild_handoff(result).iter_hops('w'))
        self.assertEqual([h['to_image'] for h in hops],
                         [[2, 1, 0], [4, 2, 0], [6, 3, 0], [8, 4, 0], [10, 5, 0]])

    def test_zero_winding_axis_violation_inside_repeat(self):
        # 每跳只沿 x 漂移；y 轴恒为 0，区域把 y 限成 [1,1]：窗口第 0 状态即越界
        data = {'operations': [IDENTITY],
                'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                           'occupancy': [1, 1]}]}
        links = [{'id': 'L', 'from': [[1, 4], [0, 1], [0, 1]],
                  'to': [[1, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]}]
        pkg = package(data, UNIT, links)
        trajectory = {'root': 'R', 'nodes': [step('s', 'L'), repeat('R', 's', 1000)]}
        win = window('w', 100, 900, (0, 0, 0), (800, 0, 0),
                     [['0', '1000'], ['1', '1'], ['0', '0']])
        result = handoff_windows(pkg, trajectory, [win])
        self.assertEqual(result['windows'][0]['status'], 'out_of_region')
        self.assertEqual(result['windows'][0]['violation']['step'], 100)
        self.assertEqual(result['windows'][0]['violation']['coordinate'], ['1/4', '0', '0'])

    def test_violation_in_last_partial_round_matches_naive(self):
        n = 6
        pkg = package()
        trajectory = cycle_m_trajectory(n)
        left_image = (2, 0, 0)
        _, states_zero = naive_states(pkg, trajectory, (0, 0, 0))
        net = tuple(states_zero[2 * n][1][i] - states_zero[1][1][i] for i in range(3))
        right_image = tuple(left_image[i] + net[i] for i in range(3))
        root_base = tuple(left_image[i] - states_zero[1][1][i] for i in range(3))
        states = [(site, tuple(root_base[i] + rel[i] for i in range(3)))
                  for site, rel in states_zero]
        region_str = [['0', '3'], ['0', '0'], ['0', '0']]
        region = tuple((Fraction(a), Fraction(b)) for a, b in region_str)
        win = window('w', 1, 2 * n, left_image, right_image, region_str)
        result = handoff_windows(pkg, trajectory, [win])
        expected = naive_violation(pkg, states, 1, 2 * n, region)
        self.assertIsNotNone(expected)
        self.assertEqual(result['windows'][0]['status'], 'out_of_region')
        violation = result['windows'][0]['violation']
        self.assertEqual(violation['step'], expected[0])
        self.assertEqual(tuple(violation['image']), expected[2])
        self.assertEqual(tuple(Fraction(v) for v in violation['coordinate']), expected[3])


class InvalidInputTests(unittest.TestCase):
    def setUp(self):
        self.pkg = package()

    def test_invalid_sequence_closure_rejected(self):
        bad = {'root': 'S', 'nodes': [step('e1', 'E'), step('e2', 'E'),
                                      sequence('S', 'e1', 'e2')]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, bad, [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                                   [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_repeat_without_site_identity_closure_rejected(self):
        bad = {'root': 'R', 'nodes': [step('e', 'E'), repeat('R', 'e', 2)]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, bad, [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                                   [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_cycle_rejected(self):
        cyclic = {'root': 'A', 'nodes': [step('e', 'E'), step('f', 'F'),
                                         sequence('A', 'f', 'B'), sequence('B', 'e', 'A')]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, cyclic,
                            [window('w', 0, 4, (0, 0, 0), (1, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])
        self_repeat = {'root': 'R', 'nodes': [repeat('R', 'R', 2)]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, self_repeat,
                            [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_unreachable_and_unknown_root_rejected(self):
        orphan = {'root': 'e', 'nodes': [step('e', 'E'), step('f', 'F')]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, orphan,
                            [window('w', 0, 1, (0, 0, 0), (2, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, {'root': 'nope', 'nodes': [step('e', 'E')]},
                            [window('w', 0, 1, (0, 0, 0), (2, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_unknown_link_instance_rejected(self):
        bad = {'root': 'x', 'nodes': [step('x', 'E', cell=(1, 0, 0))]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, bad,
                            [window('w', 0, 1, (0, 0, 0), (2, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_incomplete_package_rejected(self):
        shard = build_supercell_shard(ab_data(), [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                                      ab_links(), [[0, 0, 0]])
        opened = assemble_supercell_shards([shard])
        with self.assertRaises(ValueError):
            handoff_windows(opened, cycle_m_trajectory(2),
                            [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])
        with self.assertRaises(ValueError):
            handoff_windows({'format': 'supercell-shard/1'}, cycle_m_trajectory(2),
                            [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])

    def test_limits_enforced(self):
        base_nodes = [step('e', 'E'), step('f', 'F'), sequence('M', 'e', 'f')]
        # 结点数上限：一条可达链上的 5001 个结点
        many = [step('s0', 'E')]
        # 串接要求位点闭合，改用自连边链
        self_loop_pkg0 = package(
            data={'operations': [IDENTITY],
                  'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                             'occupancy': [1, 1]}]},
            links=[{'id': 'L', 'from': [[1, 4], [0, 1], [0, 1]],
                    'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]}])
        many = [step('s0', 'L')]
        for i in range(1, 5001):
            many.append(sequence(f's{i}', f's{i - 1}', 's0'))
        with self.assertRaises(ValueError):
            handoff_windows(self_loop_pkg0, {'root': 's5000', 'nodes': many},
                            [window('w', 0, 5001, (0, 0, 0), (1, 0, 0),
                                    [['-9', '99999'], ['0', '0'], ['0', '0']])])
        # 深度上限：49 层串接（用自连边保证每层位点身份闭合）
        self_loop_pkg = package(
            data={'operations': [IDENTITY],
                  'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                             'occupancy': [1, 1]}]},
            links=[{'id': 'L', 'from': [[1, 4], [0, 1], [0, 1]],
                    'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]}])
        deep_nodes = [step('d0', 'L')]
        for i in range(1, 49):
            deep_nodes.append(sequence(f'd{i}', f'd{i - 1}', 'd0'))
        with self.assertRaises(ValueError):
            handoff_windows(self_loop_pkg, {'root': 'd48', 'nodes': deep_nodes},
                            [window('w', 0, 49, (0, 0, 0), (98, 0, 0),
                                    [['0', '200'], ['0', '0'], ['0', '0']])])
        # 48 层（48 个结点）仍合法
        handoff_windows(self_loop_pkg, {'root': 'd47', 'nodes': deep_nodes[:-1]},
                        [window('w', 0, 48, (0, 0, 0), (96, 0, 0),
                                [['0', '200'], ['0', '0'], ['0', '0']])])
        # 总步数上限：count 本身受限且乘积受限
        huge = {'root': 'R', 'nodes': base_nodes + [repeat('R', 'M', 10 ** 18 + 1)]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, huge,
                            [window('w', 0, 2, (0, 0, 0), (1, 0, 0),
                                    [['-9', '9'], ['0', '0'], ['0', '0']])])
        # 窗口数量上限与重复 ID
        many_windows = [window(f'w{i}', 0, 2, (0, 0, 0), (1, 0, 0),
                               [['-9', '9'], ['0', '0'], ['0', '0']]) for i in range(2001)]
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, cycle_m_trajectory(2), many_windows)
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, cycle_m_trajectory(2), [])

    def test_scalar_validation(self):
        good = window('w', 0, 2, (0, 0, 0), (1, 0, 0), [['-9', '9'], ['0', '0'], ['0', '0']])
        cases = [
            dict(good, l=True), dict(good, l=2, r=2),
            dict(good, r=99), dict(good, left_image=[True, 0, 0]),
            dict(good, right_image=[0, 1.5, 0]),
            dict(good, left_image=[10 ** 60 + 1, 0, 0]),
            dict(good, region=[['0.5', '9'], ['0', '0'], ['0', '0']]),
            dict(good, region=[['9', '0'], ['0', '0'], ['0', '0']]),
        ]
        for bad in cases:
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    handoff_windows(self.pkg, cycle_m_trajectory(2), [bad])
        bad_count = {'root': 'R', 'nodes': [step('e', 'E'), step('f', 'F'),
                                            sequence('M', 'e', 'f'),
                                            repeat('R', 'M', True)]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, bad_count, [good])
        bad_count2 = {'root': 'R', 'nodes': [step('e', 'E'), step('f', 'F'),
                                             sequence('M', 'e', 'f'),
                                             repeat('R', 'M', 1)]}
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, bad_count2, [good])
        with self.assertRaises(ValueError):
            handoff_windows(self.pkg, {'root': 'R', 'nodes': [
                step('e', 'E'), step('f', 'F'), sequence('M', 'e', 'f'),
                dict(repeat('R', 'M', 2), bogus=1)]}, [good])

    def test_no_input_mutation_and_json_round_trip(self):
        pkg = copy.deepcopy(self.pkg)
        trajectory = copy.deepcopy(cycle_m_trajectory(2))
        windows = [
            window('w-out', 0, 2, (0, 0, 0), (1, 0, 0), REGION_NARROW),
            window('w-ok', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW),
            window('w-bad', 0, 2, (-1, 0, 0), (1, 0, 0), REGION_NARROW),
        ]
        snapshot = copy.deepcopy((pkg, trajectory, windows))
        result = handoff_windows(pkg, trajectory, windows)
        self.assertEqual((pkg, trajectory, windows), snapshot)
        self.assertEqual(json.loads(json.dumps(result)), result)

    def test_window_order_does_not_change_output(self):
        maker = lambda: [
            window('a', 0, 2, (0, 0, 0), (1, 0, 0), REGION_NARROW),
            window('b', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW),
            window('c', 0, 2, (-1, 0, 0), (1, 0, 0), REGION_NARROW),
        ]
        first = handoff_windows(self.pkg, cycle_m_trajectory(2), maker())
        second = handoff_windows(self.pkg, cycle_m_trajectory(2), maker()[::-1])
        self.assertEqual(json.loads(json.dumps(first)), json.loads(json.dumps(second)))


class RebuildValidationTests(unittest.TestCase):
    def setUp(self):
        self.pkg = package()
        self.result = handoff_windows(
            self.pkg, cycle_m_trajectory(2),
            [window('w-ok', 0, 2, (-1, 0, 0), (0, 0, 0), REGION_NARROW)])

    def test_reject_non_document(self):
        with self.assertRaises(ValueError):
            rebuild_handoff({'format': 'other'})
        with self.assertRaises(ValueError):
            rebuild_handoff({'format': 'orbit-handoff/1'})
        with self.assertRaises(ValueError):
            rebuild_handoff('not-a-document')

    def test_reject_broken_references(self):
        doc = copy.deepcopy(self.result)
        doc['handoff']['definitions'][0]['id'] = 'h_deadbeefdeadbeef'
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)
        doc = copy.deepcopy(self.result)
        doc['handoff']['routes'][0]['program'] = 'h_unknown000000000'
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)
        doc = copy.deepcopy(self.result)
        doc['handoff']['links'][0]['link'] = 'ZZZ'
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)

    def test_reject_unreachable_and_duplicate_definitions(self):
        doc = copy.deepcopy(self.result)
        extra = dict(doc['handoff']['definitions'][0])
        extra['id'] = 'h_extra00000000000'
        doc['handoff']['definitions'].append(extra)
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)
        doc = copy.deepcopy(self.result)
        doc['handoff']['definitions'].append(dict(doc['handoff']['definitions'][0]))
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)

    def test_state_index_bounds(self):
        rebuilt = rebuild_handoff(self.result)
        with self.assertRaises(ValueError):
            rebuilt.state_at('w-ok', 3)
        with self.assertRaises(ValueError):
            rebuilt.state_at('missing', 0)

    def test_tampered_route_fields_rejected(self):
        doc = copy.deepcopy(self.result)
        doc['handoff']['routes'][0]['hops'] = 99
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)
        doc = copy.deepcopy(self.result)
        doc['handoff']['routes'][0]['r'] = doc['handoff']['routes'][0]['l']
        with self.assertRaises(ValueError):
            rebuild_handoff(doc)
        doc = copy.deepcopy(self.result)
        doc['handoff']['routes'][0]['end_image'] = [7, 7, 7]
        # 末状态字段与程序推演不符：iter_hops 耗尽时检测
        rebuilt = rebuild_handoff(doc)
        with self.assertRaises(ValueError):
            list(rebuilt.iter_hops('w-ok'))
        doc = copy.deepcopy(self.result)
        doc['handoff']['links'][0]['winding'] = [9, 9, 9]
        with self.assertRaises(ValueError):
            list(rebuild_handoff(doc).iter_hops('w-ok'))


class MultiCellPackageTests(unittest.TestCase):
    """非单位基底、双胞片包：位点身份携带 cell，切片与重建仍然闭合。"""

    def test_basis_det_two_handoff(self):
        basis = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
        pkg = package(basis=basis, cells=((0, 0, 0), (1, 0, 0)))
        links = [
            {'id': 'E', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]],
             'shift': [2, 0, 0]},
            {'id': 'F', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]],
             'shift': [-1, 0, 0]},
        ]
        records = {(r['link'], tuple(r['cell'])): r for r in pkg['links']}
        # 基底 diag(2,1,1) 下：F 从 cell0 跨到 cell1，再由 cell1 的 E 回到 cell0 同身份
        f_rec = records[('F', (0, 0, 0))]
        self.assertEqual(tuple(f_rec['to']['cell']), (1, 0, 0))
        e = records[('E', (1, 0, 0))]
        self.assertEqual(tuple(e['to']['cell']), (1, 0, 0))
        trajectory = {'root': 'S', 'nodes': [
            step('f', 'F', (0, 0, 0)),
            step('e', 'E', (1, 0, 0)),
            sequence('S', 'f', 'e'),
        ]}
        region = [['-12', '12'], ['0', '0'], ['0', '0']]
        net = tuple(f_rec['winding'][i] + e['winding'][i] for i in range(3))
        win = window('w', 0, 2, (0, 0, 0), net, region)
        result = handoff_windows(pkg, trajectory, [win])
        self.assertEqual(result['windows'][0]['status'], 'handoff')
        rebuilt = rebuild_handoff(result)
        hops = list(rebuilt.iter_hops('w'))
        self.assertEqual([h['cell'] for h in hops], [[0, 0, 0], [1, 0, 0]])
        self.assertEqual(hops[0]['to']['cell'], [1, 0, 0])
        self.assertEqual(hops[1]['from']['cell'], [1, 0, 0])
        self.assertEqual(hops[1]['to']['cell'], [1, 0, 0])
        # 中间状态身份携带 cell1
        mid = rebuilt.state_at('w', 1)
        self.assertEqual(mid['site']['cell'], [1, 0, 0])
        self.assertEqual(mid['site']['point'], A_POINT)
        self.assertEqual(mid['site']['position'], ['5/8', '0', '0'])


if __name__ == '__main__':
    unittest.main()
