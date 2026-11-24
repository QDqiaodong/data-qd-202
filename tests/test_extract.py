import copy
import json
import unittest
from fractions import Fraction

from site_orbits import assemble_supercell_shards, build_supercell_shard, extract_supercell_groups

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]}
UNIT = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
A_POINT = ['1/4', '0', '0']
B_POINT = ['3/4', '0', '0']


def site(sid, species, x, occupancy):
    return {'id': sid, 'species': species, 'point': [x, [0, 1], [0, 1]], 'occupancy': occupancy}


def ab_data():
    return {'operations': [IDENTITY],
            'sites': [site('a', 'X', [1, 4], [1, 1]), site('b', 'Y', [3, 4], [1, 1])]}


def self_links():
    return [
        {'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 1, 0]},
        {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [8, 1, 0]},
    ]


def package(data, basis, links, cells=((0, 0, 0),)):
    shards = [build_supercell_shard(data, basis, links, [list(cell)]) for cell in cells]
    return assemble_supercell_shards(shards, finalize=True)


def entry(eid, point, replica, cell=(0, 0, 0)):
    return {'id': eid, 'cell': list(cell), 'point': list(point), 'replica': list(replica)}


class AcceptanceOneTests(unittest.TestCase):
    """恒等操作单位基底：A 带平移 (2,1,0)/(8,1,0) 两条自连边，B 无连接，重复 (6,6,1)。"""

    def setUp(self):
        self.pkg = package(ab_data(), UNIT, self_links())

    def requests(self):
        return [entry('req-a0', A_POINT, [0, 0, 0]),
                entry('req-a1', A_POINT, [4, 2, 0]),
                entry('req-b0', B_POINT, [0, 0, 0])]

    def test_delivery_grouped_and_windings(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(),
                                          [entry('keep', A_POINT, [2, 2, 0])], quota=6)
        self.assertEqual(result['format'], 'supercell-extract/1')
        self.assertEqual(result['version'], self.pkg['version'])
        self.assertEqual(result['basis'], UNIT)
        self.assertEqual(result['repeats'], [6, 6, 1])
        first, second = result['groups']
        # 两个 A 请求同组：归到同一项，交付 6 点 12 边
        self.assertEqual(first['requests'], ['req-a0', 'req-a1'])
        self.assertEqual(first['status'], 'delivered')
        self.assertEqual(first['version'], self.pkg['version'])
        self.assertEqual(first['basis'], UNIT)
        self.assertEqual(first['repeats'], [6, 6, 1])
        self.assertEqual(sorted(tuple(s['replica']) for s in first['sites']),
                         [(0, 0, 0), (0, 3, 0), (2, 1, 0), (2, 4, 0), (4, 2, 0), (4, 5, 0)])
        self.assertEqual(len(first['links']), 12)
        # B 组剩余额度不足，不截部分组
        self.assertEqual(second['requests'], ['req-b0'])
        self.assertEqual(second['status'], 'quota_exceeded')
        self.assertNotIn('sites', second)
        # 从 A 的 (4,2,0) 出发两条边都落在 A 的 (0,3,0)，绕行量不同、不合并
        selected = [link for link in first['links'] if link['replica'] == [4, 2, 0]]
        self.assertEqual([(link['link'], link['to']['replica'], link['winding']) for link in selected],
                         [('L1', [0, 3, 0], [1, 0, 0]), ('L2', [0, 3, 0], [2, 0, 0])])
        # 规范分数坐标 = (包内 position + replica) / 重复次数
        by_replica = {tuple(s['replica']): s for s in first['sites']}
        self.assertEqual(by_replica[(0, 0, 0)]['position'], ['1/24', '0', '0'])
        self.assertEqual(by_replica[(4, 2, 0)]['position'], ['17/24', '1/3', '0'])
        self.assertEqual(by_replica[(2, 4, 0)]['position'], ['3/8', '2/3', '0'])
        for s in first['sites']:
            self.assertEqual(s['cell'], [0, 0, 0])
            self.assertEqual(s['point'], A_POINT)
            self.assertEqual(s['occupancy'], '1')
            self.assertEqual(s['vacancy'], '0')
            self.assertEqual(s['contributions'],
                             [{'site': 'a', 'species': 'X', 'occupancy': '1', 'operations': ['identity']}])

    def test_winding_relation_holds(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(),
                                          [entry('keep', A_POINT, [2, 2, 0])], quota=6)
        group = result['groups'][0]
        repeats = (6, 6, 1)
        positions = {(tuple(s['cell']), tuple(s['point']), tuple(s['replica'])):
                     tuple(Fraction(v) for v in s['position']) for s in group['sites']}
        pkg_points = {(tuple(p['cell']), tuple(p['point'])):
                      tuple(Fraction(v) for v in p['position']) for p in self.pkg['points']}
        records = {r['link']: r for r in self.pkg['links']}
        self.assertEqual(len(positions), 6)
        for link in group['links']:
            source = positions[(tuple(link['from']['cell']), tuple(link['from']['point']), tuple(link['from']['replica']))]
            target = positions[(tuple(link['to']['cell']), tuple(link['to']['point']), tuple(link['to']['replica']))]
            record = records[link['link']]
            displacement = tuple(pkg_points[(tuple(record['to']['cell']), tuple(record['to']['point']))][i]
                                 + record['winding'][i]
                                 - pkg_points[(tuple(record['from']['cell']), tuple(record['from']['point']))][i]
                                 for i in range(3))
            for i in range(3):
                self.assertEqual(target[i] + link['winding'][i] - source[i],
                                 displacement[i] / repeats[i])

    def test_reserve_inside_group_blocks_and_frees_quota(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(),
                                          [entry('keep', A_POINT, [2, 1, 0])], quota=6)
        first, second = result['groups']
        # 两个 A 请求仍同组，整组保留不输出，额度让给 B
        self.assertEqual(first['requests'], ['req-a0', 'req-a1'])
        self.assertEqual(first['status'], 'reserved')
        self.assertEqual(first['reserves'], ['keep'])
        self.assertNotIn('sites', first)
        self.assertNotIn('links', first)
        self.assertEqual(second['status'], 'delivered')
        self.assertEqual(len(second['sites']), 1)
        self.assertEqual(second['sites'][0]['position'], ['1/8', '0', '0'])

    def test_request_order_controls_quota(self):
        requests = [entry('req-b0', B_POINT, [0, 0, 0]),
                    entry('req-a0', A_POINT, [0, 0, 0]),
                    entry('req-a1', A_POINT, [4, 2, 0])]
        result = extract_supercell_groups(self.pkg, [6, 6, 1], requests,
                                          [entry('keep', A_POINT, [2, 2, 0])], quota=6)
        first, second = result['groups']
        self.assertEqual(first['requests'], ['req-b0'])
        self.assertEqual(first['status'], 'delivered')
        self.assertEqual(second['requests'], ['req-a0', 'req-a1'])
        self.assertEqual(second['status'], 'quota_exceeded')
        self.assertNotIn('sites', second)

    def test_marker_on_other_replica_does_not_block(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], [entry('req', B_POINT, [0, 0, 0])],
                                          [entry('keep', B_POINT, [1, 0, 0])], quota=1)
        self.assertEqual(result['groups'][0]['status'], 'delivered')

    def test_reserves_sorted_and_order_free(self):
        reserves = [entry('m2', A_POINT, [4, 2, 0]), entry('m1', A_POINT, [2, 1, 0])]
        forward = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(), reserves, quota=6)
        self.assertEqual(forward['groups'][0]['status'], 'reserved')
        self.assertEqual(forward['groups'][0]['reserves'], ['m1', 'm2'])
        backward = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(), reserves[::-1], quota=6)
        self.assertEqual(forward, backward)

    def test_package_record_order_irrelevant(self):
        shuffled = copy.deepcopy(self.pkg)
        shuffled['points'] = shuffled['points'][::-1]
        shuffled['links'] = shuffled['links'][::-1]
        reserves = [entry('keep', A_POINT, [2, 2, 0])]
        self.assertEqual(extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(), reserves, quota=6),
                         extract_supercell_groups(shuffled, [6, 6, 1], self.requests(), reserves, quota=6))

    def test_json_round_trip(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], self.requests(),
                                          [entry('keep', A_POINT, [2, 2, 0])], quota=6)
        self.assertEqual(json.loads(json.dumps(result)), result)

    def test_same_site_two_requests_single_group(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1],
                                          [entry('r1', B_POINT, [0, 0, 0]), entry('r2', B_POINT, [0, 0, 0])],
                                          quota=1)
        (group,) = result['groups']
        self.assertEqual(group['requests'], ['r1', 'r2'])
        self.assertEqual(group['status'], 'delivered')
        self.assertEqual(len(group['sites']), 1)

    def test_request_and_marker_may_share_id(self):
        result = extract_supercell_groups(self.pkg, [6, 6, 1], [entry('x', B_POINT, [0, 0, 0])],
                                          [entry('x', B_POINT, [1, 0, 0])], quota=1)
        self.assertEqual(result['groups'][0]['status'], 'delivered')


class AcceptanceTwoTests(unittest.TestCase):
    """仅 A→B 一条边、平移 (-1,0,0)，B 混合占位，重复 (4,1,1)，额度 2。"""

    def data(self):
        return {'operations': [IDENTITY],
                'sites': [site('a', 'X', [1, 4], [1, 1]),
                          site('b1', 'Y', [3, 4], [1, 4]),
                          site('b2', 'Z', [3, 4], [1, 2])]}

    def links(self):
        return [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [-1, 0, 0]}]

    def setUp(self):
        self.pkg = package(self.data(), UNIT, self.links())

    def test_reserved_group_rejected(self):
        result = extract_supercell_groups(self.pkg, [4, 1, 1], [entry('req', B_POINT, [3, 0, 0])],
                                          [entry('keep', A_POINT, [0, 0, 0])], quota=2)
        (group,) = result['groups']
        self.assertEqual(group['status'], 'reserved')
        self.assertEqual(group['reserves'], ['keep'])
        self.assertNotIn('sites', group)

    def test_delivery_keeps_direction_sources_and_vacancy(self):
        result = extract_supercell_groups(self.pkg, [4, 1, 1], [entry('req', B_POINT, [3, 0, 0])], quota=2)
        (group,) = result['groups']
        self.assertEqual(group['status'], 'delivered')
        self.assertEqual([(s['point'], s['replica'], s['position']) for s in group['sites']],
                         [(A_POINT, [0, 0, 0], ['1/16', '0', '0']),
                          (B_POINT, [3, 0, 0], ['15/16', '0', '0'])])
        mixed = group['sites'][1]
        self.assertEqual(mixed['occupancy'], '3/4')
        self.assertEqual(mixed['vacancy'], '1/4')
        self.assertEqual([(c['site'], c['species'], c['occupancy']) for c in mixed['contributions']],
                         [('b1', 'Y', '1/4'), ('b2', 'Z', '1/2')])
        (link,) = group['links']
        self.assertEqual(link['link'], 'L1')
        self.assertEqual(link['cell'], [0, 0, 0])
        self.assertEqual(link['replica'], [0, 0, 0])
        self.assertEqual(link['from'], {'cell': [0, 0, 0], 'point': A_POINT, 'replica': [0, 0, 0]})
        self.assertEqual(link['to'], {'cell': [0, 0, 0], 'point': B_POINT, 'replica': [3, 0, 0]})
        self.assertEqual(link['winding'], [-1, 0, 0])

    def test_open_package_rejected(self):
        shard = build_supercell_shard(self.data(), [[2, 0, 0], [0, 1, 0], [0, 0, 1]], self.links(), [[0, 0, 0]])
        self.assertEqual(shard['cells'], [[0, 0, 0]])
        with self.assertRaises(ValueError):
            extract_supercell_groups(shard, [4, 1, 1], [entry('req', B_POINT, [3, 0, 0])], quota=2)
        opened = assemble_supercell_shards([shard])
        with self.assertRaises(ValueError):
            extract_supercell_groups(opened, [4, 1, 1], [entry('req', B_POINT, [3, 0, 0])], quota=2)


class AcceptanceThreeTests(unittest.TestCase):
    """单条平移 (2,1,0) 自连边，重复 (10⁹,10⁹,1)，额度 3：不展开十亿像。"""

    def setUp(self):
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 1, 0]}]
        self.pkg = package(ab_data(), UNIT, links)
        self.requests = [entry('req-a0', A_POINT, [0, 0, 0]),
                         entry('req-a1', A_POINT, [0, 500000000, 0]),
                         entry('req-b0', B_POINT, [0, 0, 0])]

    def test_reserved_by_unrequested_image(self):
        result = extract_supercell_groups(self.pkg, [10 ** 9, 10 ** 9, 1], self.requests,
                                          [entry('keep', A_POINT, [4, 500000002, 0])], quota=3)
        first, second = result['groups']
        self.assertEqual(first['requests'], ['req-a0', 'req-a1'])
        self.assertEqual(first['status'], 'reserved')
        self.assertEqual(first['reserves'], ['keep'])
        self.assertEqual(second['status'], 'delivered')
        self.assertEqual(len(second['sites']), 1)
        self.assertEqual(second['sites'][0]['point'], B_POINT)

    def test_marker_outside_group_leaves_quota_shortage(self):
        result = extract_supercell_groups(self.pkg, [10 ** 9, 10 ** 9, 1], self.requests,
                                          [entry('keep', A_POINT, [0, 1, 0])], quota=3)
        first, second = result['groups']
        self.assertEqual(first['requests'], ['req-a0', 'req-a1'])
        self.assertEqual(first['status'], 'quota_exceeded')
        self.assertNotIn('sites', first)
        self.assertEqual(second['status'], 'delivered')
        self.assertEqual(len(second['sites']), 1)


class MultiCellTests(unittest.TestCase):
    """非平凡基底两个胞片：交付同时跨胞片与重复像。"""

    def setUp(self):
        data = {'operations': [IDENTITY], 'sites': [site('a', 'X', [1, 4], [1, 1])]}
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]}]
        self.pkg = package(data, [[2, 0, 0], [0, 1, 0], [0, 0, 1]], links, ((0, 0, 0), (1, 0, 0)))

    def test_delivery_spans_cells_and_replicas(self):
        result = extract_supercell_groups(self.pkg, [3, 1, 1],
                                          [entry('req', A_POINT, [0, 0, 0])], quota=6)
        (group,) = result['groups']
        self.assertEqual(group['status'], 'delivered')
        self.assertEqual(len(group['sites']), 6)
        self.assertEqual(len(group['links']), 6)
        by_key = {(tuple(s['cell']), tuple(s['replica'])): s for s in group['sites']}
        self.assertEqual(by_key[((0, 0, 0), (0, 0, 0))]['position'], ['1/24', '0', '0'])
        self.assertEqual(by_key[((0, 0, 0), (1, 0, 0))]['position'], ['3/8', '0', '0'])
        self.assertEqual(by_key[((1, 0, 0), (0, 0, 0))]['position'], ['5/24', '0', '0'])
        self.assertEqual(by_key[((1, 0, 0), (2, 0, 0))]['position'], ['7/8', '0', '0'])
        wrapped = [link for link in group['links'] if link['winding'] != [0, 0, 0]]
        self.assertEqual(len(wrapped), 1)
        self.assertEqual(wrapped[0]['cell'], [1, 0, 0])
        self.assertEqual(wrapped[0]['replica'], [2, 0, 0])
        self.assertEqual(wrapped[0]['winding'], [1, 0, 0])
        self.assertEqual(wrapped[0]['to'], {'cell': [0, 0, 0], 'point': A_POINT, 'replica': [0, 0, 0]})

    def test_json_round_trip(self):
        result = extract_supercell_groups(self.pkg, [3, 1, 1],
                                          [entry('req', A_POINT, [0, 0, 0])], quota=6)
        self.assertEqual(json.loads(json.dumps(result)), result)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.pkg = package(ab_data(), UNIT, self_links())
        self.requests = [entry('req-a0', A_POINT, [0, 0, 0])]

    def test_repeats_validation(self):
        for bad in ([0, 1, 1], [1, 1], [1, 1, 1, 1], [10 ** 9 + 1, 1, 1], [True, 1, 1], ['6', 1, 1], 'x'):
            with self.assertRaises(ValueError):
                extract_supercell_groups(self.pkg, bad, self.requests, quota=1)

    def test_quota_validation(self):
        for bad in (0, 10001, True, '1', None, 1.5):
            with self.assertRaises(ValueError):
                extract_supercell_groups(self.pkg, [6, 6, 1], self.requests, quota=bad)

    def test_requests_validation(self):
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1],
                                     [entry('r', A_POINT, [0, 0, 0]), entry('r', A_POINT, [1, 0, 0])], quota=6)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', A_POINT, [0, 0, 0], cell=(1, 0, 0))], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', ['1/2', '0', '0'], [0, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', ['2/4', '0', '0'], [0, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', A_POINT, [6, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', A_POINT, [-1, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('r', A_POINT, [True, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [dict(entry('r', A_POINT, [0, 0, 0]), bogus=1)], quota=1)
        with self.assertRaises(ValueError):
            bad = entry('r', A_POINT, [0, 0, 0])
            del bad['replica']
            extract_supercell_groups(self.pkg, [6, 6, 1], [bad], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], [entry('bad id!', A_POINT, [0, 0, 0])], quota=1)

    def test_reserves_validation(self):
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], self.requests,
                                     [entry('m', A_POINT, [0, 0, 0]), entry('m', A_POINT, [1, 0, 0])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], self.requests,
                                     [entry('m', B_POINT, [0, 0, 0], cell=(2, 0, 0))], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], self.requests,
                                     [entry('m', A_POINT, [0, 0, 6])], quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups(self.pkg, [6, 6, 1], self.requests, reserves='not-a-list', quota=1)

    def test_open_package_rejected(self):
        shard = build_supercell_shard(ab_data(), [[2, 0, 0], [0, 1, 0], [0, 0, 1]], self_links(), [[0, 0, 0]])
        with self.assertRaises(ValueError):
            extract_supercell_groups(shard, [6, 6, 1], self.requests, quota=6)

    def test_package_field_validation(self):
        with self.assertRaises(ValueError):
            extract_supercell_groups({'format': 'supercell-shard/1'}, [6, 6, 1], self.requests, quota=1)
        with self.assertRaises(ValueError):
            extract_supercell_groups('not-a-package', [6, 6, 1], self.requests, quota=1)
        tampered = copy.deepcopy(self.pkg)
        tampered['points'][0]['position'] = ['1/3', '0', '0']
        with self.assertRaises(ValueError):
            extract_supercell_groups(tampered, [6, 6, 1], self.requests, quota=1)

    def test_no_input_mutation(self):
        requests = [entry('req-a0', A_POINT, [0, 0, 0]), entry('req-b0', B_POINT, [0, 0, 0])]
        reserves = [entry('keep', A_POINT, [2, 2, 0])]
        snapshot = copy.deepcopy((self.pkg, requests, reserves))
        extract_supercell_groups(self.pkg, [6, 6, 1], requests, reserves, quota=6)
        self.assertEqual((self.pkg, requests, reserves), snapshot)


if __name__ == '__main__':
    unittest.main()
