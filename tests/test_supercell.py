import copy
import json
import unittest

from site_orbits import assemble_supercell_shards, build_supercell_shard, process

BASIS = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]


def site(sid, species, point, occupancy):
    return {'id': sid, 'species': species, 'point': point, 'occupancy': occupancy}


def operations():
    return [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
        {'id': 'invert', 'matrix': [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
    ]


def dimer_data():
    return {'operations': operations(), 'sites': [site('a', 'X', [[1, 4], [0, 1], [0, 1]], [1, 1])]}


def dimer_links():
    return [
        {'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]},
        {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [3, 0, 0]},
    ]


def dimer_shards(cells=((0, 0, 0), (1, 0, 0))):
    return [build_supercell_shard(dimer_data(), BASIS, dimer_links(), [list(cell)]) for cell in cells]


def roundtrip(shard):
    return json.loads(json.dumps(shard))


class ShardBuildTests(unittest.TestCase):
    def test_canonical_representative_and_duplicate(self):
        direct = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[1, 0, 0]])
        wrapped = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[-1, 0, 0]])
        self.assertEqual(direct, wrapped)
        self.assertEqual(direct['cells'], [[1, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[-1, 0, 0], [1, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[0, 0, 0], [0, 0, 0]])

    def test_shard_boundary_records(self):
        shard = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[0, 0, 0]])
        self.assertEqual(shard['format'], 'supercell-shard/1')
        self.assertTrue(shard['version'].startswith('sha256:'))
        self.assertEqual(shard['cell_count'], 2)
        self.assertEqual(shard['cells'], [[0, 0, 0]])
        self.assertEqual([p['position'] for p in shard['points']], [['1/8', '0', '0'], ['3/8', '0', '0']])
        self.assertEqual(shard['links'], [])
        self.assertEqual(len(shard['out_links']), 2)
        self.assertEqual(len(shard['in_links']), 2)
        by_link = {r['link']: r for r in shard['out_links']}
        # 两条出端指向同一个 B 像点（胞片 (1,0,0) 内新 x 坐标 7/8），winding 不同
        self.assertEqual(by_link['L1']['to'], {'cell': [1, 0, 0], 'point': ['3/4', '0', '0']})
        self.assertEqual(by_link['L1']['winding'], [0, 0, 0])
        self.assertEqual(by_link['L2']['to'], {'cell': [1, 0, 0], 'point': ['3/4', '0', '0']})
        self.assertEqual(by_link['L2']['winding'], [1, 0, 0])
        other = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[1, 0, 0]])
        image = [p for p in other['points'] if p['point'] == ['3/4', '0', '0']]
        self.assertEqual(image[0]['position'], ['7/8', '0', '0'])
        # 入端按实例记录，winding 与对侧出端一一对应，不能混用
        in_windings = {r['link']: r['winding'] for r in shard['in_links']}
        self.assertEqual(in_windings, {'L1': [1, 0, 0], 'L2': [2, 0, 0]})
        out_windings = {r['link']: r['winding'] for r in other['out_links']}
        self.assertEqual(out_windings, in_windings)

    def test_open_assembly_and_finalize(self):
        shards = dimer_shards()
        opened = assemble_supercell_shards(shards[:1])
        self.assertEqual(opened, shards[0])
        with self.assertRaises(ValueError):
            assemble_supercell_shards(shards[:1], finalize=True)
        restored = [roundtrip(s) for s in shards]
        merged = assemble_supercell_shards(restored)
        self.assertEqual(len(merged['links']), 4)
        self.assertEqual(merged['out_links'], [])
        self.assertEqual(merged['in_links'], [])
        self.assertEqual(merged['cells'], [[0, 0, 0], [1, 0, 0]])
        self.assertEqual(assemble_supercell_shards(restored, finalize=True), merged)
        oneshot = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[0, 0, 0], [1, 0, 0]])
        self.assertEqual(merged, oneshot)
        windings = {(r['link'], tuple(r['cell'])): r['winding'] for r in merged['links']}
        self.assertEqual(windings, {('L1', (0, 0, 0)): [0, 0, 0], ('L1', (1, 0, 0)): [1, 0, 0],
                                    ('L2', (0, 0, 0)): [1, 0, 0], ('L2', (1, 0, 0)): [2, 0, 0]})

    def test_tampered_half_rejected(self):
        restored = [roundtrip(s) for s in dimer_shards()]
        deleted = copy.deepcopy(restored)
        del deleted[1]['in_links'][0]
        with self.assertRaises(ValueError):
            assemble_supercell_shards(deleted)
        changed = copy.deepcopy(restored)
        changed[1]['in_links'][0]['winding'] = [9, 0, 0]
        with self.assertRaises(ValueError):
            assemble_supercell_shards(changed)
        moved = copy.deepcopy(restored)
        moved[1]['in_links'][0]['to'] = {'cell': [0, 0, 0], 'point': ['1/4', '0', '0']}
        with self.assertRaises(ValueError):
            assemble_supercell_shards(moved)
        # 恢复原记录即可交付
        final = assemble_supercell_shards(restored, finalize=True)
        self.assertEqual(final['cells'], [[0, 0, 0], [1, 0, 0]])

    def test_version_mix_rejected(self):
        old = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[0, 0, 0]])
        old_other = build_supercell_shard(dimer_data(), BASIS, dimer_links(), [[1, 0, 0]])
        # 增加一条同点零平移连接：不产生跨片记录，边界记录与旧包完全吻合
        extra = dimer_links() + [{'id': 'L3', 'from': [[1, 4], [0, 1], [0, 1]],
                                  'to': [[1, 4], [0, 1], [0, 1]], 'shift': [0, 0, 0]}]
        new = build_supercell_shard(dimer_data(), BASIS, extra, [[1, 0, 0]])
        self.assertEqual(new['out_links'], old_other['out_links'])
        self.assertEqual(new['in_links'], old_other['in_links'])
        self.assertNotEqual(old['version'], new['version'])
        with self.assertRaises(ValueError):
            assemble_supercell_shards([old, new])

    def test_staged_assembly(self):
        basis = [[3, 0, 0], [0, 1, 0], [0, 0, 1]]
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]},
                 {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]}]
        shards = [build_supercell_shard(dimer_data(), basis, links, [[c, 0, 0]]) for c in range(3)]
        stage = assemble_supercell_shards(shards[:2])
        # 第一阶段保留尚未到达覆盖范围的另一端
        self.assertEqual({(r['link'], tuple(r['cell'])) for r in stage['out_links']}, {('L1', (1, 0, 0)), ('L2', (0, 0, 0))})
        self.assertEqual({(r['link'], tuple(r['cell'])) for r in stage['in_links']}, {('L1', (2, 0, 0)), ('L2', (2, 0, 0))})
        final = assemble_supercell_shards([stage, shards[2]])
        oneshot = build_supercell_shard(dimer_data(), basis, links, [[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        self.assertEqual(final, oneshot)
        self.assertEqual(len(final['links']), 6)
        self.assertEqual(assemble_supercell_shards([shards[2], shards[0], shards[1]]), oneshot)
        self.assertEqual(assemble_supercell_shards(shards, finalize=True), oneshot)

    def test_winding_carry(self):
        data = {'operations': operations()[:1],
                'sites': [site('a', 'X', [[1, 4], [0, 1], [0, 1]], [1, 1]),
                          site('b', 'Y', [[3, 4], [1, 2], [0, 1]], [1, 1])]}
        basis = [[1, 2, 0], [0, 1, 0], [0, 0, 1]]
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [1, 2], [0, 1]], 'shift': [0, 0, 0]}]
        shard = build_supercell_shard(data, basis, links, [[0, 0, 0]])
        self.assertEqual(shard['cell_count'], 1)
        self.assertEqual(len(shard['links']), 1)
        record = shard['links'][0]
        # 归一化进位：零 shift 的胞内连接 winding 仍为 (-1,0,0)
        self.assertEqual(record['winding'], [-1, 0, 0])
        self.assertEqual(record['from'], {'cell': [0, 0, 0], 'point': ['1/4', '0', '0']})
        self.assertEqual(record['to'], {'cell': [0, 0, 0], 'point': ['3/4', '1/2', '0']})
        positions = {tuple(p['point']): p['position'] for p in shard['points']}
        self.assertEqual(positions[('1/4', '0', '0')], ['1/4', '0', '0'])
        self.assertEqual(positions[('3/4', '1/2', '0')], ['3/4', '1/2', '0'])
        self.assertEqual(assemble_supercell_shards([roundtrip(shard)], finalize=True), shard)

    def test_special_and_mixed_occupancy(self):
        data = {'operations': operations(),
                'sites': [site('o', 'Y', [[0, 1], [0, 1], [0, 1]], [1, 1]),
                          site('x1', 'X', [[1, 4], [0, 1], [0, 1]], [1, 4]),
                          site('x2', 'Z', [[1, 4], [0, 1], [0, 1]], [1, 2])]}
        self.assertEqual(len(process(data)['positions']), 3)
        shard = build_supercell_shard(data, BASIS, [], [[0, 0, 0], [1, 0, 0]])
        self.assertEqual(len(shard['points']), 6)
        by_point = {(tuple(p['cell']), tuple(p['point'])): p for p in shard['points']}
        special = by_point[((0, 0, 0), ('0', '0', '0'))]
        self.assertEqual(special['occupancy'], '1')
        self.assertEqual(special['vacancy'], '0')
        self.assertEqual(len(special['contributions']), 1)
        self.assertEqual(special['contributions'][0]['operations'], ['identity', 'invert'])
        for cell in ((0, 0, 0), (1, 0, 0)):
            mixed = by_point[(cell, ('1/4', '0', '0'))]
            self.assertEqual(mixed['occupancy'], '3/4')
            self.assertEqual(mixed['vacancy'], '1/4')
            self.assertEqual([(c['site'], c['occupancy']) for c in mixed['contributions']],
                             [('x1', '1/4'), ('x2', '1/2')])
        self.assertEqual(by_point[((1, 0, 0), ('1/4', '0', '0'))]['position'], ['5/8', '0', '0'])

    def test_overoccupancy_rejected(self):
        data = {'operations': operations(),
                'sites': [site('x1', 'X', [[1, 4], [0, 1], [0, 1]], [1, 2]),
                          site('x2', 'Z', [[1, 4], [0, 1], [0, 1]], [3, 4])]}
        with self.assertRaisesRegex(ValueError, 'overoccupied'):
            build_supercell_shard(data, BASIS, [], [[0, 0, 0]])

    def test_sparse_shards(self):
        data = {'operations': operations()[:1], 'sites': [site('a', 'X', [[1, 4], [0, 1], [0, 1]], [1, 1])]}
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]}]
        basis = [[10 ** 9, 10 ** 9, 0], [0, 1, 0], [0, 0, 1]]
        shard = build_supercell_shard(data, basis, links, [[-1, 0, 0], [0, 0, 0]])
        self.assertEqual(shard['cell_count'], 10 ** 9)
        self.assertEqual(shard['cells'], [[0, 0, 0], [999999999, 0, 0]])
        self.assertEqual(len(shard['points']), 2)
        # 允许描述跨向未选胞片的连接，但不生成那些胞片的点清单
        self.assertEqual(len(shard['links']), 1)
        self.assertEqual(shard['links'][0]['winding'], [1, 0, 0])
        self.assertEqual(len(shard['out_links']), 1)
        self.assertEqual(shard['out_links'][0]['to']['cell'], [1, 0, 0])
        self.assertEqual(len(shard['in_links']), 1)
        self.assertEqual(shard['in_links'][0]['cell'], [999999998, 0, 0])
        self.assertEqual({tuple(p['cell']) for p in shard['points']}, {(0, 0, 0), (999999999, 0, 0)})
        self.assertEqual(assemble_supercell_shards([shard]), shard)
        with self.assertRaises(ValueError):
            assemble_supercell_shards([shard], finalize=True)

    def test_single_cell_finalize(self):
        data = {'operations': operations()[:1], 'sites': [site('a', 'X', [[1, 4], [0, 1], [0, 1]], [1, 1])]}
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [0, 1, 0]}]
        basis = [[1, 10 ** 9, 0], [0, 1, 0], [0, 0, 1]]
        shard = build_supercell_shard(data, basis, links, [[7, 3, 0]])
        self.assertEqual(shard['cell_count'], 1)
        self.assertEqual(shard['cells'], [[0, 0, 0]])
        self.assertEqual(shard['links'][0]['winding'], [-10 ** 9, 1, 0])
        self.assertEqual(assemble_supercell_shards([roundtrip(shard)], finalize=True), shard)

    def test_version_invariance(self):
        data, links = dimer_data(), dimer_links()
        base = build_supercell_shard(data, BASIS, links, [[0, 0, 0]])['version']
        self.assertEqual(base, build_supercell_shard(data, BASIS, links, [[1, 0, 0]])['version'])
        shuffled = copy.deepcopy(data)
        shuffled['operations'] = shuffled['operations'][::-1]
        self.assertEqual(base, build_supercell_shard(shuffled, BASIS, links, [[0, 0, 0]])['version'])
        self.assertEqual(base, build_supercell_shard(data, BASIS, links[::-1], [[0, 0, 0]])['version'])
        equivalent_data = copy.deepcopy(data)
        equivalent_data['sites'][0]['point'][0] = [5, 4]
        self.assertEqual(base, build_supercell_shard(equivalent_data, BASIS, links, [[0, 0, 0]])['version'])
        equivalent_links = copy.deepcopy(links)
        equivalent_links[0]['from'][0] = [-3, 4]
        equivalent_links[1]['to'][0] = [6, 8]
        self.assertEqual(base, build_supercell_shard(data, BASIS, equivalent_links, [[0, 0, 0]])['version'])
        changed = copy.deepcopy(links)
        changed[0]['shift'] = [2, 0, 0]
        self.assertNotEqual(base, build_supercell_shard(data, BASIS, changed, [[0, 0, 0]])['version'])
        self.assertNotEqual(base, build_supercell_shard(data, [[3, 0, 0], [0, 1, 0], [0, 0, 1]], links, [[0, 0, 0]])['version'])

    def test_build_validation(self):
        data, links = dimer_data(), dimer_links()
        for bad_basis in ([[0, 0, 0], [0, 1, 0], [0, 0, 1]],
                          [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          [[10 ** 9, 0, 0], [0, 2, 0], [0, 0, 1]],
                          [[10 ** 9 + 1, 0, 0], [0, 1, 0], [0, 0, 1]]):
            with self.assertRaises(ValueError):
                build_supercell_shard(data, bad_basis, links, [[0, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, links, [])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, links, [[i, 0, 0] for i in range(65)])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, links, [[10 ** 9 + 1, 0, 0]])
        many = [dict(link, id=f'L{i}') for i, link in enumerate(dimer_links() * 501)]
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, many, [[0, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, links + [dict(links[0])], [[0, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, [dict(links[0], to=[[1, 2], [0, 1], [0, 1]])], [[0, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, [dict(links[0], shift=[10 ** 9 + 1, 0, 0])], [[0, 0, 0]])
        with self.assertRaises(ValueError):
            build_supercell_shard(data, BASIS, [dict(links[0], id='bad id!')], [[0, 0, 0]])


class ShardAssemblyTests(unittest.TestCase):
    def test_assemble_validation(self):
        shards = dimer_shards()
        with self.assertRaises(ValueError):
            assemble_supercell_shards([])
        with self.assertRaises(ValueError):
            assemble_supercell_shards(shards[0])
        with self.assertRaises(ValueError):
            assemble_supercell_shards([shards[0], shards[0]])
        bad_version = copy.deepcopy(shards[1])
        bad_version['version'] = 'sha256:' + '0' * 64
        with self.assertRaises(ValueError):
            assemble_supercell_shards([shards[0], bad_version])
        extra = copy.deepcopy(shards[0])
        extra['bogus'] = 1
        with self.assertRaises(ValueError):
            assemble_supercell_shards([extra])
        wrong_format = copy.deepcopy(shards[0])
        wrong_format['format'] = 'supercell-shard/0'
        with self.assertRaises(ValueError):
            assemble_supercell_shards([wrong_format])
        wrong_count = copy.deepcopy(shards[0])
        wrong_count['cell_count'] = 3
        with self.assertRaises(ValueError):
            assemble_supercell_shards([wrong_count])
        not_canonical = copy.deepcopy(shards[0])
        not_canonical['cells'] = [[3, 0, 0]]
        with self.assertRaises(ValueError):
            assemble_supercell_shards([not_canonical])

    def test_local_reference_checks(self):
        shard, = dimer_shards(((0, 0, 0),))
        wrong_pos = copy.deepcopy(shard)
        wrong_pos['points'][0]['position'] = ['1/3', '0', '0']
        with self.assertRaises(ValueError):
            assemble_supercell_shards([wrong_pos])
        dropped = copy.deepcopy(shard)
        dropped['points'] = dropped['points'][:1]
        with self.assertRaises(ValueError):
            assemble_supercell_shards([dropped])
        # 内边悬空：出端记录冒充完整连接
        moved = copy.deepcopy(shard)
        moved['links'] = moved['out_links']
        moved['out_links'] = []
        with self.assertRaises(ValueError):
            assemble_supercell_shards([moved])
        # 错误方向：入端记录冒充出端
        swapped = copy.deepcopy(shard)
        swapped['out_links'] = swapped['in_links']
        swapped['in_links'] = []
        with self.assertRaises(ValueError):
            assemble_supercell_shards([swapped])
        # 两端均不属于本片
        alien = copy.deepcopy(shard)
        alien['out_links'] = [{'link': 'L9', 'cell': [1, 0, 0],
                               'from': {'cell': [1, 0, 0], 'point': ['1/4', '0', '0']},
                               'to': {'cell': [0, 0, 0], 'point': ['1/4', '0', '0']},
                               'winding': [0, 0, 0]}]
        with self.assertRaises(ValueError):
            assemble_supercell_shards([alien])
        bad_sum = copy.deepcopy(shard)
        bad_sum['points'][0]['vacancy'] = '1/2'
        with self.assertRaises(ValueError):
            assemble_supercell_shards([bad_sum])

    def test_merged_coverage_limit(self):
        data = {'operations': operations()[:1], 'sites': [site('a', 'X', [[1, 4], [0, 1], [0, 1]], [1, 1])]}
        basis = [[65, 0, 0], [0, 1, 0], [0, 0, 1]]
        first = build_supercell_shard(data, basis, [], [[i, 0, 0] for i in range(64)])
        second = build_supercell_shard(data, basis, [], [[64, 0, 0]])
        with self.assertRaises(ValueError):
            assemble_supercell_shards([first, second])

    def test_no_input_mutation(self):
        data, links, cells = dimer_data(), dimer_links(), [[0, 0, 0], [1, 0, 0]]
        snapshot = (copy.deepcopy(data), copy.deepcopy(links), copy.deepcopy(cells))
        build_supercell_shard(data, BASIS, links, cells)
        self.assertEqual((data, links, cells), snapshot)
        shards = dimer_shards()
        before = copy.deepcopy(shards)
        assemble_supercell_shards(shards, finalize=True)
        self.assertEqual(shards, before)

    def test_serialization_determinism(self):
        basis = [[3, 0, 0], [0, 1, 0], [0, 0, 1]]
        links = [{'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]},
                 {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]}]
        shards = [build_supercell_shard(dimer_data(), basis, links, [[c, 0, 0]]) for c in range(3)]
        forward = assemble_supercell_shards([roundtrip(s) for s in shards])
        shuffled = assemble_supercell_shards([roundtrip(s) for s in shards[::-1]])
        staged = assemble_supercell_shards([assemble_supercell_shards(shards[:2]), shards[2]])
        oneshot = build_supercell_shard(dimer_data(), basis, links, [[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        blob = json.dumps(oneshot, sort_keys=True)
        for result in (forward, shuffled, staged):
            self.assertEqual(json.dumps(result, sort_keys=True), blob)


if __name__ == '__main__':
    unittest.main()
