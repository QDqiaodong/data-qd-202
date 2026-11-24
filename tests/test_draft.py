import copy
import json
import unittest

from site_orbits import process, solve_draft

K_BIG = 999999937


def operations():
    return [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
        {'id': 'invert', 'matrix': [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
    ]


def case1_draft():
    return {
        'scale': 8,
        'operations': operations(),
        'sites': [
            {'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'units': [1, 7]},
            {'id': 'b', 'species': 'X', 'point': [[3, 4], [0, 1], [0, 1]], 'units': [1, 7]},
            {'id': 'c', 'species': 'Y', 'point': [[1, 4], [0, 1], [0, 1]], 'units': [1, 7]},
            {'id': 'd', 'species': 'Y', 'point': [[0, 1], [0, 1], [0, 1]], 'units': [1, 7]},
        ],
        'orbit_observations': [{'site': 'a', 'units': [16, 16]}, {'site': 'd', 'units': [4, 4]}],
        'species_observations': [{'species': 'X', 'units': [10, 14]}, {'species': 'Y', 'units': [8, 8]}],
    }


def merged_request():
    return {'groups': [{'members': ['a', 'b'], 'representative': 'a'},
                       {'members': ['c'], 'representative': 'c'}]}


def split_request():
    return {'groups': [{'members': ['a'], 'representative': 'a'},
                       {'members': ['b'], 'representative': 'b'},
                       {'members': ['c'], 'representative': 'c'}]}


def ranges(result):
    return {source['id']: source['units'] for source in result['sources']}


def check_case1_assignment(test, assignment, f_bounds=(4, 4)):
    """独立复核一份完整可行分配满足全草稿区间、容量与观测。"""
    test.assertEqual(set(assignment), {'a', 'b', 'c', 'd'})
    for value in assignment.values():
        test.assertTrue(1 <= value <= 7)
    test.assertEqual(assignment['a'] + assignment['b'] + assignment['c'], 16)   # 轨道 R 观测
    test.assertLessEqual(assignment['a'] + assignment['b'] + assignment['c'], 16)  # 轨道 R 容量 mK
    test.assertTrue(f_bounds[0] <= assignment['d'] <= f_bounds[1])              # 轨道 F 观测
    test.assertLessEqual(assignment['d'], 8)                                    # 轨道 F 容量 mK
    test.assertTrue(10 <= assignment['a'] + assignment['b'] <= 14)              # 物种 X 观测
    test.assertEqual(assignment['c'] + assignment['d'], 8)                      # 物种 Y 观测


class SolveTests(unittest.TestCase):
    def test_exact_source_ranges(self):
        result = solve_draft(case1_draft())
        self.assertEqual(result['format'], 'occupancy-draft/1')
        self.assertEqual(result['status'], 'feasible')
        self.assertEqual(result['scale'], 8)
        self.assertEqual(result['cell_points'], 3)
        self.assertEqual(result['orbits'], [{'index': 0, 'multiplicity': 2, 'sites': ['a', 'b', 'c']},
                                            {'index': 1, 'multiplicity': 1, 'sites': ['d']}])
        # 输入区间都是 [1,7]，实际范围被全晶胞约束收紧
        self.assertEqual(ranges(result), {'a': [5, 7], 'b': [5, 7], 'c': [4, 4], 'd': [4, 4]})
        self.assertNotIn('delivery', result)

    def test_infeasible_reports_no_ranges(self):
        draft = case1_draft()
        draft['species_observations'][0]['units'] = [13, 14]  # 各记录单看都在格式范围内
        for request in (None, merged_request()):
            result = solve_draft(draft, request)
            self.assertEqual(result, {'format': 'occupancy-draft/1', 'status': 'infeasible'})

    def test_zero_orbit_observation_infeasible(self):
        draft = case1_draft()
        draft['orbit_observations'][0]['units'] = [0, 0]  # 来源区间下界为 1，轨道不可能为空
        self.assertEqual(solve_draft(draft)['status'], 'infeasible')

    def test_json_round_trip(self):
        for result in (solve_draft(case1_draft()),
                       solve_draft(case1_draft(), merged_request()),
                       solve_draft(case1_draft(), split_request())):
            self.assertEqual(json.loads(json.dumps(result)), result)

    def test_reordering_invariance(self):
        base = solve_draft(case1_draft(), merged_request())
        draft = case1_draft()
        draft['operations'].reverse()
        draft['sites'] = draft['sites'][::-1]
        draft['orbit_observations'].reverse()
        draft['species_observations'].reverse()
        request = {'groups': [{'members': ['c'], 'representative': 'c'},
                              {'members': ['b', 'a'], 'representative': 'a'}]}
        self.assertEqual(base, solve_draft(draft, request))
        blocked = solve_draft(case1_draft(), split_request())
        self.assertEqual(blocked, solve_draft(draft, split_request()))

    def test_input_not_modified(self):
        draft = case1_draft()
        request = merged_request()
        draft_copy, request_copy = copy.deepcopy(draft), copy.deepcopy(request)
        solve_draft(draft, request)
        solve_draft(draft)
        self.assertEqual(draft, draft_copy)
        self.assertEqual(request, request_copy)


class DeliveryTests(unittest.TestCase):
    def test_merged_groups_delivered(self):
        result = solve_draft(case1_draft(), merged_request())
        delivery = result['delivery']
        self.assertEqual(delivery['status'], 'delivered')
        self.assertEqual(delivery['groups'], [
            {'representative': 'a', 'members': ['a', 'b'], 'orbit': 0, 'species': 'X',
             'units': 12, 'occupancy': [3, 4]},
            {'representative': 'c', 'members': ['c'], 'orbit': 0, 'species': 'Y',
             'units': 4, 'occupancy': [1, 4]}])
        structure = delivery['structure']
        self.assertEqual([site['id'] for site in structure['sites']], ['a', 'c'])
        expanded = process(structure)  # 旧入口直接接受交付结构
        self.assertEqual([p['point'] for p in expanded['positions']], [['1/4', '0', '0'], ['3/4', '0', '0']])
        for position in expanded['positions']:  # 每点物种占位与任一原始可行分配一致
            self.assertEqual(position['occupancy'], '1')
            totals = {c['species']: c['occupancy'] for c in position['contributions']}
            self.assertEqual(totals, {'X': '3/4', 'Y': '1/4'})
        self.assertEqual(expanded['cell_composition'], {'X': '3/2', 'Y': '1/2'})

    def test_split_groups_blocked_with_witnesses(self):
        result = solve_draft(case1_draft(), split_request())
        delivery = result['delivery']
        self.assertEqual(delivery['status'], 'blocked')
        self.assertEqual(ranges(result), {'a': [5, 7], 'b': [5, 7], 'c': [4, 4], 'd': [4, 4]})
        self.assertEqual([entry['representative'] for entry in delivery['blocking']], ['a', 'b'])
        for entry in delivery['blocking']:
            self.assertEqual(entry['reason'], 'not_fixed')
            self.assertEqual(entry['units'], [5, 7])
            for label, endpoint in (('min', 0), ('max', 1)):
                witness = entry['witnesses'][label]
                check_case1_assignment(self, witness)
                self.assertEqual(sum(witness[m] for m in entry['members']), entry['units'][endpoint])

    def test_unselected_orbit_observation_can_block(self):
        draft = case1_draft()
        draft['orbit_observations'][1]['units'] = [2, 4]  # 只改未选中的 F 观测
        result = solve_draft(draft, merged_request())
        self.assertEqual(result['status'], 'feasible')
        self.assertEqual(ranges(result), {'a': [3, 7], 'b': [3, 7], 'c': [4, 6], 'd': [2, 4]})
        delivery = result['delivery']
        self.assertEqual(delivery['status'], 'blocked')
        by_rep = {entry['representative']: entry for entry in delivery['blocking']}
        self.assertEqual(set(by_rep), {'a', 'c'})
        self.assertEqual(by_rep['a']['units'], [10, 12])
        self.assertEqual(by_rep['c']['units'], [4, 6])
        for entry in delivery['blocking']:
            self.assertEqual(entry['reason'], 'not_fixed')
            for label, endpoint in (('min', 0), ('max', 1)):
                witness = entry['witnesses'][label]
                check_case1_assignment(self, witness, f_bounds=(2, 4))
                self.assertEqual(sum(witness[m] for m in entry['members']), entry['units'][endpoint])

    def test_unrepresentable_fixed_occupancy(self):
        draft = {'scale': K_BIG, 'operations': operations()[:1],
                 'sites': [{'id': 's', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'units': [1, 1]}]}
        result = solve_draft(draft, {'groups': [{'members': ['s'], 'representative': 's'}]})
        self.assertEqual(result['status'], 'feasible')
        self.assertEqual(ranges(result), {'s': [1, 1]})
        blocking = result['delivery']['blocking']
        self.assertEqual(len(blocking), 1)
        entry = blocking[0]
        self.assertEqual(entry['reason'], 'occupancy_unrepresentable')  # 与 not_fixed 区分
        self.assertEqual(entry['units'], [1, 1])
        self.assertEqual(entry['occupancy'], [1, K_BIG])
        self.assertEqual(entry['witnesses']['min'], {'s': 1})
        self.assertEqual(entry['witnesses']['max'], {'s': 1})

    def test_many_orbits_delivery_scales_independently_of_k(self):
        sites, observations, groups = [], [], []
        for i in range(60):
            first, second = f's{i:02d}a', f's{i:02d}b'
            point = [[i, 60], [0, 1], [0, 1]]
            sites.append({'id': first, 'species': 'X', 'point': point, 'units': [1, K_BIG - 1]})
            sites.append({'id': second, 'species': 'X', 'point': point, 'units': [1, K_BIG - 1]})
            observations.append({'site': first, 'units': [K_BIG, K_BIG]})
            groups.append({'members': [first, second], 'representative': first})
        draft = {'scale': K_BIG, 'operations': operations()[:1], 'sites': sites,
                 'orbit_observations': observations,
                 'species_observations': [{'species': 'X', 'units': [60 * K_BIG, 60 * K_BIG]}]}
        result = solve_draft(draft, {'groups': groups})
        self.assertEqual(result['status'], 'feasible')
        self.assertEqual(result['cell_points'], 60)
        for source in result['sources']:  # 实际宽区间，不随 K 的数值铺开
            self.assertEqual(source['units'], [1, K_BIG - 1])
        delivery = result['delivery']
        self.assertEqual(delivery['status'], 'delivered')
        self.assertEqual(len(delivery['groups']), 60)
        self.assertTrue(all(group['occupancy'] == [1, 1] for group in delivery['groups']))
        expanded = process(delivery['structure'])
        self.assertEqual(len(expanded['positions']), 60)
        self.assertEqual(expanded['cell_composition'], {'X': '60'})

    def test_representative_point_and_canonicalization(self):
        request = {'groups': [{'members': ['a', 'b'], 'representative': 'b'},
                              {'members': ['c'], 'representative': 'c'}]}
        result = solve_draft(case1_draft(), request)
        site = result['delivery']['structure']['sites'][0]
        self.assertEqual(site['id'], 'b')
        self.assertEqual(site['point'], [[3, 4], [0, 1], [0, 1]])  # 用代表来源自己的坐标
        draft = case1_draft()
        draft['sites'][0]['point'] = [[2, 8], [0, 1], [0, 1]]       # 非最简坐标写法
        result = solve_draft(draft, merged_request())
        site = result['delivery']['structure']['sites'][0]
        self.assertEqual(site['point'], [[1, 4], [0, 1], [0, 1]])   # 交付结构写规范坐标

    def test_orbit_reference_via_any_member(self):
        draft = case1_draft()
        draft['orbit_observations'][0]['site'] = 'b'   # 同一轨道的另一成员，语义相同
        self.assertEqual(ranges(solve_draft(draft))['a'], [5, 7])


class ValidationTests(unittest.TestCase):
    def assert_invalid(self, mutate, request=None):
        draft = case1_draft()
        mutate(draft)
        with self.assertRaises(ValueError):
            solve_draft(draft, request)

    def test_scale_validation(self):
        for bad in (0, -1, 10 ** 9 + 1, True, 1.5, '8'):
            self.assert_invalid(lambda d, v=bad: d.update(scale=v))

    def test_source_units_validation(self):
        for bad in ([0, 7], [17, 17], [5, 3], [True, 7], [1.5, 7], [1], [1, 2, 3], '1..7', [1, None]):
            self.assert_invalid(lambda d, v=bad: d['sites'][0].update(units=v))

    def test_orbit_observation_validation(self):
        self.assert_invalid(lambda d: d['orbit_observations'].append({'site': 'b', 'units': [1, 16]}))   # 同轨道重复
        self.assert_invalid(lambda d: d['orbit_observations'].append({'site': 'a', 'units': [1, 16]}))   # 同来源重复
        self.assert_invalid(lambda d: d['orbit_observations'].append({'site': 'zz', 'units': [1, 2]}))   # 未知引用
        self.assert_invalid(lambda d: d['orbit_observations'][0].update(units=[0, 17]))                  # 超过 mK
        self.assert_invalid(lambda d: d['orbit_observations'][0].update(units=[-1, 4]))
        self.assert_invalid(lambda d: d['orbit_observations'][0].update(units=[4, 2]))
        self.assert_invalid(lambda d: d['orbit_observations'][0].update(units=[True, 4]))

    def test_species_observation_validation(self):
        self.assert_invalid(lambda d: d['species_observations'].append({'species': 'X', 'units': [1, 2]}))
        self.assert_invalid(lambda d: d['species_observations'].append({'species': 'Q', 'units': [1, 2]}))
        self.assert_invalid(lambda d: d['species_observations'][0].update(units=[0, 25]))   # 超过 K·M=24
        self.assert_invalid(lambda d: d['species_observations'][0].update(units=[-1, 4]))
        self.assert_invalid(lambda d: d['species_observations'][0].update(units=[4, 2]))
        self.assert_invalid(lambda d: d['species_observations'][0].update(units=[1.5, 4]))

    def test_draft_structure_validation(self):
        self.assert_invalid(lambda d: d.update(units=1))
        self.assert_invalid(lambda d: d.update(sites=[]))
        self.assert_invalid(lambda d: d['operations'].pop(0))           # 缺恒等
        self.assert_invalid(lambda d: d['sites'][0].update(occupancy=[1, 2]))
        self.assert_invalid(lambda d: d['sites'].append(dict(d['sites'][0])))  # 重复来源 ID

    def test_group_validation(self):
        def request(groups):
            return {'groups': groups}

        def invalid(groups):
            with self.assertRaises(ValueError):
                solve_draft(case1_draft(), request(groups))

        invalid([])                                                        # 至少一组
        invalid([{'members': [], 'representative': 'a'}])                  # 空组
        invalid([{'members': ['a', 'zz'], 'representative': 'a'},
                 {'members': ['b', 'c'], 'representative': 'b'}])          # 未知成员
        invalid([{'members': ['a', 'b'], 'representative': 'c'}])          # 代表不在组内
        invalid([{'members': ['a', 'a'], 'representative': 'a'}])          # 组内重复
        invalid([{'members': ['a', 'b'], 'representative': 'a'},
                 {'members': ['b', 'c'], 'representative': 'b'}])          # 跨组重复归组
        invalid([{'members': ['a', 'b'], 'representative': 'a'}])          # 遗漏同轨道来源 c
        invalid([{'members': ['a', 'd'], 'representative': 'a'},
                 {'members': ['b', 'c'], 'representative': 'b'}])          # 跨轨道合组
        invalid([{'members': ['a', 'c'], 'representative': 'a'},
                 {'members': ['b'], 'representative': 'b'}])               # 跨物种合组
        with self.assertRaises(ValueError):
            solve_draft(case1_draft(), {'group': []})
        with self.assertRaises(ValueError):
            solve_draft(case1_draft(), 'not-a-request')


if __name__ == '__main__':
    unittest.main()
