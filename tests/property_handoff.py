"""随机化性质测试：压缩轨迹的交接判定对朴素逐步模拟逐条等价。

在固定的 A/B 两点包（含自连边与双向边）上生成带嵌套重复与串接的随机轨迹，
保证每个串接闭合、每个重复子结点包内身份闭合；对随机窗口、定位与区域比较：
定位推得像、最早越界步/状态/坐标与获准路线的逐跳重建序列。运行：

    .venv/bin/python -I tests/property_handoff.py
"""
import json
import random
import sys
from fractions import Fraction

from site_orbits import build_supercell_shard, assemble_supercell_shards, handoff_windows, rebuild_handoff

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'translation': [[0, 1], [0, 1], [0, 1]]}
UNIT = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
A = Fraction(1, 4)
B = Fraction(3, 4)
DATA = {'operations': [IDENTITY],
        'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
                  {'id': 'b', 'species': 'Y', 'point': [[3, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]}]}
LINKS = [
    {'id': 'E', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]},
    {'id': 'H', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [0, 0, 0]},
    {'id': 'F', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [-1, 0, 0]},
    {'id': 'G', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [-2, 0, 0]},
    {'id': 'SA', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [1, 1, 0]},
    {'id': 'SB', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [1, -1, 0]},
]
A_SITE = ((0, 0, 0), (A, Fraction(0), Fraction(0)))
B_SITE = ((0, 0, 0), (B, Fraction(0), Fraction(0)))
OUTGOING = {A_SITE: ('E', 'H', 'SA'), B_SITE: ('F', 'G', 'SB')}
TARGET = {'E': B_SITE, 'H': B_SITE, 'SA': A_SITE, 'F': A_SITE, 'G': A_SITE, 'SB': B_SITE}


def package():
    shard = build_supercell_shard(DATA, UNIT, LINKS, [[0, 0, 0]])
    return assemble_supercell_shards([shard], finalize=True)


class TrajBuilder:
    """随机构造闭合 DAG：每个生成的子程序都带起止位点身份。"""

    def __init__(self, rng):
        self.rng = rng
        self.nodes = []
        self.counter = 0

    def _fresh(self, prefix):
        self.counter += 1
        return f'{prefix}{self.counter}'

    def step_from(self, site):
        lid = self.rng.choice(OUTGOING[site])
        nid = self._fresh('t')
        self.nodes.append({'id': nid, 'kind': 'step', 'link': lid, 'cell': [0, 0, 0]})
        return nid, TARGET[lid]

    def program_from(self, site, depth):
        """返回 (结点ID, 执行后位点身份)；保证起点为 site。"""
        choice = self.rng.choice(['step', 'sequence', 'repeat', 'selfloop']) if depth > 0 else 'step'
        if choice in ('step', 'selfloop'):
            return self.step_from(site)
        if choice == 'sequence':
            left_id, mid = self.program_from(site, depth - 1)
            right_id, end = self.program_from(mid, depth - 1)
            nid = self._fresh('q')
            self.nodes.append({'id': nid, 'kind': 'sequence', 'left': left_id, 'right': right_id})
            return nid, end
        # repeat：先造一个回到同身份的子程序
        child_id, child_end = self.program_from(site, depth - 1)
        if child_end != site:
            # 不满足身份闭合就退化为普通串接（仍合法）
            return child_id, child_end
        nid = self._fresh('r')
        self.nodes.append({'id': nid, 'kind': 'repeat', 'child': child_id,
                           'count': self.rng.randint(2, 6)})
        return nid, site

    def build(self):
        root_id, _ = self.program_from(A_SITE, depth=4)
        return {'root': root_id, 'nodes': self.nodes}


def total_steps(traj):
    defs = {n['id']: n for n in traj['nodes']}
    memo = {}

    def steps(nid):
        if nid in memo:
            return memo[nid]
        n = defs[nid]
        if n['kind'] == 'step':
            value = 1
        elif n['kind'] == 'sequence':
            value = steps(n['left']) + steps(n['right'])
        else:
            value = n['count'] * steps(n['child'])
        memo[nid] = value
        return value

    return steps(traj['root'])


def naive_run(pkg, traj, start_image):
    defs = {n['id']: n for n in traj['nodes']}
    records = {(r['link'], tuple(r['cell'])): r for r in pkg['links']}
    positions = {(tuple(p['cell']), tuple(Fraction(v) for v in p['point'])):
                 tuple(Fraction(v) for v in p['position']) for p in pkg['points']}
    hops, states = [], []

    def site_of(endpoint):
        return tuple(endpoint['cell']), tuple(Fraction(v) for v in endpoint['point'])

    def walk(nid):
        n = defs[nid]
        if n['kind'] == 'step':
            rec = records[(n['link'], tuple(n['cell']))]
            hops.append((n['link'], tuple(n['cell']), tuple(rec['winding']),
                         site_of(rec['from']), site_of(rec['to'])))
        elif n['kind'] == 'sequence':
            walk(n['left'])
            walk(n['right'])
        else:
            for _ in range(n['count']):
                walk(n['child'])

    walk(traj['root'])
    image = tuple(start_image)
    site = hops[0][3]
    states.append((site, image))
    for lid, cell, winding, from_site, to_site in hops:
        assert site == from_site
        image = tuple(image[i] + winding[i] for i in range(3))
        site = to_site
        states.append((site, image))
    return hops, states, positions


def main():
    rng = random.Random(20261124)
    pkg = package()
    trials = 500
    for trial in range(trials):
        traj = TrajBuilder(rng).build()
        total = total_steps(traj)
        assert total <= 6 ** 5
        hops, states, positions = naive_run(pkg, traj, (0, 0, 0))
        window_count = rng.randint(1, 4)
        windows = []
        expectations = []
        for wi in range(window_count):
            l = rng.randint(0, total - 1)
            r = rng.randint(l + 1, total)
            left_image = tuple(rng.randint(-6, 6) for _ in range(3))
            net = tuple(states[r][1][i] - states[l][1][i] for i in range(3))
            if rng.random() < 0.75:
                right_image = tuple(left_image[i] + net[i] for i in range(3))
            else:
                right_image = tuple(left_image[i] + rng.randint(-3, 3) for i in range(3))
            region = tuple((Fraction(rng.randint(-8, 2)), Fraction(rng.randint(3, 12)))
                           for _ in range(3))
            wid = f'w{trial}_{wi}'
            windows.append({'id': wid, 'l': l, 'r': r,
                            'left_image': list(left_image), 'right_image': list(right_image),
                            'region': [[str(a), str(b)] for a, b in region]})
            expectations.append((wid, l, r, left_image, right_image, region))
        result = handoff_windows(pkg, traj, windows)
        assert json.loads(json.dumps(result)) == result
        by_id = {w['id']: w for w in result['windows']}
        for wid, l, r, left_image, right_image, region in expectations:
            status = by_id[wid]
            net = tuple(states[r][1][i] - states[l][1][i] for i in range(3))
            derived = tuple(left_image[i] + net[i] for i in range(3))
            if derived != tuple(right_image):
                assert status['status'] == 'measurement_conflict'
                assert tuple(status['derived_right']) == derived
                continue
            root_base = tuple(left_image[i] - states[l][1][i] for i in range(3))
            expected = None
            for t in range(l, r + 1):
                site, rel_image = states[t]
                image = tuple(root_base[i] + rel_image[i] for i in range(3))
                coordinate = tuple(positions[site][i] + image[i] for i in range(3))
                if not all(region[i][0] <= coordinate[i] <= region[i][1] for i in range(3)):
                    expected = (t, site, image, coordinate)
                    break
            if expected is not None:
                assert status['status'] == 'out_of_region'
                v = status['violation']
                assert v['step'] == expected[0]
                assert tuple(v['image']) == expected[2]
                assert tuple(Fraction(x) for x in v['coordinate']) == expected[3]
                assert (tuple(v['site']['cell']),
                        tuple(Fraction(x) for x in v['site']['point'])) == expected[1]
                continue
            assert status['status'] == 'handoff'
            rebuilt = rebuild_handoff(result)
            got_hops = list(rebuilt.iter_hops(wid))
            assert [(h['link'], tuple(h['cell']), tuple(h['winding'])) for h in got_hops] == \
                   [(h[0], h[1], h[2]) for h in hops[l:r]]
            for t in range(r - l + 1):
                got = rebuilt.state_at(wid, t)
                site = states[l + t][0]
                rel = states[l + t][1]
                image = tuple(root_base[i] + rel[i] for i in range(3))
                assert tuple(got['image']) == image
                assert (tuple(got['site']['cell']),
                        tuple(Fraction(x) for x in got['site']['point'])) == site
            assert tuple(rebuilt.end_state(wid)['image']) == \
                   tuple(root_base[i] + states[r][1][i] for i in range(3))
        # 定义表可达闭包已由 rebuild_handoff 校验；额外核对无多余定义
        if 'handoff' in result:
            live = rebuild_handoff(result)
            reachable = set()
            stack = [route['program'] for route in result['handoff']['routes']]
            while stack:
                did = stack.pop()
                if did in reachable:
                    continue
                reachable.add(did)
                term = live._defs[did]
                if term[0] == 'sequence':
                    stack.extend([term[1], term[2]])
                elif term[0] == 'repeat':
                    stack.append(term[1])
            assert reachable == {d['id'] for d in result['handoff']['definitions']}
    print(f'{trials} randomized trials passed')


if __name__ == '__main__':
    sys.exit(main())
