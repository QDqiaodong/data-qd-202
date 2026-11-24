"""周期迁移轨迹定位校准与时间窗交接示例；from repository root run:

    .venv/bin/python -I examples/run_handoff.py

恒等操作、单位基底的完整超胞包：A 在 x=1/4、B 在 x=3/4，A→B 有平行边
E（winding (2,0,0)）与 H（零 winding），B→A 有 F（winding (-1,0,0)）。
循环段 M = E 后接 F，根为 M 的两次重复。同一批三个窗口分别演示越界、
可交接与测量冲突；再对 N=10^6 的同形轨迹给出不随轮数增长的交接定义表，
并仅凭交接文档（不使用原包与原轨迹）离线重建获准路线、逐跳复核状态。
产物写入 run/handoff.json，已存在会拒绝覆盖。
"""
import json

from site_orbits import (assemble_supercell_shards, build_supercell_shard,
                         handoff_windows, rebuild_handoff)
from site_orbits.io import publish

DATA = {
    'operations': [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
         'translation': [[0, 1], [0, 1], [0, 1]]},
    ],
    'sites': [
        {'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
        {'id': 'b', 'species': 'Y', 'point': [[3, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
    ],
}
BASIS = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
LINKS = [
    {'id': 'E', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [2, 0, 0]},
    {'id': 'H', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [0, 0, 0]},
    {'id': 'F', 'from': [[3, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [-1, 0, 0]},
]
REGION = [['-1', '2'], ['0', '0'], ['0', '0']]


def step(nid, lid):
    return {'id': nid, 'kind': 'step', 'link': lid, 'cell': [0, 0, 0]}


def trajectory(count):
    nodes = [step('e', 'E'), step('f', 'F'),
             {'id': 'M', 'kind': 'sequence', 'left': 'e', 'right': 'f'}]
    if count > 1:
        nodes.append({'id': 'R', 'kind': 'repeat', 'child': 'M', 'count': count})
        return {'root': 'R', 'nodes': nodes}
    return {'root': 'M', 'nodes': nodes}


def window(wid, l, r, left_image, right_image, region=REGION):
    return {'id': wid, 'l': l, 'r': r,
            'left_image': list(left_image), 'right_image': list(right_image),
            'region': [list(pair) for pair in region]}


if __name__ == '__main__':
    package = assemble_supercell_shards(
        [build_supercell_shard(DATA, BASIS, LINKS, [[0, 0, 0]])], finalize=True)

    windows = [
        # 两端定位相符，但第 1 步在 B（像 2，x=11/4）先越界
        window('w-out-of-region', 0, 2, (0, 0, 0), (1, 0, 0)),
        # 整体下移一格：全程在区域内，可交接 E、F
        window('w-handoff', 0, 2, (-1, 0, 0), (0, 0, 0)),
        # 左端沿用下移值、右端改回 1：与推得值 0 冲突，不查出界、不交接
        window('w-conflict', 0, 2, (-1, 0, 0), (1, 0, 0)),
    ]
    result = handoff_windows(package, trajectory(2), windows)
    publish('run/handoff.json', result)
    print(json.dumps({'format': result['format'], 'total_steps': result['total_steps'],
                      'windows': [{'id': w['id'], 'status': w['status'],
                                   **({'violation_step': w['violation']['step']}
                                      if w['status'] == 'out_of_region' else {}),
                                   **({'derived_right': w['derived_right']}
                                      if w['status'] == 'measurement_conflict' else {})}
                                  for w in result['windows']]}, ensure_ascii=False, sort_keys=True))

    # 仅凭交接文档离线重建（无原包、无原轨迹、无缓存）
    rebuilt = rebuild_handoff(result)
    for hop in rebuilt.iter_hops('w-handoff'):
        print(json.dumps({'hop': hop['index'], 'link': hop['link'],
                          'from_image': hop['from_image'], 'to_image': hop['to_image']},
                         ensure_ascii=False, sort_keys=True))
    print(json.dumps({'end_state': rebuilt.end_state('w-handoff')}, ensure_ascii=False, sort_keys=True))

    # 大重复次数：定义表规模与轮数无关，重复以计数保留
    n = 10 ** 6
    big_region = [['0', str(n + 2)], ['0', '0'], ['0', '0']]
    big = handoff_windows(
        package, trajectory(n),
        [window('w-middle', 1, 2 * n - 1, (2, 0, 0), (n + 1, 0, 0), big_region)])
    repeat_defs = [d for d in big['handoff']['definitions'] if d['kind'] == 'repeat']
    print(json.dumps({'n': n, 'definitions': len(big['handoff']['definitions']),
                      'links': len(big['handoff']['links']),
                      'repeat_count': repeat_defs[0]['count'],
                      'hops': big['handoff']['routes'][0]['hops']}, sort_keys=True))
