"""周期结构整组抽取与离线作业交付示例；run from repository root.

恒等操作、单位基底的完整包：A 在 x=1/4，带两条自连边（原胞平移 (2,1,0) 与
(8,1,0)），B 在 x=3/4 无连接。重复 (6,6,1) 次后依次请求 A 的像 (0,0,0)、
A 的像 (4,2,0) 与 B 的像 (0,0,0)，作业额度 6：两个 A 请求同组，保留标记落在
别的重复像上时整组交付 6 点 12 边，B 额度不足跳过；把标记移进 A 组则整组保留，
B 获得额度交付。产物写入 run/extract.json，已存在会拒绝覆盖。
"""
import json

from site_orbits import assemble_supercell_shards, build_supercell_shard, extract_supercell_groups
from site_orbits.io import publish

DATA = {
    'operations': [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
    ],
    'sites': [
        {'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
        {'id': 'b', 'species': 'Y', 'point': [[3, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]},
    ],
}
BASIS = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
LINKS = [
    {'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [2, 1, 0]},
    {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[1, 4], [0, 1], [0, 1]], 'shift': [8, 1, 0]},
]
A = {'cell': [0, 0, 0], 'point': ['1/4', '0', '0']}
B = {'cell': [0, 0, 0], 'point': ['3/4', '0', '0']}


def entry(eid, site, replica):
    return {'id': eid, **site, 'replica': replica}


if __name__ == '__main__':
    package = assemble_supercell_shards([build_supercell_shard(DATA, BASIS, LINKS, [[0, 0, 0]])], finalize=True)
    requests = [entry('req-a0', A, [0, 0, 0]), entry('req-a1', A, [4, 2, 0]), entry('req-b0', B, [0, 0, 0])]
    result = extract_supercell_groups(package, [6, 6, 1], requests, [entry('keep-1', A, [2, 2, 0])], quota=6)
    publish('run/extract.json', result)
    print(json.dumps({'groups': [{'requests': g['requests'], 'status': g['status']} for g in result['groups']],
                      'delivered': [{'sites': len(g.get('sites', [])), 'links': len(g.get('links', []))}
                                    for g in result['groups']]}, sort_keys=True))
    moved = extract_supercell_groups(package, [6, 6, 1], requests, [entry('keep-1', A, [2, 1, 0])], quota=6)
    print(json.dumps({'marker_moved': [g['status'] for g in moved['groups']],
                      'reserved_ids': moved['groups'][0].get('reserves')}, sort_keys=True))
