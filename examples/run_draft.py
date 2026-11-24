"""占位区间草稿求解与确值结构交付示例；run from repository root.

K=8、恒等+反演的草稿：轨道 R（x=1/4 与 x=3/4）有 X 来源 a、b 与 Y 来源 c，
固定点轨道 F（x=0）有 Y 来源 d；各源区间 [1,7]，轨道与物种总量带观测。
先只求解看每条来源的实际单位范围；再按调用方分组 {a,b}、{c} 交付确值结构
（未选中的 F 不进入输出），并演示把 {a,b} 拆成两组时被拒绝、但拿到达到两个
端点的完整可行分配见证。产物写入 run/draft.json，已存在会拒绝覆盖。
"""
import json

from site_orbits import process, solve_draft
from site_orbits.io import publish

DRAFT = {
    'scale': 8,
    'operations': [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
        {'id': 'invert', 'matrix': [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
    ],
    'sites': [
        {'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'units': [1, 7]},
        {'id': 'b', 'species': 'X', 'point': [[3, 4], [0, 1], [0, 1]], 'units': [1, 7]},
        {'id': 'c', 'species': 'Y', 'point': [[1, 4], [0, 1], [0, 1]], 'units': [1, 7]},
        {'id': 'd', 'species': 'Y', 'point': [[0, 1], [0, 1], [0, 1]], 'units': [1, 7]},
    ],
    'orbit_observations': [{'site': 'a', 'units': [16, 16]}, {'site': 'd', 'units': [4, 4]}],
    'species_observations': [{'species': 'X', 'units': [10, 14]}, {'species': 'Y', 'units': [8, 8]}],
}
MERGED = {'groups': [{'members': ['a', 'b'], 'representative': 'a'},
                     {'members': ['c'], 'representative': 'c'}]}
SPLIT = {'groups': [{'members': ['a'], 'representative': 'a'},
                    {'members': ['b'], 'representative': 'b'},
                    {'members': ['c'], 'representative': 'c'}]}

if __name__ == '__main__':
    solved = solve_draft(DRAFT)
    print(json.dumps({'status': solved['status'],
                      'ranges': {s['id']: s['units'] for s in solved['sources']}}, sort_keys=True))
    delivered = solve_draft(DRAFT, MERGED)
    publish('run/draft.json', delivered)
    expanded = process(delivered['delivery']['structure'])
    print(json.dumps({'delivery': delivered['delivery']['status'],
                      'composition': expanded['cell_composition']}, sort_keys=True))
    blocked = solve_draft(DRAFT, SPLIT)['delivery']
    print(json.dumps({'delivery': blocked['status'],
                      'blocking': [{'representative': e['representative'], 'units': e['units'],
                                    'min': e['witnesses']['min'], 'max': e['witnesses']['max']}
                                   for e in blocked['blocking']]}, sort_keys=True))
