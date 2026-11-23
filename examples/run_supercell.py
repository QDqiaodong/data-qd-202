"""离线超胞分片交付示例；run from repository root.

生成方按胞片分别产出分片包，经 JSON 写盘模拟离线交付；接收方只读回包，
不接触原胞输入，先拼出开放分片，齐全后再交付完整周期结构。
"""
import json
import tempfile
from pathlib import Path

from site_orbits import assemble_supercell_shards, build_supercell_shard
from site_orbits.io import load, publish

DATA = {
    'operations': [
        {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
        {'id': 'invert', 'matrix': [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], 'translation': [[0, 1], [0, 1], [0, 1]]},
    ],
    'sites': [{'id': 'a', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]], 'occupancy': [1, 1]}],
}
BASIS = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]
LINKS = [
    {'id': 'L1', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [1, 0, 0]},
    {'id': 'L2', 'from': [[1, 4], [0, 1], [0, 1]], 'to': [[3, 4], [0, 1], [0, 1]], 'shift': [3, 0, 0]},
]

if __name__ == '__main__':
    outbox = Path(tempfile.mkdtemp(prefix='supercell-'))
    for name, cell in (('shard-0', [0, 0, 0]), ('shard-1', [1, 0, 0])):
        publish(outbox / f'{name}.json', build_supercell_shard(DATA, BASIS, LINKS, [cell]))
    received = [load(outbox / f'shard-{i}.json') for i in range(2)]
    opened = assemble_supercell_shards(received[:1])
    final = assemble_supercell_shards(received, finalize=True)
    publish(outbox / 'structure.json', final)
    print(json.dumps({'received': len(received), 'open_cells': len(opened['cells']),
                      'open_boundary': len(opened['out_links']) + len(opened['in_links']),
                      'final_cells': len(final['cells']), 'final_links': len(final['links'])},
                     sort_keys=True))
    print(outbox)
