"""离线语料发行示例；run from repository root.

同一坐标系下的 A、B 只差整体平移，B 还改了源 ID 并把同物种占位拆写。库内把它们
视为几何等价；A、B 跨折的发行被拒绝，同折发行保留各自来源。再登记采集方证据，
撤下已经发行的跨折清单并审计首次撤下依据。
"""
import json
import tempfile
from pathlib import Path

from site_orbits import Corpus, CorpusError, FoldConflict

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'translation': [[0, 1], [0, 1], [0, 1]]}


def data(*sites):
    return {'operations': [IDENTITY], 'sites': list(sites)}


def site(sid, species, numerator, denominator, occupancy):
    return {'id': sid, 'species': species,
            'point': [[numerator, denominator], [0, 1], [0, 1]],
            'occupancy': occupancy}


A = data(site('ax', 'X', 1, 4, [1, 1]), site('ay', 'Y', 3, 4, [1, 2]))
B = data(site('bx', 'X', 3, 4, [1, 1]),
         site('by1', 'Y', 1, 4, [1, 4]), site('by2', 'Y', 1, 4, [1, 4]))
C = data(site('cx', 'X', 1, 4, [1, 1]), site('cy', 'Y', 1, 2, [1, 2]))
D = data(site('dx', 'X', 1, 4, [1, 1]), site('dy', 'Y', 1, 3, [1, 1]))

if __name__ == '__main__':
    root = Path(tempfile.mkdtemp(prefix='corpus-'))
    path = root / 'corpus.db'
    with Corpus(path) as db:
        first = db.apply_batch({'samples': [{'id': name, 'data': value}
                                            for name, value in (('A', A), ('B', B), ('C', C), ('D', D))]})
        r = db.release('R', {'A': 'train', 'C': 'test'})
        s = db.release('S', {'A': 'train', 'B': 'train'})
        try:
            db.release('crossed', {'A': 'train', 'B': 'test'})
        except FoldConflict as error:
            crossed = error.conflicts[0]
        else:
            raise AssertionError('geometry-equivalent samples must not cross folds')
        batch = db.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                                             {'id': 'E2', 'a': 'D', 'b': 'C'},
                                             {'id': 'E3', 'a': 'D', 'b': 'C'}]})
        try:
            db.fetch('R')
        except CorpusError:
            withdrawn = True
        else:
            withdrawn = False
        audit = db.audit('R')
    reopened = Corpus(path)
    try:
        summary = {
            'first_revision': first['revision'],
            'r_released_at': r['evidence_revision'],
            's_released_at': s['evidence_revision'],
            'crossed_steps': [(step['kind'], step.get('translation'), step.get('evidence'))
                              for step in crossed['steps']],
            'batch_revision': batch['revision'],
            'batch_withdrew': batch['withdrawn'],
            'r_withdrawn': withdrawn,
            'audit_chain': [step.get('evidence') for step in audit['conflicts'][0]['steps']],
            'audit_original_manifest': [entry['sample'] for entry in audit['manifest']['samples']],
            's_still_available': [entry['sample'] for entry in reopened.fetch('S')['samples']],
        }
        print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    finally:
        reopened.close()
    print(root)
