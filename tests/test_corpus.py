"""语料发行、同源证据与撤下审计的验收测试。"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from site_orbits.corpus import Corpus, CorpusError, FoldConflict

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'translation': [[0, 1], [0, 1], [0, 1]]}


def f(n, d=1):
    return [n, d]


def data(*sites):
    return {'operations': [IDENTITY], 'sites': list(sites)}


def component(value):
    return value if isinstance(value, list) else f(value)


def site(sid, species, x, occupancy, y=0, z=0):
    return {'id': sid, 'species': species,
            'point': [component(x), component(y), component(z)],
            'occupancy': component(occupancy)}


def sample_a():
    return data(site('ax', 'X', f(1, 4), f(1)),
                site('ay', 'Y', f(3, 4), f(1, 2)))


def sample_b():
    # 两点的整体平移 x+1/2：X 源 ID 改变，Y 拆成同点两条同物种 1/4
    return data(site('bx', 'X', f(3, 4), f(1)),
                site('by1', 'Y', f(1, 4), f(1, 4)),
                site('by2', 'Y', f(1, 4), f(1, 4)))


def sample_c():
    # 与 A 总组成相同，但 X–Y 间距不同
    return data(site('cx', 'X', f(1, 4), f(1)),
                site('cy', 'Y', f(1, 2), f(1, 2)))


def sample_d():
    # 与 A、B、C 均不等价
    return data(site('dx', 'X', f(1, 4), f(1)),
                site('dy', 'Y', f(1, 3), f(1)))


class CorpusCase(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / 'corpus.db'
        self.db = Corpus(self.path)
        self.db.apply_batch({'samples': [{'id': name, 'data': value}
                                        for name, value in (('A', sample_a()), ('B', sample_b()),
                                                            ('C', sample_c()), ('D', sample_d()))]})

    def tearDown(self):
        self.db.close()
        self.directory.cleanup()

    def reopen(self):
        self.db.close()
        self.db = Corpus(self.path)
        return self.db


class ReleaseTests(CorpusCase):
    def test_geometry_folds(self):
        receipt = self.db.release('R', {'A': 'train', 'C': 'test'})
        self.assertEqual(receipt['release_id'], 'R')
        self.assertEqual(receipt['evidence_revision'], 1)
        manifest = self.db.fetch('R')
        self.assertEqual(manifest['format'], 'corpus-manifest/1')
        self.assertEqual([(s['sample'], s['fold']) for s in manifest['samples']],
                         [('A', 'train'), ('C', 'test')])
        positions = {s['sample']: s['positions'] for s in manifest['samples']}
        self.assertEqual([c['site'] for c in positions['A'][0]['contributions']], ['ax'])
        with self.assertRaises(FoldConflict) as caught:
            self.db.release('bad', {'A': 'train', 'B': 'test'})
        conflict = caught.exception.conflicts[0]
        self.assertEqual((conflict['source'], conflict['target']), ('A', 'B'))
        self.assertEqual(conflict['steps'], [{'from': 'A', 'to': 'B', 'kind': 'geometry',
                                              'translation': ['1/2', '0', '0']}])
        # B 在同折时可发行，两份来源不被合并
        same = self.db.release('S', {'A': 'train', 'B': 'train'})
        self.assertEqual(same['sample_count'], 2)
        kept = {s['sample']: s['positions'] for s in self.db.fetch('S')['samples']}
        self.assertEqual(sorted(c['site'] for c in kept['A'][0]['contributions'] +
                                kept['A'][1]['contributions']), ['ax', 'ay'])
        self.assertEqual(sorted(c['site'] for c in kept['B'][0]['contributions'] +
                                kept['B'][1]['contributions']), ['bx', 'by1', 'by2'])

    def test_mapping_and_receipt_rules(self):
        receipt = self.db.release('R', [['C', 'test'], ['A', 'train']])
        shuffled = self.db.release('R', {'C': 'test', 'A': 'train'})
        self.assertEqual(receipt, shuffled)
        self.assertEqual(set(receipt), {'format', 'release_id', 'evidence_revision',
                                        'sample_count', 'fold_counts', 'manifest_digest'})
        with self.assertRaises(CorpusError):
            self.db.release('R', {'A': 'validation', 'C': 'test'})
        with self.assertRaises(CorpusError):
            self.db.release('R2', {'A': 'train', 'Z': 'test'})
        with self.assertRaises(CorpusError):
            self.db.release('R2', [['A', 'train'], ['A', 'test'], ['C', 'test']])
        with self.assertRaises(CorpusError):
            self.db.release('R2', {'A': 'training'})
        with self.assertRaises(CorpusError):
            self.db.release('R2', {})
        with self.assertRaises(CorpusError):
            self.db.fetch('R2')

    def test_failed_release_does_not_occupy_id(self):
        with self.assertRaises(FoldConflict):
            self.db.release('R', {'A': 'train', 'B': 'test'})
        receipt = self.db.release('R', {'A': 'train', 'C': 'test'})
        self.assertEqual(receipt['release_id'], 'R')
        self.assertEqual(self.db.fetch('R')['samples'][0]['sample'], 'A')
        self.assertEqual(self.db.audit('R')['released_revision'], 1)


class EvidenceTests(CorpusCase):
    def setUp(self):
        super().setUp()
        self.db.release('R', {'A': 'train', 'C': 'test'})
        self.db.release('S', {'A': 'train', 'B': 'train'})
        self.original_s = self.db.fetch('S')

    def test_batch_withdrawal_and_audit(self):
        result = self.db.apply_batch({'evidence': [
            {'id': 'E1', 'a': 'A', 'b': 'D'},
            {'id': 'E2', 'a': 'D', 'b': 'C'},
            {'id': 'E3', 'a': 'D', 'b': 'C'}]})
        self.assertEqual(result['revision'], 2)
        self.assertEqual(result['withdrawn'], ['R'])
        with self.assertRaises(CorpusError):
            self.db.fetch('R')
        self.assertEqual(self.db.fetch('S'), self.original_s)
        self.assertEqual(self.db.audit('S')['withdrawn_revision'], None)
        self.assertEqual(self.db.audit('S')['status'], 'available')
        audit = self.db.audit('R')
        self.assertEqual(audit['status'], 'withdrawn')
        self.assertEqual(audit['released_revision'], 1)
        self.assertEqual(audit['withdrawn_revision'], 2)
        steps = audit['conflicts'][0]['steps']
        self.assertEqual([(s['from'], s['to'], s['kind'], s.get('evidence')) for s in steps],
                         [('A', 'D', 'evidence', 'E1'), ('D', 'C', 'evidence', 'E2')])
        self.assertEqual(audit['manifest']['samples'][0]['sample'], 'A')
        with self.assertRaises(FoldConflict) as caught:
            self.db.release('T', {'A': 'train', 'C': 'test'})
        self.assertEqual(caught.exception.conflicts[0]['source'], 'A')

    def test_retraction_and_new_release(self):
        self.db.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                                          {'id': 'E2', 'a': 'D', 'b': 'C'},
                                          {'id': 'E3', 'a': 'D', 'b': 'C'}]})
        self.db.apply_batch({'retract': ['E2']})
        with self.assertRaises(FoldConflict):
            self.db.release('T', {'A': 'train', 'C': 'test'})
        # E1 仍有效，但唯一的 A–C 链已被切断：撤 E3 后 A 与 C 不连通
        self.db.apply_batch({'retract': ['E3']})
        receipt = self.db.release('T', {'A': 'train', 'C': 'test'})
        self.assertEqual(receipt['evidence_revision'], 4)
        with self.assertRaises(CorpusError):
            self.db.fetch('R')
        audit = self.db.audit('R')
        self.assertEqual(audit['withdrawn_revision'], 2)
        self.assertEqual([s['evidence'] for s in audit['conflicts'][0]['steps']], ['E1', 'E2'])
        # 原映射乱序重试仍返回首次回执；不同映射不能占用 R
        original = {'format': 'corpus-receipt/1', 'release_id': 'R', 'evidence_revision': 1,
                    'sample_count': 2, 'fold_counts': {'train': 1, 'test': 1, 'validation': 0},
                    'manifest_digest': self.db.audit('R')['receipt']['manifest_digest']}
        self.assertEqual(self.db.release('R', [['C', 'test'], ['A', 'train']]), original)
        with self.assertRaises(CorpusError):
            self.db.release('R', {'A': 'validation', 'C': 'test'})
        self.assertEqual(self.db.audit('R')['status'], 'withdrawn')

    def test_geometry_and_evidence_mixed_chain(self):
        # B 与 A 几何等价；A–D、D–C 同源证据使 A、B、C 同组。B 在 train、C 在 test
        # 必须跨折（链条经过未选中的 A 与 D），且证明链在几何步上带平移。
        self.db.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                                          {'id': 'E2', 'a': 'D', 'b': 'C'}]})
        with self.assertRaises(FoldConflict) as caught:
            self.db.release('T', {'B': 'train', 'C': 'test'})
        conflict = caught.exception.conflicts[0]
        self.assertEqual((conflict['source'], conflict['target']), ('B', 'C'))
        self.assertEqual([(s['kind'], s.get('translation'), s.get('evidence'))
                          for s in conflict['steps']],
                         [('geometry', ['1/2', '0', '0'], None),
                          ('evidence', None, 'E1'),
                          ('evidence', None, 'E2')])
        # 同组全在 validation 可发行
        receipt = self.db.release('U', {'B': 'validation', 'C': 'validation'})
        self.assertEqual(receipt['sample_count'], 2)

    def test_batch_rejected_atomically(self):
        before = self.db.audit('R')
        before_s = self.db.fetch('S')
        bad = [
            # 本会经 D 连接 A–C 从而撤下 R，但 E3 悬空，整批拒绝
            {'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                          {'id': 'E2', 'a': 'D', 'b': 'C'},
                          {'id': 'E3', 'a': 'C', 'b': 'Z'}]},
            {'retract': ['E0']},
            {'retract': ['E1'], 'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'}],
             },
            {'samples': [{'id': 'E1', 'data': sample_a()},
                         {'id': 'E1', 'data': sample_a()}]},
            {'samples': [{'id': 'F', 'data': sample_a()},
                         {'id': 'F', 'data': sample_b()}]},
            {'evidence': [{'id': 'E1', 'a': 'A', 'b': 'A'}]},
            {},
            # 超占位：同一展开点的两个 X 占位 3/4 相加超过 1
            {'samples': [{'id': 'Q', 'data': data(site('q', 'X', f(1, 2), f(3, 4)),
                                                  site('r', 'X', f(1, 2), f(3, 4)))}]},
            # 不闭合操作组
            {'samples': [{'id': 'Q', 'data': {'operations': [IDENTITY, {
                'id': 'half', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                'translation': [[1, 3], [0, 1], [0, 1]]}],
                'sites': [site('q', 'X', f(1, 4), f(1))]}}]},
        ]
        for batch in bad:
            with self.assertRaises((CorpusError, FoldConflict)):
                self.db.apply_batch(batch)
        self.assertEqual(self.db.evidence_revision, 1)
        self.assertEqual(self.db.audit('R'), before)
        self.assertEqual(self.db.fetch('S'), before_s)
        self.assertEqual(self.db.audit('S')['withdrawn_revision'], None)
        for eid in ('E1', 'E2', 'E3'):
            self.assertIsNone(self.db._conn.execute(
                'SELECT evidence_id FROM evidence WHERE evidence_id=?', (eid,)).fetchone())
        # 后续合法批次仍能使用被拒批次里的 ID
        self.db.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'}]})
        self.assertEqual(self.db.evidence_revision, 2)

    def test_oversized_batches_rejected(self):
        with self.assertRaises(CorpusError):
            self.db.apply_batch({'samples': [{'id': f'X{i}', 'data': sample_a()} for i in range(300)]})
        # 空批、重复 ID、悬空引用和撤回无效证据已由上一个用例覆盖；证据累计上限
        extra = {'samples': [{'id': 'Z', 'data': sample_d()}]}
        self.db.apply_batch(extra)
        self.db.apply_batch({'evidence': [{'id': f'E{i}', 'a': 'A', 'b': 'D'} for i in range(1999)]})
        with self.assertRaises(CorpusError):
            self.db.apply_batch({'evidence': [{'id': 'F1', 'a': 'A', 'b': 'D'},
                                              {'id': 'F2', 'a': 'A', 'b': 'D'}]})
        self.assertEqual(self.db.evidence_revision, 3)


class LimitTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.db = Corpus(Path(self.directory.name) / 'corpus.db')

    def tearDown(self):
        self.db.close()
        self.directory.cleanup()

    def test_position_and_release_limits(self):
        many = data(*[site(f'p{i}', 'X', f(i, 49), f(1, 49)) for i in range(49)])
        with self.assertRaises(CorpusError):
            self.db.apply_batch({'samples': [{'id': 'big', 'data': many}]})
        self.db.apply_batch({'samples': [{'id': 'A', 'data': sample_a()}]})
        for i in range(100):
            self.db.release(f'R{i}', {'A': 'train'})
        with self.assertRaises(CorpusError):
            self.db.release('R100', {'A': 'train'})
        self.assertEqual(self.db.audit('R0')['released_revision'], 1)


class InterruptTests(CorpusCase):
    def test_sample_batch_rolls_back(self):
        class Abort(BaseException):
            pass

        self.db._checkpoint = lambda stage: (_ for _ in ()).throw(Abort()) if stage == 'batch:before_commit' else None
        with self.assertRaises(Abort):
            self.db.apply_batch({'samples': [{'id': 'E', 'data': sample_a()},
                                             {'id': 'F', 'data': sample_d()}]})
        self.db._checkpoint = lambda stage: None
        self.db = self.reopen()
        self.assertEqual([r[0] for r in self.db._conn.execute('SELECT sample_id FROM samples')],
                         ['A', 'B', 'C', 'D'])
        self.assertEqual(self.db.evidence_revision, 1)
        self.db.apply_batch({'samples': [{'id': 'E', 'data': sample_a()}]})
        self.assertEqual([r[0] for r in self.db._conn.execute('SELECT sample_id FROM samples')],
                         ['A', 'B', 'C', 'D', 'E'])

    def test_interruption_boundaries(self):
        class Abort(BaseException):
            pass

        for boundary in ('batch:before_commit', 'batch:after_commit'):
            fresh = Corpus(Path(self.directory.name) / f'{boundary.replace(":", "-")}.db')
            fresh.apply_batch({'samples': [{'id': name, 'data': value}
                                           for name, value in (('A', sample_a()), ('B', sample_b()),
                                                               ('C', sample_c()), ('D', sample_d()))]})
            fresh.release('R', {'A': 'train', 'C': 'test'})
            fresh.release('S', {'A': 'train', 'B': 'train'})
            fresh._checkpoint = lambda stage, boundary=boundary: (
                (_ for _ in ()).throw(Abort()) if stage == boundary else None)
            with self.assertRaises(Abort):
                fresh.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                                                {'id': 'E2', 'a': 'D', 'b': 'C'},
                                                {'id': 'E3', 'a': 'D', 'b': 'C'},
                                                {'id': 'E4', 'a': 'C', 'b': 'B'}]})
            fresh._checkpoint = lambda stage: None
            fresh.close()
            reopened = Corpus(fresh.path)
            revision = reopened.evidence_revision
            if boundary == 'batch:before_commit':
                self.assertEqual(revision, 1)
                self.assertEqual([s['sample'] for s in reopened.fetch('R')['samples']], ['A', 'C'])
                self.assertEqual(reopened.audit('R')['status'], 'available')
                # 提交前中断后重试，只产生一个修订号
                reopened.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'},
                                                   {'id': 'E2', 'a': 'D', 'b': 'C'},
                                                   {'id': 'E3', 'a': 'D', 'b': 'C'},
                                                   {'id': 'E4', 'a': 'C', 'b': 'B'}]})
                revision = 2
            else:
                self.assertEqual(revision, 2)
                with self.assertRaises(CorpusError):
                    reopened.fetch('R')
                self.assertEqual([s['sample'] for s in reopened.fetch('S')['samples']], ['A', 'B'])
                self.assertEqual(len(reopened.audit('R')['conflicts'][0]['steps']), 2)
                # 提交后中断不重复提交证据
                with self.assertRaises(CorpusError):
                    reopened.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'D'}]})
            self.assertEqual(reopened.evidence_revision, revision)
            reopened.close()

    def test_release_response_loss(self):
        class Abort(BaseException):
            pass

        self.db._checkpoint = lambda stage: (_ for _ in ()).throw(Abort()) if stage == 'release:after_commit' else None
        with self.assertRaises(Abort):
            self.db.release('R', {'A': 'train', 'C': 'test'})
        self.db._checkpoint = lambda stage: None
        self.db = self.reopen()
        self.assertEqual(self.db.evidence_revision, 1)
        self.assertEqual(self.db.audit('R')['status'], 'available')
        receipt = self.db.release('R', {'C': 'test', 'A': 'train'})
        self.assertEqual(receipt['release_id'], 'R')
        manifests = [r for r in self.db._conn.execute('SELECT release_id FROM releases')]
        self.assertEqual(manifests, [('R',)])
        self.assertEqual(self.db.audit('R')['manifest']['format'], 'corpus-manifest/1')

    def test_release_before_commit(self):
        class Abort(BaseException):
            pass

        self.db._checkpoint = lambda stage: (_ for _ in ()).throw(Abort()) if stage == 'release:before_commit' else None
        with self.assertRaises(Abort):
            self.db.release('R', {'A': 'train', 'C': 'test'})
        self.db._checkpoint = lambda stage: None
        self.db = self.reopen()
        with self.assertRaises(CorpusError):
            self.db.audit('R')
        receipt = self.db.release('R', {'A': 'train', 'C': 'test'})
        self.assertEqual(receipt['evidence_revision'], 1)


DATA_A = {'operations': [IDENTITY],
          'sites': [{'id': 'ax', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 1]},
                    {'id': 'ay', 'species': 'Y', 'point': [[3, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 2]}]}
DATA_B = {'operations': [IDENTITY],
          'sites': [{'id': 'bx', 'species': 'X', 'point': [[3, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 1]},
                    {'id': 'by1', 'species': 'Y', 'point': [[1, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 4]},
                    {'id': 'by2', 'species': 'Y', 'point': [[1, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 4]}]}
DATA_C = {'operations': [IDENTITY],
          'sites': [{'id': 'cx', 'species': 'X', 'point': [[1, 4], [0, 1], [0, 1]],
                     'occupancy': [1, 1]},
                    {'id': 'cy', 'species': 'Y', 'point': [[1, 2], [0, 1], [0, 1]],
                     'occupancy': [1, 2]}]}


class ProcessKillTests(unittest.TestCase):
    """用真实子进程 os._exit 验证提交前后中断，可重复执行。"""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / 'corpus.db'

    def tearDown(self):
        self.directory.cleanup()

    def run_child(self, stage):
        package = Path(__import__('site_orbits').__file__).resolve().parent.parent
        code = f"""
import os
import sys
sys.path.insert(0, {str(package)!r})
from site_orbits import Corpus

DATA_A = {DATA_A!r}
DATA_B = {DATA_B!r}

with Corpus({str(self.path)!r}) as db:
    db.apply_batch({{'samples': [{{'id': 'A', 'data': DATA_A}},
                                  {{'id': 'B', 'data': DATA_B}}]}})
    db.release('R', {{'A': 'train'}})
    db._checkpoint = lambda stage: os._exit(9) if stage == {stage!r} else None
    db.release('S', {{'A': 'train', 'B': 'train'}})
"""
        return __import__('subprocess').run([__import__('sys').executable, '-c', code],
                                            capture_output=True, text=True)

    def test_kill_after_release_commit(self):
        result = self.run_child('release:after_commit')
        self.assertEqual(result.returncode, 9, result.stderr)
        reopened = Corpus(self.path)
        try:
            self.assertEqual(reopened.evidence_revision, 1)
            self.assertEqual(reopened.audit('R')['status'], 'available')
            self.assertEqual(reopened.audit('S')['status'], 'available')
            first = reopened.audit('S')['receipt']
            self.assertEqual(reopened.release('S', {'B': 'train', 'A': 'train'}), first)
            self.assertEqual([row[0] for row in reopened._conn.execute(
                'SELECT release_id FROM releases')], ['R', 'S'])
            self.assertEqual(reopened.audit('S')['manifest']['format'], 'corpus-manifest/1')
        finally:
            reopened.close()

    def test_kill_after_batch_commit_withdraws(self):
        package = Path(__import__('site_orbits').__file__).resolve().parent.parent
        code = f"""
import os
import sys
sys.path.insert(0, {str(package)!r})
from site_orbits import Corpus

DATA_A = {DATA_A!r}
DATA_B = {DATA_B!r}
DATA_C = {DATA_C!r}

with Corpus({str(self.path)!r}) as db:
    db.apply_batch({{'samples': [{{'id': 'A', 'data': DATA_A}},
                                  {{'id': 'B', 'data': DATA_B}},
                                  {{'id': 'C', 'data': DATA_C}}]}})
    db.release('R', {{'A': 'train', 'C': 'test'}})
    db.release('S', {{'A': 'train', 'B': 'train'}})
    db._checkpoint = lambda stage: os._exit(9) if stage == 'batch:after_commit' else None
    db.apply_batch({{'evidence': [{{'id': 'E1', 'a': 'A', 'b': 'C'}}]}})
"""
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 9, result.stderr)
        reopened = Corpus(self.path)
        try:
            self.assertEqual(reopened.evidence_revision, 2)
            with self.assertRaises(CorpusError):
                reopened.fetch('R')
            self.assertEqual([entry['sample'] for entry in reopened.fetch('S')['samples']], ['A', 'B'])
            audit = reopened.audit('R')
            self.assertEqual(audit['withdrawn_revision'], 2)
            self.assertEqual(audit['conflicts'][0]['steps'][0]['evidence'], 'E1')
        finally:
            reopened.close()

    def test_kill_before_batch_commit(self):
        package = Path(__import__('site_orbits').__file__).resolve().parent.parent
        code = f"""
import os
import sys
sys.path.insert(0, {str(package)!r})
from site_orbits import Corpus

DATA_A = {DATA_A!r}
DATA_B = {DATA_B!r}
DATA_C = {DATA_C!r}

with Corpus({str(self.path)!r}) as db:
    db.apply_batch({{'samples': [{{'id': 'A', 'data': DATA_A}},
                                  {{'id': 'B', 'data': DATA_B}},
                                  {{'id': 'C', 'data': DATA_C}}]}})
    db.release('R', {{'A': 'train', 'C': 'test'}})
    db._checkpoint = lambda stage: os._exit(9) if stage == 'batch:before_commit' else None
    db.apply_batch({{'evidence': [{{'id': 'E1', 'a': 'A', 'b': 'C'}}]}})
"""
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 9, result.stderr)
        reopened = Corpus(self.path)
        try:
            self.assertEqual(reopened.evidence_revision, 1)
            self.assertEqual([entry['sample'] for entry in reopened.fetch('R')['samples']], ['A', 'C'])
            self.assertEqual(reopened.audit('R')['status'], 'available')
            reopened.apply_batch({'evidence': [{'id': 'E1', 'a': 'A', 'b': 'C'}]})
            self.assertEqual(reopened.evidence_revision, 2)
            with self.assertRaises(CorpusError):
                reopened.fetch('R')
        finally:
            reopened.close()

    def test_kill_before_release_commit(self):
        result = self.run_child('release:before_commit')
        self.assertEqual(result.returncode, 9, result.stderr)
        fresh = Corpus(self.path)
        try:
            self.assertEqual([row[0] for row in fresh._conn.execute(
                'SELECT release_id FROM releases')], ['R'])
            with self.assertRaises(CorpusError):
                fresh.audit('S')
            receipt = fresh.release('S', {'A': 'train', 'B': 'train'})
            self.assertEqual(receipt['evidence_revision'], 1)
            self.assertEqual([row[0] for row in fresh._conn.execute(
                'SELECT release_id FROM releases')], ['R', 'S'])
        finally:
            fresh.close()


if __name__ == '__main__':
    unittest.main()
