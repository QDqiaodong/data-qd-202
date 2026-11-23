"""本地语料发行：样本登记、同源证据批次、分折发行与审计。

本模块只做“把已经由 :func:`site_orbits.process` 校验并展开的样本收进本地
语料库、按证据连通性约束候选分折、发行不可变清单”这一件事。几何等价使用
精确有理数平移判定，不读取晶格度量、不做结构预测、也不自动调整调用方给出
的折映射。

持久化使用单个标准库 SQLite 数据库文件。证据修订号、证据行、发行回执、清单
和撤下记录都在同一事务中提交：进程在提交前退出时，重新打开只看到旧状态；在
提交后退出时，重新打开看到完整的新状态。单写者即可，不需要网络或并发服务。
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
from collections import deque
from pathlib import Path

from .checks import fields, ident
from .geometry import structure_key, translation
from .service import process

FORMAT_MANIFEST = 'corpus-manifest/1'
FORMAT_RECEIPT = 'corpus-receipt/1'
FOLDS = ('train', 'validation', 'test')
MAX_SAMPLES = 300
MAX_POSITIONS_PER_SAMPLE = 48
MAX_EVIDENCE_IDS = 2000
MAX_RELEASES = 100
SCHEMA_VERSION = '1'

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,
    canonical TEXT NOT NULL,
    positions TEXT NOT NULL,
    added_revision INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS evidence (
    evidence_id TEXT PRIMARY KEY,
    sample_a TEXT NOT NULL,
    sample_b TEXT NOT NULL,
    active INTEGER NOT NULL CHECK (active IN (0, 1)),
    created_revision INTEGER NOT NULL,
    retracted_revision INTEGER,
    FOREIGN KEY (sample_a) REFERENCES samples(sample_id),
    FOREIGN KEY (sample_b) REFERENCES samples(sample_id)
);
CREATE TABLE IF NOT EXISTS releases (
    release_id TEXT PRIMARY KEY,
    mapping TEXT NOT NULL,
    receipt TEXT NOT NULL,
    released_revision INTEGER NOT NULL,
    manifest TEXT NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('available', 'withdrawn')),
    withdrawn_revision INTEGER,
    conflicts TEXT
);
"""


class CorpusError(ValueError):
    """语料库输入、状态或持久化错误。"""


class FoldConflict(CorpusError):
    """候选分折在同一连通组内跨折。

    ``conflicts`` 为至少一组内容：``source``/``target`` 两端样本及一条
    ``steps`` 证据链。每一步含 ``from``、``to``、``kind``（``geometry`` 或
    ``evidence``），几何步还带 ``translation``，同源步还带 ``evidence``。
    """

    def __init__(self, conflicts):
        self.conflicts = conflicts
        first = conflicts[0]
        super().__init__(
            f'release would place connected samples {first["source"]!r} and '
            f'{first["target"]!r} in different folds')


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def _read_mapping(mapping):
    """规范化样本到折的映射，忽略输入顺序，返回 (canonical_text, dict)。"""
    if isinstance(mapping, dict):
        items = list(mapping.items())
    elif isinstance(mapping, (list, tuple)):
        items = []
        for entry in mapping:
            if isinstance(entry, dict):
                fields(entry, ('sample', 'fold'))
                items.append((entry['sample'], entry['fold']))
            elif isinstance(entry, (list, tuple)) and len(entry) == 2:
                items.append((entry[0], entry[1]))
            else:
                raise CorpusError('mapping entries must be sample/fold pairs')
    else:
        raise CorpusError('mapping must be a dict or a list of sample/fold pairs')
    if not items:
        raise CorpusError('non-empty fold mapping required')
    result, seen = {}, set()
    for sample, fold in items:
        sid = ident(sample)
        if not isinstance(fold, str) or fold not in FOLDS:
            raise CorpusError('fold must be one of train, validation, test')
        if sid in seen:
            raise CorpusError('duplicate sample in fold mapping')
        seen.add(sid)
        result[sid] = fold
    return _json([[sid, result[sid]] for sid in sorted(result)]), result


class _UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}
        self.rank = {value: 0 for value in values}

    def find(self, value):
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != root:
            self.parent[value], value = root, self.parent[value]
        return root

    def union(self, left, right):
        a, b = self.find(left), self.find(right)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def _edge_graph(samples, evidence):
    """构建全库连通图：几何等价类为完全边，当前有效同源证据为边。"""
    adjacency = {sid: [] for sid in samples}
    groups = {}
    for sid, canonical in samples.items():
        groups.setdefault(canonical['key'], []).append(sid)
    union = _UnionFind(samples)
    for group in groups.values():
        ordered = sorted(group)
        for position, left in enumerate(ordered):
            for right in ordered[position + 1:]:
                union.union(left, right)
                forward = translation(samples[left], samples[right])
                backward = translation(samples[right], samples[left])
                adjacency[left].append({'to': right, 'kind': 'geometry', 'translation': forward})
                adjacency[right].append({'to': left, 'kind': 'geometry', 'translation': backward})
    for row in evidence:
        left, right, evidence_id = row
        union.union(left, right)
        adjacency[left].append({'to': right, 'kind': 'evidence', 'evidence': evidence_id})
        adjacency[right].append({'to': left, 'kind': 'evidence', 'evidence': evidence_id})
    for sid in adjacency:
        adjacency[sid].sort(key=lambda step: (step['to'], 0 if step['kind'] == 'evidence' else 1,
                                              step.get('evidence') or ''))
    return adjacency, union


def _path(adjacency, source, target):
    parent = {source: (None, None)}
    queue = deque([source])
    while queue:
        node = queue.popleft()
        if node == target:
            break
        for step in adjacency[node]:
            other = step['to']
            if other not in parent:
                parent[other] = (node, step)
                queue.append(other)
    if target not in parent:
        raise CorpusError('no evidence path between conflicting samples')
    steps = []
    node = target
    while node != source:
        previous, step = parent[node]
        steps.append({'from': previous, 'to': node,
                      'kind': step['kind'],
                      **({'translation': step['translation']} if step['kind'] == 'geometry'
                         else {'evidence': step['evidence']})})
        node = previous
    steps.reverse()
    return steps


def _conflicts_for(mapping, adjacency, union):
    """返回候选分折的冲突（可能为空列表）；每组含端点和完整证据链。"""
    folds_by_root, members_by_root = {}, {}
    for sid, fold in mapping.items():
        root = union.find(sid)
        folds_by_root.setdefault(root, set()).add(fold)
        members_by_root.setdefault(root, []).append(sid)
    conflicting = [root for root in folds_by_root if len(folds_by_root[root]) > 1]
    if not conflicting:
        return []
    root = min(conflicting, key=lambda item: min(members_by_root[item]))
    members = sorted(members_by_root[root])
    source = members[0]
    target = next(sid for sid in members if mapping[sid] != mapping[source])
    return [{'source': source, 'target': target, 'steps': _path(adjacency, source, target)}]


class Corpus:
    """本地语料库；调用方持有路径，可在进程退出后重新打开同一路径。"""

    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.path), isolation_level=None, timeout=30)
        self._conn.execute('PRAGMA foreign_keys=ON')
        self._conn.execute('PRAGMA synchronous=FULL')
        self._conn.executescript(_SCHEMA)
        self._conn.execute("INSERT OR IGNORE INTO meta(key, value) VALUES('evidence_revision', '0')")
        self._conn.execute("INSERT OR IGNORE INTO meta(key, value) VALUES('schema_version', ?)",
                           (SCHEMA_VERSION,))
        version = self._conn.execute(
            "SELECT value FROM meta WHERE key='schema_version'").fetchone()[0]
        if version != SCHEMA_VERSION:
            self.close()
            raise CorpusError('unsupported corpus schema version')

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def _checkpoint(self, stage):
        """测试用中断边界钩子；生产路径为空操作。"""

    @property
    def evidence_revision(self):
        row = self._conn.execute("SELECT value FROM meta WHERE key='evidence_revision'").fetchone()
        return int(row[0])

    def _load_samples(self):
        samples = {}
        for sid, canonical, positions in self._conn.execute(
                'SELECT sample_id, canonical, positions FROM samples'):
            fingerprint = json.loads(canonical)
            samples[sid] = {'key': fingerprint['key'], 'anchor': fingerprint['anchor'],
                            'positions': json.loads(positions)}
        return samples

    def _active_evidence(self):
        return [row for row in self._conn.execute(
            'SELECT sample_a, sample_b, evidence_id FROM evidence WHERE active=1 '
            'ORDER BY evidence_id')]

    def apply_batch(self, batch):
        """一起登记新样本、新同源证据并撤回批前有效证据。

        成功时生成下一个证据修订号；任一项错误则整个批次不生效。返回含
        ``revision`` 与受影响撤下发行 ID 的字典。
        """
        try:
            return self._apply_batch(batch)
        except CorpusError:
            raise
        except ValueError as error:
            raise CorpusError(str(error)) from error

    def _apply_batch(self, batch):
        if not isinstance(batch, dict):
            raise CorpusError('batch must be an object')
        unknown = set(batch) - {'samples', 'evidence', 'retract'}
        if unknown:
            raise CorpusError(f'unknown batch fields: {sorted(unknown)}')
        new_samples = batch.get('samples', [])
        new_evidence = batch.get('evidence', [])
        retract = batch.get('retract', [])
        for name, value in (('samples', new_samples), ('evidence', new_evidence), ('retract', retract)):
            if not isinstance(value, list):
                raise CorpusError(f'{name} must be a list')
        if not new_samples and not new_evidence and not retract:
            raise CorpusError('empty batch')
        prepared, sample_ids = {}, set()
        for raw in new_samples:
            fields(raw, ('id', 'data'))
            sid = ident(raw['id'])
            if sid in sample_ids:
                raise CorpusError('duplicate sample id in batch')
            sample_ids.add(sid)
            try:
                result = process(raw['data'])
            except ValueError as error:
                raise CorpusError(f'sample {sid}: {error}') from error
            positions = result['positions']
            if len(positions) > MAX_POSITIONS_PER_SAMPLE:
                raise CorpusError('sample exceeds 48 expanded positions')
            prepared[sid] = {'canonical': structure_key(positions),
                             'positions': positions}
        evidence_rows, evidence_ids = {}, set()
        for raw in new_evidence:
            fields(raw, ('id', 'a', 'b'))
            eid = ident(raw['id'])
            a, b = ident(raw['a']), ident(raw['b'])
            if a == b:
                raise CorpusError('evidence endpoints must be different samples')
            if eid in evidence_ids:
                raise CorpusError('duplicate evidence id in batch')
            evidence_ids.add(eid)
            evidence_rows[eid] = (a, b) if a < b else (b, a)
        retract_ids, retract_list = set(), []
        for raw in retract:
            eid = ident(raw)
            if eid in retract_ids:
                raise CorpusError('duplicate evidence retraction in batch')
            retract_ids.add(eid)
            retract_list.append(eid)
        if retract_ids & evidence_ids:
            raise CorpusError('cannot retract evidence created in the same batch')

        self._conn.execute('BEGIN IMMEDIATE')
        try:
            revision = int(self._conn.execute(
                "SELECT value FROM meta WHERE key='evidence_revision'").fetchone()[0]) + 1
            existing_count = self._conn.execute('SELECT COUNT(*) FROM samples').fetchone()[0]
            if existing_count + len(prepared) > MAX_SAMPLES:
                raise CorpusError('corpus exceeds 300 samples')
            existing_ids = {row[0] for row in self._conn.execute('SELECT sample_id FROM samples')}
            if existing_ids & sample_ids:
                raise CorpusError('sample id already registered')
            evidence_count = self._conn.execute('SELECT COUNT(*) FROM evidence').fetchone()[0]
            if evidence_count + len(evidence_rows) > MAX_EVIDENCE_IDS:
                raise CorpusError('corpus exceeds 2000 evidence ids')
            used_ids = {row[0] for row in self._conn.execute('SELECT evidence_id FROM evidence')}
            if used_ids & evidence_ids:
                raise CorpusError('evidence id already used')
            known = existing_ids | sample_ids
            for eid, (a, b) in evidence_rows.items():
                if a not in known or b not in known:
                    raise CorpusError('evidence references unknown sample')
            for eid in retract_list:
                row = self._conn.execute(
                    'SELECT sample_a, sample_b, active FROM evidence WHERE evidence_id=?',
                    (eid,)).fetchone()
                if row is None or row[2] != 1:
                    raise CorpusError('retraction references inactive or unknown evidence')
            for sid in sorted(prepared):
                item = prepared[sid]
                self._conn.execute(
                    'INSERT INTO samples(sample_id, canonical, positions, added_revision) '
                    'VALUES(?,?,?,?)',
                    (sid, _json(item['canonical']), _json(item['positions']), revision))
            for eid in sorted(evidence_rows):
                a, b = evidence_rows[eid]
                self._conn.execute(
                    'INSERT INTO evidence(evidence_id, sample_a, sample_b, active, '
                    'created_revision, retracted_revision) VALUES(?,?,?,1,?,NULL)',
                    (eid, a, b, revision))
            for eid in sorted(retract_ids):
                self._conn.execute(
                    'UPDATE evidence SET active=0, retracted_revision=? WHERE evidence_id=?',
                    (revision, eid))
            self._conn.execute(
                "UPDATE meta SET value=? WHERE key='evidence_revision'", (str(revision),))
            withdrawn = self._withdraw_conflicts(revision)
            self._checkpoint('batch:before_commit')
            self._conn.execute('COMMIT')
        except BaseException:
            if self._conn.in_transaction:
                self._conn.execute('ROLLBACK')
            raise
        self._checkpoint('batch:after_commit')
        return {'revision': revision,
                'samples': sorted(prepared),
                'evidence': sorted(evidence_rows),
                'retracted': sorted(retract_ids),
                'withdrawn': withdrawn}

    def _withdraw_conflicts(self, revision):
        """在同一事务内检查仍可取用的清单，首次跨折即撤下。"""
        available = list(self._conn.execute(
            "SELECT release_id, mapping FROM releases WHERE status='available' "
            'ORDER BY release_id'))
        if not available:
            return []
        samples = self._load_samples()
        adjacency, union = _edge_graph(samples, self._active_evidence())
        withdrawn = []
        for release_id, mapping_text in available:
            mapping = {sid: fold for sid, fold in json.loads(mapping_text)}
            conflicts = _conflicts_for(mapping, adjacency, union)
            if conflicts:
                self._conn.execute(
                    'UPDATE releases SET status=\'withdrawn\', withdrawn_revision=?, conflicts=? '
                    'WHERE release_id=?',
                    (revision, _json(conflicts), release_id))
                withdrawn.append(release_id)
        return sorted(withdrawn)

    def release(self, release_id, mapping):
        """发行固定清单与回执；同 ID 同映射重试直接返回原回执。"""
        try:
            return self._release(release_id, mapping)
        except CorpusError:
            raise
        except ValueError as error:
            raise CorpusError(str(error)) from error

    def _release(self, release_id, mapping):
        rid = ident(release_id)
        mapping_text, normalized = _read_mapping(mapping)
        row = self._conn.execute(
            'SELECT mapping, receipt, status FROM releases WHERE release_id=?', (rid,)).fetchone()
        if row is not None:
            if row[0] != mapping_text:
                raise CorpusError('release id already used with a different mapping')
            return json.loads(row[1])
        self._conn.execute('BEGIN IMMEDIATE')
        try:
            row = self._conn.execute(
                'SELECT mapping, receipt FROM releases WHERE release_id=?', (rid,)).fetchone()
            if row is not None:
                if row[0] != mapping_text:
                    raise CorpusError('release id already used with a different mapping')
                self._conn.execute('COMMIT')
                return json.loads(row[1])
            count = self._conn.execute('SELECT COUNT(*) FROM releases').fetchone()[0]
            if count >= MAX_RELEASES:
                raise CorpusError('corpus exceeds 100 successful releases')
            samples = self._load_samples()
            unknown = sorted(set(normalized) - set(samples))
            if unknown:
                raise CorpusError(f'release references unknown sample: {unknown[0]}')
            adjacency, union = _edge_graph(samples, self._active_evidence())
            conflicts = _conflicts_for(normalized, adjacency, union)
            if conflicts:
                raise FoldConflict(conflicts)
            revision = int(self._conn.execute(
                "SELECT value FROM meta WHERE key='evidence_revision'").fetchone()[0])
            manifest = {
                'format': FORMAT_MANIFEST,
                'release_id': rid,
                'evidence_revision': revision,
                'samples': [{'sample': sid, 'fold': normalized[sid],
                             'positions': samples[sid]['positions']} for sid in sorted(normalized)],
            }
            manifest_text = _json(manifest)
            digest = 'sha256:' + hashlib.sha256(manifest_text.encode('utf-8')).hexdigest()
            receipt = {
                'format': FORMAT_RECEIPT,
                'release_id': rid,
                'evidence_revision': revision,
                'sample_count': len(normalized),
                'fold_counts': {fold: sum(1 for value in normalized.values() if value == fold)
                                for fold in FOLDS},
                'manifest_digest': digest,
            }
            receipt_text = _json(receipt)
            self._conn.execute(
                'INSERT INTO releases(release_id, mapping, receipt, released_revision, '
                'manifest, status, withdrawn_revision, conflicts) '
                "VALUES(?,?,?,?,?,'available',NULL,NULL)",
                (rid, mapping_text, receipt_text, revision, manifest_text))
            self._checkpoint('release:before_commit')
            self._conn.execute('COMMIT')
        except BaseException:
            if self._conn.in_transaction:
                self._conn.execute('ROLLBACK')
            raise
        self._checkpoint('release:after_commit')
        return json.loads(receipt_text)

    def fetch(self, release_id):
        """按发行 ID 取用当前清单；已撤下的清单不可用。"""
        try:
            return self._fetch(release_id)
        except CorpusError:
            raise
        except ValueError as error:
            raise CorpusError(str(error)) from error

    def _fetch(self, release_id):
        rid = ident(release_id)
        row = self._conn.execute(
            'SELECT manifest, status FROM releases WHERE release_id=?', (rid,)).fetchone()
        if row is None:
            raise CorpusError('unknown release id')
        if row[1] != 'available':
            raise CorpusError(f'release {rid} is withdrawn')
        return json.loads(row[0])

    def audit(self, release_id):
        """审计入口：始终返回原清单、发行修订号和首次撤下依据。"""
        try:
            return self._audit(release_id)
        except CorpusError:
            raise
        except ValueError as error:
            raise CorpusError(str(error)) from error

    def _audit(self, release_id):
        rid = ident(release_id)
        row = self._conn.execute(
            'SELECT receipt, manifest, released_revision, status, withdrawn_revision, conflicts '
            'FROM releases WHERE release_id=?', (rid,)).fetchone()
        if row is None:
            raise CorpusError('unknown release id')
        return {'release_id': rid,
                'status': row[3],
                'receipt': json.loads(row[0]),
                'manifest': json.loads(row[1]),
                'released_revision': row[2],
                'withdrawn_revision': row[4],
                'conflicts': json.loads(row[5]) if row[5] is not None else None}
