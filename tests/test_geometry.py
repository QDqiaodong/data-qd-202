"""精确几何等价判定测试：不复用矩阵、不吞样本，也不接受近邻。"""
import unittest
from fractions import Fraction

from site_orbits.geometry import structure_key, translation
from site_orbits import process

IDENTITY = {'id': 'identity', 'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'translation': [[0, 1], [0, 1], [0, 1]]}


def f(n, d=1):
    return [n, d]


def positions(*sites):
    return process({'operations': [IDENTITY], 'sites': list(sites)})['positions']


def component(value):
    return value if isinstance(value, list) else f(value)


def point(x, y, z, species='X', occupancy=f(1), sid='p'):
    return {'id': sid, 'species': species,
            'point': [component(x), component(y), component(z)], 'occupancy': component(occupancy)}


class GeometryTests(unittest.TestCase):
    def test_pure_translation_equivalence(self):
        a = positions(point(f(1, 4), 0, 0, sid='a'), point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'))
        b = positions(point(f(3, 4), 0, 0, sid='a2'), point(f(1, 4), 0, 0, 'Y', f(1, 4), sid='b1'),
                      point(f(1, 4), 0, 0, 'Y', f(1, 4), sid='b2'))
        c = positions(point(f(1, 4), f(1, 2), 0, sid='c'), point(f(3, 4), f(1, 2), 0, 'Y', f(1, 2), sid='d'))
        ka, kb, kc = structure_key(a), structure_key(b), structure_key(c)
        self.assertEqual(ka['key'], kb['key'])
        self.assertEqual(ka['key'], kc['key'])
        self.assertEqual(translation(ka, kb), ['1/2', '0', '0'])
        self.assertEqual(translation(ka, kc), ['0', '1/2', '0'])

    def test_rotation_scale_and_composition_are_not_equivalence(self):
        a = positions(point(f(1, 4), 0, 0, sid='a'), point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'))
        # 同一组点但旋转到 y 轴不是三维分数平移
        rotated = positions(point(0, f(1, 4), 0, sid='a'), point(0, f(3, 4), 0, 'Y', f(1, 2), sid='b'))
        # 等间距缩放
        scaled = positions(point(f(1, 8), 0, 0, sid='a'), point(f(3, 8), 0, 0, 'Y', f(1, 2), sid='b'))
        # 总组成相同但间距不同
        composition = positions(point(f(1, 4), 0, 0, sid='a'), point(f(1, 2), 0, 0, 'Y', f(1, 2), sid='b'))
        # 同坐标不同物种
        species = positions(point(f(1, 4), 0, 0, 'Y', sid='a'), point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'))
        # 同坐标同物种不同占位
        occupancy = positions(point(f(1, 4), 0, 0, 'X', f(1, 2), sid='a'),
                              point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'))
        keys = {structure_key(a)['key']}
        for other in (rotated, scaled, composition, species, occupancy):
            self.assertNotIn(structure_key(other)['key'], keys)

    def test_translation_reproduces_target(self):
        from site_orbits.geometry import _labeled_positions

        def shifted(records, vector):
            return {(tuple((point[i] + vector[i]) % 1 for i in range(3)), label)
                    for point, label in _labeled_positions(records)}
        a = positions(point(f(1, 4), 0, 0, sid='a'), point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'))
        b = positions(point(f(3, 4), 0, 0, sid='a2'), point(f(1, 4), 0, 0, 'Y', f(1, 4), sid='b1'),
                      point(f(1, 4), 0, 0, 'Y', f(1, 4), sid='b2'))
        # 周期结构：同一份周期内多组平移都能复原
        periodic = positions(point(f(1, 4), 0, 0, sid='p1'), point(f(1, 4), f(1, 2), 0, sid='p2'))
        for source, target in ((a, b), (b, a), (periodic, periodic)):
            t = [Fraction(v) for v in translation(structure_key(source), structure_key(target))]
            self.assertEqual(shifted(source, t), shifted(target, [Fraction(0)] * 3))

    def test_permutation_invariance(self):
        sites = [point(f(1, 4), 0, 0, sid='a'), point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='b'),
                 point(f(1, 2), f(1, 2), 0, 'Y', f(1, 2), sid='c')]
        base = structure_key(positions(*sites))['key']
        for order in (sites[::-1], sites[1:] + sites[:1], sites[2:] + sites[:2]):
            self.assertEqual(structure_key(positions(*order))['key'], base)

    def test_float_tolerance_not_used(self):
        a = positions(point(0, 0, 0, sid='a'), point(f(1, 3), 0, 0, sid='b'))
        b = positions(point(0, 0, 0, sid='a'), point(f(33333, 100000), 0, 0, sid='b'))
        self.assertNotEqual(structure_key(a)['key'], structure_key(b)['key'])

    def test_pointwise_species_vector_must_match(self):
        # 同一点 X+Y 混合与另一结构 X、Y 分点不可等价
        mixed = positions(point(f(1, 4), 0, 0, 'X', f(1, 2), sid='m1'),
                          point(f(1, 4), 0, 0, 'Y', f(1, 2), sid='m2'), point(f(3, 4), 0, 0, sid='m3'))
        separate = positions(point(f(1, 4), 0, 0, 'X', f(1, 2), sid='s1'),
                             point(f(3, 4), 0, 0, 'Y', f(1, 2), sid='s2'),
                             point(f(3, 4), 0, 0, 'X', f(1, 2), sid='s3'))
        self.assertNotEqual(structure_key(mixed)['key'], structure_key(separate)['key'])


if __name__ == '__main__':
    unittest.main()
