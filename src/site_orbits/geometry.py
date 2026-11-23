"""展开位点的精确几何指纹。

本模块只回答“两份已经展开的位点集合是否只差一个三维分数平移”。按物种把同
一点的贡献占位相加得到可比较的物种占位向量；取全部精确成对差坐标的分母公倍
数，把相对坐标化成整数，再对每个候选锚点平移、排序并取字典序最小者作为指纹。
成对差对整体平移不变，且最小表示对参考点选择不变，所以指纹不依赖输入顺序。
全部算术语精确有理数，不使用浮点容差、旋转、换轴、缩放、晶格度量或连接关系。
"""
import json
from fractions import Fraction
from math import lcm


def _labeled_positions(positions):
    """把展开点整理为 (精确坐标三元组, 可比较的物种占位向量) 的列表。"""
    labeled = []
    for record in positions:
        point = tuple(Fraction(v) for v in record['point'])
        totals = {}
        for contribution in record['contributions']:
            species = contribution['species']
            totals[species] = totals.get(species, Fraction(0)) + Fraction(contribution['occupancy'])
        labeled.append((point, tuple((name, str(value)) for name, value in sorted(totals.items()))))
    return labeled


def _common_scale(relative):
    """相对坐标分母的最小公倍数；只差平移的位点集合得到同一个值。"""
    value = 1
    for point in relative:
        for component in point:
            value = lcm(value, component.denominator)
    return value


def structure_key(positions):
    """返回平移不变的结构指纹与参考点。

    ``key`` 为规范化字符串（公分母与整数余量），两份位点集合逐点相同（模 1、
    带物种占位标签）当且仅当 key 相等，且与输入顺序无关。``anchor`` 是某个样本
    点；若两样本 key 相等，则后者坐标等于前者加 ``anchor_b - anchor_a``（模 1）。
    """
    labeled = _labeled_positions(positions)
    points = [point for point, _ in labeled]
    labels = [label for _, label in labeled]
    reference = points[0]
    relative = [tuple((point[i] - reference[i]) % 1 for i in range(3)) for point in points]
    scale = _common_scale(relative)
    grid = [tuple(int(component * scale) % scale for component in point) for point in relative]
    candidates = []
    for index, origin in enumerate(grid):
        shifted = sorted((tuple((point[i] - origin[i]) % scale for i in range(3)), label)
                         for point, label in zip(grid, labels))
        candidates.append((json.dumps([scale, [[list(point), [list(item) for item in label]]
                                               for point, label in shifted]],
                                       ensure_ascii=False, separators=(',', ':')),
                           [str(value) for value in points[index]]))
    key, anchor = min(candidates, key=lambda entry: (entry[0], entry[1]))
    return {'key': key, 'anchor': anchor}


def translation(source, target):
    """返回精确三元平移 t，使 target 的位点集合等于 source 加上 t（模 1）。

    调用者只需保证两份指纹 key 相同；本函数直接使用各自参考点，不做近似匹配。
    """
    a = [Fraction(v) for v in source['anchor']]
    b = [Fraction(v) for v in target['anchor']]
    return [str((b[i] - a[i]) % 1) for i in range(3)]
