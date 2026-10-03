"""Independent high-precision scalar stencil; no solver/runtime imports."""
from decimal import Decimal as D, localcontext
import json
import sys


def add(a, b):
    return [a[k] + b[k] for k in range(2)]


def sub(a, b):
    return [a[k] - b[k] for k in range(2)]


def scale(a, s):
    return [v * s for v in a]


def mean(a, b):
    return scale(add(a, b), D('.5'))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def speed(cell, t, gamma, h0):
    nodes = [[add(p, scale(v, t)) for p, v in zip(row, vel)]
             for row, vel in zip(cell['nodes'], cell['nodeTangents'])]
    lower, upper = nodes
    left = [sub(lower[i + 1], lower[i]) for i in range(2)]
    right = [sub(upper[i + 1], upper[i]) for i in range(2)]
    tangents = [mean(a, b) for a, b in zip(left, right)]
    normals = [mean(sub(upper[i], lower[i]), sub(upper[i + 1], lower[i + 1])) for i in range(2)]
    area = cross(mean(tangents[0], tangents[1]), mean(normals[0], normals[1]))
    assert area > 0
    # Direct signed bank bends, independently of the runtime's relative-gap
    # rearrangement intended to reduce cancellation in double precision.
    curvature = (cross(left[0], left[1]) - cross(right[0], right[1])) / (2 * area)
    mass = cell['mass'] + t * cell['massTangent']
    qs, ms = [], []
    for i in range(2):
        length = sum(v * v for v in tangents[i]).sqrt()
        normal_area = cross(tangents[i], normals[i]) / length
        density = cell['densities'][i] * (t * cell['logDensityTangents'][i]).exp()
        q = mass / (density * normal_area)
        m = q * q / ((gamma - 1) * (h0 - q * q / 2))
        assert normal_area > 0 and q > 0 and 0 < m < 1
        qs.append(q)
        ms.append(m)
    q, m = sum(qs) / 2, sum(ms) / 2
    return q + D('.2') * q * m * (m - 1) * curvature


data = json.load(sys.stdin, parse_float=D, parse_int=D)
with localcontext() as ctx:
    ctx.prec = 80
    def residual(t):
        return data['ue'] + t * data['ueTangent'] - sum(
            speed(c, t, data['gamma'], data['h0']) for c in data['cells']) / len(data['cells'])
    base = residual(D(0))
    checks = []
    for h in [D('1e-8'), D('1e-10'), D('1e-12')]:
        a = (residual(h) - residual(-h)) / (2 * h)
        b = (residual(h / 2) - residual(-h / 2)) / h
        checks.append({'h': str(h), 'derivative': str((4 * b - a) / 3)})
    print(json.dumps({'precision': ctx.prec, 'residual': str(base), 'checks': checks}))
