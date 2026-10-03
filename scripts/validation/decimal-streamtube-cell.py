# SPDX-License-Identifier: GPL-2.0-or-later
"""Offline, independent 80-digit intrinsic-cell reference; no solver imports.

JSON numbers are first parsed as binary64 and then converted *exactly* to
Decimal. Thus this evaluates the actual supplied JS inputs, not their short
decimal print approximations. Derivatives use two fixed tiny central steps,
not the JS derivative formulas or an adaptively selected matching step.
"""
import json
import sys
from decimal import Decimal as D, localcontext


def exact(value):
    if isinstance(value, dict):
        return {k: exact(v) for k, v in value.items()}
    if isinstance(value, list):
        return [exact(v) for v in value]
    if isinstance(value, float):
        return D.from_float(value)
    if isinstance(value, int):
        return D(value)
    return value


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def cell(p):
    banks = [[(q['x'], q['y']) for q in p[side]] for side in ['lower', 'upper']]
    segments = [[sub(row[i + 1], row[i]) for i in range(2)] for row in banks]
    # Difference vectors first. Midpoint distances never subtract large
    # absolute coordinates; this is algebraically the published geometry.
    gaps = [sub(banks[1][i], banks[0][i]) for i in range(3)]
    areas = [mul(add(gaps[i], gaps[i + 1]), D('.5')) for i in range(2)]
    tangents = [mul(add(segments[0][i], segments[1][i]), D('.5')) for i in range(2)]
    tangents = [mul(t, 1 / dot(t, t).sqrt()) for t in tangents]
    sides = [mul(add(*row), D('.5')) for row in segments]
    s, n = mul(add(*sides), D('.5')), mul(add(*areas), D('.5'))
    volume = cross(s, n)
    assert volume > 0
    curvature = (cross(*segments[0]) - cross(*segments[1])) / (2 * volume)
    mass, h0, gamma = p['massFlow'], p['stagnationEnthalpy'], p.get('gamma', D('1.4'))
    states = []
    for rho, a, t in zip(p['densities'], areas, tangents):
        normal_area = cross(t, a)
        assert normal_area > 0
        q = mass / (rho * normal_area)
        h = h0 - q*q / 2
        assert h > 0
        states.append({'q': q, 'p': (gamma - 1) / gamma * rho * h,
                       'h': h, 'm2': q*q / ((gamma - 1) * h), 'rho': rho})
    a, b = states
    pm, m2 = (a['p'] + b['p']) / 2, (a['m2'] + b['m2']) / 2
    factor = p.get('pressureCorrectionFactor', D('.1'))
    pc = factor * gamma * pm * m2 * (1 - m2) * curvature if m2 < 1 else D(0)
    # Projected vector momentum and the auxiliary pressure sum determine
    # both streamline pressures. This code supplies values only.
    delta = (mass * (a['q'] * dot(tangents[0], n) - b['q'] * dot(tangents[1], n))
             + pc * cross(*areas)) / volume
    pl, pu = pm + pc - delta / 2, pm + pc + delta / 2
    momentum = (mass * (b['q'] * dot(tangents[1], s) - a['q'] * dot(tangents[0], s))
                - pc * cross(*sides)) / volume + b['p'] - a['p']
    entropy = (b['h'] / a['h']).ln() / (gamma - 1) - (b['rho'] / a['rho']).ln()
    qbar = (a['q'] + b['q']) / 2
    ue = qbar + D('.2') * qbar * m2 * (m2 - 1) * curvature if m2 < 1 else qbar
    return {'isentropicResidual': -pm * entropy, 'streamwiseResidual': momentum,
            'lowerPressure': pl, 'upperPressure': pu, 'pressureCorrection': pc,
            'edgeVelocity': ue, 'pressureCurvature': curvature, 'area': volume,
            'q0': a['q'], 'q1': b['q']}


def shifted(base, tangent, h):
    if isinstance(base, dict):
        return {k: shifted(v, tangent.get(k, {} if isinstance(v, dict) else [] if isinstance(v, list) else D(0)), h)
                for k, v in base.items()}
    if isinstance(base, list):
        return [shifted(v, tangent[i] if i < len(tangent) else {} if isinstance(v, dict) else D(0), h)
                for i, v in enumerate(base)]
    return base + h * tangent if isinstance(base, D) else base


def central(p, t, h):
    plus, minus = cell(shifted(p, t, h)), cell(shifted(p, t, -h))
    return {key: (plus[key] - minus[key]) / (2*h) for key in plus}


def serial(value):
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serial(v) for v in value]
    return str(value) if isinstance(value, D) else value


def run(request):
    with localcontext() as ctx:
        ctx.prec = 80
        cases = []
        for raw in request['cases']:
            p = exact(raw['parameters'])
            case = {'id': raw['id'], 'value': cell(p), 'directions': []}
            for direction in raw.get('directions', []):
                t = exact(direction['tangent'])
                case['directions'].append({'direction': direction['direction'],
                    'derivatives': [{'h': h, 'value': central(p, t, D(h))} for h in ['1e-20', '1e-25']],
                    'samples': [{'step': sample['step'], 'value': cell(exact(sample['parameters']))}
                                for sample in direction.get('samples', [])]})
            cases.append(case)
        return serial({'precision': ctx.prec, 'inputConversion': 'exact binary64', 'cases': cases})


if __name__ == '__main__':
    json.dump(run(json.load(sys.stdin)), sys.stdout, allow_nan=False)
    sys.stdout.write('\n')
