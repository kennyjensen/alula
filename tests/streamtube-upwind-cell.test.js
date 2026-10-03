// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';

const parameters = {
  lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
  upper: [{ x: .1, y: .25 }, { x: .58, y: .3 }, { x: 1.18, y: .4 }],
  densities: [1.05, 1.02], massFlow: .27, stagnationEnthalpy: 8,
  gamma: 1.4, pressureCorrectionFactor: .1, transportSpeeds: [.88, 1.15],
};
const close = (a, b, tolerance = 2e-12, label = '') => assert.ok(
  Number.isFinite(a) && Number.isFinite(b)
    && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${label}: ${a} != ${b}`);
const add = (a, b) => ({ x: a.x + b.x, y: a.y + b.y });
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const scale = (a, v) => ({ x: a.x * v, y: a.y * v });
const midpoint = (a, b) => scale(add(a, b), .5);
const dot = (a, b) => a.x * b.x + a.y * b.y;
const cross = (a, b) => a.x * b.y - a.y * b.x;
const clockwise = a => ({ x: a.y, y: -a.x });
const unit = a => scale(a, 1 / Math.hypot(a.x, a.y));

// Independent outward face integration on the midpoint control polygon.
// Inertia transports m*qtilde*t; pressure uses the physical thermodynamics.
function vectorFlux(input, value) {
  const lower = [midpoint(input.lower[0], input.lower[1]), midpoint(input.lower[1], input.lower[2])];
  const upper = [midpoint(input.upper[0], input.upper[1]), midpoint(input.upper[1], input.upper[2])];
  const centers = input.lower.map((p, i) => midpoint(p, input.upper[i]));
  const tangents = [unit(sub(centers[1], centers[0])), unit(sub(centers[2], centers[1]))];
  const sections = lower.map((p, i) => sub(upper[i], p));
  const cuts = sections.map((section, i) => scale(add(
    scale(clockwise(section), value.states[i].p),
    scale(tangents[i], input.massFlow * input.transportSpeeds[i])), i === 0 ? -1 : 1));
  const banks = [scale(clockwise(sub(lower[1], lower[0])), value.interfacePressure.lower),
    scale(clockwise(sub(upper[1], upper[0])), -value.interfacePressure.upper)];
  return { cuts, banks, flux: [...cuts, ...banks].reduce(add, { x: 0, y: 0 }),
    transverse: midpoint(...sections), streamwise: midpoint(sub(lower[1], lower[0]), sub(upper[1], upper[0])) };
}

test('explicit transport speeds preserve independent vector momentum balance on curved skew cells', () => {
  for (const pressureCorrectionFactor of [0, .1, .2]) {
    const p = { ...parameters, pressureCorrectionFactor }, value = evaluateStreamtubeCell(p);
    const { flux, transverse, streamwise } = vectorFlux(p, value);
    const area = cross(streamwise, transverse);
    close(flux.x, value.streamwiseResidual * transverse.y);
    close(flux.y, -value.streamwiseResidual * transverse.x);
    close(dot(flux, transverse) / area, 0);
    close(dot(flux, streamwise) / area, value.streamwiseResidual);
    close(value.interfacePressure.lower + value.interfacePressure.upper,
      value.states[0].p + value.states[1].p + 2 * value.pressureCorrection);
    const physical = evaluateStreamtubeCell({ ...p, transportSpeeds: undefined });
    // Unequal transport corrections must affect normal as well as S momentum.
    assert.ok(Math.abs(value.interfacePressure.lower - physical.interfacePressure.lower) > 1e-4);
    assert.ok(Math.abs(value.streamwiseResidual - physical.streamwiseResidual) > 1e-3);
  }
});

test('artificial entropy uses transport enthalpy and physical pressure without changing physical states', () => {
  const p = parameters, value = evaluateStreamtubeCell(p);
  const physical = evaluateStreamtubeCell({ ...p, transportSpeeds: undefined });
  for (const key of ['states', 'geometry', 'pressureCorrection', 'entropyJump'])
    assert.deepEqual(value[key], physical[key], key);
  assert.deepEqual(value.transportSpeeds, p.transportSpeeds);
  const h = p.transportSpeeds.map(q => p.stagnationEnthalpy - .5 * q * q);
  assert.deepEqual(value.artificialEnthalpies, h);
  const [a, b] = value.states;
  const entropy = p.gamma / (p.gamma - 1) * Math.log(h[1] / h[0]) - Math.log(b.p / a.p);
  close(value.artificialEntropyJump, entropy);
  close(value.isentropicResidual, -.5 * (a.p + b.p) * entropy);
  const invalidDensitySubstitution = Math.log(h[1] / h[0]) / (p.gamma - 1) - Math.log(b.rho / a.rho);
  assert.ok(Math.abs(entropy - invalidDensitySubstitution) > .01,
    'This fixture must distinguish physical pressure from the invalid artificial ideal-gas substitution.');
});

test('the same shared transport speed cancels both momentum flux components at an adjacent curved-cell face', () => {
  const lower = [...parameters.lower, { x: 1.7, y: .19 }];
  const upper = [...parameters.upper, { x: 1.77, y: .59 }];
  const density = [1.05, 1.02, .99], transport = [.88, 1.15, 1.09];
  const pieces = [0, 1].map(i => {
    const input = { ...parameters, lower: lower.slice(i, i + 3), upper: upper.slice(i, i + 3),
      densities: density.slice(i, i + 2), transportSpeeds: transport.slice(i, i + 2) };
    const value = evaluateStreamtubeCell(input);
    return { value, ...vectorFlux(input, value) };
  });
  const [a, b] = pieces;
  assert.deepEqual(a.value.states[1], b.value.states[0]);
  for (const key of ['x', 'y']) close(a.cuts[1][key] + b.cuts[0][key], 0);
  const exterior = [a.cuts[0], b.cuts[1], ...a.banks, ...b.banks].reduce(add, { x: 0, y: 0 });
  for (const key of ['x', 'y']) close(a.flux[key] + b.flux[key], exterior[key]);
  close(exterior.x, a.value.streamwiseResidual * a.transverse.y + b.value.streamwiseResidual * b.transverse.y);
  close(exterior.y, -a.value.streamwiseResidual * a.transverse.x - b.value.streamwiseResidual * b.transverse.x);
});

const outputs = value => {
  const fields = { streamwiseResidual: value.streamwiseResidual, isentropicResidual: value.isentropicResidual,
    entropyJump: value.entropyJump, pressureCorrection: value.pressureCorrection, interfacePressure: value.interfacePressure,
    states: value.states.map(({ rho, q, p, enthalpy, machSquared }) => ({ rho, q, p, enthalpy, machSquared })),
    geometry: value.geometry };
  for (const key of ['transportSpeeds', 'artificialEnthalpies', 'artificialEntropyJump'])
    if (key in value) fields[key] = value[key];
  const result = new Map();
  const visit = (v, name) => {
    if (typeof v === 'number') result.set(name, v);
    else for (const key of Object.keys(v).sort()) visit(v[key], `${name}.${key}`);
  };
  visit(fields, '');
  return result;
};
const shift = (base, tangent, step) => {
  const p = structuredClone(base);
  for (const side of ['lower', 'upper']) if (tangent[side]) for (let i = 0; i < 3; i++)
    for (const key of ['x', 'y']) p[side][i][key] += step * tangent[side][i][key];
  for (const key of ['densities', 'transportSpeeds']) if (tangent[key])
    p[key] = p[key].map((v, i) => v + step * tangent[key][i]);
  for (const key of ['massFlow', 'stagnationEnthalpy', 'pressureCorrectionFactor'])
    if (tangent[key]) p[key] += step * tangent[key];
  return p;
};
const tangents = () => {
  const list = [{ transportSpeeds: [1, 0] }, { transportSpeeds: [0, 1] }];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const row = Array.from({ length: 3 }, () => ({ x: 0, y: 0 })); row[i][key] = 1;
    list.push({ [side]: row });
  }
  list.push({ densities: [1, 0] }, { densities: [0, 1] }, { massFlow: 1 },
    { stagnationEnthalpy: 1 }, { pressureCorrectionFactor: 1 });
  list.push({ transportSpeeds: [.17, -.12], densities: [.06, -.08], massFlow: .025,
    stagnationEnthalpy: .3, pressureCorrectionFactor: .02,
    lower: [{ x: .1, y: .07 }, { x: -.12, y: .02 }, { x: .03, y: -.05 }] });
  return list;
};

test('upwind cell derivatives match fourth-order differences in both smooth pressure-correction branches', () => {
  for (const massFlow of [.27, .65]) {
    const p = { ...parameters, massFlow }, linear = linearizeStreamtubeCell(p);
    assert.deepEqual(linear.value, evaluateStreamtubeCell(p));
    assert.equal(linear.value.states.reduce((sum, v) => sum + v.machSquared, 0) < 2, massFlow === .27);
    for (const tangent of tangents()) {
      const expected = outputs(linear.apply(tangent));
      for (const h of [2e-4, 1e-4]) {
        const samples = [-2, -1, 1, 2].map(n => outputs(evaluateStreamtubeCell(shift(p, tangent, n * h))));
        assert.deepEqual([...expected.keys()], [...samples[0].keys()]);
        for (const [key, value] of expected) {
          const fd = (samples[0].get(key) - 8 * samples[1].get(key) + 8 * samples[2].get(key) - samples[3].get(key)) / (12 * h);
          close(value, fd, 3e-9, `${key}, h=${h}, direction=${JSON.stringify(tangent)}`);
        }
      }
    }
  }
});

test('unfiltered transport recovers physical residuals and its chained physical-speed derivative', () => {
  const { transportSpeeds: ignored, ...p } = parameters;
  const ordinary = linearizeStreamtubeCell(p), q = ordinary.value.states.map(s => s.q);
  const explicit = linearizeStreamtubeCell({ ...p, transportSpeeds: q });
  assert.deepEqual(evaluateStreamtubeCell({ ...p, transportSpeeds: undefined }), ordinary.value);
  for (const key of ['transportSpeeds', 'artificialEnthalpies', 'artificialEntropyJump'])
    assert.equal(Object.hasOwn(ordinary.value, key), false, 'The default result shape must remain unchanged.');
  for (const key of ['states', 'geometry', 'pressureCorrection', 'interfacePressure', 'streamwiseResidual', 'entropyJump'])
    assert.deepEqual(explicit.value[key], ordinary.value[key], key);
  close(explicit.value.isentropicResidual, ordinary.value.isentropicResidual);
  for (const tangent of tangents().filter(t => !t.transportSpeeds)) {
    const physical = ordinary.apply(tangent);
    const filtered = explicit.apply({ ...tangent, transportSpeeds: physical.states.map(s => s.q) });
    const originalFields = outputs(physical), filteredFields = outputs(filtered);
    for (const [key, value] of originalFields) close(value, filteredFields.get(key), 2e-12, key);
  }
});

test('supplied transport speeds stay fixed under omitted tangents and absent options reject a transport perturbation', () => {
  const explicit = linearizeStreamtubeCell(parameters);
  const mass = explicit.apply({ massFlow: 1 });
  assert.deepEqual(mass.transportSpeeds, [0, 0]);
  assert.ok(mass.states.every(state => state.q > 0));
  assert.deepEqual(mass.artificialEnthalpies, [0, 0]);
  const enthalpy = explicit.apply({ stagnationEnthalpy: 1 });
  assert.deepEqual(enthalpy.artificialEnthalpies, [1, 1]);
  for (const direction of [[1, 0], [0, 1]]) {
    const value = explicit.apply({ transportSpeeds: direction });
    value.states.forEach(state => Object.values(state).forEach(v => close(v, 0)));
    close(value.entropyJump, 0);
  }
  const { transportSpeeds: ignored, ...p } = parameters;
  const ordinary = linearizeStreamtubeCell(p);
  assert.throws(() => ordinary.apply({ transportSpeeds: [1, 0] }));
  assert.deepEqual(ordinary.apply({ transportSpeeds: [0, 0] }), ordinary.apply());
  for (const value of [[NaN, 0], [1], [1, 2, 3]]) assert.throws(() => explicit.apply({ transportSpeeds: value }));
});
