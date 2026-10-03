import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';

const parameters = { lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
  upper: [{ x: .1, y: .25 }, { x: .58, y: .3 }, { x: 1.18, y: .4 }],
  densities: [1.05, 1.02], massFlow: .27, stagnationEnthalpy: 8, gamma: 1.4, pressureCorrectionFactor: .1 };
const close = (a, b, tolerance = 2e-8) => assert.ok(Math.abs(a - b) < tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const flatten = value => {
  const result = [];
  const visit = v => { if (typeof v === 'number') result.push(v); else for (const key of Object.keys(v).sort()) if (key !== 'derivatives') visit(v[key]); };
  // Fix a common named-key order recursively, excluding section partials.
  // Geometry value and tangent objects need not share insertion order.
  visit({ streamwiseResidual: value.streamwiseResidual, entropyJump: value.entropyJump, isentropicResidual: value.isentropicResidual, interfacePressure: value.interfacePressure,
    pressureCorrection: value.pressureCorrection, states: value.states.map(s => ({ rho: s.rho, q: s.q, p: s.p, enthalpy: s.enthalpy, machSquared: s.machSquared })), geometry: value.geometry });
  return result;
};
const shift = (base, tangent, step) => {
  const p = structuredClone(base);
  for (const side of ['lower', 'upper']) if (tangent[side]) for (let i = 0; i < 3; i++)
    for (const key of ['x', 'y']) p[side][i][key] += step * tangent[side][i][key];
  if (tangent.densities) p.densities = p.densities.map((v, i) => v + step * tangent.densities[i]);
  for (const key of ['massFlow', 'stagnationEnthalpy', 'pressureCorrectionFactor']) if (tangent[key]) p[key] += step * tangent[key];
  return p;
};

test('intrinsic cell chain rule differentiates every coordinate, density, mass, enthalpy and auxiliary coefficient', () => {
  const tangents = [];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const row = Array.from({ length: 3 }, () => ({ x: 0, y: 0 })); row[i][key] = 1;
    tangents.push({ [side]: row });
  }
  tangents.push({ densities: [1, 0] }, { densities: [0, 1] }, { massFlow: 1 }, { stagnationEnthalpy: 1 }, { pressureCorrectionFactor: 1 });
  // Subsonic and supersonic local states exercise both smooth Pc branches.
  // This remains a derivative check, not a shock-capturing demonstration.
  for (const m of [.27, .65]) {
    const p = { ...parameters, massFlow: m }, linearization = linearizeStreamtubeCell(p);
    assert.deepEqual(linearization.value, evaluateStreamtubeCell(p));
    for (const tangent of tangents) {
      const step = 1e-6, analytic = flatten(linearization.apply(tangent));
      const plus = flatten(evaluateStreamtubeCell(shift(p, tangent, step)));
      const minus = flatten(evaluateStreamtubeCell(shift(p, tangent, -step)));
      assert.equal(analytic.length, plus.length);
      analytic.forEach((v, i) => close(v, (plus[i] - minus[i]) / (2 * step)));
    }
  }
});

test('intrinsic geometry derivatives preserve translation, rotation and similarity identities', () => {
  const { apply, value } = linearizeStreamtubeCell(parameters);
  const translation = apply({ lower: parameters.lower.map(() => ({ x: .3, y: -.2 })), upper: parameters.upper.map(() => ({ x: .3, y: -.2 })) });
  flatten(translation).forEach(v => close(v, 0, 1e-12));
  const rotate = p => ({ x: -p.y, y: p.x });
  const rotation = apply({ lower: parameters.lower.map(rotate), upper: parameters.upper.map(rotate) });
  close(rotation.streamwiseResidual, 0, 1e-12); close(rotation.pressureCorrection, 0, 1e-12);
  Object.values(rotation.interfacePressure).forEach(v => close(v, 0, 1e-12));
  const similarity = apply({ lower: parameters.lower, upper: parameters.upper, massFlow: parameters.massFlow });
  close(similarity.streamwiseResidual, 0, 1e-12); close(similarity.pressureCorrection, 0, 1e-12);
  Object.values(similarity.interfacePressure).forEach(v => close(v, 0, 1e-12));
  similarity.states.forEach(s => Object.values(s).forEach(v => close(v, 0, 1e-12)));
  close(similarity.geometry.area, 2 * value.geometry.area, 1e-12);
  similarity.geometry.normalAreas.forEach((v, i) => close(v, value.geometry.normalAreas[i], 1e-12));
});

test('cell linearization is linear in simultaneous perturbations and rejects malformed tangents', () => {
  const { apply } = linearizeStreamtubeCell(parameters);
  const a = { densities: [.1, -.07], massFlow: .03 }, b = { stagnationEnthalpy: .2, pressureCorrectionFactor: .01 };
  const da = flatten(apply(a)), db = flatten(apply(b)), sum = flatten(apply({ ...a, ...b }));
  sum.forEach((v, i) => close(v, da[i] + db[i], 1e-12));
  assert.throws(() => apply({ densities: [NaN, 0] }), /Invalid/);
  assert.throws(() => apply({ lower: [{ x: 0, y: 0 }] }), /Invalid/);
});

test('reused geometry scratch never aliases retained derivatives or caller buffers', () => {
  const p = { ...parameters, transportSpeeds: [.9, .95] };
  const { apply } = linearizeStreamtubeCell(p);
  const tangent = { lower: parameters.lower.map(p => ({ ...p })), upper: parameters.upper.map(p => ({ ...p })),
    densities: [.1, -.07], transportSpeeds: [.02, -.01] };
  const first = apply(tangent), snapshot = structuredClone(first);
  tangent.lower[0].x += .1; tangent.densities[0] += .3; tangent.transportSpeeds[0] += .2;
  const next = apply(tangent);
  assert.deepEqual(first, snapshot);
  assert.deepEqual(next, linearizeStreamtubeCell(p).apply(tangent));
  assert.throws(() => apply({ densities: [Infinity, 0] }), /Invalid/);
  assert.deepEqual(apply(tangent), next);
});
