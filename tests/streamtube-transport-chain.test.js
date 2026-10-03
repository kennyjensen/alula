// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { prepareStreamtubeTransportChain as prepare } from '../src/euler/streamtube-transport-chain.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const x = [0, .2, .6, 1.1, 1.5, 2];
const base = { lower: x.map(x => ({ x, y: .02 * x * x })),
  upper: x.map(x => ({ x: x + .08, y: .35 + .03 * x + .02 * x * x })),
  densities: [.9, 1, .85, .82, .91], massFlow: .42, stagnationEnthalpy: 2.7, gamma: 1.4,
  upwind: { mucon: .75, mcrit: .9, boundary: { kind: 'unfiltered-first-two' } } };
const sectionKeys = ['rho', 'q', 'p', 'enthalpy', 'machSquared'];
const observe = p => {
  const v = prepare(p);
  return v.sections.flatMap((s, k) => [...sectionKeys.map(key => s[key]), v.sectionGeometry[k].length,
    v.sectionGeometry[k].normalArea, v.transportSpeeds[k]]);
};
const perturb = (p, d, h) => ({ ...p,
  ...Object.fromEntries(['lower', 'upper'].map(side => [side, p[side].map((v, i) => ({ x: v.x + h * d[side][i].x, y: v.y + h * d[side][i].y }))])),
  densities: p.densities.map((v, i) => v + h * d.densities[i]), massFlow: p.massFlow + h * d.massFlow,
  stagnationEnthalpy: p.stagnationEnthalpy + h * d.stagnationEnthalpy });
function tangents(chain, d) {
  return chain.sections.map((_, k) => chain.sectionTangent(k, { lower: d.lower.slice(k, k + 2), upper: d.upper.slice(k, k + 2),
    density: d.densities[k], massFlow: d.massFlow, stagnationEnthalpy: d.stagnationEnthalpy }));
}

test('nonuniform affine speeds retain exact second-order cancellation across the full tube', () => {
  const arc = x.slice(1).map((v, i) => .5 * (x[i] + v)), q = arc.map(s => 1 + .2 * s), width = .4;
  const p = { ...base, lower: x.map(x => ({ x, y: 0 })), upper: x.map(x => ({ x, y: width })),
    densities: q.map(v => base.massFlow / (width * v)) };
  const value = prepare(p);
  value.sections.forEach((s, k) => { close(s.q, q[k]); close(value.transportSpeeds[k], q[k]); close(value.sectionArc[k], arc[k]); });
  assert.deepEqual(value.upwind.filtered, [false, false, true, true, true]);
  assert.ok(value.upwind.coefficients.slice(2).every(v => v > .01), 'Cancellation must occur with an active filter.');
});

test('complete local 2D section/filter derivatives match fourth-order differences for both upwind orders', t => {
  let maximum = 0, comparisons = 0;
  for (const mucon of [-.75, .75]) {
    const p = { ...base, upwind: { ...base.upwind, mucon } }, chain = prepare(p, { linearize: true });
    assert.ok(chain.sections.some(s => s.machSquared > 1), 'Include physical supersonic sections.');
    for (let seed = 0; seed < 6; seed++) {
      const d = { lower: x.map((_, i) => ({ x: .03 * Math.sin(i + seed), y: .025 * Math.cos(2 * i + seed) })),
        upper: x.map((_, i) => ({ x: .02 * Math.cos(3 * i + seed), y: .03 * Math.sin(2 * i - seed) })),
        densities: p.densities.map((_, i) => .04 * Math.sin(i * 2 + seed)),
        massFlow: .03 * Math.cos(seed), stagnationEnthalpy: .12 * Math.sin(seed + .4) };
      const ds = tangents(chain, d), exact = ds.flatMap((s, k) => [...sectionKeys.map(key => s[key]), s.length, s.normalArea,
        chain.transportTangent(k, i => ds[i])]);
      for (const h of [1e-4, 5e-5]) {
        const values = [2, 1, -1, -2].map(a => observe(perturb(p, d, a * h)));
        exact.forEach((v, i) => {
          const difference = (-values[0][i] + 8 * values[1][i] - 8 * values[2][i] + values[3][i]) / (12 * h);
          maximum = Math.max(maximum, Math.abs(v - difference) / Math.max(1, Math.abs(v), Math.abs(difference)));
          close(v, difference, 2e-9); comparisons++;
        });
      }
    }
  }
  t.diagnostic(JSON.stringify({ comparisons, maximumNormalizedError: maximum }));
});

test('rigid rotation and length-unit changes preserve physical gas and shared momentum speeds', () => {
  const original = prepare(base), angle = .73, c = Math.cos(angle), s = Math.sin(angle);
  for (const scale of [.001, 3, 1000]) {
    const p = { ...base, massFlow: base.massFlow * scale };
    for (const side of ['lower', 'upper']) p[side] = base[side].map(v => ({
      x: scale * (c * v.x - s * v.y + .7), y: scale * (s * v.x + c * v.y - .3) }));
    const value = prepare(p, { linearize: true });
    value.sections.forEach((v, k) => {
      for (const key of sectionKeys) close(v[key], original.sections[k][key]);
      close(value.transportSpeeds[k], original.transportSpeeds[k]);
      close(value.sectionArc[k], scale * original.sectionArc[k]);
    });
    const d = { lower: p.lower, upper: p.upper, massFlow: p.massFlow,
      densities: p.densities.map(() => 0), stagnationEnthalpy: 0 };
    const ds = tangents(value, d);
    ds.forEach((v, k) => { close(v.q, 0); close(v.machSquared, 0); close(value.transportTangent(k, i => ds[i]), 0); });
  }
});

test('local queries use at most three sections, retain snapshots, and reject invalid geometry/controls', () => {
  const p = structuredClone(base), original = structuredClone(p), chain = prepare(p, { linearize: true });
  assert.deepEqual(p, original);
  const before = chain.sectionTangent(2, { density: 1 });
  p.lower[2].x += 10; p.densities[2] = 99; p.upwind.mucon = 99;
  assert.deepEqual(chain.sectionTangent(2, { density: 1 }), before);
  for (let k = 0; k < 5; k++) {
    const visited = [];
    close(chain.transportTangent(k, i => { visited.push(i); return { q: 0, machSquared: 0, length: 0 }; }), 0);
    assert.deepEqual(visited, k < 2 ? [k] : [k - 2, k - 1, k]);
  }
  assert.throws(() => prepare({ ...base, lower: base.upper, upper: base.lower }), /Folded/);
  assert.throws(() => prepare({ ...base, densities: base.densities.map(() => .01) }), /enthalpy/);
  assert.throws(() => prepare({ ...base, upwind: undefined }), /explicit/);
  assert.throws(() => prepare({ ...base, upwind: { ...base.upwind, mcrit: 1.01 } }), /controls/);
  assert.throws(() => chain.sectionTangent(5), /index/);
  assert.throws(() => chain.sectionTangent(0, { massFlow: NaN }), /tangent/);
});

test('thermal failure identifies the physical section before any speed filtering', () => {
  const p = structuredClone(base); p.densities[3] = .01;
  const before = structuredClone(p);
  assert.throws(() => prepare(p), error => {
    assert.equal(error.code, 'streamtube-static-enthalpy');
    const d = error.diagnostics;
    assert.deepEqual(d.section, { i: 3, bankStations: [3, 4] });
    assert.equal(d.rho, .01); assert.equal(d.massFlow, p.massFlow);
    close(d.q, d.massFlow / (d.rho * d.normalArea));
    close(d.enthalpy, p.stagnationEnthalpy - .5 * d.q ** 2);
    assert.ok(d.enthalpy < 0);
    return true;
  });
  assert.deepEqual(p, before);
});
