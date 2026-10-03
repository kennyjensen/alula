import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directChannelConservation } from './oracles/streamtube.js';

const maximum = values => Math.max(...Array.from(values, Math.abs));
const build = file => {
  const f = JSON.parse(fs.readFileSync(new URL(`fixtures/${file}`, import.meta.url)));
  return createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
};

test('inactive coupled constraints recover the ordinary dogleg step at the retained isentropic gradient stall', () => {
  const solve = projectedSteps => solveCoupledStreamtubeBody(build('streamtube-isentropic-gradient-stall.json'), {
    stepMethod: 'dogleg', maxIterations: 1, initialTrustRadius: .25, projectedSteps });
  const plain = solve(false), constrained = solve(true);
  assert.equal(constrained.history.length, 2); assert.equal(plain.history.length, 2);
  const h = constrained.history[1]; assert.ok(h.trialRadius < .25); assert.equal(h.stepKind, 'dogleg');
  assert.ok(h.actualReduction > 1e-8); assert.ok(h.reductionRatio > .5);
  assert.deepEqual(constrained.x, plain.x); assert.deepEqual(constrained.history, plain.history);
  assert.equal(constrained.converged, false); assert.equal(constrained.mesh.quality.valid, true);
});

for (const [name, file, n] of [['coarse', 'streamtube-isentropic-coupled-root.json', 3185], ['refined', 'streamtube-isentropic-refined-root.json', 7919],
  ['nested', 'streamtube-isentropic-nested-root.json', 11929],
  ['stagnation-scaled', 'streamtube-isentropic-stagnation-root.json', 23097],
  ['section-speed', 'streamtube-section-velocity-root.json', 3185],
  ['section-speed refined', 'streamtube-section-velocity-refined-root.json', 23097],
  ['selective wall refinement', 'streamtube-selective-wall-root.json', 25191]])
test(`the ${name} isentropic coupled root closes mass, enthalpy and total pressure with reciprocal Jacobian coupling`, t => {
  const system = build(file), x = system.initial, v = system.evaluate(x), j = system.jacobian(x);
  assert.equal(system.n, n); assert.ok(system.admissible(x)); assert.ok(maximum(v.residual) < 1e-10);
  const { gamma, pInf, mach } = system.euler.conditions;
  const p0Inf = pInf * (1 + .5 * (gamma - 1) * mach ** 2) ** (gamma / (gamma - 1));
  for (const s of v.outer.sections.flat(2)) {
    const m2 = s.q ** 2 / (gamma * s.p / s.rho);
    assert.ok(Math.abs(s.p * (1 + .5 * (gamma - 1) * m2) ** (gamma / (gamma - 1)) / p0Inf - 1) < 1e-10);
  }
  const momentum = [];
  v.outer.nodes.forEach((nodes, g) => {
    const c = directChannelConservation({ nodes, sections: v.outer.sections.map(row => row[g]), cells: v.outer.cells.map(row => row[g]) }, gamma);
    for (const key of ['maxLocal', 'total']) for (const index of [0, 3]) assert.ok(Math.abs(c[key][index]) < 2e-9);
    assert.ok(maximum(c.internalCancellation) < 2e-9);
    momentum.push({ local: c.maxLocal.slice(1, 3), total: c.total.slice(1, 3) });
  });
  // Isentropy replaces a discrete momentum row; explicitly retain this
  // defect rather than claiming exact vector-momentum conservation.
  assert.ok(maximum(momentum.flatMap(m => m.local)) > 1e-6);
  const errors = [];
  for (const mode of ['euler', 'bl', 'both']) {
    const d = x.map((value, i) => i < system.ne ? (mode === 'bl' ? 0 : .001 * Math.sin(i * .73 + .2))
      : mode === 'euler' ? 0 : Math.max(.01, Math.abs(value)) * Math.sin(i * .43 + .4));
    const h = 2e-6, exact = sparseProduct(j, d), p = system.residual(x.map((value, i) => value + h * d[i])), m = system.residual(x.map((value, i) => value - h * d[i]));
    const error = maximum(exact.map((value, i) => Math.abs((p[i] - m[i]) / (2 * h) - value) / Math.max(1, Math.abs(value), Math.abs((p[i] - m[i]) / (2 * h)))));
    assert.ok(error < 5e-6); errors.push({ mode, error });
    if (mode === 'bl') assert.ok(maximum(exact.slice(0, system.ne)) > 1e-6);
    if (mode === 'euler') assert.ok(maximum(exact.slice(system.ne)) > 1e-6);
  }
  t.diagnostic(JSON.stringify({ residual: maximum(v.residual), errors, momentum }));
});
