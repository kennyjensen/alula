// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { directPolygonGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { assertConvexStreamtubeGrid } from '../src/geometry/streamtube-convex-step.js';
import { prepareConvexWakeGridUpdate } from '../src/euler/streamtube-grid-update.js';

function retained(name) {
  const f = JSON.parse(fs.readFileSync(new URL(`fixtures/${name}.json`, import.meta.url)));
  return createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
}
function calls(system) {
  const counters = { euler: 0, layers: 0 }, euler = system.euler.evaluate, layers = system.bl.evaluate;
  system.euler.evaluate = (...args) => { counters.euler++; return euler(...args); };
  system.bl.evaluate = (...args) => { counters.layers++; return layers(...args); };
  return { counters, reset() { counters.euler = counters.layers = 0; } };
}
function compare(system, x, options = {}) {
  const count = calls(system), before = x.slice(), active = system.bl.snapshotActive();
  const baseline = system.admissible(x, options), expectedCounts = { ...count.counters }; count.reset();
  let failure, failures = 0;
  const reported = system.admissible(x, { ...options, onFailure: d => { failure = d; failures++; } });
  assert.equal(reported, baseline); assert.deepEqual(count.counters, expectedCounts);
  assert.deepEqual(x, before); assert.deepEqual(system.bl.snapshotActive(), active);
  assert.equal(failures, baseline ? 0 : 1);
  return { accepted: baseline, failure, counters: expectedCounts, count };
}

test('retained roundoff shape failure reports the exact controlling halfspace before any Euler evaluation', () => {
  const system = retained('default-coupled-roundoff-shape'), x = system.initial;
  const r = compare(system, x); assert.equal(r.accepted, false);
  assert.deepEqual(r.counters, { euler: 0, layers: 0 });
  const d = r.failure, id = 51, k = system.ne + 4 * id;
  assert.equal(d.kind, 'bl-domain'); assert.equal(d.station.id, id);
  assert.deepEqual(d.station, system.bl.stations[id]);
  assert.deepEqual(d.failedConstraints, ['kinematic-shape']);
  const { mach, gamma } = system.bl.kernel.parameters;
  const mf = mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2), h = 1 - .5 * (gamma - 1) * mf * x[k + 3] ** 2;
  const shape = (x[k + 2] - x[k + 1]) * h - x[k + 1] * (.29 + .113) * mf * x[k + 3] ** 2;
  // This witness is deliberately at cancellation roundoff: an independent
  // multiplication association need not reproduce its last bits. The
  // reported value must be the exact original polynomial gate's value.
  assert.equal(d.shape, -5.095750210681871e-18);
  assert.ok(shape < 0 && d.shape <= 0);
  assert.ok(Math.abs(d.shape - shape) <= 8 * Number.EPSILON * Math.max(Math.abs((x[k + 2] - x[k + 1]) * h), Math.abs(x[k + 1] * (.29 + .113) * mf * x[k + 3] ** 2)));
  assert.equal(d.enthalpy, h);
  assert.equal(d.packed.wakeGap, 0);
  assert.equal(d.physicalLengthPerPackedThickness, system.bl.scale * system.euler.conditions.lengthScale);
});

test('finite wake diagnostics subtract the current physical dead-air gap and retain both failed thermal margins', () => {
  for (const mode of ['shape', 'thermal']) {
    const { system, x } = twoActiveFiniteBaseWakes(), id = system.bl.wakes[1].ids[0], k = system.ne + 4 * id;
    const geometry = system.bl.geometry(x.subarray(0, system.ne));
    const gap = geometry.coordinates[id].wakeGap / system.bl.scale;
    if (mode === 'shape') x[k + 2] = .5 * gap + x[k + 1];
    else {
      const { mach, gamma } = system.bl.kernel.parameters;
      const hstinv = (gamma - 1) * mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
      x[k + 3] = 1.01 * Math.sqrt(2 / hstinv);
    }
    assert.equal(system.bl.admissible(x.subarray(system.ne)), true, 'The old total-H precheck alone does not resolve this failure.');
    const r = compare(system, x); assert.equal(r.accepted, false); assert.deepEqual(r.counters, { euler: 0, layers: 0 });
    const d = r.failure; assert.equal(d.kind, 'bl-domain'); assert.equal(d.station.id, id); assert.equal(d.station.body, 1);
    assert.equal(d.station.kind, 'wake'); assert.equal(d.packed.wakeGap, gap);
    assert.ok(d.shape < 0);
    if (mode === 'shape') { assert.deepEqual(d.failedConstraints, ['kinematic-shape']); assert.ok(d.enthalpy > 0 && d.rawHk < 1); }
    else { assert.deepEqual(d.failedConstraints, ['kinematic-shape', 'static-enthalpy']); assert.ok(d.enthalpy < 0); assert.equal(d.rawHk, null); assert.equal(d.edgeMachSquared, null); }
  }
});

test('actual native amplification exception is retained without inventing an interval location', () => {
  const { system, x } = twoActiveFiniteBaseWakes(), s = system.bl.stations.find(s => s.regime === 'similarity');
  x[system.ne + 4 * s.id] = system.bl.kernel.parameters.ncrit;
  const r = compare(system, x); assert.equal(r.accepted, false);
  assert.deepEqual(r.counters, { euler: 1, layers: 1 });
  assert.equal(r.failure.kind, 'evaluation-error'); assert.equal(r.failure.stage, 'evaluation');
  assert.equal(r.failure.message, 'Natural transition is outside the fixed-trip research model.');
  assert.equal(r.failure.station, undefined); assert.equal(r.failure.firstCell, undefined);
});

test('retained convexity failure reports the actual first cell and still allows the original positive-simple iteration mode', () => {
  const system = retained('default-coupled-ises-wake-fold'), r = compare(system, system.initial);
  assert.equal(r.accepted, false); assert.deepEqual(r.counters, { euler: 1, layers: 1 });
  const d = r.failure; assert.equal(d.kind, 'grid-convexity'); assert.ok(d.invalidCellCount > 0);
  assert.ok(d.minCornerSine < d.minimumCornerSine);
  const { group: g, interval: i, tube: j, cell } = d.firstCell, { nx, tubes } = system.euler.layout;
  assert.equal(cell, tubes.slice(0, g).reduce((sum, n) => sum + nx * n, 0) + i * tubes[g] + j);
  const nodes = system.euler.decode(system.initial.subarray(0, system.ne)).nodes;
  assert.deepEqual(d.firstCell.points, [nodes[g][i][j], nodes[g][i + 1][j], nodes[g][i + 1][j + 1], nodes[g][i][j + 1]]);
  assert.ok(directPolygonGeometry(d.firstCell.points).minCornerSine <= d.minimumCornerSine);
  const permissive = compare(system, system.initial, { requireConvex: false });
  assert.equal(permissive.accepted, true); assert.equal(permissive.failure, undefined);
  assert.deepEqual(permissive.counters, { euler: 1, layers: 1 });
});

test('valid observers add no evaluation and failed observers cannot mutate state or change the boolean', () => {
  const { system, x } = twoActiveFiniteBaseWakes(), baseline = system.evaluate(x).residual;
  const valid = compare(system, x); assert.equal(valid.accepted, true); assert.deepEqual(valid.counters, { euler: 1, layers: 1 });
  assert.deepEqual(system.evaluate(x).residual, baseline);
  const s = system.bl.stations[0], k = system.ne + 4 * s.id; x[k + 2] = x[k + 1];
  const before = x.slice(), metadata = { ...s };
  assert.equal(system.admissible(x, { onFailure: d => { d.station.id = -1; d.packed.theta = -100; throw new Error('observer'); } }), false);
  assert.deepEqual(x, before); assert.deepEqual(system.bl.stations[0], metadata);
});

// Real ISES control source with explicit fake numerical boundaries, matching
// the repository's continuation-controller tests. This tests propagation,
// not Euler/BL physics or convergence, and performs no numerical solve.
let serial = 0;
async function routing(failInitial) {
  const nodes = Array.from({ length: 2 }, (_, g) => Array.from({ length: 3 }, (_, i) => [{ x: i, y: g }, { x: i, y: g + 1 }]));
  const outer = { nodes, undisplacedNodes: nodes, stagnation: [0], diagnostics: { maxMach: 0 } }, families = { euler: .1, boundaryLayer: 0, edgeMatching: 0 };
  const payload = { kind: 'bl-domain', stage: 'bl-domain', station: { id: 7, body: 0, side: 'lower', i: 2 },
    failedConstraints: ['kinematic-shape'], shape: -.0003, enthalpy: .9 };
  let constructors = 0, fakeLinearCalls = 0;
  const bindings = {
    assertConvexStreamtubeGrid, prepareConvexWakeGridUpdate,
    createStreamtubeStationOrdering: () => assert.fail('This diagnostic fixture does not request station ordering.'),
    solveCoupledLinearSystem: () => assert.fail('This diagnostic fixture retains explicit automatic ordering.'),
    requireCoupledResidualDecrease: () => assert.fail('This diagnostic fixture does not request Armijo acceptance.'),
    respondToCoupledProjectionGeometry: () => assert.fail('This diagnostic fixture does not request projection geometry.'),
    createCoupledStreamtubeBody() {
      const ordinal = ++constructors;
      return { ne: 1, n: 5, initial: new Float64Array(5),
        euler: { layout: { densityCount: 1, bodies: [{ leadingIndex: 1 }] }, decode: () => outer,
          setDisplacement() {}, adoptGeometry: x => x.slice() },
        bl: { thicknesses: () => ({}), snapshotActive: () => [] }, jacobian: () => null,
        admissible(_, { onFailure } = {}) { if (ordinal === (failInitial ? 2 : 3)) { onFailure?.(structuredClone(payload)); return false; } return true; },
        admissibleValue(x, options) { return this.admissible(x, options) ? this.evaluate(x) : null; },
        evaluate: () => ({ residual: Float64Array.of(.1, 0, 0, 0, 0), families, outer }) };
    },
    coupledStreamtubeTripEvents: () => ({ snapshot: () => [], restore() {}, prepare: () => ({ changed: false }) }),
    coupledStreamtubeResult: (_, x, details) => ({ x, ...details }),
    proposeCoupledDensityNewton: (_, x) => ({ x: x.slice(), step: 1 }),
    proposeCoupledXfoilBLUpdate: (_, x) => ({ x: x.slice(), step: 1 }),
    solveSparseDirect: () => { fakeLinearCalls++; return { x: new Float64Array(5), relativeResidual: 0, refinements: 0, attempts: [] }; },
    redistributeStreamtubeTangentially: group => ({ nodes: group, solution: { pairs: 5, relativeResidual: 0 }, maxDisplacement: 0 }),
    captureStreamtubeInletFractions: () => [[0, 1]],
    adjustStreamtubeInlets: nodes => ({ nodes }), dekinkStreamtubeInteriors: nodes => ({ nodes, repairs: [] }),
  };
  const key = `__admissibilityRouting${++serial}`; globalThis[key] = bindings;
  const proxy = 'data:text/javascript;base64,' + Buffer.from(`const h=globalThis[${JSON.stringify(key)}];\n`
    + Object.keys(bindings).map(k => `export const ${k}=h.${k};`).join('\n')).toString('base64');
  try {
    const sourceURL = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
    const source = fs.readFileSync(sourceURL, 'utf8').replace(/import \{([^}]+)\} from '([^']+)';/g, (_, names, specifier) => {
      const imported = names.split(',').map(name => name.trim());
      const mocked = imported.filter(name => Object.hasOwn(bindings, name));
      const actual = imported.filter(name => !Object.hasOwn(bindings, name));
      return [mocked.length ? `import { ${mocked.join(', ')} } from '${proxy}';` : '',
        actual.length ? `import { ${actual.join(', ')} } from '${new URL(specifier, sourceURL)}';` : ''].join('\n');
    });
    const module = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
    const result = module.solveCoupledStreamtubeIses({ bodies: [{ leadingIndex: 1 }] }, {
      linearOrdering: 'auto', maxIterations: failInitial ? 0 : 1 });
    return { result, payload, fakeLinearCalls };
  } finally { delete globalThis[key]; }
}

test('ISES retains the factual prefix/suffix and structured failure in initial-maintenance and trial rejection records', async () => {
  const first = await routing(true);
  assert.equal(first.fakeLinearCalls, 0);
  assert.deepEqual(first.result.initialRedistribution.admissibility, first.payload);
  assert.match(first.result.reason, /Maintained coupled state fails physical\/grid admissibility\. BL station 7/);
  const trial = await routing(false);
  assert.equal(trial.fakeLinearCalls, 1, 'Mock boundary, not an actual linear solve.');
  assert.deepEqual(trial.result.lastRejectedStep.admissibility, trial.payload);
  assert.deepEqual(trial.result.lastRejectedStep.rejections[0].admissibility, trial.payload);
  assert.match(trial.result.reason, /shape=-0.0003, enthalpy=0.9/);
  assert.equal(trial.result.history.length, 1, 'Rejected trial did not publish an accepted iterate.');
});
