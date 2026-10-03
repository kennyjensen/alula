// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { directUpwindChannelConservation } from './oracles/upwind-streamtube.js';

const upwind = { mucon: 1, mcrit: .5, boundary: { kind: 'unfiltered-first-two' } };
const fixture = (elements = 1) => ({ ...intrinsicBodyFixture({ elements, bodySegments: 4, tubes: 2,
  mach: .5, alpha: .25 }), streamwiseMode: 'momentum', upwind: structuredClone(upwind) });
const close = (a, b, tolerance = 2e-12, label = '') => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${label}: ${a} != ${b}`);
const fourth = (values, h) => (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h);
const max = a => Math.max(...Array.from(a, Math.abs));
function stateFor(system) {
  const { layout } = system;
  const state = system.initial.map((v, col) => v + (col < layout.densityCount ? .002 : .000005) * Math.sin(.7 * (col + 1)));
  for (const [key, value] of [['circulation', .03], ['source', .005], ['doubletX', .001], ['doubletY', -.002]])
    state[layout.globals[key]] = value;
  layout.globals.stagnation.forEach((col, b) => { if (col !== null) state[col] = .00004 * (b + 1); });
  layout.globals.capture.forEach(col => { if (col !== null) state[col] = .0001; });
  return state;
}
function compareDirection(system, state, jacobian, direction, h, statistics) {
  const samples = [-2, -1, 1, 2].map(k => system.residual(state.map((v, col) => v + k * h * direction[col])));
  const product = sparseProduct(jacobian, direction);
  for (let row = 0; row < system.layout.n; row++) {
    const numerical = fourth(samples.map(r => r[row]), h), exact = product[row];
    const error = Math.abs(exact - numerical) / Math.max(1, Math.abs(exact), Math.abs(numerical));
    statistics.comparisons++;
    if (error > statistics.maximumError) Object.assign(statistics,
      { maximumError: error, row, family: system.layout.rows[row].kind, exact, numerical, h });
    close(exact, numerical, 4e-7, `row ${row}, ${system.layout.rows[row].kind}, h=${h}`);
  }
}

test('upwind body Jacobian matches every fourth-order column including global constraints', t => {
  const system = createStreamtubeBodySystem(fixture()), state = stateFor(system), { n } = system.layout;
  assert.ok(n < 200, 'Keep the complete-column fixture small.');
  const dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense);
  const statistics = { comparisons: 0, maximumError: 0, unknowns: n };
  for (let col = 0; col < n; col++) {
    const direction = new Float64Array(n); direction[col] = 1;
    compareDirection(system, state, sparse, direction, 1e-6, statistics);
  }
  system.jacobian(state.map(v => .99 * v), { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense, 'Later assembly must not overwrite the retained sparse matrix.');
  t.diagnostic(JSON.stringify(statistics));
});

test('two-body capture and stagnation derivatives survive normal-chart rebase', t => {
  const system = createStreamtubeBodySystem(fixture(2));
  let state = stateFor(system);
  const { layout } = system, statistics = { comparisons: 0, maximumError: 0, unknowns: layout.n };
  const globals = [...new Set(Object.values(layout.globals).flat().filter(col => col !== null))];
  assert.ok(layout.globals.capture.some(col => col !== null));
  for (let chart = 0; chart < 2; chart++) {
    const sparse = system.jacobian(state, { sparse: true });
    assert.deepEqual(sparseDense(sparse), system.jacobian(state));
    const directions = globals.map(col => { const v = new Float64Array(layout.n); v[col] = 1; return v; });
    for (const phase of [.2, .9]) directions.push(state.map((_, col) => .1 * Math.cos(.43 * col + phase)));
    for (const direction of directions) for (const h of [1e-6, 5e-7])
      compareDirection(system, state, sparse, direction, h, statistics);
    if (chart === 0) {
      const before = system.evaluate(state);
      state = system.rebase(state);
      const after = system.evaluate(state);
      assert.deepEqual(after.nodes, before.nodes);
      assert.deepEqual(after.allocation, before.allocation);
      assert.deepEqual(after.transportSpeeds, before.transportSpeeds);
      after.residual.forEach((r, i) => close(r, before.residual[i], 2e-12));
      layout.positions.forEach(p => assert.equal(state[p.column], 0));
    }
  }
  t.diagnostic(JSON.stringify(statistics));
});

test('filtered body momentum uses the full upstream stencil across both leading and trailing edges', () => {
  const input = fixture(2), system = createStreamtubeBodySystem(input), state = stateFor(system), value = system.evaluate(state);
  const touchedEdges = new Set();
  for (let g = 0; g <= system.layout.elements; g++) for (let j = 0; j < system.layout.tubes[g]; j++) {
    const centers = value.nodes[g].map(row => ({ x: .5 * (row[j].x + row[j + 1].x), y: .5 * (row[j].y + row[j + 1].y) }));
    const lengths = centers.slice(1).map((p, i) => Math.hypot(p.x - centers[i].x, p.y - centers[i].y));
    const q = value.sections.map(row => row[g][j].q), m2 = value.sections.map(row => row[g][j].machSquared);
    for (let k = 0; k < system.layout.nx; k++) {
      let expected = q[k];
      if (k >= 2) {
        const d0 = .5 * (lengths[k - 2] + lengths[k - 1]), d1 = .5 * (lengths[k - 1] + lengths[k]);
        const activation = .5 * (m2[k - 1] + m2[k]), e = 1 - upwind.mcrit;
        const coefficient = upwind.mucon / system.conditions.gamma * e * Math.log1p(Math.exp((1 - 1 / activation) / e));
        expected += coefficient * (-(q[k] - q[k - 1]) + d1 / d0 * (q[k - 1] - q[k - 2]));
      }
      close(value.transportSpeeds[k][g][j], expected);
      if (k > 0) assert.equal(value.cells[k - 1][g][j].transportSpeeds[1], value.transportSpeeds[k][g][j]);
      if (k < system.layout.nx - 1) assert.equal(value.cells[k][g][j].transportSpeeds[0], value.transportSpeeds[k][g][j]);
      input.bodies.forEach((b, body) => { for (const edge of ['leadingIndex', 'trailingIndex'])
        if (k === b[edge] && Math.abs(expected - q[k]) > 1e-9) touchedEdges.add(`${body}/${edge}`); });
    }
  }
  assert.equal(touchedEdges.size, 2 * input.bodies.length, 'Both edges of every element retain a nonzero upstream correction.');
  const i = input.bodies[0].trailingIndex + 1, g = 0, j = 0, { layout } = system;
  const row = layout.rows.find(r => r.kind === 'streamwise' && r.i === i && r.group === g && r.tube === j).index;
  const remote = [layout.densityIndex(i - 3, g, j), layout.nodes[g][i - 3][j + 1].column];
  const biased = system.jacobian(state);
  const plain = createStreamtubeBodySystem({ ...input, upwind: { ...upwind, mucon: 0 } }).jacobian(state);
  for (const col of remote) {
    assert.notEqual(col, null);
    assert.ok(Math.abs(biased[row * layout.n + col]) > 1e-10, `Missing wider upstream column ${col}`);
    assert.equal(plain[row * layout.n + col], 0);
  }
});

test('independent body passage flux integration keeps physical mass and enthalpy separate from biased momentum', () => {
  const system = createStreamtubeBodySystem(fixture(2)), value = system.evaluate(stateFor(system));
  let largestCorrection = 0;
  for (let g = 0; g <= system.layout.elements; g++) {
    const check = directUpwindChannelConservation({ nodes: value.nodes[g], sections: value.sections.map(row => row[g]),
      cells: value.cells.map(row => row[g]), massFlows: value.allocation.groups[g].map(tube => tube.massFlow) },
    system.conditions.h0, system.conditions.gamma);
    close(check.maximumMassMismatch, 0); close(check.maximumEnthalpyMismatch, 0);
    largestCorrection = Math.max(largestCorrection, check.maximumTransportMomentumDifference);
    for (const kind of ['biased', 'physical']) for (const field of ['total', 'external', 'maxLocal'])
      for (const k of [0, 3]) close(check[kind][field][k], 0);
    for (const c of check.local) {
      const row = system.layout.rows.find(r => r.kind === 'streamwise' && r.i === c.i && r.group === g && r.tube === c.j);
      const residual = value.residual[row.index] * system.conditions.pressureScale;
      close(c.biased[1], residual * c.transverse.y);
      close(c.biased[2], -residual * c.transverse.x);
    }
    for (let k = 0; k < 4; k++) close(check.physicalMinusBiased.total[k], check.physicalMinusBiased.external[k]);
  }
  assert.ok(largestCorrection > 1e-5);
  // Pressure matching is not solved in this frozen fixture: no claim of a
  // converged body balance or physical shock accuracy follows from the check.
});

test('upwind Euler-to-displacement block differentiates wall and wake thickness on a fixed chart', t => {
  const input = fixture(), nx = input.outerLower.length - 1;
  input.displacement = { surfaces: input.bodies.map(b => Object.fromEntries(['upper', 'lower'].map((side, s) => [side,
    Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0001 * (1 + .02 * (s + 1) * i))]))),
  wakes: input.bodies.map(b => Array.from({ length: nx - b.trailingIndex }, (_, i) => .00024 * (1 + .01 * i))) };
  const system = createStreamtubeBodySystem(input), state = stateFor(system);
  const block = system.jacobian(state, { includeDisplacement: true });
  assert.deepEqual(block.state, system.jacobian(state));
  assert.deepEqual(sparseDense(system.jacobian(state, { sparse: true, includeDisplacement: true }).state), block.state);
  const original = system.evaluate(state), statistics = { comparisons: 0, maximumError: 0, columns: block.parameters.length };
  assert.ok(block.parameters.some(p => p.kind === 'wake'));
  for (const [col, parameter] of block.parameters.entries()) {
    const h = 1e-7, samples = [];
    try {
      for (const factor of [-2, -1, 1, 2]) {
        const displacement = structuredClone(input.displacement);
        if (parameter.kind === 'wake') displacement.wakes[parameter.body][parameter.index] += factor * h;
        else for (const side of parameter.side === 'both' ? ['upper', 'lower'] : [parameter.side])
          displacement.surfaces[parameter.body][side][parameter.index] += factor * h;
        system.setDisplacement(displacement);
        samples.push(system.residual(state));
      }
    } finally { system.setDisplacement(input.displacement); }
    let effect = 0;
    for (let row = 0; row < system.layout.n; row++) {
      const exact = block.displacement[row].get(col) ?? 0, fd = fourth(samples.map(r => r[row]), h);
      effect = Math.max(effect, Math.abs(exact));
      statistics.maximumError = Math.max(statistics.maximumError, Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)));
      statistics.comparisons++; close(exact, fd, 3e-6, `displacement ${JSON.stringify(parameter)}, row ${row}`);
    }
    assert.ok(effect > 1e-5, `Missing physical displacement effect ${JSON.stringify(parameter)}`);
  }
  assert.deepEqual(system.evaluate(state).residual, original.residual);
  assert.deepEqual(system.evaluate(state).nodes, original.nodes);
  t.diagnostic(JSON.stringify(statistics));
});

test('explicit upwinding admits a supersonic interior and rejects unsupported incompressible or isentropic use', () => {
  const input = fixture(), system = createStreamtubeBodySystem(input), state = system.initial.slice();
  const i = input.bodies[0].trailingIndex + 1;
  state[system.layout.densityIndex(i, 0, 0)] = -.8;
  const value = system.evaluate(state);
  assert.ok(value.sections[i][0][0].machSquared > 1);
  assert.ok(value.residual.every(Number.isFinite));
  for (const endpoint of [0, system.layout.nx - 1]) {
    const bad = state.slice(); bad[system.layout.densityIndex(endpoint, 0, 0)] = -.8;
    assert.throws(() => system.evaluate(bad), /inlet\/outlet section must remain subsonic/);
  }
  const { upwind: ignored, ...plain } = input;
  assert.throws(() => createStreamtubeBodySystem(plain).evaluate(state), /sonic/);
  assert.throws(() => createStreamtubeBodySystem({ ...input, flowModel: 'incompressible', mach: 0 }));
  assert.throws(() => createStreamtubeBodySystem({ ...input, streamwiseMode: 'isentropic' }));
  const zero = createStreamtubeBodySystem({ ...input, upwind: { ...upwind, mucon: 0 } });
  assert.deepEqual(zero.residual(zero.initial), createStreamtubeBodySystem(plain).residual(zero.initial));
});
