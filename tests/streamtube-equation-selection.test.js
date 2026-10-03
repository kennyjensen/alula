// SPDX-License-Identifier: GPL-2.0-or-later
// Tiny residual/Jacobian checks only; no Euler/BL Newton solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { normalizeStreamtubeEquationSelection, streamtubeEquationAt, ISMOM3_LEADING_REGION }
  from '../src/euler/streamtube-equation-selection.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const serial = v => JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x);
const parameters = (ismom, elements = 1, tubes = 2) => ({
  ...intrinsicBodyFixture({ elements, bodySegments: 4, tubes, mach: .5, alpha: .25 }),
  streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-3, ...(ismom === undefined ? {} : { ismom }) },
  upwind: { mucon: 1, mcrit: .5, boundary: { kind: 'unfiltered-first-two' } },
});
const seeded = s => s.initial.map((_, k) => (k < s.layout.densityCount ? 1e-3 : 1e-5) * Math.sin(k + .3));
const relative = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));

test('ISMOM controls preserve the default schema, copy explicit selection, and reject invalid values', () => {
  const h = { epsilonP: .003 }, old = normalizeStreamtubeEquationSelection(h);
  assert.deepEqual(old, h); assert.ok(Object.isFrozen(old)); h.epsilonP = .4; assert.equal(old.epsilonP, .003);
  for (const ismom of [1, 2, 3, 4]) {
    const p = { epsilonP: .001, ismom }, q = normalizeStreamtubeEquationSelection(p);
    assert.deepEqual(q, p); p.ismom = 7; assert.equal(q.ismom, ismom);
    assert.equal(createStreamtubeBodySystem(parameters(ismom)).conditions.hybrid.ismom, ismom);
  }
  for (const ismom of [0, 5, -1, 1.5, '3', null, NaN, Infinity])
    assert.throws(() => normalizeStreamtubeEquationSelection({ epsilonP: .01, ismom }), /ISMOM/);
  for (const epsilonP of [0, -1, NaN, Infinity])
    assert.throws(() => normalizeStreamtubeEquationSelection({ epsilonP, ismom: 1 }), /epsilonP/);
});

test('ISMOM3 uses an explicit inlet/LE logical region with exact endpoints and multielement union', () => {
  assert.deepEqual(ISMOM3_LEADING_REGION, { downstreamCells: 10, transverseCells: 4 });
  const input = { hybrid: { epsilonP: .01, ismom: 3 }, bodies: [{ leadingIndex: 2 }, { leadingIndex: 8 }], tubes: [6, 3, 7] };
  const saved = serial(input), at = (i, globalTube) => {
    let group = 0, tube = globalTube;
    while (tube >= input.tubes[group]) tube -= input.tubes[group++];
    return streamtubeEquationAt({ ...input, i, group, tube });
  };
  // Cuts at global streamline6 and9. Their four-tube neighborhoods overlap.
  for (const i of [1, 2, 12, 13, 18, 19]) for (let tube = 0; tube < 16; tube++) {
    const expected = (i <= 12 && tube >= 2 && tube <= 9) || (i <= 18 && tube >= 5 && tube <= 12);
    assert.equal(at(i, tube), expected ? 'isentropic' : 'momentum', `i=${i},tube=${tube}`);
  }
  assert.equal(serial(input), saved);
  for (const bad of [{ i: 0 }, { group: 3 }, { tube: 6 }, { bodies: [{ leadingIndex: 0 }, { leadingIndex: 8 }] }])
    assert.throws(() => streamtubeEquationAt({ ...input, i: 1, group: 0, tube: 0, ...bad }), /logical/);
  for (const [ismom, expected] of [[1, 'momentum'], [2, 'isentropic'], [4, 'hybrid']])
    assert.equal(streamtubeEquationAt({ hybrid: { ismom } }), expected, 'Fixed modes do not require irrelevant region data.');
});

test('default and explicit ISMOM4 retain archived default residual and Jacobian arithmetic', async () => {
  const root = path.resolve('src/euler');
  const movedModules = {
    './streamtube-base-geometry.js': './streamtube-geometry.js',
    './streamtube-motion-directions.js': './streamtube-geometry.js',
    './streamtube-geometry-transfer.js': './streamtube-geometry.js',
    './streamtube-force-coefficients.js': './streamtube-forces.js',
    './streamtube-solid-pressure-forces.js': './streamtube-forces.js',
    './streamtube-wake-gap.js': './streamtube-wake-geometry.js',
    './streamtube-wake-width-increment.js': './streamtube-wake-geometry.js',
    './streamtube-displacement-initializer.js': './streamtube-displacement.js',
    './streamtube-equation-controls.js': './streamtube-equation-selection.js',
  };
  const resolve = (text, replacement) => text.replace(/from '(\.\.?\/[^']+)'/g, (_, specifier) =>
    `from '${specifier === './streamtube-body-jacobian.js' && replacement ? replacement
      : pathToFileURL(path.resolve(root, movedModules[specifier] ?? specifier)).href}'`);
  const dataURL = text => `data:text/javascript;base64,${Buffer.from(text).toString('base64')}`;
  const oldJac = dataURL(resolve(fs.readFileSync('docs/ismom-equation-selection/before/streamtube-body-jacobian.js.txt', 'utf8')));
  const oldBody = dataURL(resolve(fs.readFileSync('docs/ismom-equation-selection/before/streamtube-body.js.txt', 'utf8'), oldJac));
  const { createStreamtubeBodySystem: createOld } = await import(oldBody);
  const old = createOld(parameters()), now = createStreamtubeBodySystem(parameters()), explicit = createStreamtubeBodySystem(parameters(4));
  const x = seeded(old), before = old.evaluate(x), current = now.evaluate(x);
  assert.ok(before.diagnostics.hybrid.blendedCells > 0);
  assert.equal(serial(current), serial(before)); assert.deepEqual(now.conditions, old.conditions);
  assert.deepEqual(now.jacobian(x), old.jacobian(x));
  assert.deepEqual(explicit.evaluate(x).residual, before.residual); assert.deepEqual(explicit.jacobian(x), old.jacobian(x));
});

test('ISMOM1 selects exact conservative rows/J; ISMOM2 selects exact biased entropy rows', () => {
  const aInput = parameters(1), rawInput = { ...aInput, streamwiseMode: 'momentum' }; delete rawInput.hybrid;
  const a = createStreamtubeBodySystem(aInput), raw = createStreamtubeBodySystem(rawInput), x = seeded(a);
  assert.deepEqual(a.evaluate(x).residual, raw.evaluate(x).residual);
  assert.deepEqual(a.jacobian(x), raw.jacobian(x));
  const b = createStreamtubeBodySystem(parameters(2)), flow = b.evaluate(x);
  let biasedEntropyDifferences = 0;
  for (const row of b.layout.rows) if (row.kind === 'streamwise') {
    const cell = flow.cells[row.i - 1][row.group][row.tube];
    assert.equal(flow.residual[row.index], cell.isentropicResidual / b.conditions.pressureScale);
    assert.equal(flow.hybridCells[row.i - 1][row.group][row.tube].fraction, 0);
    if (cell.artificialEntropyJump !== cell.entropyJump) biasedEntropyDifferences++;
  }
  assert.ok(biasedEntropyDifferences > 0, 'Mode2 must use the shared biased entropy, not silently reset physical isentropy.');
  const zero = parameters(2); zero.upwind.mucon = 0;
  const unfiltered = { ...zero, streamwiseMode: 'isentropic' }; delete unfiltered.hybrid; delete unfiltered.upwind;
  const z = createStreamtubeBodySystem(zero), old = createStreamtubeBodySystem(unfiltered);
  assert.deepEqual(z.evaluate(x).residual, old.evaluate(x).residual);
  const jz = z.jacobian(x), jo = old.jacobian(x);
  assert.ok(jz.every((v, k) => relative(v, jo[k]) < 3e-12));
});

test('ISMOM1/2/3 analytic Jacobians include density, moving geometry and globals on two elements', t => {
  let worst = 0, comparisons = 0;
  const evidence = [];
  for (const ismom of [1, 2, 3]) {
    const s = createStreamtubeBodySystem(parameters(ismom, 2, 5)), x = seeded(s), flow = s.evaluate(x), j = s.jacobian(x, { sparse: true });
    if (ismom === 3) {
      assert.ok(flow.diagnostics.hybrid.entropyCells > 0 && flow.diagnostics.hybrid.momentumCells > 0);
      let entropy = 0, momentum = 0;
      for (const row of s.layout.rows) if (row.kind === 'streamwise') {
        const globalTube = s.layout.tubes.slice(0, row.group).reduce((a, n) => a + n, 0) + row.tube;
        const entropyRegion = s.layout.bodies.some((body, b) => {
          const cut = s.layout.tubes.slice(0, b + 1).reduce((a, n) => a + n, 0);
          return row.i <= body.leadingIndex + 10 && globalTube >= cut - 4 && globalTube < cut + 4;
        });
        const cell = flow.cells[row.i - 1][row.group][row.tube];
        assert.equal(flow.residual[row.index], (entropyRegion ? cell.isentropicResidual : cell.streamwiseResidual) / s.conditions.pressureScale);
        entropyRegion ? entropy++ : momentum++;
      }
      evidence.push({ ismom, entropy, momentum });
    }
    const globals = Object.values(s.layout.globals).flat().filter(c => c !== null);
    const directions = [Float64Array.from(x, (_, k) => .02 * Math.cos(k + .4)),
      Float64Array.from(x, (_, k) => globals.includes(k) ? .02 * Math.sin(k + .7) : 0)];
    for (const d of directions) {
      const h = 7e-7, values = [-2, -1, 1, 2].map(m => s.residual(x.map((v, k) => v + m * h * d[k]))), exact = sparseProduct(j, d);
      for (let k = 0; k < s.layout.n; k++) {
        const fd = (values[0][k] - 8 * values[1][k] + 8 * values[2][k] - values[3][k]) / (12 * h);
        worst = Math.max(worst, relative(exact[k], fd)); comparisons++;
      }
    }
  }
  assert.ok(worst < 2e-7, String(worst)); t.diagnostic(JSON.stringify({ comparisons, worst, evidence }));
});

test('each fixed equation admits a thermal supersonic interior while retaining the remote subsonic gate', () => {
  for (const ismom of [1, 2, 3]) {
    const s = createStreamtubeBodySystem(parameters(ismom)), x = s.initial.slice(), i = s.layout.bodies[0].trailingIndex + 1;
    x[s.layout.densityIndex(i, 0, 0)] = -.8;
    const flow = s.evaluate(x); assert.ok(flow.sections[i][0][0].machSquared > 1);
    assert.ok(flow.residual.every(Number.isFinite));
    for (const boundary of [0, s.layout.nx - 1]) {
      const bad = x.slice(); bad[s.layout.densityIndex(boundary, 0, 0)] = -.8;
      assert.throws(() => s.evaluate(bad), /inlet\/outlet section must remain subsonic/);
    }
  }
});

test('fixed equation derivatives retain surface and wake displacement terms', t => {
  let worst = 0, comparisons = 0;
  for (const ismom of [1, 2, 3]) {
    const input = parameters(ismom), nx = input.outerLower.length - 1;
    input.displacement = { surfaces: input.bodies.map(b => Object.fromEntries(['upper', 'lower'].map(side =>
      [side, Array(b.trailingIndex - b.leadingIndex + 1).fill(.0001)]))),
      wakes: input.bodies.map(b => Array(nx - b.trailingIndex).fill(.0002)) };
    const s = createStreamtubeBodySystem(input), x = seeded(s), block = s.jacobian(x, { includeDisplacement: true });
    const direction = block.parameters.map((_, k) => .01 * Math.cos(k + .6)), samples = [], h = 2e-6;
    try {
      for (const m of [-2, -1, 1, 2]) {
        const displacement = structuredClone(input.displacement);
        block.parameters.forEach((p, k) => {
          if (p.kind === 'wake') displacement.wakes[p.body][p.index] += m * h * direction[k];
          else for (const side of p.side === 'both' ? ['upper', 'lower'] : [p.side])
            displacement.surfaces[p.body][side][p.index] += m * h * direction[k];
        });
        s.setDisplacement(displacement); samples.push(s.residual(x));
      }
    } finally { s.setDisplacement(input.displacement); }
    for (let row = 0; row < s.layout.n; row++) {
      const exact = [...block.displacement[row]].reduce((sum, [col, value]) => sum + value * direction[col], 0);
      const fd = (samples[0][row] - 8 * samples[1][row] + 8 * samples[2][row] - samples[3][row]) / (12 * h);
      worst = Math.max(worst, relative(exact, fd)); comparisons++;
    }
  }
  assert.ok(worst < 3e-7, String(worst)); t.diagnostic(JSON.stringify({ comparisons, worst }));
});
