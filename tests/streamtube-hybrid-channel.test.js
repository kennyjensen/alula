// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-upwind-channel.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { channelGas, nozzleChannel } from './oracles/streamtube.js';

const input = () => {
  const gas = channelGas(.75), { exact: unusedExact, ...geometry } = nozzleChannel(7, 3);
  return { ...geometry, ...gas, outletPressure: gas.referencePressure,
    upwind: { mucon: 1, mcrit: .6, boundary: { kind: 'unfiltered-first-two' } } };
};
const hybridInput = () => ({ ...input(), streamwiseMode: 'hybrid', hybrid: { epsilonP: .001 } });
const skew = system => system.initial.map((v, i) => v + .003 * Math.sin(.7 * (i + 1)));
const fourth = (a, h) => (a[0] - 8 * a[1] + 8 * a[2] - a[3]) / (12 * h);
const close = (a, b, tolerance = 3e-9, label = '') => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${label}: ${a} != ${b}`);

test('omitted and explicit momentum retain the archived complete channel values and Jacobians', async () => {
  const archive = JSON.parse(fs.readFileSync(new URL('../docs/hybrid-flow/channel-before.json', import.meta.url)));
  const bytes = fs.readFileSync(new URL('../' + archive.archive, import.meta.url));
  assert.equal(createHash('sha256').update(bytes).digest('hex'), archive.sha256);
  const source = new URL('../' + archive.path, import.meta.url);
  const code = bytes.toString().replace(/from '([^']+)'/g, (_, specifier) => {
    const original = new URL(specifier, source);
    const resolved = fs.existsSync(original) ? original
      : new URL(specifier, new URL('../src/euler/tests/streamtube-upwind-channel.js', import.meta.url));
    return `from '${resolved.href}'`;
  });
  const previous = await import(`data:text/javascript;base64,${Buffer.from(code).toString('base64')}`);
  for (const mode of [undefined, 'momentum']) {
    const p = { ...input(), ...(mode ? { streamwiseMode: mode } : {}) };
    const current = createUpwindStreamtubeChannel(p), old = previous.createUpwindStreamtubeChannel(p), state = skew(current);
    assert.deepEqual(current.conditions, old.conditions);
    assert.deepEqual(current.initial, old.initial);
    assert.deepEqual(current.evaluate(state), old.evaluate(state));
    assert.deepEqual(current.jacobian(state), old.jacobian(state));
    assert.deepEqual(current.jacobian(state, { sparse: true }), old.jacobian(state, { sparse: true }));
    assert.equal(Object.hasOwn(current.evaluate(state), 'hybrid'), false);
  }
});

test('complete hybrid channel derivatives include the changing blend fraction within the existing sparse stencil', t => {
  const p = hybridInput(), system = createUpwindStreamtubeChannel(p), state = skew(system), { n, nx, nt } = system;
  assert.ok(n < 100);
  const value = system.evaluate(state), dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
  const momentum = createUpwindStreamtubeChannel({ ...p, streamwiseMode: 'momentum' });
  const ordinary = momentum.evaluate(state), ordinaryJac = momentum.jacobian(state, { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense);
  assert.deepEqual(sparse.rowPtr, ordinaryJac.rowPtr);
  assert.deepEqual(sparse.colIndex, ordinaryJac.colIndex);
  assert.deepEqual(value.cells, ordinary.cells, 'The physical S/N kernel must be unchanged.');
  assert.deepEqual(value.sections, ordinary.sections);
  assert.deepEqual(value.transportSpeeds, ordinary.transportSpeeds);
  const firstBoundary = (nx - 1) * nt;
  assert.deepEqual(value.residual.slice(firstBoundary), ordinary.residual.slice(firstBoundary));
  const ordinaryDense = sparseDense(ordinaryJac);
  assert.deepEqual(dense.slice(firstBoundary * n), ordinaryDense.slice(firstBoundary * n));
  assert.ok(value.hybrid.flat().some(b => b.fraction > 0 && b.fraction < 1), 'Exercise a mixed streamwise row.');
  // Independent transcription of the declared research switch. This is a
  // numerical assembly check, not an independent shock-accuracy reference.
  const step = x => x <= 0 ? 0 : x >= 1 ? 1 : 6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3;
  for (let i = 0; i < nx - 1; i++) for (let j = 0; j < nt; j++) {
    const cell = value.cells[i][j], [a, b] = cell.states;
    const change = b.rho * b.q / b.p - a.rho * a.q / a.p;
    const losses = cell.transportSpeeds.map((q, k) => -(q - cell.states[k].q) * change);
    const [wa, wb] = losses.map(l => step(l / p.hybrid.epsilonP));
    const fraction = step(Math.log(a.q / b.q) / Math.sqrt(p.hybrid.epsilonP)) * (wa + (1 - wa) * wb);
    close(value.hybrid[i][j].fraction, fraction, 2e-12);
    losses.forEach((l, k) => close(value.hybrid[i][j].lossIndicators[k], l, 2e-12));
    close(value.residual[i * nt + j],
      (fraction * cell.streamwiseResidual + (1 - fraction) * cell.isentropicResidual) / p.referencePressure, 2e-12);
  }
  const stats = { unknowns: n, comparisons: 0, maximumNormalizedError: 0, maximumBlendDerivativeEffect: 0 };
  for (let col = 0; col < n; col++) for (const h of [2e-5, 1e-5]) {
    const samples = [-2, -1, 1, 2].map(k => { const x = state.slice(); x[col] += k * h; return system.evaluate(x); });
    for (let row = 0; row < n; row++) {
      const fd = fourth(samples.map(v => v.residual[row]), h), exact = dense[row * n + col];
      const error = Math.abs(fd - exact) / Math.max(1, Math.abs(fd), Math.abs(exact));
      stats.maximumNormalizedError = Math.max(stats.maximumNormalizedError, error); stats.comparisons++;
      close(exact, fd, 3e-9, `row ${row}, column ${col}, h=${h}`);
      if (row < firstBoundary) {
        const i = Math.floor(row / nt), j = row % nt, f = value.hybrid[i][j].fraction;
        const dR1 = fourth(samples.map(v => v.cells[i][j].streamwiseResidual), h);
        const dR2 = fourth(samples.map(v => v.cells[i][j].isentropicResidual), h);
        const frozenFraction = (f * dR1 + (1 - f) * dR2) / p.referencePressure;
        stats.maximumBlendDerivativeEffect = Math.max(stats.maximumBlendDerivativeEffect, Math.abs(exact - frozenFraction));
      }
    }
  }
  assert.ok(stats.maximumBlendDerivativeEffect > 1e-7, 'The test must detect omission of df*(R1-R2).');
  const direction = state.map((_, i) => .3 * Math.cos(i + .2)), h = 1e-5;
  const samples = [-2, -1, 1, 2].map(k => system.residual(state.map((v, i) => v + k * h * direction[i])));
  sparseProduct(sparse, direction).forEach((v, row) => close(v, fourth(samples.map(r => r[row]), h)));
  const row = 3 * nt + 1;
  for (const col of [nt + 1, system.positionIndex(1, 1)])
    assert.ok(Math.abs(dense[row * n + col]) > 1e-8, 'The full upstream speed-filter footprint remains present.');
  system.jacobian(state.map(v => .9 * v), { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense);
  t.diagnostic(JSON.stringify(stats));
});

test('zero artificial-speed bias selects the physical isentropic row without changing normal momentum or boundaries', () => {
  const p = hybridInput(); p.upwind.mucon = 0;
  const system = createUpwindStreamtubeChannel(p), state = skew(system), value = system.evaluate(state);
  const ordinary = createUpwindStreamtubeChannel({ ...p, streamwiseMode: 'momentum' }).evaluate(state);
  for (let i = 0; i < system.nx - 1; i++) for (let j = 0; j < system.nt; j++) {
    const cell = value.cells[i][j], blend = value.hybrid[i][j];
    assert.equal(blend.fraction, 0);
    assert.equal(value.residual[i * system.nt + j], cell.isentropicResidual / p.referencePressure);
    assert.deepEqual(cell.transportSpeeds, cell.states.map(s => s.q));
  }
  assert.deepEqual(value.cells, ordinary.cells);
  assert.deepEqual(value.residual.slice((system.nx - 1) * system.nt), ordinary.residual.slice((system.nx - 1) * system.nt));
  const matrix = system.jacobian(state), row = 3 * system.nt + 1;
  assert.equal(matrix[row * system.n + system.nt + 1], 0);
  assert.equal(matrix[row * system.n + system.positionIndex(1, 1)], 0);
});

test('hybrid pressure-loss scale is explicit, positive, and isolated from caller edits', () => {
  for (const hybrid of [undefined, {}, { epsilonP: 0 }, { epsilonP: -1 }, { epsilonP: NaN }, { epsilonP: Infinity }])
    assert.throws(() => createUpwindStreamtubeChannel({ ...input(), streamwiseMode: 'hybrid', hybrid }), /positive finite epsilonP/);
  const p = hybridInput(), before = structuredClone(p), system = createUpwindStreamtubeChannel(p), state = skew(system);
  const value = system.evaluate(state);
  assert.deepEqual(p, before);
  assert.equal(system.conditions.streamwiseMode, 'hybrid');
  assert.deepEqual(system.conditions.hybrid, before.hybrid);
  p.hybrid.epsilonP = NaN;
  assert.deepEqual(system.evaluate(state), value);
  const rescaled = createUpwindStreamtubeChannel({ ...before, referencePressure: 7 * before.referencePressure }).evaluate(state);
  assert.deepEqual(rescaled.cells, value.cells);
  assert.deepEqual(rescaled.hybrid, value.hybrid, 'Reference-pressure normalization is applied once, after the raw-pressure blend.');
  for (let row = 0; row < (system.nx - 1) * system.nt; row++)
    close(rescaled.residual[row], value.residual[row] / 7, 2e-13);
});
