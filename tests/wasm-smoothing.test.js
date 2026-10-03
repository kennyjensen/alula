import test from 'node:test';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createWasmSmoother } from '../src/geometry/wasm-smoothing.js';
import { smoothPairedBoundaryGrid } from '../src/geometry/paired-boundary-slor.js';

function problem(nt, sourceForm, scheme, sides) {
  const nx = 16, xi = Array.from({ length: nx + 1 }, (_, i) => (i / nx) ** 1.1);
  const eta = Array.from({ length: nt + 1 }, (_, j) => (j / nt) ** 1.3);
  const nodes = xi.map(u => eta.map(e => ({ x: u + .02 * Math.sin(Math.PI * u) * (1 + e),
    y: e + .03 * Math.sin(Math.PI * u) * (1 + .2 * e) + .01 * Math.sin(2 * Math.PI * u) * Math.sin(Math.PI * e) })));
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: eta.slice(1).map((v, j) => v - eta[j]), streamwiseCoordinates: xi,
    discretization: 'giles-1985', lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
    streamwiseSourceDiscretization: scheme, orthogonalBoundaryControl: { sourceForm, sides,
      background: nodes.map(row => row.map(() => .02)), decay: { lower: 3, upper: 4 },
      corners: Object.fromEntries(sides.map(side => [side, [6]])),
      activeStations: Object.fromEntries(sides.map(side => [side, xi.map((_, i) => i > 2 && i < nx - 2)])) } });
  return { nodes, system };
}
function close(a, b, tolerance = 2e-10) {
  assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
}
for (const nt of [3, 4, 11]) for (const form of ['poisson', 'metric-stretch']) for (const scheme of ['centered', 'grape-1980'])
  test(`WASM complete sweeps match JS: tubes=${nt}, ${form}, ${scheme}`, () => {
    const { system, nodes } = problem(nt, form, scheme, ['lower', 'upper']);
    const kernel = createWasmSmoother(system);
    let js = nodes, wasm = structuredClone(nodes), rejected = false;
    for (let iteration = 0; iteration < 8; iteration++) {
      const state = kernel.evaluate(wasm), reference = system.residuals(js);
      close(state.residual, reference.residual);
      close(state.merit, .5 * reference.rows.reduce((s, p) => s + p.x ** 2 + p.y ** 2, 0));
      if (reference.residual < 1e-10) break;
      const omega = iteration === 0 ? 1.9999 : 1;
      const expected = system.sweep(js, omega), actual = kernel.sweep(wasm, omega, { lineTolerance: 0 });
      actual.nodes.forEach((row, i) => row.forEach((p, j) => { close(p.x, expected.nodes[i][j].x); close(p.y, expected.nodes[i][j].y); }));
      assert.deepEqual(actual.rowSteps.map(r => r.fraction), expected.rowSteps.map(r => r.fraction));
      close(actual.maxUpdate, expected.maxUpdate);
      actual.rowSteps.forEach((r, i) => { close(r.baseMerit, expected.rowSteps[i].baseMerit); close(r.merit, expected.rowSteps[i].merit); });
      rejected ||= actual.rowSteps.some(r => r.trials.length > 1);
      js = expected.nodes; wasm = actual.nodes;
    }
    assert.equal(kernel.replays, 0, 'parity must be achieved in WASM, not via JS replay');
    if (nt === 11) assert.ok(rejected, 'exercise rejected steps');
  });

test('WASM instances remain independent across nested observer solves', () => {
  const first = problem(7, 'metric-stretch', 'centered', ['lower']);
  const second = problem(4, 'poisson', 'centered', ['lower', 'upper']);
  const reference = smoothPairedBoundaryGrid(first.system, { maxSweeps: 3, backend: 'javascript' });
  const result = smoothPairedBoundaryGrid(first.system, { maxSweeps: 3, onSweep: () => {
    const nested = smoothPairedBoundaryGrid(second.system, { maxSweeps: 2 });
    assert.equal(nested.backend, 'wasm');
  } });
  assert.equal(result.backend, 'wasm'); assert.equal(result.referenceReplays, 0);
  result.nodes.forEach((row, i) => row.forEach((p, j) => { close(p.x, reference.nodes[i][j].x); close(p.y, reference.nodes[i][j].y); }));
});

test('invalid starting grids retain original diagnostics and coordinates', () => {
  const { system, nodes } = problem(7, 'metric-stretch', 'centered', ['lower', 'upper']);
  const invalid = structuredClone(nodes); invalid[4][3].y = -2;
  const actual = smoothPairedBoundaryGrid(system, { initial: invalid });
  const expected = smoothPairedBoundaryGrid(system, { initial: invalid, backend: 'javascript' });
  assert.deepEqual(actual, expected);
});

for (const mode of ['harmonic', 'stretch', 'poisson']) for (const scheme of ['centered', 'grape-1980'])
  test(`WASM scalar SLOR matches JS: ${mode}/${scheme}`, () => {
    const { nodes, system: source } = problem(11, 'metric-stretch', scheme, ['lower', 'upper']);
    const controls = nodes.map(row => row.map(() => .04));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: source.eta.slice(1).map((v, j) => v - source.eta[j]),
      streamwiseCoordinates: source.xi, discretization: 'giles-1985', streamwiseSourceDiscretization: scheme,
      ...(mode === 'stretch' ? { streamwiseStretch: controls } : mode === 'poisson' ? { streamwiseSource: controls } : {}) });
    const expected = smoothEllipticStreamtubeGrid(system, { backend: 'javascript', maxSweeps: 20, requireConvex: true });
    const actual = smoothEllipticStreamtubeGrid(system, { maxSweeps: 20, requireConvex: true });
    assert.equal(actual.backend, 'wasm'); assert.equal(actual.referenceReplays, 0);
    assert.deepEqual(actual.nodes, expected.nodes);
    assert.deepEqual(actual.history, expected.history);
  });

test('unsupported moving-boundary formulation keeps the JS solver', () => {
  const nodes = Array.from({ length: 5 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({ x: i / 4, y: j / 3 })));
  const system = createEllipticStreamtubeGrid({ nodes, massFlows: [1, 1, 1], discretization: 'giles-1985', boundaryConditions: { lower: 'giles-vertical' } });
  assert.equal(smoothEllipticStreamtubeGrid(system).backend, 'javascript');
});

test('shipped WASM artifact matches its source and build manifest', () => {
  const directory = new URL('../src/geometry/wasm/', import.meta.url);
  const build = JSON.parse(readFileSync(new URL('BUILD.json', directory)));
  const hash = name => createHash('sha256').update(readFileSync(new URL(name, directory))).digest('hex');
  assert.equal(hash('smoothing.ts'), build.sourceSha256, 'run npm run build:smoothing after kernel edits');
  assert.equal(hash('smoothing.wasm'), build.wasmSha256);
});
