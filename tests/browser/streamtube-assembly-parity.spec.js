// SPDX-License-Identifier: GPL-2.0-or-later
// Tiny complete matrix and one accepted update in a real browser worker.
// The same current Newton driver uses current versus archived assembly modules.
import fs from 'node:fs';
import { test, expect } from '@playwright/test';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

const fixturePath = 'docs/coupled-current-profile-preparation/two-element-six-update/initial.json';
const fixture = JSON.parse(fs.readFileSync(fixturePath)).checkpoint;
const paths = [
  ['src/numerics/sparse.js', 'docs/performance-assembly/before/sparse.js.txt'],
  ['src/euler/streamtube-body-jacobian.js', 'docs/performance-assembly/before/streamtube-body-jacobian.js.txt'],
  ['src/euler/streamtube-body.js', 'docs/performance-assembly/before/streamtube-body.js.txt'],
  ['src/euler/streamtube-coupled.js', 'docs/performance-assembly/before/streamtube-coupled.js.txt'],
  ['src/euler/streamtube-coupled-ises.js', 'src/euler/streamtube-coupled-ises.js'],
];

test('optimized and archived coupled assembly give identical matrices and one Newton step in a browser worker', async ({ page }) => {
  const receipt = process.env.MSES_ASSEMBLY_BROWSER_RECEIPT ?? 'docs/performance-assembly/browser-matrix-step.json';
  expect(fs.existsSync(receipt), 'Keep prior receipts; select a new MSES_ASSEMBLY_BROWSER_RECEIPT path.').toBe(false);
  const sources = numericalSourceHashes(['tests/browser/streamtube-assembly-parity.spec.js', fixturePath,
    'third_party/klu/BUILD.json', ...paths.map(([, archive]) => archive)]);
  const entries = paths.map(([original, archive]) => ({ original, source: fs.readFileSync(archive, 'utf8') }));
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.route('**/assembly-parity-check', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Assembly check</title>' }));
  await page.goto('/assembly-parity-check');
  const result = await page.evaluate(payload => new Promise((resolve, reject) => {
    const run = async ({ checkpoint, entries, origin }) => {
      const blobs = [], urls = new Map(), moduleURL = relative => new URL('/' + relative, origin).href;
      for (const { original, source } of entries) {
        const originalURL = moduleURL(original);
        const rewritten = source.replace(/from\s+(['"])(\.[^'"]+)\1/g, (_, quote, path) => {
          const absolute = new URL(path, originalURL).href;
          return `from ${quote}${urls.get(absolute) ?? absolute}${quote}`;
        });
        const blob = URL.createObjectURL(new Blob([rewritten + '\n//# sourceURL=assembly-before/' + original + '\n'], { type: 'text/javascript' }));
        blobs.push(blob); urls.set(originalURL, blob);
      }
      const currentBody = await import(moduleURL('src/euler/streamtube-coupled.js'));
      const oldBody = await import(urls.get(moduleURL('src/euler/streamtube-coupled.js')));
      const currentDriver = await import(moduleURL('src/euler/streamtube-coupled-ises.js'));
      const oldDriver = await import(urls.get(moduleURL('src/euler/streamtube-coupled-ises.js')));
      const before = JSON.stringify(checkpoint), f = structuredClone(checkpoint.restart), c = checkpoint.continuation;
      f.input.streamwiseMode = 'hybrid'; f.input.hybrid = { epsilonP: 1e-5, ismom: 4 };
      f.input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
      f.options.blThermodynamics = 'historical-common-isentrope';
      const selected = JSON.stringify(f), build = create => create(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
      const a = build(currentBody.createCoupledStreamtubeBody), b = build(oldBody.createCoupledStreamtubeBody);
      const equal = (x, y) => x.length === y.length && x.every((v, i) => v === y[i]);
      const difference = (x, y) => {
        if (x.length !== y.length || !x.every(Number.isFinite) || !y.every(Number.isFinite)) throw new Error('Nonfinite or incompatible result comparison.');
        return Math.max(0, ...x.map((v, i) => Math.abs(v - y[i])));
      };
      const started = performance.now(), ja = a.jacobian(a.initial), jb = b.jacobian(b.initial);
      const matrix = { unknowns: a.n, entries: ja.values.length, sameInitial: equal(a.initial, b.initial),
        sameRowPointers: equal(ja.rowPtr, jb.rowPtr), sameColumnIndices: equal(ja.colIndex, jb.colIndex),
        maximumValueDifference: difference(ja.values, jb.values),
        sameResidual: equal(a.residual(a.initial), b.residual(b.initial)) };
      const options = { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL,
        iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance, stagnationLimiter: c.stagnationLimiter,
        ...(c.blUpdate === undefined ? {} : { blUpdate: c.blUpdate }),
        ...(c.projectionGeometry === undefined ? {} : { projectionGeometry: c.projectionGeometry }),
        maxIterations: 1, tolerance: 1e-10, linearOrdering: 'station-auto' };
      const current = currentDriver.solveCoupledStreamtubeIses(f.input, options);
      const archived = oldDriver.solveCoupledStreamtubeIses(f.input, options);
      const nodes = r => r.flow.nodes.flat(2).flatMap(p => [p.x, p.y]);
      const cp = r => r.flow.surfaces.flatMap(s => s.points.map(p => p.cp));
      const summary = r => ({ unknowns: r.x.length, families: r.families, conditions: r.conditions, phase: r.boundaryLayer.transitionState,
        ismom: r.solverInput.hybrid.ismom, surfaceCount: r.boundaryLayer.surfaces.length, wakeCount: r.boundaryLayer.wakes.length,
        pressureCount: cp(r).length, quality: r.mesh.quality, iterations: r.history.length - 1,
        step: r.history[1]?.step, backtracks: r.history[1]?.backtracks, linearOrdering: r.linearOrdering,
        linearDiagnostics: r.linearDiagnostics, initialRedistribution: r.initialRedistribution });
      const answer = { matrix, current: summary(current), archived: summary(archived),
        maximumStateDifference: difference(current.x, archived.x), maximumResidualDifference: difference(current.residual, archived.residual),
        maximumNodeDifference: difference(nodes(current), nodes(archived)), maximumPressureDifference: difference(cp(current), cp(archived)),
        samePhysics: JSON.stringify(current.solverInput) === JSON.stringify(archived.solverInput)
          && JSON.stringify(current.coupledOptions) === JSON.stringify(archived.coupledOptions),
        unchangedSelectedInput: JSON.stringify(f) === selected, unchangedFixture: before === JSON.stringify(checkpoint),
        fixtureConvertedExplicitlyToIsmom4: JSON.stringify(checkpoint.restart) !== selected,
        elapsedMilliseconds: performance.now() - started,
        wasmResources: performance.getEntriesByType('resource').filter(r => r.name.endsWith('/klu.wasm')).map(r => r.name) };
      blobs.forEach(URL.revokeObjectURL); return answer;
    };
    const source = `self.onmessage=async({data})=>{try{self.postMessage(await (${run.toString()})(data));}catch(error){self.postMessage({error:error.stack??error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({ ...payload, origin: location.href });
  }), { checkpoint: fixture, entries });
  expect(errors).toEqual([]); expect(result.wasmResources).toHaveLength(1);
  expect(result.matrix).toMatchObject({ unknowns: 349, sameInitial: true, sameRowPointers: true,
    sameColumnIndices: true, maximumValueDifference: 0, sameResidual: true });
  for (const r of [result.current, result.archived]) {
    expect(r.unknowns).toBe(349); expect(r.iterations).toBe(1); expect(r.ismom).toBe(4);
    expect(r.surfaceCount).toBe(4); expect(r.wakeCount).toBe(2); expect(r.pressureCount).toBe(20);
    expect(r.quality.valid).toBe(true); expect(r.linearOrdering).toBe('station-auto');
    expect(r.linearDiagnostics.solves).toBe(1); expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
    expect(r.initialRedistribution.accepted).toBe(true); expect(r.initialRedistribution.resumed).not.toBe(true);
  }
  expect(result.current.phase).toEqual(result.archived.phase); expect(result.current.step).toBe(result.archived.step);
  expect(result.current.backtracks).toBe(result.archived.backtracks);
  for (const key of ['maximumStateDifference', 'maximumResidualDifference', 'maximumNodeDifference', 'maximumPressureDifference']) expect(result[key]).toBe(0);
  expect(result.samePhysics).toBe(true); expect(result.unchangedSelectedInput).toBe(true); expect(result.unchangedFixture).toBe(true);
  expect(result.fixtureConvertedExplicitlyToIsmom4).toBe(true);
  expect(changedSources(sources)).toEqual([]);
  fs.writeFileSync(receipt, JSON.stringify({ passed: true, date: new Date().toISOString(), sourceHashes: sources,
    fixture: { path: fixturePath, sha256: sha256(fixturePath) }, physicalAcceptance: false,
    scope: 'Real Chromium worker and shipped KLU WASM: explicit ISMOM4 supplied349 state, one complete matrix and one Newton update using current versus archived assembly. Current Newton policy is identical. No mesh startup, full convergence, transonic or experimental-accuracy claim.',
    ...result }, null, 2) + '\n', { flag: 'wx' });
  console.log(JSON.stringify({ receipt, matrix: result.matrix, elapsedMilliseconds: result.elapsedMilliseconds,
    maximumStateDifference: result.maximumStateDifference, maximumPressureDifference: result.maximumPressureDifference }));
});
