import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

test('browser public initializer recovers the attached natural case and preserves the default seed', async ({ page }) => {
  test.setTimeout(150000);
  const sourceHashes = numericalSourceHashes(['tests/browser/quad-iset-wake-startup.spec.js']);
  const parents = ['docs/current-natural-attached-cold.json', 'docs/current-quad-automatic-mrchue-public.json', 'docs/current-iset-wake-integration.json']
    .map(path => ({ path, sha256: sha256(path) }));
  await page.route('**/quad-iset-wake-harness', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Wake startup check</title>' }));
  await page.goto('/quad-iset-wake-harness');
  const report = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `
      import { createStreamtubeBodySystem } from '${location.origin}/src/euler/streamtube-body.js';
      import { transferStreamtubeGeometry } from '${location.origin}/src/euler/streamtube-geometry.js';
      import { initializeCoupledStreamtubeBody } from '${location.origin}/src/euler/streamtube-coupled-initializer.js';
      import { solveCoupledStreamtubeIses } from '${location.origin}/src/euler/streamtube-coupled-ises.js';
      try {
        const expected = await (await fetch('${location.origin}/docs/current-iset-wake-integration.json')).json(), cases = [];
        let attached;
        for (const [name, path] of [['attached','current-natural-attached-cold.json'], ['default','current-quad-automatic-mrchue-public.json']]) {
          const saved = await (await fetch('${location.origin}/docs/' + path)).json();
          const p = saved.stages.euler.result, f = saved.stages['boundary-layer-initialization'].restart;
          const source = createStreamtubeBodySystem(p.solverInput), sx = source.adoptGeometry(Float64Array.from(p.flow.x), p.flow.nodes);
          const input = f.input, options = { ...f.options }; delete options.transitionState;
          const displacement = { surfaces: input.bodies.map(b => ({ upper: Array(b.trailingIndex-b.leadingIndex+1).fill(0), lower: Array(b.trailingIndex-b.leadingIndex+1).fill(0) })),
            wakes: input.bodies.map(b => Array(input.outerLower.length-1-b.trailingIndex).fill(0)) };
          const target = createStreamtubeBodySystem({ ...input, displacement }), tx = transferStreamtubeGeometry(source, sx, target);
          const prepared = initializeCoupledStreamtubeBody(input, { ...options, initialEuler: { x: tx, ...target.decode(tx) } });
          const s = prepared.system, v = s.evaluate(s.initial), reference = expected.cases.find(c => c.name === name).seed;
          const ref = [...reference.initialEuler.x, ...reference.initialBL];
          const stateDifference = Math.max(...s.initial.map((v,i) => Math.abs(v-ref[i])/Math.max(1,Math.abs(ref[i]))));
          cases.push({ name, initialization: prepared.initialization, families: v.families, stateDifference, quality: prepared.mesh.quality,
            surfaces: s.bl.surfaces.length, wakes: s.bl.wakes.length });
          if (name === 'attached') attached = { input, options: { ...options, transitionState: s.bl.snapshotActive() },
            initialEuler: { x: s.initial.slice(0,s.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes }, initialBL: s.initial.slice(s.ne) };
        }
        const r = solveCoupledStreamtubeIses(attached.input, { ...attached.options, initialEuler: attached.initialEuler, initialBL: attached.initialBL,
          maxIterations: 24, tolerance: 1e-10, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' });
        self.postMessage({ cases, result: { converged: r.converged, reason: r.reason, families: r.families,
          quality: r.mesh.quality, iterations: r.history.length-1, history: r.history, surfaces: r.boundaryLayer.surfaces.length,
          wakes: r.boundaryLayer.wakes.length, transitions: r.boundaryLayer.transitions, linearDiagnostics: r.linearDiagnostics },
          restart: r.checkpoint.restart, checkpoint: r.checkpoint });
      } catch (error) { self.postMessage({ error: error.stack }); }
    `;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const finish = () => { clearTimeout(timer); worker.terminate(); URL.revokeObjectURL(url); };
    const timer = setTimeout(() => { finish(); reject(new Error('Wake startup worker timed out')); }, 120000);
    worker.onerror = event => { finish(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { finish(); data.error ? reject(new Error(data.error)) : resolve(data); };
  }));
  expect(changedSources(sourceHashes)).toEqual([]);
  for (const p of parents) expect(sha256(p.path)).toBe(p.sha256);
  // Retain the endpoint before assertions, including failed physical/numerical checks.
  fs.writeFileSync('docs/current-iset-wake-browser.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes, parents,
    physicalAcceptance: false, inProgress: false,
    scope: 'Chromium worker recomputes public attached/default BL initialization from saved converged Euler precursors, then solves the full attached natural-transition system. No new inviscid solve or independent physical acceptance.',
    ...report }, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v) + '\n');
  const [attached, defaults] = report.cases;
  expect(attached.initialization.wakeInitialization).toBe('iset-linear-shape');
  expect(attached.initialization.thicknessFactor).toBe(1);
  expect(defaults.initialization.history.length).toBe(1); expect(defaults.surfaces).toBe(4); expect(defaults.wakes).toBe(2);
  for (const c of report.cases) { expect(c.quality.valid).toBe(true); expect(c.stateDifference).toBeLessThan(1e-9); }
  expect(report.result.converged, report.result.reason).toBe(true); expect(report.result.quality.valid).toBe(true);
  expect(Math.max(...Object.values(report.result.families))).toBeLessThan(1e-10);
  expect(report.result.linearDiagnostics.maxRelativeResidual).toBeLessThan(1e-10);
  expect(report.result.transitions.length).toBe(2); expect(report.result.transitions.every(t => t.kind === 'natural')).toBe(true);
  console.log(JSON.stringify({ iterations: report.result.iterations, families: report.result.families, transitions: report.result.transitions,
    seedDifferences: report.cases.map(c => ({ name: c.name, difference: c.stateDifference })) }));
});
