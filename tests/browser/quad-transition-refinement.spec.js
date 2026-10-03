import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('Chromium worker refines all four natural-transition intervals and reconverges the coupled system', async ({ page }) => {
  const sourceHashes = numericalSourceHashes(['tests/browser/quad-transition-refinement.spec.js', 'tests/fixtures/intrinsic-body.js']);
  await page.route('**/transition-refinement-harness', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Transition refinement</title>' }));
  await page.goto('/transition-refinement-harness');
  const result = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `
      import { intrinsicBodyFixture } from '${location.origin}/tests/fixtures/intrinsic-body.js';
      import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '${location.origin}/src/euler/streamtube-coupled.js';
      import { refineCoupledStreamtubeBody } from '${location.origin}/src/euler/streamtube-coupled-refinement.js';
      try {
        const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
        const parent = createCoupledStreamtubeBody(input, { transitionMode: 'automatic', edgeMatching: 'section-velocity' });
        const root = solveCoupledStreamtubeBody(parent, { tolerance: 1e-10, maxIterations: 12 });
        if (!root.converged) throw new Error(root.reason);
        const before = parent.evaluate(root.x), counts = Array(parent.euler.layout.nx).fill(1);
        for (const q of before.layers.transitions) counts[parent.bl.stations[q.id].i - 1] = 4;
        const mapped = refineCoupledStreamtubeBody(input, parent, { initial: root.x, streamwiseSubdivisions: counts, normalFactor: 1 });
        const child = mapped.system, seed = child.evaluate(child.initial);
        const solved = solveCoupledStreamtubeBody(child, { tolerance: 1e-10, maxIterations: 16 });
        const replay = createCoupledStreamtubeBody(mapped.input, { ...mapped.options, initialEuler: solved.flow,
          initialBL: solved.x.slice(child.ne), transitionState: solved.boundaryLayer.transitionState });
        self.postMessage({ converged: solved.converged, reason: solved.reason, families: solved.families,
          replayFamilies: replay.evaluate(replay.initial).families, iterations: solved.history.length - 1,
          unknowns: child.n, surfaces: solved.boundaryLayer.surfaces.length, wakes: solved.boundaryLayer.wakes.length,
          transitions: solved.boundaryLayer.transitions, qualityValid: solved.mesh.quality.valid,
          allocationExact: JSON.stringify(before.outer.allocation) === JSON.stringify(seed.outer.allocation) });
      } catch (error) { self.postMessage({ error: error.stack }); }
    `;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const finish = () => { clearTimeout(timer); worker.terminate(); URL.revokeObjectURL(url); };
    const timer = setTimeout(() => { finish(); reject(new Error('Transition refinement worker timed out')); }, 35000);
    worker.onerror = e => { finish(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { finish(); data.error ? reject(new Error(data.error)) : resolve(data); };
  }));
  expect(result.converged, result.reason).toBe(true); expect(result.qualityValid).toBe(true);
  expect(result.allocationExact).toBe(true); expect(result.surfaces).toBe(4); expect(result.wakes).toBe(2);
  expect(result.transitions).toHaveLength(4);
  expect(Math.max(...Object.values(result.families))).toBeLessThan(1e-10);
  expect(Math.max(...Object.values(result.replayFamilies))).toBeLessThan(1e-10);
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-transition-refinement-browser.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes,
    physicalAcceptance: false, scope: 'Small multielement automatic-transition streamwise refinement, full coupled solve and restart in Chromium; no independent physical acceptance.', result }) + '\n');
  console.log(JSON.stringify(result));
});
