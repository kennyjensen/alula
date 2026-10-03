import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('browser worker solves automatic multielement transition and replays its active intervals', async ({ page }) => {
  const sourceHashes = numericalSourceHashes(['tests/browser/quad-automatic-transition.spec.js', 'tests/fixtures/intrinsic-body.js']);
  await page.route('**/quad-transition-harness', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Quad transition test</title>' }));
  await page.goto('/quad-transition-harness');
  const results = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `
      import { intrinsicBodyFixture } from '${location.origin}/tests/fixtures/intrinsic-body.js';
      import { solveCoupledStreamtubeIses } from '${location.origin}/src/euler/streamtube-coupled-ises.js';
      import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '${location.origin}/src/euler/streamtube-coupled.js';
      try {
        const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
        const options = { edgeMatching: 'section-velocity', transitionMode: 'automatic', reynolds: 1e6, ncrit: 9 };
        const controls = { ...options, tolerance: 1e-10, maxIterations: 16, stepAcceptance: 'admissible' };
        const full = solveCoupledStreamtubeIses(input, controls);
        let checkpoint;
        try {
          solveCoupledStreamtubeIses(input, { ...controls, onCheckpoint: (c, details) => {
            if (details.history.at(-1).iteration === 3) {
              checkpoint = JSON.parse(JSON.stringify(c, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
              throw new Error('deliberate checkpoint interruption');
            }
          } });
        } catch (error) { if (error.message !== 'deliberate checkpoint interruption') throw error; }
        if (!checkpoint) throw new Error('Missing accepted automatic checkpoint');
        const zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: checkpoint, maxIterations: 0 });
        const resumed = solveCoupledStreamtubeIses(undefined, { ...controls, resume: checkpoint, maxIterations: 13 });
        const mixed = solveCoupledStreamtubeBody(createCoupledStreamtubeBody(input, { ...options, reynolds: 3e5, ncrit: 14 }),
          { maxIterations: 16, tolerance: 1e-10 });
        const compact = r => ({ converged: r.converged, reason: r.reason, families: r.families,
          surfaces: r.boundaryLayer.surfaces.length, wakes: r.boundaryLayer.wakes.length,
          transitions: r.boundaryLayer.transitions, transitionState: r.boundaryLayer.transitionState,
          qualityValid: r.mesh.quality.valid, iterations: r.history.length - 1 });
        self.postMessage({ full: compact(full), mixed: compact(mixed), resumed: compact(resumed),
          replayFamiliesExact: JSON.stringify(zero.families) === JSON.stringify(checkpoint.families),
          replayPhasesExact: JSON.stringify(zero.boundaryLayer.transitionState) === JSON.stringify(checkpoint.restart.options.transitionState),
          maximumResumeDifference: Math.max(...resumed.x.map((v, i) => Math.abs(v - full.x[i]))) });
      } catch (error) { self.postMessage({ error: error.stack }); }
    `;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const finish = () => { clearTimeout(timer); worker.terminate(); URL.revokeObjectURL(url); };
    const timer = setTimeout(() => { finish(); reject(new Error('Automatic transition worker timed out')); }, 35000);
    worker.onerror = event => { finish(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { finish(); data.error ? reject(new Error(data.error)) : resolve(data); };
  }));
  for (const r of [results.full, results.mixed, results.resumed]) {
    expect(r.converged, r.reason).toBe(true); expect(r.qualityValid).toBe(true);
    expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
    expect(r.surfaces).toBe(4); expect(r.wakes).toBe(2); expect(r.transitions.length).toBe(4);
  }
  expect(results.full.transitions.every(t => t.kind === 'natural')).toBe(true);
  expect(results.mixed.transitions.filter(t => t.kind === 'trailing-edge').length).toBe(3);
  expect(results.replayFamiliesExact).toBe(true); expect(results.replayPhasesExact).toBe(true);
  expect(results.maximumResumeDifference).toBeLessThan(1e-11);
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-quad-automatic-transition-browser.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes,
    scope: 'Small two-element coupled API solves in a real Chromium worker: four natural transitions, mixed natural/laminar-to-TE surfaces and accepted-checkpoint resume. This does not enable the GUI or establish independent physical accuracy.',
    physicalAcceptance: false, results }, null, 2) + '\n');
  console.log(JSON.stringify(results));
});
