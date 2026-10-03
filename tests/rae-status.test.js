import test from 'node:test';
import assert from 'node:assert/strict';
import { verifiedRaeReceipt, renderRaeStatus } from '../scripts/validation/rae-status.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { resolveStreamtubeStartup } from '../src/euler/streamtube-result.js';

const { caseData: input } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
  changes: { mach: .74, alpha: 2.68, gridIntervals: 64, gridTubes: 9, eulerStartup: 'standard' } });

test('standard startup chooses the fine transonic method without changing conditions', () => {
  const before = structuredClone(input);
  for (const eulerStartup of [undefined, 'standard', 'harmonic-shock'])
    assert.equal(resolveStreamtubeStartup({ ...input, eulerStartup }), 'harmonic-shock');
  for (const changes of [{ mach: .2 }, { gridIntervals: 32 },
    { eulerIsmom: undefined }, { eulerIsmom: 2 }, { gridEllipticSmoothing: false }])
    assert.equal(resolveStreamtubeStartup({ ...input, ...changes }), 'standard');
  assert.equal(resolveStreamtubeStartup({ ...input, quadBoundaryLayers: true }), 'harmonic-shock');
  assert.throws(() => resolveStreamtubeStartup({ ...input, eulerStartup: 'typo' }), /Unknown/);
  assert.deepEqual(input, before);
});

test('status cannot promote an experiment, changed settings, stale source or provisional equations', () => {
  const receipt = { input, sourceHash: 'current', route: 'browser-form', browser: 'test', recordedAt: 'test',
    result: { status: 'research-converged', mach: .74, alpha: 2.68, gridValid: true, residual: 1e-11,
      upwind: { mucon: 1, mcrit: .99 } } };
  assert.equal(verifiedRaeReceipt(receipt, 'current'), true);
  for (const change of [r => { r.sourceHash = 'old'; }, r => { r.route = 'experiment'; },
    r => { r.input = { ...r.input, eulerStartup: 'harmonic-shock' }; },
    r => { r.input = { ...r.input, maxIterations: 1000 }; },
    r => { r.result = { ...r.result, residual: .1 }; },
    r => { r.result = { ...r.result, upwind: { mucon: -2, mcrit: .75 } }; },
    r => { r.result = { ...r.result, gridValid: false }; }]) {
    const altered = structuredClone(receipt); change(altered);
    assert.equal(verifiedRaeReceipt(altered, 'current'), false);
    assert.doesNotMatch(renderRaeStatus([altered], 'current'), /\| works in app \|/);
  }
});

test('status without local receipts makes no verification claims', () => {
  const report = renderRaeStatus([], 'current');
  assert.match(report, /unresolved/);
  assert.doesNotMatch(report, /\| works in app \|/);
});
