import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledNcrit, quadCoupledNcritLabel, quadCoupledNcritResult } from '../src/ui/quad-coupled-ncrit.js';

test('actual Ncrit comes from the displayed checkpoint; requested settings cannot relabel it', () => {
  const frame = { checkpoint: { restart: { options: { ncrit: 7.5 } } },
    actualNcrit: 9, targetNcrit: 8, conditions: { ncrit: 8 }, ncritContinuation: { actualNcrit: 4 } };
  const before = structuredClone(frame), point = quadCoupledNcrit(frame, 9);
  assert.deepEqual(point, { actualNcrit: 7.5, targetNcrit: 9 });
  assert.equal(quadCoupledNcritLabel(point), ' · Ncrit 7.5 → target 9');
  assert.deepEqual(frame, before);
  assert.deepEqual(quadCoupledNcrit({ targetNcrit: 9 }), { targetNcrit: 9 });
  assert.equal(quadCoupledNcritLabel({ targetNcrit: 9 }), '');
});

test('legacy same-condition results remain identical and intermediate roots remain unsuccessful', () => {
  const request = { ncrit: 9, mach: .185, reynolds: 2.7e6 }, legacy = { conditions: { ncrit: 9, mach: .185 }, converged: true };
  assert.equal(quadCoupledNcritResult(legacy, request), legacy);
  const checkpoint = { restart: { options: { ncrit: 8 } } }, raw = { ...legacy, checkpoint,
    conditions: { ncrit: 8, mach: .185 }, mesh: { quality: { valid: true } },
    ncritContinuation: { actualNcrit: 8, targetNcrit: 9, reachedTarget: false },
    sourceCase: { ...request }, cl: .9, cd: .01, cm: -.2 };
  const before = structuredClone(raw), result = quadCoupledNcritResult(raw, request);
  assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
  assert.equal(result.actualNcrit, 8); assert.equal(result.targetNcrit, 9);
  assert.equal(result.sourceCase.ncrit, 8); assert.deepEqual(result.requestedCase, request);
  assert.equal(result.checkpoint, checkpoint); assert.equal(result.cl, raw.cl); assert.equal(result.cd, raw.cd);
  assert.deepEqual(raw, before);
  const final = quadCoupledNcritResult({ ...raw, checkpoint: { restart: { options: { ncrit: 9 } } },
    conditions: { ncrit: 9, mach: .185 }, ncritContinuation: { actualNcrit: 9, targetNcrit: 9, reachedTarget: true } }, request);
  assert.equal(final.converged, true); assert.equal(final.sourceCase.ncrit, 9);
  assert.equal(quadCoupledNcritLabel(final), ' · Ncrit 9');
});

test('same numerical threshold does not promote a failed target initialization to success', () => {
  const result = quadCoupledNcritResult({ converged: false, actualNcrit: 9, targetNcrit: 9,
    ncritContinuation: { actualNcrit: 9, targetNcrit: 9, reachedTarget: false },
    mesh: { quality: { valid: true } } }, { ncrit: 9 });
  assert.equal(result.converged, false); assert.equal(result.stateConverged, false);
  assert.equal(result.status, 'research-coupled-target-not-reached');
});
