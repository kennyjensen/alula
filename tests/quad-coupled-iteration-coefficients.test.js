import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledIterationCoefficients } from '../src/ui/quad-coupled-iteration-coefficients.js';
import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';

// Frozen independently captured GUI coefficients include the nontrivial
// flap/main passage ordering and the kernel/reference length conversion.
const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
const r = saved.result;
const input = { checkpoint: saved.checkpoint, flow: r.flow,
  bl: { ...r.numericalBoundaryLayer, scale: 1 / Math.sqrt(r.kernelReynolds) },
  bodies: saved.checkpoint.restart.input.bodies, solverLength: r.solverLength,
  referenceChord: r.referenceChord, momentReference: r.momentReference, alpha: r.alpha, mach: r.mach };

test('live coefficients match the captured public result without changing the accepted state', () => {
  const before = JSON.stringify(input.checkpoint);
  assert.deepEqual(quadCoupledIterationCoefficients(input), r.coefficients);
  assert.equal(JSON.stringify(input.checkpoint), before);
});

test('compact startup and recovery frames retain the captured coefficients exactly', () => {
  const compact = structuredClone({ ...input, flow: observableFlow(input.flow), bl: observableBL(input.bl) });
  assert.deepEqual(quadCoupledIterationCoefficients(compact), r.coefficients);
});

test('live coefficients are independent of the physical length unit', () => {
  const factor = 3;
  const actual = quadCoupledIterationCoefficients({ ...input,
    solverLength: factor * input.solverLength, referenceChord: factor * input.referenceChord,
    momentReference: { x: factor * input.momentReference.x, y: factor * input.momentReference.y },
    flow: { undisplacedNodes: input.flow.undisplacedNodes.map(g => g.map(row => row.map(p => ({ x: factor * p.x, y: factor * p.y })))) } });
  for (const key of ['cl', 'cd', 'cm', 'pressureIntegralDrag']) assert.ok(Math.abs(actual[key] - r.coefficients[key]) < 1e-12, key);
});

test('live coefficients require a complete versioned checkpoint', () => {
  assert.throws(() => quadCoupledIterationCoefficients({ ...input, checkpoint: { ...input.checkpoint, version: 2 } }));
  assert.throws(() => quadCoupledIterationCoefficients({ ...input, checkpoint: { ...input.checkpoint,
    restart: { ...input.checkpoint.restart, initialBL: [] } } }));
});
