import test from 'node:test';
import assert from 'node:assert/strict';
import { normalRefinementCheckpoint } from '../scripts/validation/normal-refinement-checkpoint.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
const serialize = x => JSON.parse(JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));

test('normal refinement preserves complete two-element maintenance history and reconverges all BLs and wakes', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const controls = { edgeMatching: 'section-velocity', transitionMode: 'automatic', maxIterations: 12,
    tolerance: 1e-10, stepAcceptance: 'admissible' };
  const original = solveCoupledStreamtubeIses(input, controls); assert.equal(original.converged, true, original.reason);
  const c = serialize(original.checkpoint), f = c.restart;
  const parent = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const mapped = refineCoupledStreamtubeBody(f.input, parent, { streamwiseFactor: 1, normalFactor: 2 });
  const s = mapped.system, v = s.evaluate(s.initial), restart = { input: mapped.input, options: mapped.options,
    initialEuler: { x: s.initial.slice(0, s.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes }, initialBL: s.initial.slice(s.ne) };
  const r = normalRefinementCheckpoint(c, restart, mapped.diagnostics.normalSubdivisions);
  const frozen = serialize(r.checkpoint);
  const zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: r.checkpoint, maxIterations: 0 });
  assert.equal(zero.initialRedistribution.resumed, true); assert.deepEqual(zero.initialRedistribution.passages, []);
  assert.deepEqual(zero.families, v.families); assert.deepEqual(zero.flow.nodes, v.outer.nodes);
  assert.deepEqual(serialize(zero.checkpoint.continuation), c.continuation);
  const result = solveCoupledStreamtubeIses(undefined, { ...controls, resume: r.checkpoint });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.boundaryLayer.surfaces.length, 4); assert.equal(result.boundaryLayer.wakes.length, 2);
  assert.ok(Math.max(...Object.values(result.families)) < 1e-10);
  assert.deepEqual(serialize(r.checkpoint), frozen); assert.deepEqual(serialize(original.checkpoint), c);
  const corrupt = structuredClone(frozen.restart); corrupt.input.alpha += 1;
  assert.throws(() => normalRefinementCheckpoint(c, corrupt, mapped.diagnostics.normalSubdivisions));
  t.diagnostic(JSON.stringify({ parentUnknowns: parent.n, unknowns: s.n, iterations: result.history.length - 1, families: result.families, transfer: r.diagnostics }));
});
