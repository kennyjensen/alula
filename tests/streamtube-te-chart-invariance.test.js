import test from 'node:test';
import assert from 'node:assert/strict';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody, coupledStreamtubeResult } from '../src/euler/streamtube-coupled.js';
import { transferStreamtubeGeometry } from '../src/euler/streamtube-geometry.js';

const closeArrays = (a, b, tolerance = 1e-11) => {
  assert.equal(a.length, b.length);
  a.forEach((v, i) => assert.ok(Math.abs(v - b[i]) <= tolerance * Math.max(1, Math.abs(v), Math.abs(b[i])), `${i}: ${v} != ${b[i]}`));
};
const positions = nodes => nodes.flat(2).flatMap(p => [p.x, p.y]);

test('changing only the TE coordinate chart preserves physical Euler, Cp, BL arcs and both finite-base gap models', () => {
  const f = twoActiveFiniteBaseWakes(), original = f.system.evaluate(f.x);
  let source = f.system, x = f.x;
  for (const mode of ['te-center', 'fixed']) {
    const input = { ...f.input, wakeDisplacementMotion: mode };
    const euler = createStreamtubeBodySystem({ ...input, displacement: source.euler.displacement });
    const state = transferStreamtubeGeometry(source.euler, x.subarray(0, source.ne), euler);
    const transferred = euler.decode(state);
    const target = createCoupledStreamtubeBody(input, { reynolds: source.conditions.reynolds,
      ncrit: source.conditions.ncrit, edgeMatching: source.conditions.edgeMatching,
      initialBL: x.subarray(source.ne), initialEuler: { x: state, nodes: transferred.nodes,
        undisplacedNodes: transferred.undisplacedNodes } });
    const value = target.evaluate(target.initial);
    closeArrays(positions(value.outer.nodes), positions(original.outer.nodes), 2e-14);
    closeArrays(value.residual, original.residual);
    assert.deepEqual(value.outer.allocation, original.outer.allocation);
    closeArrays(value.layers.states.flatMap(s => [s.s, s.theta, s.deltaStar, s.ue, s.wakeGap ?? 0]),
      original.layers.states.flatMap(s => [s.s, s.theta, s.deltaStar, s.ue, s.wakeGap ?? 0]));
    const before = coupledStreamtubeResult(f.system, f.x), after = coupledStreamtubeResult(target, target.initial);
    closeArrays(before.flow.surfaces.flatMap(s => s.points.map(p => p.cp)), after.flow.surfaces.flatMap(s => s.points.map(p => p.cp)));
    source = target; x = target.initial;
  }
});
