import test from 'node:test';
import assert from 'node:assert/strict';
import { directDisplacedAssemblyLedger, directLinearPressureChain } from './oracles/multielement-displaced-ledger.js';

const close = (a, b, message) => assert.ok(Math.abs(a - b) < 1e-12, `${message}: ${a} versus ${b}`);
const fixture = () => {
  const p = 2, gamma = 1.4, nx = 3, groups = 3, tubes = 2;
  return { nodes: Array.from({ length: groups }, (_, g) => Array.from({ length: nx + 1 }, (_, i) =>
    Array.from({ length: tubes + 1 }, (_, j) => ({ x: i, y: g + j / tubes })))),
  sections: Array.from({ length: nx }, () => Array.from({ length: groups }, () =>
    Array.from({ length: tubes }, () => ({ rho: 1, q: 1, p })))),
  cells: Array.from({ length: nx - 1 }, () => Array.from({ length: groups }, () =>
    Array.from({ length: tubes }, () => ({ interfacePressure: { lower: p, upper: p }, transportSpeeds: [1, 1] })))),
  massFlows: Array.from({ length: groups }, () => Array(tubes).fill(.5)),
  bodies: [{ element: 1, leadingIndex: 1, trailingIndex: 2 }, { element: 0, leadingIndex: 1, trailingIndex: 2 }],
  stagnationEnthalpy: gamma / (gamma - 1) * p + .5, gamma };
};

test('three passages telescope without double counting either shared cut', () => {
  const input = fixture(), before = structuredClone(input), value = directDisplacedAssemblyLedger(input);
  assert.equal(value.cutPairs.paired, 4); assert.equal(value.cutPairs.coincident, 4);
  assert.equal(value.cutPairs.noncoincident, 0); assert.equal(value.counts.volumes, 12);
  assert.ok(value.maximumAlgebraError < 1e-12);
  for (const v of value.whole.localPhysical) close(v, 0, 'constant flow local total');
  for (const v of value.cutPairs.retainedCoincidentFlux) close(v, 0, 'coincident face pair');
  assert.deepEqual(input, before);
});

test('internal pressure mismatch is retained exactly once, independently of exterior faces', () => {
  const input = fixture();
  for (const row of input.cells) { row[1][0].interfacePressure.upper += .1; row[1][1].interfacePressure.lower -= .2; }
  const value = directDisplacedAssemblyLedger(input);
  close(value.whole.internalPressureMismatch[2], .6, 'known outward upper-face traction mismatch over length2');
  close(value.whole.externalPhysical[2], 0, 'unchanged exterior traction');
  close(value.whole.localPhysical[2], .6, 'local defects include mismatch once');
  assert.ok(value.maximumAlgebraError < 1e-12);
});

test('noncoincident displaced wake banks remain a force term even with equal pressure', () => {
  const input = fixture(); input.nodes[1][3][0].y += .1;
  const value = directDisplacedAssemblyLedger(input);
  assert.equal(value.cutPairs.noncoincident, 1); assert.equal(value.cutPairs.coincident, 3);
  close(value.cutPairs.maximumPressureMismatch, 0, 'equal pressure');
  close(value.cutPairs.retainedNoncoincidentFlux[1], .1, 'P times unmatched dy=.05');
  close(value.whole.cuts[1], .1, 'noncoincident force retained in whole assembly');
  assert.ok(value.maximumAlgebraError < 1e-12);
});

test('closed finite-base solid pressure includes each retained base segment once', () => {
  const wetted = [{ x: 1, y: .1 }, { x: 0, y: 0 }, { x: 1, y: -.1 }];
  const base = [wetted.at(-1), { x: 1.02, y: -.04 }, { x: 1.02, y: .04 }, wetted[0]];
  const a = directLinearPressureChain(wetted, wetted.map(() => 2));
  const b = directLinearPressureChain(base, base.map(() => 2));
  const closed = [...wetted, ...base.slice(1)];
  const full = directLinearPressureChain(closed, closed.map(() => 2));
  assert.equal(b.segments, 3); assert.equal(full.segments, a.segments + b.segments);
  for (const key of ['x', 'y', 'noseUpMoment']) {
    close(a[key] + b[key], full[key], 'wetted + one base = closed contour');
    close(full[key], 0, 'uniform pressure on closed solid');
  }
  assert.ok(Math.abs(a.x + 2 * b.x) > .1, 'duplicate base would break force closure');
});
