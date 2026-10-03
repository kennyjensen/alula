import test from 'node:test';
import assert from 'node:assert/strict';
import { reconcileBodyEndpointStations } from '../src/geometry/body-endpoint-stations.js';

const component = (n, length, firstIndex, split = false) => {
  const positions = Array.from({ length: n + 1 }, (_, i) => .5 * length * (1 - Math.cos(Math.PI * i / n)));
  return { positions, anchors: [0, ...(split ? [n / 2] : []), n].map(i => ({ index: firstIndex + i, position: positions[i] })) };
};
const input = () => ({ upper: component(24, 1, 16, true), lower: component(24, .8, 16),
  upstream: component(16, 4, 0), wake: component(32, 5, 40) });

test('the three incident lengths match at each body end while actual hits and internal intervals remain fixed', () => {
  const source = input(), before = structuredClone(source), r = reconcileBodyEndpointStations(source);
  assert.deepEqual(source, before);
  const le = [r.upper.blocks[0].firstSpacing, r.lower.blocks[0].firstSpacing, r.upstream.blocks.at(-1).lastSpacing];
  const te = [r.upper.blocks.at(-1).lastSpacing, r.lower.blocks.at(-1).lastSpacing, r.wake.blocks[0].firstSpacing];
  for (const [j, achieved] of [le, te].entries()) {
    const join = r.joins[j], logs = join.natural.map(Math.log), optimum = logs.reduce((a, b) => a + b, 0) / 3;
    assert.ok(Math.abs(Math.log(join.target) - optimum) < 1e-14);
    achieved.forEach(h => assert.ok(Math.abs(h / join.target - 1) < 1e-12));
    const cost = z => logs.reduce((sum, log) => sum + (z - log) ** 2, 0);
    assert.ok(cost(optimum) < cost(optimum + .01) && cost(optimum) < cost(optimum - .01));
  }
  for (const name of Object.keys(source)) {
    const c = source[name], f = r[name];
    assert.equal(f.positions.length, c.positions.length);
    c.anchors.forEach(a => assert.equal(f.positions[a.index - c.anchors[0].index], a.position));
  }
  assert.ok(Math.abs(r.upper.blocks[0].lastSpacing - (source.upper.positions[12] - source.upper.positions[11])) < 1e-14);
  assert.ok(Math.abs(r.upper.blocks[1].firstSpacing - (source.upper.positions[13] - source.upper.positions[12])) < 1e-14);
  assert.equal(r.facingReconciled, false); assert.equal(r.joinDerivativeMatched, false);
});

test('body endpoint reconciliation is invariant under a change of physical units', () => {
  const source = input(), a = reconcileBodyEndpointStations(source), scale = 3.7;
  const b = reconcileBodyEndpointStations(Object.fromEntries(Object.entries(source).map(([name, c]) => [name,
    { positions: c.positions.map(s => scale * s), anchors: c.anchors.map(p => ({ ...p, position: scale * p.position })) }])));
  for (const name of Object.keys(source)) a[name].positions.forEach((s, i) => assert.ok(Math.abs(scale * s - b[name].positions[i]) < 1e-12));
});

test('unattainable endpoint targets and invalid anchors fail rather than moving hits or adding silent refinement', () => {
  const source = input();
  assert.throws(() => reconcileBodyEndpointStations({ ...source, upper: component(3, .000001, 16) }), /upper endpoint fit/);
  const invalid = input(); invalid.upper.anchors[1].position += .001;
  assert.throws(() => reconcileBodyEndpointStations(invalid), /exactly/);
});
