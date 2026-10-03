import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createOuterStreamtubeSpacing } from '../src/geometry/outer-streamtube-spacing.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/outer-streamtube-spacing.json', import.meta.url)));
const maxRatio = values => {
  const h = values.slice(1).map((v, i) => v - values[i]);
  return Math.max(...h.slice(1).map((v, i) => Math.max(v / h[i], h[i] / v)));
};

test('PCHIP improves the saved outer transition and preserves anchors, while the isolated old skeleton still fails its screen', () => {
  const before = structuredClone(fixture);
  for (const [side, targets] of fixture.targets.entries()) {
    const args = { ...fixture, targets }, current = createOuterStreamtubeSpacing(args);
    const previous = createOuterStreamtubeSpacing({ ...args, derivatives: 'harmonic' });
    if (!side) assert.ok(maxRatio(previous.potentials) > 1.65, 'The fixture must retain the reported boundary defect');
    assert.ok(maxRatio(current.potentials) < maxRatio(previous.potentials));
    if (!side) assert.ok(maxRatio(current.potentials) > 1.5, 'Interpolation alone must not be reported as a completed spacing correction');
    for (const rank of fixture.anchors) {
      const i = fixture.ranks.indexOf(rank);
      assert.equal(current.potentials[i], previous.potentials[i]);
    }
    for (const { rank, potential } of targets) assert.equal(current.potentials[fixture.ranks.indexOf(rank)], potential);
    assert.equal(current.diagnostics.indexInterpolation, 'pchip');
    assert.equal(current.diagnostics.exactModernMsetX, false);
    const shifted = createOuterStreamtubeSpacing({ ...args,
      targets: targets.map(p => ({ rank: p.rank, potential: 7 + 3 * p.potential })) });
    current.potentials.forEach((p, i) => assert.ok(Math.abs(shifted.potentials[i] - (7 + 3 * p)) < 2e-14));
  }
  assert.deepEqual(fixture, before);
});

test('outer spacing retains the chosen blend and rejects missing or nonmonotone block data', () => {
  const args = { ranks: [0, .25, .75, 1, 1.5, 2], anchors: [0, 1, 2],
    targets: [{ rank: 0, potential: 2 }, { rank: 1, potential: 3 }, { rank: 2, potential: 5 }] };
  for (const spread of [0, .4, 1]) {
    const result = createOuterStreamtubeSpacing({ ...args, spread });
    result.potentials.forEach((p, i) => assert.ok(Math.abs(p - ((1 - spread) * result.rankMap.value(args.ranks[i])
      + spread * result.indexMap.value(i))) < 1e-14));
    assert.ok(result.potentials.every((v, i, row) => !i || v > row[i - 1]));
  }
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, ranks: [0, .25, .75, 1.5, 2] }), /anchor/);
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, anchors: [0, 2] }), /target/);
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, spread: 1.01 }), /Invalid/);
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, derivatives: 'unknown' }), /Invalid/);
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, endSpacing: 'unknown' }), /Invalid/);
  assert.throws(() => createOuterStreamtubeSpacing({ ...args, targets: [{ rank: 0, potential: 2 }, { rank: 2, potential: 1 }] }), /Invalid/);
});

test('Vinokur end segments preserve all anchors and interior correspondence with C1 joins and zero far-end curvature', () => {
  for (const targets of fixture.targets) {
    const args = { ...fixture, targets }, cubic = createOuterStreamtubeSpacing(args);
    const fitted = createOuterStreamtubeSpacing({ ...args, endSpacing: 'vinokur' });
    const { indices } = fitted.diagnostics;
    assert.equal(fitted.diagnostics.endSpacing, 'vinokur');
    assert.equal(fitted.diagnostics.exactModernMsetX, false);
    for (const i of indices) assert.equal(fitted.potentials[i], cubic.potentials[i]);
    for (let i = indices[1]; i <= indices.at(-2); i++)
      assert.equal(fitted.potentials[i], cubic.potentials[i]);
    for (const [side, segment] of fitted.indexMap.endSegments.entries()) {
      const i = side === 0 ? indices[1] : indices.at(-2);
      assert.equal((side === 0 ? segment.right : segment.left).derivative, cubic.indexMap.slopes[indices.indexOf(i)]);
      assert.ok((side === 0 ? segment.left : segment.right).secondDerivative === 0);
    }
    const shifted = createOuterStreamtubeSpacing({ ...args, endSpacing: 'vinokur',
      targets: targets.map(p => ({ rank: p.rank, potential: 7 + 3 * p.potential })) });
    fitted.potentials.forEach((p, i) => assert.ok(Math.abs(shifted.potentials[i] - (7 + 3 * p)) < 3e-14));
    assert.ok(maxRatio(fitted.potentials) <= 1.5, 'The saved scalar spacing defect must pass its unchanged screen');
  }
});

test('Vinokur outer ends support a shared middle anchor, a single affine interval and all blend values', () => {
  for (const anchors of [[0, 1], [0, .5, 1]]) for (const spread of [0, .4, 1]) {
    const args = { ranks: [0, .25, .5, .75, 1], anchors,
      targets: anchors.map(rank => ({ rank, potential: 2 + rank + rank * rank })), spread, endSpacing: 'vinokur' };
    const result = createOuterStreamtubeSpacing(args);
    result.potentials.forEach((p, i) => assert.ok(Math.abs(p - ((1 - spread) * result.rankMap.value(args.ranks[i])
      + spread * result.indexMap.value(i))) < 1e-14));
    assert.ok(result.potentials.every((p, i, row) => !i || p > row[i - 1]));
    for (const bad of [-1, 5, NaN]) assert.throws(() => result.indexMap.evaluate(bad), /outside/);
  }
});
