import test from 'node:test';
import assert from 'node:assert/strict';
import { distributeSurfaceByDensity } from '../src/geometry/surface-density-stations.js';

const cosine = n => Array.from({ length: n + 1 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / n)));
test('unchanged count reproduces every requested station exactly', () => {
  const requested = cosine(16), anchors = [{ index: 4, position: 0 }, { index: 20, position: 1 }];
  const r = distributeSurfaceByDensity({ requested, anchors });
  assert.deepEqual(r.positions, requested); assert.deepEqual(r.coordinates, requested.map((_, i) => i));
  assert.equal(r.blocks[0].densityCoordinatePerInterval, 1);
});

test('a larger count refines the local distribution without borrowing any other surface monitor', () => {
  const requested = cosine(16), anchors = [{ index: 32, position: 0 }, { index: 78, position: 1 }];
  const r = distributeSurfaceByDensity({ requested, anchors });
  assert.equal(r.positions.length, 47);
  r.coordinates.forEach((q, i) => assert.ok(Math.abs(q - 16 * i / 46) < 1e-14));
  const ds = r.positions.slice(1).flatMap((s, i) => r.positions[i] >= .2 && s <= .8 ? [s - r.positions[i]] : []);
  assert.ok(Math.max(...ds) / Math.min(...ds) < 1.3);
});

test('physical hits are exact even between requested stations and each block reports its density scale', () => {
  const requested = cosine(16), anchors = [{ index: 32, position: 0 }, { index: 58, position: .887 }, { index: 78, position: 1 }];
  const r = distributeSurfaceByDensity({ requested, anchors });
  assert.equal(r.positions[26], .887); assert.equal(r.positions[0], 0); assert.equal(r.positions.at(-1), 1);
  for (const block of r.blocks) {
    const start = block.fromIndex - r.firstIndex;
    for (let k = 1; k <= block.intervals; k++) assert.ok(Math.abs(r.coordinates[start + k] - r.coordinates[start + k - 1]
      - block.densityCoordinatePerInterval) < 1e-13);
  }
  assert.notEqual(r.blocks[0].densityCoordinatePerInterval, r.blocks[1].densityCoordinatePerInterval);
  assert.equal(r.joinDerivativeMatched, false); assert.equal(r.facingReconciled, false);
});

test('physical unit/origin changes and logical index translation preserve the distribution', () => {
  const requested = cosine(16), anchors = [{ index: 3, position: 0 }, { index: 19, position: .81 }, { index: 35, position: 1 }];
  const a = distributeSurfaceByDensity({ requested, anchors });
  const b = distributeSurfaceByDensity({ requested: requested.map(s => 3 + 2.7 * s),
    anchors: anchors.map(p => ({ index: p.index + 17, position: 3 + 2.7 * p.position })) });
  a.positions.forEach((s, i) => assert.ok(Math.abs(b.positions[i] - (3 + 2.7 * s)) < 1e-13));
});

test('invalid hit constraints are rejected rather than sorted or moved', () => {
  const requested = cosine(16), anchors = [{ index: 0, position: 0 }, { index: 16, position: 1 }];
  for (const other of [[...anchors].reverse(), [{ index: 0, position: .1 }, anchors[1]],
    [anchors[0], { index: 4, position: .8 }, { index: 7, position: .7 }, anchors[1]],
    [anchors[0], { index: 0, position: 1 }]]) assert.throws(() => distributeSurfaceByDensity({ requested, anchors: other }), /Ordered/);
  assert.throws(() => distributeSurfaceByDensity({ requested: [0, .5, .5, 1], anchors }), /Ordered/);
});

test('endpoint fitting retains physical hits and matches physical intervals on a nonuniform requested density', () => {
  const requested = cosine(32), anchors = [{ index: 0, position: 0 }, { index: 29, position: .74 }, { index: 46, position: 1 }];
  const endpointSpacings = [{ firstSpacing: .002, lastSpacing: .013 }, { firstSpacing: .013, lastSpacing: .003 }];
  const r = distributeSurfaceByDensity({ requested, anchors, endpointSpacings });
  assert.equal(r.positions[29], .74); assert.equal(r.positions.at(-1), 1);
  for (let k = 0; k < 2; k++) {
    const block = r.blocks[k];
    assert.ok(Math.abs(block.firstSpacing - endpointSpacings[k].firstSpacing) < 3e-15);
    assert.ok(Math.abs(block.lastSpacing - endpointSpacings[k].lastSpacing) < 3e-15);
    assert.equal(block.endpointFit.uniformDensityIncrements, false);
  }
  assert.equal(r.joinDerivativeMatched, false); // Cell lengths are not derivatives.
  assert.equal(r.facingReconciled, false);
  assert.throws(() => distributeSurfaceByDensity({ requested, anchors, endpointSpacings: [] }), /pair/);
  assert.throws(() => distributeSurfaceByDensity({ requested, anchors,
    endpointSpacings: [{ firstSpacing: .5, lastSpacing: .5 }, endpointSpacings[1]] }), /pair/);
});
