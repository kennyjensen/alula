import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { fitConstrainedStationPositions } from '../src/numerics/constrained-station-positions.js';
import { projectPositionKkt } from '../src/numerics/position-kkt-projection.js';

const read = name => JSON.parse(fs.readFileSync(new URL(`./fixtures/${name}.json`, import.meta.url)));
const check = result => {
  const h = result.positions.slice(1).map((p, i) => p - result.positions[i]);
  assert.ok(h.every(v => v > 0));
  assert.ok(Math.max(...h.slice(1).map((v, i) => Math.max(v / h[i], h[i] / v))) <= 1.5 + 1e-10);
  assert.ok(result.forwardCertificate.physicalPositionErrorBound < 1e-10);
  assert.equal(result.forwardCertificate.maximumInactiveViolationBound, 0);
  assert.ok(result.forwardCertificate.minimumActiveMultiplierBound >= 0);
  assert.ok(Math.max(result.primalResidual, result.stationarityResidual, result.complementarity) < 1e-10);
  assert.ok(Math.max(result.rawStationarityBackwardError, result.rawComplementarityBackwardError) <= 256 * Number.EPSILON);
};

test('cove station positions agree with independent dense reorthogonalized QR', () => {
  const input = read('30p-cove-position-stations'), before = structuredClone(input);
  const r = fitConstrainedStationPositions(input), reference = read('30p-cove-position-dense-qr'); check(r);
  assert.ok(Math.max(...r.positions.map((v, i) => Math.abs(v - reference.positions[i]))) < 3e-14);
  assert.deepEqual(input, before);
  assert.equal(r.positions[1], input.firstSpacing);
  assert.equal(r.positions.at(-2), 1 - input.lastSpacing);
});

test('reflection exchanges endpoint requests without changing the physical optimum', () => {
  const input = read('30p-cove-position-stations'), r = fitConstrainedStationPositions(input);
  const reflected = fitConstrainedStationPositions({ ...input, positions: input.positions.toReversed().map(p => 1 - p),
    firstSpacing: input.lastSpacing, lastSpacing: input.firstSpacing }); check(reflected);
  assert.ok(Math.max(...r.positions.map((v, i) => Math.abs(1 - v - reflected.positions.at(-1 - i)))) < 3e-13);
});

test('physical arc translation and scaling preserve the normalized cove solution', () => {
  const input = read('30p-cove-position-stations'), r = fitConstrainedStationPositions(input);
  const length = 7, origin = 5, physical = input.positions.map(p => origin + length * p);
  const transformed = fitConstrainedStationPositions({ ...input,
    positions: physical.map(p => (p - origin) / length),
    firstSpacing: input.firstSpacing * length / length, lastSpacing: input.lastSpacing * length / length }); check(transformed);
  assert.ok(Math.max(...r.positions.map((p, i) => Math.abs(origin + length * p - (origin + length * transformed.positions[i])))) < 2e-12);
});

test('the coarse main cove has an independent physical position certificate despite large raw multipliers', () => {
  const r = fitConstrainedStationPositions(read('30p-main-position-stations')); check(r);
  assert.equal(r.certificate, 'normalized-objective KKT; original primal rows and physical growth');
  assert.ok(Number.isFinite(r.rawStationarityResidual)); assert.ok(Number.isFinite(r.rawComplementarity));
});

test('ordered requests and finite objective data are mandatory, including the three-interval edge case', () => {
  const input = read('30p-cove-position-stations');
  for (const positions of [null, [0, .5, NaN, 1], [0, .5, .5, 1], [0, .5, .6, Infinity]]) {
    assert.throws(() => fitConstrainedStationPositions({ ...input, positions }), /finite ordered/);
  }
  assert.throws(() => fitConstrainedStationPositions({ ...input, tolerance: NaN }), /finite ordered/);
  assert.throws(() => fitConstrainedStationPositions({ ...input, maxSweeps: 1 }), /did not converge/);
  assert.throws(() => fitConstrainedStationPositions({ positions: [0, 1e-170, 2e-170, 3e-170, 1], firstSpacing: .2, lastSpacing: .2 }), /Unrepresentable/);
  const r = fitConstrainedStationPositions({ positions: [0, 1 / 3, 2 / 3, 1], firstSpacing: 1 / 3, lastSpacing: 1 / 3 });
  r.positions.forEach((p, i) => assert.ok(Math.abs(p - i / 3) < 1e-12));
});

test('banded KKT preserves a known projection and rejects invalid or dependent active rows', () => {
  const r = projectPositionKkt({ scales: [1, 1], rows: [{ terms: [[0, 1], [1, 1]], rhs: 2 }], active: [0], certify: true });
  assert.deepEqual(r.target, [1, 1]); assert.deepEqual(r.multipliers, [-1]);
  assert.ok(r.forwardCertificate.physicalPositionErrorBound < 1e-12);
  for (const terms of [[[0, NaN]], [[0, 1], [0, 2]], [[-1, 1]], [[3, 1]]]) {
    assert.throws(() => projectPositionKkt({ scales: [1, 1], rows: [{ terms, rhs: 1 }], active: [0] }), /Invalid/);
  }
  assert.throws(() => projectPositionKkt({ scales: [1], rows: [{ terms: [[0, 1]], rhs: 1 }], active: [0, 0] }), /Invalid/);
  assert.throws(() => projectPositionKkt({ scales: [1], rows: [1, 2].map(rhs => ({ terms: [[0, 1]], rhs })), active: [0, 1] }), /Unresolved/);
});
