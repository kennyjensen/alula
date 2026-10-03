import test from 'node:test';
import assert from 'node:assert/strict';
import { certifyBernsteinCell, nextUp, nextDown, intervalPoint, intervalAdd, intervalSub, intervalMul, intervalDiv } from '../src/geometry/bernstein-cell-certificate.js';

const directions = [[{ x: 0, y: 1 }, { x: 0, y: 1 }], [{ x: 0, y: 1 }, { x: 0, y: 1 }]];
const conformal = () => Array.from({ length: 3 }, (_, i) => Array.from({ length: 3 }, (_, j) => ({
  x: i / 2 + .2 * (i === 2 ? 1 : 0) - .2 * (j === 2 ? 1 : 0), y: j / 2 + .1 * i * j,
})));
const integralQuadratic = (center, offset) => {
  const a = center * center + offset, x = [0, a / 3, (2 * a - center) / 3, a - center + 1 / 3];
  return x.map(value => [{ x: value, y: 0 }, { x: value, y: 1 }]);
};

// Exact rational arithmetic from the IEEE-754 bit pattern, independent of
// the production interval operations and floating-point expected results.
const rational = value => {
  const view = new DataView(new ArrayBuffer(8)); view.setFloat64(0, value);
  const b = view.getBigUint64(0), exponent = Number((b >> 52n) & 2047n), fraction = b & ((1n << 52n) - 1n);
  let numerator = exponent ? (1n << 52n) + fraction : fraction;
  const shift = exponent ? exponent - 1023 - 52 : -1074;
  if (b >> 63n) numerator = -numerator;
  return shift >= 0 ? { n: numerator << BigInt(shift), d: 1n } : { n: numerator, d: 1n << BigInt(-shift) };
};
const operations = {
  add: (a, b) => ({ n: a.n * b.d + b.n * a.d, d: a.d * b.d }),
  sub: (a, b) => ({ n: a.n * b.d - b.n * a.d, d: a.d * b.d }),
  mul: (a, b) => ({ n: a.n * b.n, d: a.d * b.d }),
  div: (a, b) => {
    const sign = b.n < 0 ? -1n : 1n;
    return { n: sign * a.n * b.d, d: sign * a.d * b.n };
  },
};
const compare = (a, b) => a.n * b.d - b.n * a.d;
const encloses = (bounds, exact) => {
  assert.ok(bounds[0] === -Infinity || compare(rational(bounds[0]), exact) <= 0n, `Lower bound ${bounds[0]} excludes exact value`);
  assert.ok(bounds[1] === Infinity || compare(rational(bounds[1]), exact) >= 0n, `Upper bound ${bounds[1]} excludes exact value`);
};

test('outward interval arithmetic encloses exact rational results, including underflow and overflow', () => {
  assert.equal(nextUp(0), Number.MIN_VALUE); assert.equal(nextDown(0), -Number.MIN_VALUE);
  assert.equal(nextUp(-Number.MIN_VALUE), -0); assert.equal(nextDown(Number.MIN_VALUE), 0);
  assert.equal(nextUp(Number.MAX_VALUE), Infinity); assert.equal(nextDown(-Number.MAX_VALUE), -Infinity);
  assert.equal(nextDown(Infinity), Number.MAX_VALUE); assert.equal(nextUp(-Infinity), -Number.MAX_VALUE);
  const cases = [[1, 2 ** -53], [.1, .2], [-.3, .1], [1e-200, 1e-200], [1e200, 1e200], [Number.MIN_VALUE, .25]];
  for (const [name, implementation] of Object.entries({ add: intervalAdd, sub: intervalSub, mul: intervalMul, div: intervalDiv }))
    for (const [a, b] of cases) encloses(implementation(intervalPoint(a), intervalPoint(b)), operations[name](rational(a), rational(b)));
  const a = [nextDown(-.3), nextUp(.1)], b = [nextDown(.2), nextUp(.7)];
  for (const [name, implementation] of Object.entries({ add: intervalAdd, sub: intervalSub, mul: intervalMul, div: intervalDiv })) {
    const bounds = implementation(a, b);
    for (const x of a) for (const y of b) encloses(bounds, operations[name](rational(x), rational(y)));
  }
  assert.deepEqual(intervalDiv([1, 1], [-1, 1]), [-Infinity, Infinity]);
  assert.deepEqual(intervalAdd([1, 1], [2, 2]), [nextDown(3), nextUp(3)]);
});

test('a nonaffine conformal polynomial receives whole-cell positive and transverse bounds', () => {
  const input = { controlPoints: conformal(), directions }, before = structuredClone(input);
  const result = certifyBernsteinCell(input, { minimumJacobian: .9, minimumTransversality: .6 });
  assert.equal(result.status, 'certified'); assert.equal(result.valid, true);
  assert.ok(result.lowerBounds.jacobian > .99 && result.lowerBounds.jacobian <= 1);
  assert.ok(result.lowerBounds.transversalityNumerator > .99 && result.lowerBounds.transversalityNumerator <= 1);
  assert.ok(result.lowerBounds.normalizedTransversality > .6);
  assert.equal(result.globallyInjective, false); assert.deepEqual(input, before);
});

test('fully bicubic geometry uses the complete degree-five Jacobian polynomial', () => {
  // x=s+.1s^3+.05st^2, y=t+.1t^3+.025s^2t. Its analytic J>=1.
  const controlPoints = Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({
    x: i / 3 + .1 * (i === 3 ? 1 : 0) + .05 * (i / 3) * j * (j - 1) / 6,
    y: j / 3 + .1 * (j === 3 ? 1 : 0) + .025 * i * (i - 1) / 6 * (j / 3),
  })));
  const result = certifyBernsteinCell({ controlPoints, directions });
  assert.equal(result.valid, true); assert.equal(result.jacobianStatus, 'certified');
  assert.ok(result.lowerBounds.jacobian > .99 && result.lowerBounds.jacobian <= 1);
});

test('subdivision certifies a positive polynomial whose initial Bernstein hull is inconclusive', () => {
  // x_s=(s-.5)^2+.01 > 0; its degree-2 Bernstein middle coefficient is negative.
  const input = { controlPoints: integralQuadratic(.5, .01), directions };
  const coarse = certifyBernsteinCell(input, { maxDepth: 0 });
  assert.equal(coarse.status, 'unresolved'); assert.equal(coarse.rejection, null);
  assert.ok(coarse.lowerBounds.jacobian < 0);
  const result = certifyBernsteinCell(input);
  assert.equal(result.valid, true); assert.ok(result.deepestSubdivision > 0);
  assert.ok(result.lowerBounds.jacobian > 0 && result.lowerBounds.jacobian < .010000001);
  const limited = certifyBernsteinCell(input, { maxPatches: 1 });
  assert.equal(limited.status, 'unresolved'); assert.equal(limited.patches, 1);
  assert.equal(limited.budgets.patchLimited, true); assert.equal(limited.valid, false);
});

test('negative interior Jacobian is detected despite positive corner and all 3x3 Gauss samples', () => {
  // x_s=(s-.25)^2-(1/16)^2. Its narrow negative interval misses the
  // usual Gauss abscissae and the four corners of an otherwise rectangular cell.
  const controlPoints = integralQuadratic(.25, -((1 / 16) ** 2));
  const derivative = controlPoints.slice(1).map((row, i) => 3 * (row[0].x - controlPoints[i][0].x));
  const at = s => derivative[0] * (1 - s) ** 2 + 2 * derivative[1] * s * (1 - s) + derivative[2] * s * s;
  for (const s of [0, 1, .5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2]) assert.ok(at(s) > 0);
  assert.ok(at(.25) < 0);
  const result = certifyBernsteinCell({ controlPoints, directions });
  assert.equal(result.status, 'rejected'); assert.equal(result.valid, false);
  assert.ok(result.rejection.upperBound < 0); assert.ok(result.deepestSubdivision >= 2);
});

test('uncertain coefficients remain unresolved and reversed guides are proved nontransverse', () => {
  const controlPoints = [[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: [-1, 1], y: 0 }, { x: [-1, 1], y: 1 }]];
  const uncertain = certifyBernsteinCell({ controlPoints, directions }, { maxDepth: 1 });
  assert.equal(uncertain.valid, false); assert.equal(uncertain.status, 'unresolved');
  const reversed = directions.map(row => row.map(d => ({ x: -d.x, y: -d.y })));
  const result = certifyBernsteinCell({ controlPoints: conformal(), directions: reversed });
  assert.equal(result.positive, true); assert.equal(result.transverse, false);
  assert.equal(result.status, 'rejected'); assert.equal(result.rejection.condition, 'transverse');
});

test('certificate bounds retain physical scaling and accept rigorously enclosed control coefficients', () => {
  const base = conformal(), original = certifyBernsteinCell({ controlPoints: base, directions });
  for (const scale of [1e-50, 1e50]) {
    const controlPoints = base.map(row => row.map(p => ({ x: intervalMul([p.x, p.x], [scale, scale]), y: intervalMul([p.y, p.y], [scale, scale]) })));
    const result = certifyBernsteinCell({ controlPoints, directions });
    assert.equal(result.valid, true);
    assert.ok(Math.abs(result.lowerBounds.jacobian / (scale * scale) - original.lowerBounds.jacobian) < 3e-13);
    assert.ok(Math.abs(result.lowerBounds.normalizedTransversality - original.lowerBounds.normalizedTransversality) < 3e-13);
  }
  assert.throws(() => certifyBernsteinCell({ controlPoints: base, directions }, { maxPatches: 0 }), /controls/);
  const bad = structuredClone(base); bad[0][0].x = [1, -1];
  assert.throws(() => certifyBernsteinCell({ controlPoints: bad, directions }), /coefficient/);
  assert.throws(() => intervalPoint(Infinity), /finite/);
});
