// SPDX-License-Identifier: GPL-2.0-or-later
// Outward-rounded tensor Bernstein bounds for a single polynomial cell.
// This certifies local Jacobian/transversality conditions, not global
// injectivity, shared-edge conformity, body clearance, or a line-search path.
const bits = new DataView(new ArrayBuffer(8));
export function nextUp(x) {
  if (Number.isNaN(x) || x === Infinity) return x;
  if (x === -Infinity) return -Number.MAX_VALUE;
  if (x === 0) return Number.MIN_VALUE;
  bits.setFloat64(0, x);
  bits.setBigUint64(0, bits.getBigUint64(0) + (x > 0 ? 1n : -1n));
  return bits.getFloat64(0);
}
export const nextDown = x => -nextUp(-x);
const entire = () => [-Infinity, Infinity];
const zero = a => a[0] === 0 && a[1] === 0;
const one = a => a[0] === 1 && a[1] === 1;
const rounded = (lo, hi) => Number.isNaN(lo) || Number.isNaN(hi) ? entire() : [nextDown(lo), nextUp(hi)];
export function intervalPoint(value) {
  if (!Number.isFinite(value)) throw new Error('An exact interval input must be finite.');
  return [value, value];
}
export function intervalAdd(a, b) {
  if (zero(a)) return b.slice(); if (zero(b)) return a.slice();
  return rounded(a[0] + b[0], a[1] + b[1]);
}
export function intervalSub(a, b) {
  if (zero(b)) return a.slice();
  if (zero(a)) return [-b[1], -b[0]];
  return rounded(a[0] - b[1], a[1] - b[0]);
}
export function intervalMul(a, b) {
  if (zero(a) || zero(b)) return [0, 0];
  if (one(a)) return b.slice(); if (one(b)) return a.slice();
  const products = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]];
  if (products.some(Number.isNaN)) return entire();
  return rounded(Math.min(...products), Math.max(...products));
}
export function intervalDiv(a, b) {
  if (b[0] <= 0 && b[1] >= 0) return entire();
  if (zero(a)) return [0, 0]; if (one(b)) return a.slice();
  const quotients = [a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1]];
  if (quotients.some(Number.isNaN)) return entire();
  return rounded(Math.min(...quotients), Math.max(...quotients));
}
const coefficient = value => {
  if (Number.isFinite(value)) return [value, value];
  if (Array.isArray(value) && value.length === 2 && value.every(Number.isFinite) && value[0] <= value[1]) return value.slice();
  throw new Error('A Bernstein coefficient must be a finite number or enclosing [lo,hi] interval.');
};

const choose = (n, k) => {
  let result = 1;
  for (let i = 1; i <= k; i++) result = result * (n - i + 1) / i;
  return result; // Exact small integer: polynomial degrees never exceed six.
};
const derivative = (p, axis) => {
  const ns = p.length - 1, nt = p[0].length - 1, degree = axis ? nt : ns;
  return Array.from({ length: axis ? ns + 1 : Math.max(1, ns) }, (_, i) =>
    Array.from({ length: axis ? Math.max(1, nt) : nt + 1 }, (_, j) => degree
      ? intervalMul([degree, degree], intervalSub(p[i + (axis ? 0 : 1)][j + (axis ? 1 : 0)], p[i][j])) : [0, 0]));
};
const subtract = (a, b) => a.map((row, i) => row.map((value, j) => intervalSub(value, b[i][j])));
const product = (a, b) => {
  const as = a.length - 1, at = a[0].length - 1, bs = b.length - 1, bt = b[0].length - 1;
  const result = Array.from({ length: as + bs + 1 }, () => Array.from({ length: at + bt + 1 }, () => [0, 0]));
  for (let i = 0; i <= as; i++) for (let j = 0; j <= at; j++) for (let k = 0; k <= bs; k++) for (let l = 0; l <= bt; l++) {
    const ws = intervalDiv([choose(as, i) * choose(bs, k), choose(as, i) * choose(bs, k)], [choose(as + bs, i + k), choose(as + bs, i + k)]);
    const wt = intervalDiv([choose(at, j) * choose(bt, l), choose(at, j) * choose(bt, l)], [choose(at + bt, j + l), choose(at + bt, j + l)]);
    result[i + k][j + l] = intervalAdd(result[i + k][j + l], intervalMul(intervalMul(a[i][j], b[k][l]), intervalMul(ws, wt)));
  }
  return result;
};
const splitCurve = points => {
  let row = points.map(p => p.slice()); const left = [row[0]], right = [row.at(-1)];
  while (row.length > 1) {
    row = row.slice(1).map((p, i) => intervalMul([.5, .5], intervalAdd(row[i], p)));
    left.push(row[0]); right.unshift(row.at(-1));
  }
  return [left, right];
};
const splitAxis = (p, axis) => {
  if (axis) {
    const rows = p.map(splitCurve);
    return [rows.map(row => row[0]), rows.map(row => row[1])];
  }
  const columns = p[0].map((_, j) => splitCurve(p.map(row => row[j])));
  return [0, 1].map(side => p.map((_, i) => columns.map(column => column[side][i])));
};
const quarters = p => splitAxis(p, 0).flatMap(half => splitAxis(half, 1));
const range = p => [Math.min(...p.flat().map(a => a[0])), Math.max(...p.flat().map(a => a[1]))];
const maxAbs = a => Math.max(Math.abs(a[0]), Math.abs(a[1]));
const normUpper = (x, y) => {
  const a = maxAbs(x), b = maxAbs(y), fallback = intervalAdd([a, a], [b, b])[1];
  const square = intervalAdd(intervalMul([a, a], [a, a]), intervalMul([b, b], [b, b]))[1];
  if (!Number.isFinite(square)) return fallback;
  if (square === 0) return 0;
  // sqrt is only a starting estimate. The directed square verifies the
  // returned upper bound; an L1 bound remains valid if range is unresolved.
  let value = Math.sqrt(square);
  for (let k = 0; k < 4 && Number.isFinite(value); k++, value = nextUp(value))
    if (nextDown(value * value) >= square) return value;
  return fallback;
};
const vectorUpper = (x, y) => {
  const a = x.flat(), b = y.flat(); let upper = 0;
  for (let k = 0; k < a.length; k++) upper = Math.max(upper, normUpper(a[k], b[k]));
  return upper;
};

export function certifyBernsteinCell({ controlPoints, directions }, {
  maxDepth = 10, maxPatches = 4096, minimumJacobian = 0, minimumTransversality = 0,
} = {}) {
  const ns = controlPoints?.length - 1, nt = controlPoints?.[0]?.length - 1;
  if (!Array.isArray(controlPoints) || !(ns >= 1 && ns <= 3 && nt >= 1 && nt <= 3)
    || controlPoints.some(row => !Array.isArray(row) || row.length !== nt + 1)
    || !Array.isArray(directions) || directions.length !== 2 || directions.some(row => !Array.isArray(row) || row.length !== 2
      || row.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))) throw new Error('Invalid polynomial cell or Q1 guide dimensions.');
  if (!Number.isInteger(maxDepth) || maxDepth < 0 || maxDepth > 30 || !Number.isInteger(maxPatches) || maxPatches < 1
    || !Number.isFinite(minimumJacobian) || minimumJacobian < 0 || !Number.isFinite(minimumTransversality)
    || minimumTransversality < 0 || minimumTransversality >= 1) throw new Error('Invalid Bernstein certificate controls.');
  const x = controlPoints.map(row => row.map(p => coefficient(p?.x))), y = controlPoints.map(row => row.map(p => coefficient(p?.y)));
  const ax = derivative(x, 0), ay = derivative(y, 0), bx = derivative(x, 1), by = derivative(y, 1);
  const dx = directions.map(row => row.map(p => intervalPoint(p.x))), dy = directions.map(row => row.map(p => intervalPoint(p.y)));
  const polynomial = { jacobian: subtract(product(ax, by), product(ay, bx)),
    transverse: subtract(product(ax, dy), product(ay, dx)), ax, ay, dx, dy };
  const summarize = (poly, box, depth) => {
    const j = range(poly.jacobian), p = range(poly.transverse);
    const derivativeUpper = vectorUpper(poly.ax, poly.ay), guideUpper = vectorUpper(poly.dx, poly.dy);
    const denominator = intervalMul([derivativeUpper, derivativeUpper], [guideUpper, guideUpper])[1];
    const normalized = p[0] > 0 && denominator > 0 ? intervalDiv([p[0], p[0]], [denominator, denominator])[0] : -Infinity;
    const positive = j[0] > minimumJacobian;
    const transverse = p[0] > 0 && (minimumTransversality === 0 || normalized > minimumTransversality);
    let rejected = null; const rejectedConditions = new Set();
    for (const [key, margin] of [['jacobian', minimumJacobian], ['transverse', 0]]) {
      const q = poly[key];
      for (const [i, k] of [[0, 0], [q.length - 1, 0], [q.length - 1, q[0].length - 1], [0, q[0].length - 1]])
        if (q[i][k][1] <= margin) {
          rejectedConditions.add(key);
          rejected ??= { condition: key, upperBound: q[i][k][1], requiredGreaterThan: margin,
            point: { s: i === 0 ? box[0] : box[1], t: k === 0 ? box[2] : box[3] } };
        }
    }
    return { poly, box, depth, j, p, normalized, positive, transverse, rejected, rejectedConditions };
  };
  const first = summarize(polynomial, [0, 1, 0, 1], 0), queue = [first], leaves = new Set([first]);
  let patches = 1, cursor = 0, deepestSubdivision = 0, depthLimited = false, patchLimited = false, rejection = null;
  while (cursor < queue.length) {
    const patch = queue[cursor++];
    if (patch.rejected) { rejection = patch.rejected; break; }
    if (patch.positive && patch.transverse) continue;
    if (patch.depth >= maxDepth) { depthLimited = true; continue; }
    if (patches + 4 > maxPatches) { patchLimited = true; continue; }
    const pieces = Object.fromEntries(Object.entries(patch.poly).map(([key, p]) => [key, quarters(p)]));
    const [s0, s1, t0, t1] = patch.box, sm = .5 * (s0 + s1), tm = .5 * (t0 + t1);
    const boxes = [[s0, sm, t0, tm], [s0, sm, tm, t1], [sm, s1, t0, tm], [sm, s1, tm, t1]];
    const children = boxes.map((box, k) => summarize(Object.fromEntries(Object.entries(pieces).map(([key, p]) => [key, p[k]])), box, patch.depth + 1));
    leaves.delete(patch);
    for (const child of children) { leaves.add(child); queue.push(child); }
    patches += 4; deepestSubdivision = Math.max(deepestSubdivision, patch.depth + 1);
  }
  const frontier = [...leaves], positive = frontier.every(p => p.positive), transverse = frontier.every(p => p.transverse);
  const valid = positive && transverse;
  const lowerBounds = { jacobian: Infinity, transversalityNumerator: Infinity, normalizedTransversality: Infinity };
  for (const p of frontier) {
    lowerBounds.jacobian = Math.min(lowerBounds.jacobian, p.j[0]);
    lowerBounds.transversalityNumerator = Math.min(lowerBounds.transversalityNumerator, p.p[0]);
    lowerBounds.normalizedTransversality = Math.min(lowerBounds.normalizedTransversality, p.normalized);
  }
  return { valid, positive, transverse, status: valid ? 'certified' : rejection ? 'rejected' : 'unresolved',
    jacobianStatus: positive ? 'certified' : frontier.some(p => p.rejectedConditions.has('jacobian')) ? 'rejected' : 'unresolved',
    transversalityStatus: transverse ? 'certified' : frontier.some(p => p.rejectedConditions.has('transverse')) ? 'rejected' : 'unresolved',
    lowerBounds, rejection, patches, leafPatches: frontier.length, deepestSubdivision,
    unresolvedPatches: frontier.filter(p => !p.rejected && (!p.positive || !p.transverse)).length,
    rejectedPatches: frontier.filter(p => p.rejected).length,
    budgets: { maxDepth, maxPatches, depthLimited, patchLimited },
    margins: { minimumJacobian, minimumTransversality }, globallyInjective: false,
    scope: 'Whole-cell polynomial Jacobian and interpolated-guide transversality bounds, using outward-rounded interval Bernstein arithmetic. Endpoint geometry only; no global injectivity or path-validity claim.' };
}
