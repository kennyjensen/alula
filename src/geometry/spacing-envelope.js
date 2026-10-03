// SPDX-License-Identifier: GPL-2.0-or-later
// Positive piecewise-linear spacing requests and their lower envelope.
// The density integral is analytic on each segment and on each crossing.
// Full requests can extend outside the integration window. Clipping the
// window must not insert fictitious request nodes or change local density.
// This is the same construction used by the common station scheduler.
export function createSpacingEnvelope({ profiles, start, end, growth = .25 }) {
  if (!Array.isArray(profiles) || !profiles.length || ![start, end, growth].every(Number.isFinite)
    || !(start < end) || !(growth > 0)
    || profiles.some(row => !Array.isArray(row) || row.length < 3 || !row.every(Number.isFinite)
      || row.some((v, i) => i && v <= row[i - 1])))
    throw new Error('Spacing profiles must increase and integration bounds must be finite and ordered.');
  const span = end - start;
  const monitors = profiles.map(row => {
    const centers = row.slice(1).map((v, i) => ({ x: .5 * (row[i] + v), h: v - row[i] }));
    const points = [{ x: row[0], h: centers[0].h }, ...centers, { x: row.at(-1), h: centers.at(-1).h }];
    const spacing = x => {
      if (x <= points[0].x) return points[0].h + growth * (points[0].x - x);
      if (x >= points.at(-1).x) return points.at(-1).h + growth * (x - points.at(-1).x);
      let lo = 0, hi = points.length - 1;
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (points[mid].x < x) lo = mid; else hi = mid; }
      const a = points[lo], b = points[hi]; return a.h + (b.h - a.h) * (x - a.x) / (b.x - a.x);
    };
    return { points, spacing };
  });
  const anchors = [...new Set([start, end, ...profiles.flatMap(row => [row[0], row.at(-1)]).filter(x => x > start && x < end)])].sort((a, b) => a - b);
  if (anchors.some((v, i) => i && v - anchors[i - 1] < 1e-10 * span)) throw new Error('Unresolved nearly coincident potential block corners.');
  const knots = [...new Set([...anchors, ...monitors.flatMap(m => m.points.map(p => p.x)).filter(x => x >= start && x <= end)])].sort((a, b) => a - b);
  const integral = (h, slope, width) => slope === 0 ? width / h : Math.log1p(slope * width / h) / slope;
  const pieces = []; let total = 0;
  for (let k = 1; k < knots.length; k++) {
    const left = knots[k - 1], right = knots[k], width = right - left;
    const lines = monitors.map(m => ({ a: m.spacing(left), b: m.spacing(right) })), splits = [left, right];
    // The minimum of affine segments can switch at pairwise intersections.
    for (let i = 0; i < lines.length; i++) for (let j = i + 1; j < lines.length; j++) {
      const a = lines[i].a - lines[j].a, b = lines[i].b - lines[j].b;
      if (a * b < 0) { const x = left + width * a / (a - b); if (x > left && x < right) splits.push(x); }
    }
    splits.sort((a, b) => a - b);
    for (let j = 1; j < splits.length; j++) {
      const a = splits[j - 1], b = splits[j]; if (!(b > a)) continue;
      const midpoint = .5 * (a + b); let chosen = 0;
      for (let m = 1; m < monitors.length; m++) if (monitors[m].spacing(midpoint) < monitors[chosen].spacing(midpoint)) chosen = m;
      const h = monitors[chosen].spacing(a), slope = (lines[chosen].b - lines[chosen].a) / width;
      const mass = integral(h, slope, b - a);
      if (!(h > 0 && mass > 0) || !Number.isFinite(mass)) throw new Error('Invalid potential spacing metric.');
      pieces.push({ left: a, right: b, h, slope, mass, total }); total += mass;
    }
  }
  const metric = x => {
    if (!Number.isFinite(x) || x < start || x > end) throw new Error('Spacing metric argument is outside its bounds.');
    if (x === end) return total;
    let lo = 0, hi = pieces.length - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (pieces[mid].right <= x) lo = mid + 1; else hi = mid; }
    const p = pieces[lo]; return p.total + integral(p.h, p.slope, x - p.left);
  };
  const inverse = mass => {
    if (!Number.isFinite(mass) || mass < 0 || mass > total) throw new Error('Spacing integral argument is outside its bounds.');
    let lo = 0, hi = pieces.length - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (pieces[mid].total + pieces[mid].mass < mass) lo = mid + 1; else hi = mid; }
    const p = pieces[lo], dm = mass - p.total;
    return p.left + (p.slope === 0 ? dm * p.h : p.h * Math.expm1(p.slope * dm) / p.slope);
  };
  return { total, anchors, metric, inverse };
}
