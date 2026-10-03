// Controlled body-grid fixtures, not an MSET initializer or an experimental
// geometry. Every supplied cut is an initial guess and moves in the solve.
import { naca4, transform } from '../../src/geometry/airfoil.js';
import { createContourCurve } from '../../src/geometry/contour-curve.js';

export function intrinsicBodyFixture({ bodySegments = 8, tubes = 3, elements = 1, alpha = 0, mach = .2, height = 2,
  tubeGrowth = 2, surfaceSpacing = 'uniform', cutSpacing = 'legacy', contourPanels = 160 } = {}) {
  if (![1, 2].includes(elements)) throw new Error('This fixture defines one or two elements.');
  if (!Number.isInteger(bodySegments) || bodySegments < 4 || !Number.isInteger(tubes) || tubes < 1
    || !Number.isFinite(tubeGrowth) || tubeGrowth < 1 || !['uniform', 'cosine'].includes(surfaceSpacing)
    || !Number.isFinite(height) || height <= .3 || !['legacy', 'surface-matched'].includes(cutSpacing)
    || !Number.isInteger(contourPanels) || contourPanels < 40 || contourPanels % 2) throw new Error('Invalid intrinsic body fixture resolution.');
  const contours = [naca4('0012', contourPanels)];
  if (elements === 2) contours.push(transform(naca4('0012', contourPanels), { chord: .3, x: -.4, y: .2 }));
  const curves = contours.map(createContourCurve), stagnation = curves.map(c => c.length / 2);
  const round = x => Math.round(x * 1e12) / 1e12;
  const spacing = t => surfaceSpacing === 'cosine' ? .5 * (1 - Math.cos(Math.PI * t)) : t;
  const xs = [-4, -1, -.6, ...curves.flatMap((curve, body) => Array.from({ length: bodySegments + 1 }, (_, i) => round(curve.branch('upper', spacing(i / bodySegments), stagnation[body]).point.x))), 1.1, 1.4, 2, 4];
  // Controlled study of the interval-size jump at LE/TE. This is not the
  // general MSET block/gap-spacing algorithm: retain legacy fixtures too.
  if (cutSpacing === 'surface-matched') curves.forEach((curve, body) => {
    const le = curve.evaluate(stagnation[body]).point, te = contours[body][0];
    const ds = .5 * curve.length * spacing(1 / bodySegments);
    for (const [origin, sign] of [[le.x, -1], [te.x, 1]]) {
      let distance = ds, step = ds;
      while (distance < .2 * (te.x - le.x)) { xs.push(round(origin + sign * distance)); step *= 1.4; distance += step; }
    }
  });
  const x = [...new Set(xs)].sort((a, b) => a - b);
  // Keep an upstream cross-line close to the first leading edge.
  const firstLE = Math.min(...curves.map((c, b) => c.evaluate(stagnation[b]).point.x));
  if (!x.includes(round(firstLE - .08))) x.push(round(firstLE - .08)); x.sort((a, b) => a - b);
  const bodies = contours.map((points, body) => {
    const curve = curves[body], leadingIndex = x.indexOf(round(curve.evaluate(stagnation[body]).point.x)), trailingIndex = x.indexOf(round(points[0].x));
    const surfaceFractions = x.slice(leadingIndex, trailingIndex + 1).map((xi, i, row) => {
      if (i === 0) return 0; if (i === row.length - 1) return 1;
      let lo = 0, hi = 1;
      for (let k = 0; k < 50; k++) { const mid = (lo + hi) / 2; if (curve.branch('upper', mid, stagnation[body]).point.x < xi) lo = mid; else hi = mid; }
      return (lo + hi) / 2;
    });
    return { points, leadingIndex, trailingIndex, surfaceFractions, stagnationParameter: stagnation[body] };
  });
  const cutPaths = contours.map(points => x.map(x => ({ x, y: points[0].y })));
  const weights = Array.from({ length: elements + 1 }, (_, group) => Array.from({ length: tubes }, (_, j) =>
    group === 0 ? tubeGrowth ** (tubes - 1 - j) : group === elements ? tubeGrowth ** j : Math.sin(Math.PI * (j + .5) / tubes)));
  return { bodies, outerLower: x.map(x => ({ x, y: -height })), outerUpper: x.map(x => ({ x, y: height })), cutPaths, weights, alpha, mach };
}
