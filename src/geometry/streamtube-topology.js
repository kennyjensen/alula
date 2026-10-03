// SPDX-License-Identifier: GPL-2.0-or-later
// Initial cross-line topology for downstream-facing, x-monotone
// bodies. Panel tracing subsequently moves the cuts and free nodes. General
// MSET block selection remains incomplete. The optional panel/potential
// matching and elliptic stages subsequently reconcile these initial anchors.
import { prepareContour, validateAssembly } from './airfoil.js';
import { createContourCurve } from './contour-curve.js';
import { createContourTopology, createSurfaceContourCurve } from './contour-topology.js';
import { gradedIntervals } from './graded-intervals.js';
import { prepareAirfoilElement } from './airfoil-element.js';

export function createInitialStreamtubeTopology({ elements, alpha = 0, mach = .2 }, {
  // Sample a fixed exponential mass-coordinate map as the tube count grows.
  // At the original seven-tube resolution this is exactly the old ratio 3.
  // Holding that *adjacent* ratio at 3 instead made the first tube shrink
  // exponentially under refinement, below the panel wall-field accuracy.
  // An explicitly supplied tubeGrowth still means the adjacent mass ratio.
  surfaceIntervals = 16, surfaceChordExponent = 0, tubes = 7, lowerTubes = tubes, upperTubes = tubes, gapTubes = tubes,
  tubeGrowth, padding = 4, height = 2,
  inletIntervals = Math.ceil(surfaceIntervals / 2), outletIntervals = Math.ceil(surfaceIntervals / 2), parametricLowerBranches = false } = {}) {
  if (typeof parametricLowerBranches !== 'boolean' || !Array.isArray(elements) || elements.length < 1 || elements.length > 6
    || !Number.isInteger(surfaceIntervals) || surfaceIntervals < 8 || surfaceIntervals > 128
    || ![tubes, lowerTubes, upperTubes, gapTubes].every(n => Number.isInteger(n) && n >= 3 && n <= 31)
    || ![tubeGrowth ?? 1, padding, height, surfaceChordExponent].every(Number.isFinite)
    || (tubeGrowth !== undefined && tubeGrowth < 1) || padding < 1 || height <= 0 || surfaceChordExponent < 0 || surfaceChordExponent > 1
    || ![inletIntervals, outletIntervals].every(n => Number.isInteger(n) && n >= 4 && n <= 256))
    throw new Error('Invalid initial streamtube topology controls.');
  elements = elements.map(prepareAirfoilElement);
  const topologies = elements.map(e => e.trailingEdge?.kind === 'finite-base' ? createContourTopology(e.points, e) : null);
  const contours = elements.map((e, k) => topologies[k]?.points ?? prepareContour(e.points)); validateAssembly(contours);
  // Explicit base panels have constrained strengths rather than independent
  // surface gamma unknowns. Preserve them without spending the surface budget.
  if (contours.reduce((sum, p, k) => sum + (topologies[k]?.surface.panels.length ?? p.length - 1), 0) > 700)
    throw new Error('The panel streamtube initializer supports at most 700 total surface panels.');
  const reentrantElements = new Set();
  const profiles = contours.map((solidPoints, element) => {
    const points = topologies[element]?.surface.points ?? solidPoints;
    const curve = topologies[element] ? createSurfaceContourCurve(solidPoints, elements[element]) : createContourCurve(points); let j = 0;
    for (let i = 1; i < points.length - 1; i++) if (points[i].x < points[j].x) j = i;
    if (j === 0 || j === points.length - 2) throw new Error('Initial streamtube topology requires a downstream-facing trailing edge.');
    let lo = curve.knots[j - 1], hi = curve.knots[j + 1];
    for (let k = 0; k < 60; k++) {
      const a = (2 * lo + hi) / 3, b = (lo + 2 * hi) / 3;
      if (curve.evaluate(a).point.x < curve.evaluate(b).point.x) hi = b; else lo = a;
    }
    const stagnationParameter = .5 * (lo + hi), nose = curve.evaluate(stagnationParameter).point;
    // Check extrema of each quadratic spline derivative, not just a set of
    // plotted samples, before using monotone x inversion on either branch.
    for (let i = 0; i < curve.knots.length - 1; i++) {
      const a = curve.knots[i], b = curve.knots[i + 1], va = curve.evaluate(a), vb = curve.evaluate(b), locations = [a, b];
      if (va.secondDerivative.x * vb.secondDerivative.x < 0)
        locations.push(a - va.secondDerivative.x * (b - a) / (vb.secondDerivative.x - va.secondDerivative.x));
      for (const s of locations) {
        const dx = curve.evaluate(s).derivative.x;
        if ((s < stagnationParameter && dx > 1e-8) || (s > stagnationParameter && dx < -1e-8)) {
          const side = s < stagnationParameter ? 'upper' : 'lower';
          if (parametricLowerBranches && side === 'lower') reentrantElements.add(element);
          else throw Object.assign(new Error(`Element ${element + 1} needs a nonmonotone body topology.`),
            { code: 'STREAMTUBE_NONMONOTONE_BODY', diagnostics: { element, side, parameter: s, derivativeX: dx } });
        }
      }
    }
    return { points, solidPoints, curve, element, stagnationParameter, nose, chord: points[0].x - nose.x,
      surfaceChord: Math.max(...points.map(p => Math.hypot(p.x - points[0].x, p.y - points[0].y))),
      meanY: points.reduce((sum, p) => sum + p.y, 0) / points.length };
  }).sort((a, b) => a.meanY - b.meanY);
  const chord = Math.max(...profiles.map(p => p.chord)), origin = Math.min(...profiles.map(p => p.nose.x));
  const largestSurfaceChord = Math.max(...profiles.map(p => p.surfaceChord));
  profiles.forEach(p => { p.intervals = Math.max(4, Math.round(surfaceIntervals * (p.surfaceChord / largestSurfaceChord) ** surfaceChordExponent)); });
  const xmin = origin - padding * chord, xmax = Math.max(...profiles.map(p => p.points[0].x)) + padding * chord;
  const snap = x => origin + Math.round((x - origin) / (1e-12 * chord)) * 1e-12 * chord;
  const firstBody = profiles.reduce((a, b) => a.nose.x < b.nose.x ? a : b);
  const lastBody = profiles.reduce((a, b) => a.points[0].x > b.points[0].x ? a : b);
  const firstFraction = .5 * (1 - Math.cos(Math.PI / firstBody.intervals));
  const lastFraction = .5 * (1 - Math.cos(Math.PI / lastBody.intervals));
  const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
  const inlet = gradedIntervals({ extent: firstBody.nose.x - xmin, intervals: inletIntervals,
    firstSpacing: distance(firstBody.nose, firstBody.curve.branch('upper', firstFraction, firstBody.stagnationParameter).point) });
  const outlet = gradedIntervals({ extent: xmax - lastBody.points[0].x, intervals: outletIntervals,
    firstSpacing: distance(lastBody.points[0], lastBody.curve.branch('upper', 1 - lastFraction, lastBody.stagnationParameter).point) });
  const xs = [...inlet.x.map(d => firstBody.nose.x - d), ...outlet.x.map(d => lastBody.points[0].x + d),
    ...profiles.flatMap(p => [
      // Local anchors inside the assembly still connect adjacent body
      // blocks. Exterior anchors belong to the graded distribution above.
      ...[p.nose.x - .08 * p.chord, p.points[0].x + .1 * p.chord]
        .filter(x => x > firstBody.nose.x && x < lastBody.points[0].x),
      ...Array.from({ length: p.intervals + 1 }, (_, i) => p.curve.branch('upper',
        .5 * (1 - Math.cos(Math.PI * i / p.intervals)), p.stagnationParameter).point.x)])];
  const x = [...new Set(xs.map(snap))].sort((a, b) => a - b);
  const bodies = profiles.map(p => {
    const leadingIndex = x.indexOf(snap(p.nose.x)), trailingIndex = x.indexOf(snap(p.points[0].x));
    const surfaceFractions = x.slice(leadingIndex, trailingIndex + 1).map((xi, i, row) => {
      if (i === 0) return 0; if (i === row.length - 1) return 1;
      let lo = 0, hi = 1;
      for (let k = 0; k < 55; k++) {
        const f = .5 * (lo + hi);
        if (p.curve.branch('upper', f, p.stagnationParameter).point.x < xi) lo = f; else hi = f;
      }
      return .5 * (lo + hi);
    });
    return { points: p.solidPoints, leadingIndex, trailingIndex, stagnationParameter: p.stagnationParameter, surfaceFractions, element: p.element,
      ...(topologies[p.element] ? { trailingEdge: structuredClone(elements[p.element].trailingEdge) } : {}) };
  });
  const y = contours.flatMap(points => points.map(p => p.y)), lower = Math.min(...y) - height * chord, upper = Math.max(...y) + height * chord;
  return { ...(reentrantElements.size ? { potentialTopologyRequired: true, reentrantElements: [...reentrantElements] } : {}), bodies, primaryBody: profiles.findIndex(p => p.chord === chord), alpha, mach,
    gridSpacing: { surfaceIntervals, surfaceIntervalsByElement: bodies.map(b => ({ element: b.element, intervals: b.trailingIndex - b.leadingIndex })),
      surfaceChordExponent, requestedSurfaceIntervalsByElement: profiles.map(p => ({ element: p.element, intervals: p.intervals })),
      inlet: { intervals: inlet.intervals, growth: inlet.growth, firstSpacing: inlet.firstSpacing },
      outlet: { intervals: outlet.intervals, growth: outlet.growth, firstSpacing: outlet.firstSpacing } },
    outerLower: x.map(x => ({ x, y: lower })), outerUpper: x.map(x => ({ x, y: upper })),
    cutPaths: profiles.map(p => x.map(x => ({ x, y: p.meanY }))),
    weights: Array.from({ length: bodies.length + 1 }, (_, g) => {
      const n = g === 0 ? lowerTubes : g === bodies.length ? upperTubes : gapTubes, growth = tubeGrowth ?? 3 ** (7 / n);
      return Array.from({ length: n }, (_, j) => g === 0 ? growth ** (n - 1 - j) : g === bodies.length ? growth ** j : Math.sin(Math.PI * (j + .5) / n));
    }) };
}
