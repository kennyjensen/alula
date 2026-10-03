// SPDX-License-Identifier: GPL-2.0-or-later
// An acceptance inventory, not a list of known successes.
import { getBenchmarkAirfoil } from '../../src/geometry/benchmark-airfoils.js';

export const gridRobustnessControls = Object.freeze({
  surfaceIntervals: [8, 16, 32, 64, 128], tubes: [7, 9, 11, 24],
  inletIntervals: ['auto', 16, 32, 64, 128], outletIntervals: ['auto', 16, 32, 64, 128],
  overrideStreamlines: { min: 4, max: 32, step: 1 },
  chordExponent: { min: 0, max: 1 }, endSpacingRatio: { min: .01, max: 2 },
  curvatureExponent: { min: .05, max: 4 }, aspectRatio: { min: .1, max: 20 },
});
export const gridRobustnessPresets = ['rae2822-mses', 'nlr7301'];
const curvature = { exponent: .5, leadingSpacingRatio: .2, trailingSpacingRatio: .4 };
const advanced = [
  ...[16, 32, 64, 128].flatMap(n => [
    { suffix: `inlet-${n}`, changes: { gridInletIntervals: n } },
    { suffix: `wake-${n}`, changes: { gridOutletIntervals: n } },
  ]),
  { suffix: 'ends-16-128', changes: { gridInletIntervals: 16, gridOutletIntervals: 128 } },
  { suffix: 'ends-128-16', changes: { gridInletIntervals: 128, gridOutletIntervals: 16 } },
  ...['Upper', 'Lower', 'Gap'].flatMap(region => [4, 32].map(n => ({
    suffix: `${region.toLowerCase()}-streamlines-${n}`, changes: { [`grid${region}Tubes`]: n - 1 },
  }))),
  { suffix: 'asymmetric-low-gap', changes: { gridUpperTubes: 31, gridLowerTubes: 3, gridGapTubes: 3 } },
  { suffix: 'asymmetric-high-gap', changes: { gridUpperTubes: 3, gridLowerTubes: 31, gridGapTubes: 31 } },
  ...[.4, 1].map(n => ({ suffix: `chord-${n}`, changes: { gridChordExponent: n } })),
  { suffix: 'curvature', changes: { gridSurfaceSpacing: 'curvature', gridCurvatureSpacing: curvature } },
  ...['leadingSpacingRatio', 'trailingSpacingRatio'].flatMap(key => [.01, 2].map(n => ({
    suffix: `${key}-${n}`, changes: { gridSurfaceSpacing: 'curvature', gridCurvatureSpacing: { ...curvature, [key]: n } },
  }))),
  ...[.05, 4].map(n => ({ suffix: `curvature-exponent-${n}`, changes: {
    gridSurfaceSpacing: 'curvature', gridCurvatureSpacing: { ...curvature, exponent: n },
  } })),
  ...[.1, 2.5, 20].map(n => ({ suffix: `aspect-${n}`, changes: { gridStagnationAspectRatio: n } })),
];

export const gridRobustnessCases = gridRobustnessPresets.flatMap(preset => [
  ...gridRobustnessControls.surfaceIntervals.flatMap(gridIntervals =>
    gridRobustnessControls.tubes.flatMap(gridTubes => [false, true].map(gridEllipticSmoothing => ({
      id: `${preset}-${gridIntervals}x${gridTubes}-slor-${gridEllipticSmoothing ? 'on' : 'off'}`,
      preset, tier: 'surface-tube-matrix', changes: { gridIntervals, gridTubes, gridEllipticSmoothing },
    })))),
  ...advanced.map(({ suffix, changes }) => ({ id: `${preset}-advanced-${suffix}`, preset, tier: 'advanced', changes })),
]);

export function buildGridRobustnessCase(spec, { mode = 'coupled' } = {}) {
  if (!gridRobustnessPresets.includes(spec?.preset) || !['coupled', 'euler'].includes(mode))
    throw new Error('Unknown grid robustness preset or mode.');
  const preset = getBenchmarkAirfoil(spec.preset), rae = spec.preset === 'rae2822-mses';
  return { elements: preset.elements, alpha: rae ? 2.68 : 6, mach: rae ? .74 : .185,
    reynolds: rae ? 2510000 : 2700000, ncrit: rae ? 4 : 9, referenceChord: 1,
    momentReference: { x: .25, y: 0 }, flowModel: 'streamtube-grid', quadBoundaryLayers: mode === 'coupled',
    transitionMode: 'automatic', materialTrips: preset.elements.map(() => [1, 1]), eulerIsmom: 4,
    gridIntervals: 16, gridTubes: 7, gridCrosslinePlacement: 'potential', gridChordExponent: 0,
    gridSurfaceSpacing: 'supplied', gridEllipticSmoothing: true, gridSmoothingMethod: 'elliptic',
    geometrySource: { id: preset.id, referenceChord: preset.referenceChord, ...preset.provenance },
    ...structuredClone(spec.changes),
  };
}

const bounded = (x, min, max, integer = false) => Number.isFinite(x) && x >= min && x <= max
  && (!integer || Number.isInteger(x));

// Resolve the requested logical counts only. Tracing, matching, aspect
// allocation and turn refinement may increase counts; no actual mesh is inferred.
export function planGridRobustnessCase(input) {
  if (!Array.isArray(input?.elements) || !input.elements.length || input.elements.length > 6
    || !bounded(input.gridIntervals, 8, 128, true) || !bounded(input.gridTubes, 3, 31, true)
    || !bounded(input.gridChordExponent, 0, 1) || typeof input.gridEllipticSmoothing !== 'boolean'
    || input.gridSmoothingMethod !== 'elliptic' || input.gridCrosslinePlacement !== 'potential'
    || !bounded(input.mach, Number.MIN_VALUE, 1 - Number.EPSILON)
    || !Number.isFinite(input.alpha) || !Number.isFinite(input.referenceChord) || !(input.referenceChord > 0))
    throw new Error('Invalid grid robustness case controls.');
  const inletIntervals = input.gridInletIntervals ?? 2 * input.gridIntervals;
  const outletIntervals = input.gridOutletIntervals ?? 2 * input.gridIntervals;
  if (![inletIntervals, outletIntervals].every(n => bounded(n, 4, 256, true)))
    throw new Error('Invalid inlet/wake interval count.');
  const baseTubes = Array.from({ length: input.elements.length + 1 }, (_, g) =>
    g === 0 ? input.gridLowerTubes ?? input.gridTubes
      : g === input.elements.length ? input.gridUpperTubes ?? input.gridTubes : input.gridGapTubes ?? input.gridTubes);
  // Validate even an unused gap override on a single-element case.
  if (![...baseTubes, input.gridGapTubes ?? input.gridTubes].every(n => bounded(n, 3, 31, true)))
    throw new Error('Invalid region tube count.');
  if (!['supplied', 'curvature'].includes(input.gridSurfaceSpacing)) throw new Error('Invalid surface distribution.');
  if (input.gridSurfaceSpacing === 'curvature') {
    const c = input.gridCurvatureSpacing;
    if (!c || !bounded(c.exponent, .05, 4)
      || ![c.leadingSpacingRatio, c.trailingSpacingRatio].every(n => bounded(n, .01, 2)))
      throw new Error('Invalid curvature controls.');
  }
  if (input.gridStagnationAspectRatio !== undefined && !bounded(input.gridStagnationAspectRatio, .1, 20))
    throw new Error('Invalid stagnation aspect ratio.');
  if (input.quadBoundaryLayers && (!Number.isFinite(input.reynolds) || !Number.isFinite(input.ncrit)
    || !(input.reynolds > 0) || !(input.ncrit > 0)
    || input.transitionMode !== 'automatic' || input.materialTrips?.length !== input.elements.length
    || input.materialTrips.some(p => !Array.isArray(p) || p.length !== 2 || p.some(n => !bounded(n, Number.MIN_VALUE, 1)))))
    throw new Error('Invalid coupled physical controls.');
  return { kind: 'requested-logical-counts', performed: 'parameter validation and topology planning only',
    inletIntervals, outletIntervals, requestedSideIntervals: input.gridIntervals,
    chordExponent: input.gridChordExponent, baseTubes, baseStreamlines: baseTubes.map(n => n + 1),
    regions: baseTubes.length, smoothingRequested: input.gridEllipticSmoothing,
    actualMesh: null, geometricAdmissibility: 'not-run', gasAdmissibility: 'not-run', convergence: 'not-run',
    countCaveat: 'Actual passage station counts and normal allocations require geometry construction; requested counts are not an actual cell estimate.' };
}

// Fail when visible values change without extending this independent inventory.
export function assertGridRobustnessGuiCoverage(html) {
  const expected = { 'grid-intervals': gridRobustnessControls.surfaceIntervals,
    'grid-tubes': gridRobustnessControls.tubes, 'grid-inlet': gridRobustnessControls.inletIntervals,
    'grid-outlet': gridRobustnessControls.outletIntervals };
  for (const [id, values] of Object.entries(expected)) {
    const select = html.match(new RegExp(`<select\\b[^>]*\\bid="${id}"[^>]*>([\\s\\S]*?)<\\/select>`));
    if (!select) throw new Error(`Missing GUI select ${id}.`);
    const actual = [...select[1].matchAll(/<option\b[^>]*\bvalue="([^"]+)"/g)].map(m => m[1]);
    if (JSON.stringify(actual) !== JSON.stringify(values.map(String))) throw new Error(`Grid robustness coverage differs from GUI ${id}.`);
  }
  const ranges = { 'grid-upper-streamlines': [4, 32], 'grid-lower-streamlines': [4, 32],
    'grid-gap-streamlines': [4, 32], 'grid-chord-exponent': [0, 1],
    'grid-le-ratio': [.01, 2], 'grid-te-ratio': [.01, 2],
    'grid-curvature-exponent': [.05, 4], 'grid-aspect-ratio': [.1, 20] };
  for (const [id, limits] of Object.entries(ranges)) {
    const tag = html.match(new RegExp(`<input\\b[^>]*\\bid="${id}"[^>]*>`))?.[0];
    if (!tag || limits.some((n, k) => Number(tag.match(new RegExp(`\\b${k ? 'max' : 'min'}="([^"]+)"`))?.[1]) !== n))
      throw new Error(`Grid robustness coverage differs from GUI ${id}.`);
  }
  if (!/<input\b[^>]*\bid="grid-elliptic"[^>]*\btype="checkbox"|<input\b[^>]*\btype="checkbox"[^>]*\bid="grid-elliptic"/.test(html))
    throw new Error('Missing smoothing checkbox.');
  return true;
}
