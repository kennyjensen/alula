// SPDX-License-Identifier: GPL-2.0-or-later
// Acceptance inventory, not a list of known successes. Do not suppress cases
// because a previous solver run failed.
import { builtinAirfoils } from '../../src/geometry/builtin-airfoils.js';
import { naca4Standard as naca4, transform } from '../../src/geometry/airfoil.js';
import { getBenchmarkAirfoil } from '../../src/geometry/benchmark-airfoils.js';
import { prepareAirfoilElement } from '../../src/geometry/airfoil-element.js';
import { flowModelRoute } from '../../src/ui/flow-model.js';
import { benchmarkRestriction } from '../../src/ui/tests/benchmark-permissions.js';

export const solverModes = ['inviscid', 'coupled', 'streamtube-grid', 'streamtube-bl'];
export const builtinDefinitions = builtinAirfoils;
const defaults = ['single', 'flap', 'three', 'rae2822', 'rae2822-mses', 'nlr7301', '30p30n'];
export const solverReliabilityCases = [
  ...defaults.flatMap(preset => solverModes.map(mode => ({ id: `${preset}-${mode}`, preset, mode,
    tier: preset === 'flap' ? 'smoke' : 'standard', changes: {} }))),
  // Every built-in preset, with automatic inlet/outlet counts. Pair these
  // with its ordinary streamtube-bl case to cover three default-grid points.
  ...defaults.flatMap(preset => [
    ['alpha0', { alpha: 0 }],
    ['re3m', { alpha: 2, mach: .3, reynolds: 3e6, ncrit: 4 }],
  ].map(([variant, changes]) => ({ id: `default-robustness-${preset}-${variant}`,
    preset, mode: 'streamtube-bl', tier: 'nearby', changes }))),
  // Additional cold worker qualifications at ordinary app settings. Retain
  // failed slatted cases in the same inventory as successful single elements.
  ...[
    ['single', 'negative2', { alpha: -2 }],
    ['single', 'mach04', { alpha: 2, mach: .4, reynolds: 2e6, ncrit: 4 }],
    ['flap', 'alpha2', { alpha: 2 }],
    ['flap', 'mach03', { alpha: 2, mach: .3 }],
    ['nlr7301', 'alpha2', { alpha: 2 }],
    ['nlr7301', 're2m', { alpha: 4, reynolds: 2e6, ncrit: 4 }],
    ...['rae2822', 'rae2822-mses'].flatMap(preset => [
      [preset, 'negative2', { alpha: -2 }],
      [preset, 'mach05', { alpha: 2, mach: .5, reynolds: 2e6, ncrit: 4 }],
    ]),
    ...['three', '30p30n'].flatMap(preset => [
      [preset, 'alpha2', { alpha: 2 }],
      [preset, 're2m', { alpha: 4, reynolds: 2e6, ncrit: 4 }],
    ]),
  ].map(([preset, variant, changes]) => ({ id: `app-nearby-${preset}-${variant}`,
    preset, mode: 'streamtube-bl', tier: 'app-nearby', changes })),
  ...['streamtube-grid', 'streamtube-bl'].flatMap(mode => [
    { id: `flap-${mode}-slor-off`, preset: 'flap', mode, tier: 'nearby', changes: { gridEllipticSmoothing: false } },
    { id: `flap-${mode}-32x11`, preset: 'flap', mode, tier: 'nearby', changes: { gridIntervals: 32, gridTubes: 11 } },
    { id: `flap-${mode}-32x11-slor-off`, preset: 'flap', mode, tier: 'nearby', changes: { gridIntervals: 32, gridTubes: 11, gridEllipticSmoothing: false } },
    { id: `nlr7301-${mode}-slor-off`, preset: 'nlr7301', mode, tier: 'nearby', changes: { gridEllipticSmoothing: false } },
  ]),
  ...[false, true].map(slor => ({ id: `flap-streamtube-bl-fixed-trip-slor-${slor ? 'on' : 'off'}`, preset: 'flap',
    mode: 'streamtube-bl', tier: 'nearby', changes: { transitionMode: 'fixed-trip', tripLimits: [.05, .05], gridEllipticSmoothing: slor } })),
  ...['coupled', 'streamtube-bl'].flatMap(mode => [
    { id: `flap-${mode}-re500k`, preset: 'flap', mode, tier: 'nearby', changes: { reynolds: 5e5 } },
    { id: `flap-${mode}-re3m`, preset: 'flap', mode, tier: 'nearby', changes: { reynolds: 3e6 } },
  ]),
  { id: 'rae2822-user-mach074', preset: 'rae2822', mode: 'streamtube-bl', tier: 'reported',
    changes: { alpha: 2.68, mach: .74, reynolds: 2.7e6, ncrit: 4, transitionMode: 'automatic', tripLimits: [1, 1], gridEllipticSmoothing: false } },
  { id: 'rae2822-mses-user-mach074-64x11-automatic-slor', preset: 'rae2822-mses', mode: 'streamtube-bl', tier: 'reported',
    changes: { alpha: 2.68, mach: .74, reynolds: 2.7e6, ncrit: 4, transitionMode: 'automatic',
      tripLimits: [1, 1], gridIntervals: 64, gridTubes: 11, gridEllipticSmoothing: true } },
  ...[24, 7].map(gridTubes => ({ id: `rae2822-mses-early-trips-64x${gridTubes}`,
    preset: 'rae2822-mses', mode: 'streamtube-bl', tier: gridTubes === 24 ? 'reported' : 'nearby',
    changes: { alpha: 2.68, mach: .74, reynolds: 2.71e6, ncrit: 4, transitionMode: 'automatic',
      tripLimits: [.03, .07], gridIntervals: 64, gridTubes, gridInletIntervals: 16, gridOutletIntervals: 16 } })),
  // Reported coarse transonic case and nearby controls, including known failures.
  ...[
    ['16x7', 16, 7, .74, 1e6, 9], ['16x9', 16, 9, .74, 1e6, 9],
    ['32x7', 32, 7, .74, 1e6, 9], ['32x9', 32, 9, .74, 1e6, 9],
    ['8x7', 8, 7, .74, 1e6, 9], ['16x7-mach070', 16, 7, .70, 1e6, 9],
    ['16x7-re271n4', 16, 7, .74, 2.71e6, 4],
    ['32x7-re271n4', 32, 7, .74, 2.71e6, 4],
    ['64x7-re271n4', 64, 7, .74, 2.71e6, 4],
    ['64x9-re271n4', 64, 9, .74, 2.71e6, 4],
    ['16x11', 16, 11, .74, 1e6, 9], ['32x11', 32, 11, .74, 1e6, 9],
    ['8x9', 8, 9, .74, 1e6, 9], ['16x7-mach072', 16, 7, .72, 1e6, 9],
    ['16x7-mach076', 16, 7, .76, 1e6, 9],
  ].map(([id, gridIntervals, gridTubes, mach, reynolds, ncrit]) => ({
    id: `rae2822-mses-robustness-${id}`, preset: 'rae2822-mses', mode: 'streamtube-bl',
    tier: ['16x7', '64x7-re271n4'].includes(id) ? 'reported' : 'nearby',
    changes: { alpha: 2.68, mach, reynolds, ncrit, gridIntervals, gridTubes,
      gridInletIntervals: 16, gridOutletIntervals: 16, transitionMode: 'automatic', tripLimits: [1, 1] },
  })),
  ...[
    ['rae-alpha15', 'rae2822-mses', 16, 7, .74, 1.5],
    ['rae-alpha35', 'rae2822-mses', 16, 7, .74, 3.5],
    ['naca-alpha6', 'single', 32, 9, .2, 6],
    ['flap-alpha4', 'flap', 32, 9, .2, 4],
  ].map(([id, preset, gridIntervals, gridTubes, mach, alpha]) => ({
    id: `wider-robustness-${id}`, preset, mode: 'streamtube-bl', tier: 'nearby',
    changes: { alpha, mach, reynolds: 1e6, ncrit: 9, gridIntervals, gridTubes,
      gridInletIntervals: 16, gridOutletIntervals: 16, transitionMode: 'automatic', tripLimits: [1, 1] },
  })),
  // Generalization sweep: incidence sign, Reynolds number, compressibility,
  // and transition sensitivity. These are tests, never solver routing rules.
  ...[
    ['naca-zero', 'single', 0, .2, 1e6, 9],
    ['naca-negative', 'single', -4, .2, 1e6, 9],
    ['naca-re300k', 'single', 4, .2, 3e5, 9],
    ['naca-re100k', 'single', 4, .2, 1e5, 9],
    ['naca-m05', 'single', 4, .5, 1e6, 9],
    ['naca-m065', 'single', 2, .65, 1e6, 4],
    ['naca-alpha8', 'single', 8, .2, 1e6, 9],
    ['rae-subsonic', 'rae2822-mses', 2.68, .3, 1e6, 9],
    ['rae-ncrit4', 'rae2822-mses', 2.68, .72, 1e6, 4],
    ['rae-m078', 'rae2822-mses', 2.68, .78, 1e6, 9],
  ].map(([id, preset, alpha, mach, reynolds, ncrit]) => ({
    id: `general-robustness-${id}`, preset, mode: 'streamtube-bl', tier: 'nearby',
    changes: { alpha, mach, reynolds, ncrit, gridIntervals: 16, gridTubes: 7,
      gridInletIntervals: 16, gridOutletIntervals: 16, transitionMode: 'automatic', tripLimits: [1, 1] },
  })),
  ...[
    ['naca-re3m', 'single', { alpha: 2, mach: .3, reynolds: 3e6 }],
    ['naca-negative6', 'single', { alpha: -6, mach: .2, ncrit: 4 }],
    ['naca-trip8', 'single', { alpha: 8, mach: .2, transitionMode: 'fixed-trip', tripLimits: [.05, .05] }],
    ['rae-zero', 'rae2822-mses', { alpha: 0, mach: .74 }],
    ['rae-re500k', 'rae2822-mses', { alpha: 2, mach: .7, reynolds: 5e5, ncrit: 4 }],
    ['naca-m06', 'single', { alpha: 4, mach: .6, ncrit: 4 }],
  ].map(([id, preset, changes]) => ({
    id: `progress-robustness-${id}`, preset, mode: 'streamtube-bl', tier: 'nearby',
    changes: { reynolds: 1e6, ncrit: 9, gridIntervals: 16, gridTubes: 7,
      gridInletIntervals: 16, gridOutletIntervals: 16, ...changes },
  })),
  // Preserve the three initially recorded unsmoothed controls under their
  // original IDs; the confirmed automatic + SLOR case is separate below.
  ...['fixed', 'automatic'].map(transitionMode => ({ id: `nlr7301-user-mach0185-${transitionMode}`,
    preset: 'nlr7301', mode: 'streamtube-bl', tier: 'reported',
    changes: { alpha: 6, mach: .185, reynolds: 2.51e6, ncrit: 9, transitionMode: transitionMode === 'fixed' ? 'fixed-trip' : transitionMode,
      tripLimits: transitionMode === 'automatic' ? [1, 1] : [.05, .05], gridEllipticSmoothing: false } })),
  { id: 'nlr7301-user-mach0185-automatic-slor', preset: 'nlr7301', mode: 'streamtube-bl', tier: 'reported',
    changes: { alpha: 6, mach: .185, reynolds: 2.51e6, ncrit: 9, transitionMode: 'automatic',
      tripLimits: [1, 1], gridIntervals: 16, gridTubes: 7, gridEllipticSmoothing: true } },
  ...[16, 64].map(gridIntervals => ({ id: `nlr7301-user-mach0185-${gridIntervals}x11-automatic-slor`,
    preset: 'nlr7301', mode: 'streamtube-bl', tier: 'reported',
    changes: { alpha: 6, mach: .185, reynolds: 2.51e6, ncrit: 9, transitionMode: 'automatic',
      tripLimits: [1, 1], gridIntervals, gridTubes: 11, gridEllipticSmoothing: true } })),
  { id: 'nlr7301-user-mach0185-64x11', preset: 'nlr7301', mode: 'streamtube-grid', tier: 'reported',
    changes: { alpha: 6, mach: .185, gridIntervals: 64, gridTubes: 11, gridEllipticSmoothing: false } },
];

export function buildReliabilityCase(spec) {
  if (!solverModes.includes(spec.mode)) throw new Error('Unknown reliability solver mode.');
  if (spec.legacyEulerEquations !== undefined && typeof spec.legacyEulerEquations !== 'boolean')
    throw new Error('Legacy Euler equation omission must be explicitly designated.');
  const benchmark = getBenchmarkAirfoil(spec.preset);
  const definitions = benchmark?.elements ?? builtinDefinitions[spec.preset];
  if (!definitions) throw new Error(`Unknown reliability airfoil: ${spec.preset}`);
  const { tripLimits, panels = 160, ...changes } = spec.changes ?? {};
  const elements = definitions.map(d => prepareAirfoilElement(d.points ? structuredClone(d) : {
    name: d.name, points: transform(naca4(d.code, panels), { chord: d.chord, x: d.x, y: d.y, angle: -d.deflection }),
  }));
  const quad = spec.mode.startsWith('streamtube-');
  if (spec.legacyEulerEquations && (!quad || Object.hasOwn(changes, 'eulerIsmom')))
    throw new Error('Legacy Euler equation omission requires a quad case without explicit ISMOM.');
  if (quad && Object.hasOwn(changes, 'eulerIsmom') && changes.eulerIsmom === undefined)
    throw new Error('Use legacyEulerEquations:true to reproduce historical equation omission.');
  const viscous = spec.mode === 'coupled' || spec.mode === 'streamtube-bl';
  const caseData = { elements, alpha: 4, referenceChord: 1, momentReference: { x: .25, y: 0 },
    ...flowModelRoute(spec.mode, changes.mach ?? .2),
    ...(quad ? { ...(!spec.legacyEulerEquations ? { eulerIsmom: 4 } : {}),
      gridIntervals: 16, gridTubes: 7, gridCrosslinePlacement: 'potential', gridChordExponent: 0,
      gridSurfaceSpacing: 'automatic', gridEllipticSmoothing: true, gridSmoothingMethod: 'elliptic' } : {}),
    ...(viscous ? { reynolds: 1e6, ncrit: 9 } : {}),
    ...(spec.mode === 'coupled' ? { trips: tripLimits ?? [1, 1] } : {}),
    ...(spec.mode === 'streamtube-bl' ? { transitionMode: changes.transitionMode ?? 'automatic',
      materialTrips: elements.map(() => [...(tripLimits ?? (changes.transitionMode === 'fixed-trip' ? [.05, .05] : [1, 1]))]) } : {}),
    ...changes,
    ...(benchmark ? { geometrySource: { id: benchmark.id, referenceChord: benchmark.referenceChord, ...benchmark.provenance } } : {}),
  };
  // All presets offer all four modes. Keep the compatibility field in
  // receipts so older inventories remain readable; numerical failures are
  // reported by the solver, not by preset-specific availability lists.
  return { caseData, guiRestriction: benchmarkRestriction(benchmark, spec.mode, { quadBoundaryLayers: spec.mode === 'streamtube-bl' }) };
}
