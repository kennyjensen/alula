// SPDX-License-Identifier: GPL-2.0-or-later
// Fast specification/acceptance checks. Flow qualification is performed by
// check-coupled-gui-default.js and must actually pass for every catalog entry.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { buildGuiDefaultCoupledCase, guiDefaultCoupledCases, guiDefaultCoupledControls,
  guiDefaultGeometrySummary, compactGuiDefaultIteration, assessGuiDefaultCoupledResult } from '../scripts/validation/gui-default-coupled-cases.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { prepareAirfoilElement } from '../src/geometry/airfoil-element.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';

test('exact reported GUI defaults preserve point sets, finite-base corners and explicit ISMOM4', () => {
  assert.deepEqual(guiDefaultCoupledCases.map(c => c.id), ['flap', 'single', 'three', 'rae2822', 'nlr7301', '30p30n']);
  const expected = { flap: [162, 162], single: [162], three: [162, 162, 162], rae2822: [129], nlr7301: [457, 249], '30p30n': [201, 221, 243] };
  for (const { id, preset } of guiDefaultCoupledCases) {
    const c = buildGuiDefaultCoupledCase(id);
    assert.deepEqual(c.elements.map(e => e.points.length), expected[id]);
    assert.equal(c.alpha, 4); assert.equal(c.mach, .2); assert.equal(c.reynolds, 1e6); assert.equal(c.ncrit, 9);
    assert.equal(c.referenceChord, 1); assert.equal(c.eulerIsmom, 4); assert.equal(c.flowModel, 'streamtube-grid'); assert.equal(c.quadBoundaryLayers, true);
    assert.equal(c.gridIntervals, 16); assert.equal(c.gridTubes, 7); assert.equal(c.gridCrosslinePlacement, 'potential');
    assert.equal(c.gridSurfaceSpacing, 'automatic'); assert.equal(c.gridChordExponent, 0); assert.equal(c.gridSmoothingMethod, 'elliptic');
    assert.equal(c.gridEllipticSmoothing, true); assert.equal(c.transitionMode, 'automatic');
    assert.deepEqual(c.materialTrips, c.elements.map(() => [1, 1]));
    for (const key of ['gridInletIntervals', 'gridOutletIntervals', 'gridUpperTubes', 'gridLowerTubes', 'gridGapTubes', 'gridStagnationAspectRatio', 'coupledNativeHk'])
      assert.equal(Object.hasOwn(c, key), false, `GUI Auto/off must stay omitted: ${key}`);
    if (['flap', 'single', 'three'].includes(id)) for (const element of c.elements) {
      assert.equal(element.sourcePoints.length, 161);
      assert.deepEqual(element.points.at(-1), element.points[0]);
      assert.equal(element.trailingEdge.lowerIndex, 160);
    }
    const benchmark = getBenchmarkAirfoil(preset);
    if (benchmark) assert.deepEqual(c.elements, benchmark.elements.map(prepareAirfoilElement));
  }
  const nlr = guiDefaultGeometrySummary(buildGuiDefaultCoupledCase('nlr7301'));
  assert.deepEqual(nlr.map(e => e.distinctPoints), [456, 248]);
  assert.deepEqual(nlr.map(e => e.trailingEdge), [{ kind: 'finite-base', upperIndex: 0, lowerIndex: 424 }, { kind: 'finite-base', upperIndex: 0, lowerIndex: 216 }]);
  assert.deepEqual(guiDefaultCoupledControls, { maxIterations: 40, eulerMaxIterations: 20, tolerance: 1e-10, maxStartupAttempts: 2 });
  const changed = buildGuiDefaultCoupledCase('nlr7301'); changed.elements[0].points[0].x = 99; changed.materialTrips[0][0] = .1;
  assert.notEqual(buildGuiDefaultCoupledCase('nlr7301').elements[0].points[0].x, 99);
  assert.deepEqual(buildGuiDefaultCoupledCase('nlr7301').materialTrips, [[1, 1], [1, 1]]);
  assert.throws(() => buildGuiDefaultCoupledCase('rae2822-mses'), /Unknown GUI default/);
});

test('catalog matches selected HTML default controls instead of a historical solver default', () => {
  const html = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8');
  const presets = html.match(/<select\b[^>]*id="preset"[^>]*>([\s\S]*?)<\/select>/)[1];
  assert.deepEqual([...presets.matchAll(/<option\b[^>]*value="([^"]+)"/g)].map(m => m[1]).filter(id => id !== 'custom'),
    guiDefaultCoupledCases.map(c => c.preset));
  const tag = id => html.match(new RegExp(`<input\\b[^>]*id="${id}"[^>]*>`))?.[0];
  for (const [id, value] of [['alpha-number', 4], ['quad-mach', .2], ['quad-reynolds', 1e6], ['quad-ncrit', 9], ['quad-trip-upper', 1], ['quad-trip-lower', 1], ['grid-chord-exponent', 0]])
    assert.equal(Number(tag(id)?.match(/\bvalue="([^"]+)"/)?.[1]), value, id);
  const selected = id => {
    const select = html.match(new RegExp(`<select\\b[^>]*id="${id}"[^>]*>([\\s\\S]*?)<\\/select>`))?.[1];
    const options = [...select.matchAll(/<option\b([^>]*)>/g)];
    return (options.find(m => /\bselected\b/.test(m[1])) ?? options[0])[1].match(/\bvalue="([^"]+)"/)[1];
  };
  for (const [id, value] of [['resolution', '160'], ['euler-ismom', '4'], ['quad-transition', 'automatic'], ['grid-intervals', '16'], ['grid-tubes', '7'], ['grid-surface-spacing', 'automatic'], ['grid-inlet', 'auto'], ['grid-outlet', 'auto']])
    assert.equal(selected(id), value, id);
  assert.match(tag('grid-elliptic'), /\bchecked\b/); assert.doesNotMatch(tag('grid-match-aspect'), /\bchecked\b/);
});

test('general reliability catalog defaults to GUI ISMOM4 and requires an explicit legacy designation for omission', () => {
  const spec = { preset: 'flap', mode: 'streamtube-bl', changes: {} };
  assert.equal(buildReliabilityCase(spec).caseData.eulerIsmom, 4);
  assert.equal(Object.hasOwn(buildReliabilityCase({ ...spec, legacyEulerEquations: true }).caseData, 'eulerIsmom'), false);
  assert.equal(buildReliabilityCase({ ...spec, changes: { eulerIsmom: 2 } }).caseData.eulerIsmom, 2);
  assert.throws(() => buildReliabilityCase({ ...spec, changes: { eulerIsmom: undefined } }), /legacyEulerEquations/);
  assert.throws(() => buildReliabilityCase({ ...spec, legacyEulerEquations: true, changes: { eulerIsmom: 4 } }), /without explicit ISMOM/);
});

function accepted() {
  const input = buildGuiDefaultCoupledCase('flap');
  const station = { theta: .001, deltaStar: .002, ue: 1 };
  const displayed = { converged: true, status: 'research-converged', mach: .2, alpha: 4, cl: .8, cd: .02, cm: -.1,
    materialTrips: structuredClone(input.materialTrips),
    families: { euler: 1e-11, boundaryLayer: 2e-11, edgeMatching: 1e-12 }, mesh: { quality: { valid: true } },
    boundaryLayer: { surfaces: input.elements.flatMap((_, element) => ['upper', 'lower'].map(side => ({ element, side, stations: [station, station] }))),
      wakes: input.elements.map((_, element) => ({ element, stations: [station, station] })) } };
  const raw = { ...displayed, referenceReynolds: 1e6, referenceChord: 1, solverLength: .95, kernelReynolds: 950000,
    conditions: { reynolds: 950000, ncrit: 9, transitionMode: 'automatic' }, solverInput: { hybrid: { ismom: 4 } },
    materialTrips: [[1, 1], [1, 1]], checkpoint: { restart: { options: { reynolds: 950000, ncrit: 9 } } } };
  return { input, displayed: structuredClone(displayed), raw: structuredClone(raw) };
}
test('automatic default acceptance rejects stalled, retained-condition, wrong-equation and incomplete-BL results', () => {
  let { input, raw, displayed } = accepted(); assert.equal(assessGuiDefaultCoupledResult(input, raw, displayed).passed, true);
  const changes = [
    [p => { p.displayed.converged = false; p.displayed.status = 'unconverged'; p.displayed.families.boundaryLayer = 2.73; }, 'equations'],
    [p => { p.displayed.mach = .1; }, 'requestedMach'],
    [p => { p.raw.checkpoint.restart.options.ncrit = 7; }, 'requestedNcrit'],
    [p => { p.raw.referenceReynolds = 2e6; }, 'requestedReynolds'],
    [p => { p.raw.referenceChord = 2; }, 'requestedReynolds'],
    [p => { p.raw.kernelReynolds = 1e6; }, 'reynoldsNormalization'],
    [p => { p.raw.checkpoint.restart.options.reynolds = 1e6; }, 'reynoldsNormalization'],
    [p => { p.raw.solverLength = NaN; }, 'reynoldsNormalization'],
    [p => { p.raw.solverInput.hybrid.ismom = 2; }, 'selectedEquations'],
    [p => { p.raw.gridSequence = { reachedTarget: false }; }, 'requestedGrid'],
    [p => { p.displayed.mesh.quality.valid = false; }, 'convexGrid'],
    [p => { p.displayed.boundaryLayer.wakes.pop(); }, 'completeBoundaryLayers'],
    [p => { p.raw.materialTrips[1][0] = .1; }, 'materialTrips'],
  ];
  for (const [change, failed] of changes) {
    const packet = accepted(); change(packet);
    const report = assessGuiDefaultCoupledResult(packet.input, packet.raw, packet.displayed);
    assert.equal(report.passed, false); assert.ok(report.failures.includes(failed), `${failed}: ${report.failures}`);
  }
});

test('compact observations retain the limiter and distinguish a near-zero accepted step from initialization', () => {
  const h = { stage: 'coupled', iteration: 40, residual: 2.73, step: 1e-16, actualNcrit: 5, targetNcrit: 9,
    limiter: { kind: 'density', row: 123 },
    changes: [{ body: 0, from: 5, to: 6 }], projection: { active: true, displacementChanges: [{ id: 3, correction: .1 }] },
    maintenance: { geometryRedistribution: true, triggeredBodies: [],
      passages: [{ group: 1, correctionScale: .25, maxDisplacement: .001, referenceBank: 0, fixedBanks: [true, true] }] } };
  const before = structuredClone(h), compact = compactGuiDefaultIteration(h);
  assert.equal(compact.tinyAcceptedStep, true); assert.equal(compact.residual, 2.73);
  assert.equal(compact.actualNcrit, 5); assert.equal(compact.targetNcrit, 9);
  assert.deepEqual(compact.limiter, h.limiter); assert.deepEqual(compact.projection.displacementChanges, h.projection.displacementChanges);
  assert.equal(compact.maintenance.geometryRedistribution, true);
  assert.deepEqual(compact.maintenance.passages, h.maintenance.passages);
  assert.deepEqual(h, before); assert.equal(compactGuiDefaultIteration({ iteration: 0, step: 0 }).tinyAcceptedStep, false);
});
