// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { solverModes, solverReliabilityCases, buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { assessSolverResult } from '../scripts/validation/solver-reliability-acceptance.js';
import { runBoundedWorker } from '../scripts/validation/bounded-worker.js';

test('inventory covers every preset in all four modes, including every solver choice and reported conditions', () => {
  const ids = solverReliabilityCases.map(c => c.id); assert.equal(new Set(ids).size, ids.length);
  for (const preset of ['single', 'flap', 'three', 'rae2822', 'rae2822-mses', 'nlr7301', '30p30n'])
    for (const mode of solverModes) assert.ok(solverReliabilityCases.some(c => c.id === `${preset}-${mode}`));
  const request = buildReliabilityCase(solverReliabilityCases.find(c => c.id === 'rae2822-user-mach074')).caseData;
  assert.equal(request.alpha, 2.68); assert.equal(request.mach, .74); assert.equal(request.reynolds, 2.7e6);
  assert.equal(request.ncrit, 4); assert.equal(request.transitionMode, 'automatic'); assert.deepEqual(request.materialTrips, [[1, 1]]);
  const fineRAE = buildReliabilityCase(solverReliabilityCases.find(c =>
    c.id === 'rae2822-mses-user-mach074-64x11-automatic-slor')).caseData;
  assert.equal(fineRAE.geometrySource.id, 'rae2822-mses');
  assert.deepEqual([fineRAE.alpha, fineRAE.mach, fineRAE.reynolds, fineRAE.ncrit], [2.68, .74, 2.7e6, 4]);
  assert.deepEqual([fineRAE.gridIntervals, fineRAE.gridTubes, fineRAE.gridEllipticSmoothing], [64, 11, true]);
  assert.equal(fineRAE.transitionMode, 'automatic'); assert.deepEqual(fineRAE.materialTrips, [[1, 1]]);
  assert.equal(fineRAE.gridInletIntervals, undefined); assert.equal(fineRAE.gridOutletIntervals, undefined);
  for (const mode of ['fixed', 'automatic']) {
    const nlr = buildReliabilityCase(solverReliabilityCases.find(c => c.id === `nlr7301-user-mach0185-${mode}`)).caseData;
    assert.equal(nlr.alpha, 6); assert.equal(nlr.mach, .185); assert.equal(nlr.reynolds, 2.51e6);
    assert.equal(nlr.transitionMode, mode === 'fixed' ? 'fixed-trip' : mode);
    assert.equal(nlr.gridEllipticSmoothing, false);
    assert.deepEqual(nlr.materialTrips, Array.from({ length: 2 }, () => mode === 'automatic' ? [1, 1] : [.05, .05]));
  }
  for (const spec of solverReliabilityCases) assert.equal(buildReliabilityCase(spec).guiRestriction, '', spec.id);
});

test('small residual cannot hide an invalid grid, retained wrong Mach or missing BL', () => {
  const input = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .74, alpha: 4,
    reynolds: 1e6, ncrit: 9, materialTrips: [[1, 1], [1, 1]], elements: [{}, {}] };
  const stations = [{ theta: .001, deltaStar: .003, ue: 1 }, { theta: .002, deltaStar: .004, ue: .9 }];
  const result = { converged: true, mach: .74, alpha: 4, cl: 1, cm: -.1, cd: .01, families: { euler: 1e-12, bl: 1e-12, edge: 1e-12 },
    referenceReynolds: 1e6, conditions: { ncrit: 9 }, materialTrips: [[1, 1], [1, 1]],
    mesh: { quality: { valid: true } }, boundaryLayer: { surfaces: Array.from({ length: 4 }, (_, i) => ({ stations, element: Math.floor(i / 2), side: i % 2 ? 'lower' : 'upper' })), wakes: [{ stations, element: 0 }, { stations, element: 1 }] } };
  assert.equal(assessSolverResult(input, result).passed, true);
  for (const [change, expected] of [[{ mach: .2 }, 'requestedMach'], [{ alpha: 0 }, 'requestedIncidence'],
    [{ referenceReynolds: 3e6 }, 'requestedReynolds'], [{ referenceReynolds: undefined }, 'requestedReynolds'],
    [{ conditions: { ncrit: 4 } }, 'requestedNcrit'], [{ ncritContinuation: { reachedTarget: false } }, 'requestedNcrit'],
    [{ materialTrips: [[.05, .05], [1, 1]] }, 'requestedTrips'], [{ materialTrips: undefined }, 'requestedTrips'],
    [{ diagnostics: { equationResidual: 1e-13 }, families: { boundaryLayer: .2 } }, 'equations'],
    [{ machContinuation: { reachedTarget: false } }, 'requestedMach'], [{ mesh: { quality: { valid: false } } }, 'convexGrid'],
    [{ boundaryLayer: { surfaces: result.boundaryLayer.surfaces.slice(0, 2), wakes: [{ stations }] } }, 'completeBoundaryLayers'],
    [{ converged: false, coefficientStatus: 'unconverged' }, 'terminalConvergence'], [{ cd: null }, 'coefficientsFinite'],
    [{ status: 'solved', converged: false }, 'terminalConvergence'],
    [{ boundaryLayer: { surfaces: result.boundaryLayer.surfaces.map(s => ({ ...s, element: 0 })), wakes: result.boundaryLayer.wakes } }, 'completeBoundaryLayers']]) {
    const assessment = assessSolverResult(input, { ...result, ...change });
    assert.equal(assessment.passed, false); assert.ok(assessment.failures.includes(expected));
  }
});

test('panel/BL cannot pass until both coupled equations and wake directions converge', () => {
  const input = { flowModel: 'coupled', mach: 0, alpha: 4, elements: [{}] }, stations = [{ theta: .001, deltaStar: .003, ue: 1 }, { theta: .002, deltaStar: .004, ue: .9 }];
  const result = { status: 'solved', mach: 0, alpha: 4, cl: 1, cm: -.1, cd: .01,
    diagnostics: { equationResidual: 1e-10, wakeResidual: 2e-3 },
    boundaryLayer: { surfaces: [{ stations, element: 0, side: 'upper' }, { stations, element: 0, side: 'lower' }], wakes: [{ stations, element: 0 }] } };
  assert.deepEqual(assessSolverResult(input, result).failures, ['wakeConvergence']);
});

test('watchdog interrupts a synchronous infinite loop and preserves its last iteration', async () => {
  const outcome = await runBoundedWorker(new URL('./fixtures/reliability/stuck-worker.js', import.meta.url), {}, { timeoutMs: 1000 });
  assert.equal(outcome.type, 'timeout'); assert.equal(outcome.history.length, 1); assert.equal(outcome.history[0].residual, 3);
  assert.ok(outcome.seconds < 5);
});

test('invalid BL thickness cannot pass a finite-gap station check', () => {
  const input = { flowModel: 'coupled', mach: 0, alpha: 4, elements: [{}] };
  const stations = [{ theta: .01, deltaStar: .005, wakeGap: .006, ue: 1 }, { theta: .01, deltaStar: .03, ue: 1 }];
  const result = { status: 'solved', mach: 0, alpha: 4, cl: 1, cm: 0, cd: .01,
    diagnostics: { equationResidual: 1e-12, wakeResidual: 1e-9 },
    boundaryLayer: { surfaces: [{ element: 0, side: 'upper', stations }, { element: 0, side: 'lower', stations }], wakes: [{ element: 0, stations }] } };
  assert.ok(assessSolverResult(input, result).failures.includes('completeBoundaryLayers'));
});

test('progress persistence failure returns retained diagnostics and stops the worker', async () => {
  const outcome = await runBoundedWorker(new URL('./fixtures/reliability/stuck-worker.js', import.meta.url), {}, {
    timeoutMs: 1000, onProgress: () => { throw new Error('Simulated disk full'); },
  });
  assert.equal(outcome.code, 'PROGRESS_CAPTURE_FAILED'); assert.equal(outcome.history[0].residual, 3);
});

test('acceptance recognizes the retained complete automatic multielement root without solving it again', () => {
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
  const result = assessSolverResult(saved.input, saved.result);
  assert.equal(result.passed, true, JSON.stringify(result));
  assert.equal(result.physicalValidation, 'not evaluated by this numerical gate');
});

test('current defaults use automatic transition and smoothing while explicit controls stay independent', () => {
  for (const preset of ['single', 'flap', 'three', 'rae2822', 'rae2822-mses', 'nlr7301', '30p30n']) {
    for (const mode of ['streamtube-grid', 'streamtube-bl']) {
      const input = buildReliabilityCase({ preset, mode }).caseData;
      assert.equal(input.gridEllipticSmoothing, true); assert.equal(input.gridIntervals, 16); assert.equal(input.gridTubes, 7);
      if (mode === 'streamtube-bl') {
        assert.equal(input.transitionMode, 'automatic'); assert.equal(input.ncrit, 9);
        assert.deepEqual(input.materialTrips, input.elements.map(() => [1, 1]));
      } else assert.equal(input.transitionMode, undefined);
    }
  }
  const explicit = buildReliabilityCase({ preset: 'flap', mode: 'streamtube-bl',
    changes: { transitionMode: 'fixed-trip', gridEllipticSmoothing: false } }).caseData;
  assert.equal(explicit.transitionMode, 'fixed-trip'); assert.equal(explicit.gridEllipticSmoothing, false);
  assert.deepEqual(explicit.materialTrips, [[.05, .05], [.05, .05]]);
  const confirmed = buildReliabilityCase(solverReliabilityCases.find(c => c.id === 'nlr7301-user-mach0185-automatic-slor')).caseData;
  assert.deepEqual([confirmed.alpha, confirmed.mach, confirmed.reynolds, confirmed.ncrit], [6, .185, 2.51e6, 9]);
  assert.equal(confirmed.gridEllipticSmoothing, true); assert.equal(confirmed.transitionMode, 'automatic');
  assert.deepEqual(confirmed.materialTrips, [[1, 1], [1, 1]]);
  for (const intervals of [16, 64]) {
    const reported = buildReliabilityCase(solverReliabilityCases.find(c =>
      c.id === `nlr7301-user-mach0185-${intervals}x11-automatic-slor`)).caseData;
    assert.deepEqual([reported.alpha, reported.mach, reported.reynolds, reported.ncrit], [6, .185, 2.51e6, 9]);
    assert.deepEqual([reported.gridIntervals, reported.gridTubes, reported.gridEllipticSmoothing], [intervals, 11, true]);
    assert.equal(reported.quadBoundaryLayers, true); assert.equal(reported.transitionMode, 'automatic');
    assert.deepEqual(reported.materialTrips, [[1, 1], [1, 1]]);
    assert.equal(reported.gridInletIntervals, undefined); assert.equal(reported.gridOutletIntervals, undefined);
  }
  const fine = buildReliabilityCase(solverReliabilityCases.find(c => c.id === 'nlr7301-user-mach0185-64x11')).caseData;
  assert.deepEqual([fine.gridIntervals, fine.gridTubes, fine.gridEllipticSmoothing], [64, 11, false]);
  for (const mode of ['streamtube-grid', 'streamtube-bl']) {
    const ids = [`flap-${mode}-slor-off`, `flap-${mode}-32x11-slor-off`, `nlr7301-${mode}-slor-off`];
    for (const id of ids) assert.equal(buildReliabilityCase(solverReliabilityCases.find(c => c.id === id)).caseData.gridEllipticSmoothing, false);
  }
});

test('generated reliability presets retain the app finite trailing-edge bases in every mode', () => {
  for (const preset of ['single', 'flap', 'three']) for (const mode of solverModes) {
    const { caseData } = buildReliabilityCase({ preset, mode });
    for (const element of caseData.elements) {
      assert.equal(element.trailingEdge.kind, 'finite-base');
      const upper = element.points[element.trailingEdge.upperIndex];
      const lower = element.points[element.trailingEdge.lowerIndex];
      assert(Math.hypot(upper.x - lower.x, upper.y - lower.y) > 0);
      assert.deepEqual(element.points.at(-1), element.points[0]);
    }
    if (preset === 'single') assert(Math.abs(caseData.elements[0].points[0].y - .00126) < 1e-14);
  }
});

test('the transonic robustness inventory retains reported cases and failing nearby grids', () => {
  const cases = solverReliabilityCases.filter(c => c.id.startsWith('rae2822-mses-robustness-'));
  assert.equal(cases.length, 15);
  const reported = buildReliabilityCase(cases.find(c => c.id.endsWith('-16x7'))).caseData;
  for (const [key, value] of Object.entries({ alpha: 2.68, mach: .74, reynolds: 1e6, ncrit: 9,
    gridIntervals: 16, gridTubes: 7, gridInletIntervals: 16, gridOutletIntervals: 16, eulerIsmom: 4 }))
    assert.equal(reported[key], value, key);
  assert.deepEqual(reported.materialTrips, [[1, 1]]);
  assert.equal(reported.elements[0].points.length, 129);
  for (const id of ['32x7', '32x9', '8x7', '16x11', '32x11', '8x9', '16x7-mach072', '16x7-mach076']) assert.ok(cases.some(c => c.id.endsWith(`-${id}`)));
  for (const [id, gridIntervals, gridTubes] of [['32x7-re271n4', 32, 7], ['64x7-re271n4', 64, 7], ['64x9-re271n4', 64, 9]]) {
    const input = buildReliabilityCase(cases.find(c => c.id.endsWith(`-${id}`))).caseData;
    for (const [key, value] of Object.entries({ alpha: 2.68, mach: .74, reynolds: 2.71e6, ncrit: 4,
      gridIntervals, gridTubes, gridInletIntervals: 16, gridOutletIntervals: 16, eulerIsmom: 4 }))
      assert.equal(input[key], value, `${id}: ${key}`);
    assert.deepEqual(input.elements, reported.elements);
    assert.deepEqual(input.materialTrips, [[1, 1]]);
  }
});

test('every preset has three default-grid viscous operating points with automatic inlet and outlet counts', () => {
  for (const preset of ['single', 'flap', 'three', 'rae2822', 'rae2822-mses', 'nlr7301', '30p30n']) {
    const ids = [`${preset}-streamtube-bl`, `default-robustness-${preset}-alpha0`, `default-robustness-${preset}-re3m`];
    const inputs = ids.map(id => {
      const spec = solverReliabilityCases.find(c => c.id === id);
      assert(spec, id);
      return buildReliabilityCase(spec).caseData;
    });
    assert.deepEqual(inputs.map(i => [i.alpha, i.mach, i.reynolds, i.ncrit]),
      [[4, .2, 1e6, 9], [0, .2, 1e6, 9], [2, .3, 3e6, 4]]);
    for (const input of inputs) {
      assert.equal(input.gridIntervals, 16);
      assert.equal(input.gridTubes, 7);
      assert.equal(input.gridEllipticSmoothing, true);
      assert.equal(input.gridSurfaceSpacing, 'automatic');
      assert.equal(input.gridInletIntervals, undefined);
      assert.equal(input.gridOutletIntervals, undefined);
      assert.equal(input.transitionMode, 'automatic');
      assert.deepEqual(input.materialTrips, input.elements.map(() => [1, 1]));
    }
  }
});
