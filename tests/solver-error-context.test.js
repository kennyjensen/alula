import test from 'node:test';
import assert from 'node:assert/strict';
import { buildSolverErrorContext } from '../src/ui/solver-error-context.js';

function rae() {
  return { geometrySource: { id: 'rae2822-mses', sourceHash: 'exact-geometry' }, elements: [{ name: 'RAE 2822', points: [{ x: 1, y: 0 }, { x: 0, y: .02 }] }],
    flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .74, alpha: 2.68, reynolds: 2700000, ncrit: 4,
    transitionMode: 'automatic', materialTrips: [[1, 1]], gridIntervals: 64, gridTubes: 11,
    gridEllipticSmoothing: true, gridSmoothingMethod: 'elliptic', gridSurfaceSpacing: 'supplied', gridChordExponent: 0 };
}

test('error context separates requested Ncrit from the retained checkpoint and drops stale-stage guesses', () => {
  const input = { ...rae(), ncrit: 9 }, event = { stage: 'coupled-ncrit-startup', actualNcrit: 8.5, targetNcrit: 10,
    checkpoint: { restart: { options: { ncrit: 8 } } }, conditions: { ncrit: 8 },
    ncritContinuation: { actualNcrit: 8, targetNcrit: 9, reachedTarget: false } };
  const { report, summary } = buildSolverErrorContext({ input, event, message: 'target stopped',
    progress: { stage: 'coupled-ncrit-startup', actualNcrit: 9, targetNcrit: 9, iteration: 3 } });
  assert.equal(report.actualNcrit, 8); assert.equal(report.targetNcrit, 9);
  assert.match(summary, /Requested Ncrit 9 · Actual Ncrit 8/);
  const unknown = buildSolverErrorContext({ input, event: { stage: 'boundary-layer-initialization' }, message: 'no state',
    progress: { stage: 'coupled-ncrit-startup', actualNcrit: 8, targetNcrit: 9 } });
  assert.equal(unknown.report.actualNcrit, undefined);
  assert.match(unknown.summary, /Requested Ncrit 9 · Actual Ncrit unknown/);
});

test('failed actual Mach and stage beat stale progress while requested Mach remains the submitted value', () => {
  const input = rae(), message = 'Body cell i=124, group=0, tube=2: Nonpositive interface pressure.';
  const { summary, report } = buildSolverErrorContext({ message, input,
    event: { actualMach: .2, stage: 'cold-euler', code: 'STARTUP', diagnostics: { reason: message, section: { i: 124, group: 0, tube: 2 } } },
    progress: { actualMach: .74, stage: 'coupled', iteration: 18, residual: .3 } });
  assert.equal(report.actualMach, .2); assert.equal(report.requestedMach, .74); assert.equal(report.stage, 'cold-euler');
  assert.equal(report.iteration, undefined); assert.equal(report.progress.iteration, 18);
  assert.equal(report.error.message, message); assert.equal(report.error.code, 'STARTUP');
  assert.match(summary, /Requested Mach 0\.74 · Actual Mach 0\.2/);
  assert.match(summary, /Stage: cold-euler · Iteration: unknown/);
  assert.match(summary, /Inlet Auto · Wake Auto/); assert.match(summary, /ISMOM: Auto/);
});

test('whole failed results supply precise failure context and only compact mesh/stop diagnostics', () => {
  const mesh = { cells: Array.from({ length: 5 }, () => [0, 1, 2, 3]), vertices: [{ x: 0, y: 1 }],
    initialization: { streamwiseSegments: 7, tubes: [10, 13, 10], fullHistory: ['not captured'] },
    quality: { valid: false, invalidCells: [2470], minCornerSine: -.01 } };
  const event = { mach: .74, conditions: { mach: .185 }, failureDiagnostic: { actualMach: .185, stage: 'cold-coupled', invalidCells: [2470] },
    failure: { stage: 'cold-coupled' }, iteration: 0, mesh, solverStopReason: 'Invalid final displacement grid',
    diagnostics: { reason: 'invalid grid', residual: new Float64Array([1, 2]), flow: { giant: true }, mesh, section: { i: 90 } },
    lastRejectedStep: { stage: 'admissibility', step: .01, message: 'bad polygon', physicalProposal: { huge: true }, x: [1, 2] },
    flow: { huge: true }, checkpoint: { huge: true } };
  const { report, summary } = buildSolverErrorContext({ input: rae(), message: 'failed', event,
    progress: { mach: .4, stage: 'coupled', iteration: 9, flow: { huge: true } } });
  assert.equal(report.actualMach, .185); assert.equal(report.stage, 'cold-coupled'); assert.equal(report.iteration, 0);
  assert.deepEqual(report.mesh, { cellCount: 5, streamwiseIntervals: 7, tubes: [10, 13, 10], quality: mesh.quality });
  assert.deepEqual(report.error.diagnostics, { reason: 'invalid grid', section: { i: 90 } });
  assert.deepEqual(report.error.lastRejectedStep, { stage: 'admissibility', message: 'bad polygon', step: .01 });
  assert.equal(report.error.solverStopReason, event.solverStopReason);
  assert.doesNotMatch(JSON.stringify(report), /giant|huge|fullHistory/);
  assert.match(summary, /invalid cells \[2470\]/);
});

test('NLR per-element trips and advanced tube counts, spacing and overrides remain explicit', () => {
  const input = { ...rae(), geometrySource: { id: 'nlr7301' }, elements: [{ name: 'Main', points: [] }, { name: 'Flap', points: [] }],
    mach: .185, alpha: 6, reynolds: 2510000, ncrit: 9, materialTrips: [[.05, 1], [.12, .7]], eulerIsmom: 3,
    gridIntervals: 16, gridTubes: 7, gridUpperTubes: 9, gridLowerTubes: 11, gridGapTubes: 14,
    gridInletIntervals: 32, gridSurfaceSpacing: 'curvature', gridCurvatureSpacing: { exponent: .5, leadingSpacingRatio: .2, trailingSpacingRatio: .4 },
    gridStagnationAspectRatio: 2.5 };
  const { summary, report } = buildSolverErrorContext({ message: 'failed', input, progress: { mach: .185, stage: 'coupled', iteration: 7 } });
  assert.equal(report.actualMach, .185); assert.equal(report.iteration, 7);
  assert.match(summary, /Main: upper 0\.05, lower 1; Flap: upper 0\.12, lower 0\.7/);
  assert.match(summary, /upper 9, lower 11, gap 14 \(tube counts\)/);
  assert.match(summary, /Inlet 32 · Wake Auto/); assert.match(summary, /ISMOM: 3/);
  assert.match(summary, /Stagnation aspect ratio: 2\.5/); assert.match(summary, /LE ratio 0\.2, TE ratio 0\.4/);
  assert.deepEqual(report.input, input);
});

test('unknown actual Mach is never filled from the requested Mach and zero values are retained', () => {
  const a = buildSolverErrorContext({ message: 'mesh failure', input: rae() });
  assert.equal(a.report.actualMach, undefined); assert.match(a.summary, /Actual Mach unknown/);
  const b = buildSolverErrorContext({ message: 'panel failure', input: { flowModel: 'inviscid', mach: .74, alpha: 0 }, event: { mach: 0, iteration: 0, stage: 'panel' } });
  assert.equal(b.report.actualMach, 0); assert.equal(b.report.iteration, 0); assert.match(b.summary, /α 0°/);
  const c = buildSolverErrorContext({ input: null, event: null, progress: null });
  assert.equal(c.report.schemaVersion, 1); assert.doesNotThrow(() => JSON.stringify(c.report));
});

test('submitted geometry, all unknown controls and error causes are detached and JSON serializable', () => {
  const input = rae(); input.futureControl = { enabled: false, level: 0, values: new Float64Array([.1, .2]) };
  const cause = new Error('original cause'); cause.code = 'CAUSE'; cause.diagnostics = { location: { i: 3 } };
  const event = { cause, diagnostics: { values: [1, 2] } }, progress = { stage: 'euler', iteration: 1 };
  const a = buildSolverErrorContext({ input, message: 'failed', event, progress });
  input.elements[0].points[0].x = 99; input.materialTrips[0][0] = .5; cause.diagnostics.location.i = 9; event.diagnostics.values.push(3); progress.iteration = 9;
  assert.equal(a.report.input.elements[0].points[0].x, 1); assert.equal(a.report.input.materialTrips[0][0], 1);
  assert.deepEqual(a.report.input.futureControl, { enabled: false, level: 0, values: [.1, .2] });
  assert.equal(a.report.error.cause.message, 'original cause'); assert.equal(a.report.error.cause.diagnostics.location.i, 3);
  assert.deepEqual(a.report.error.diagnostics.values, [1, 2]); assert.equal(a.report.progress.iteration, 1);
  assert.deepEqual(JSON.parse(JSON.stringify(a.report)), a.report);
  a.report.input.elements[0].points[0].x = -3; assert.equal(input.elements[0].points[0].x, 99);
});

test('messages and names remain literal text for the caller textContent renderer', () => {
  const message = '<img src=x onerror="bad()"> & original\nsecond line', input = rae(); input.elements[0].name = '<b>RAE & test</b>';
  const { summary, report } = buildSolverErrorContext({ input, message, label: '<script>airfoil</script>' });
  assert.equal(report.error.message, message);
  assert.match(summary, /<script>airfoil<\/script>/); assert.match(summary, /<b>RAE & test<\/b>/);
  assert.doesNotMatch(summary, /&lt;|&amp;|<div|<br/);
});

test('explicit error payload wins over conflicting derived diagnostics and keeps a compact precursor with stack', () => {
  const event = { actualMach: .1, stage: 'euler', failure: { stage: 'failed-baseline' },
    failureDiagnostic: { actualMach: .2, stage: 'cold-coupled' },
    stack: 'Error: original\n    at example.js:12',
    diagnostics: { residual: .3, precursor: { lastRejectedStep: { message: 'nested rejection', step: .01, x: [1, 2] } } },
    precursor: { reason: 'pressure failure', quality: { valid: false, invalidCells: [12] },
      lastRejectedStep: { stage: 'admissibility', message: 'pressure', step: .5,
        code: 'streamtube-interface-pressure', diagnostics: { interfacePressure: -1, cell: { i: 124, group: 0, tube: 2 } },
        physicalProposal: { giant: true } },
      flow: { giant: true }, mesh: { vertices: [] } } };
  const a = buildSolverErrorContext({ input: rae(), message: 'failed', event, progress: { actualMach: .74, stage: 'coupled', iteration: 9 } });
  assert.equal(a.report.actualMach, .1); assert.equal(a.report.stage, 'euler');
  assert.equal(a.report.error.stack, event.stack);
  assert.equal(a.report.error.diagnostics.residual, .3);
  assert.deepEqual(a.report.error.diagnostics.precursor.lastRejectedStep, { message: 'nested rejection', step: .01 });
  assert.deepEqual(a.report.error.precursor, { reason: 'pressure failure', quality: event.precursor.quality,
    lastRejectedStep: { stage: 'admissibility', message: 'pressure', step: .5,
      code: 'streamtube-interface-pressure', diagnostics: { interfacePressure: -1, cell: { i: 124, group: 0, tube: 2 } } } });
  assert.doesNotMatch(JSON.stringify(a.report), /giant|physicalProposal/);
  delete event.stage;
  assert.equal(buildSolverErrorContext({ input: rae(), event }).report.stage, 'failed-baseline');
});

test('pressure failures are readable without copying and nonfinite diagnostics survive JSON export', () => {
  for (const wrap of [diagnostics => ({ code: 'streamtube-interface-pressure', diagnostics }),
    diagnostics => ({ code: 'coupled-euler-precursor', precursor: { lastRejectedStep: { diagnostics } } }),
    diagnostics => ({ code: 'coupled-euler-precursor', diagnostics: { precursor: { lastRejectedStep: { diagnostics } } } }),
    diagnostics => ({ code: 'sonic-capacity', diagnostics: { stagnationDensityFallback: { diagnostics } } })]) {
    const diagnostics = { cell: { i: 124, group: 0, tube: 2 }, interfacePressure: { lower: -.002, upper: NaN },
      pressureDifference: Infinity, normalInertia: -Infinity };
    const { report, summary } = buildSolverErrorContext({ input: rae(), message: 'Nonpositive or nonfinite streamline interface pressure.', event: wrap(diagnostics) });
    assert.match(summary, /Failure code:/);
    assert.match(summary, /Failing cell: i=124, group=0, tube=2 \(solver indices\)/);
    assert.match(summary, /Interface pressure \(solver units\): lower -0\.002, upper NaN/);
    const json = JSON.stringify(report);
    assert.match(json, /"upper":"NaN"/); assert.match(json, /"pressureDifference":"Infinity"/);
    assert.match(json, /"normalInertia":"-Infinity"/);
    assert.deepEqual(JSON.parse(json), report);
    assert.ok(Number.isNaN(diagnostics.interfacePressure.upper));
  }
});

test('alpha diagnostics separate the retained root from the failed trial and retain recovery history', () => {
  const alphaContinuation = { actualAlpha: 2.1155, targetAlpha: 2.68, reachedTarget: false,
    stopReason: 'minimum alpha step', lastAttempt: { alpha: 2.1156, families: { boundaryLayer: 5.1e-9 },
      progress: { cause: 'two-state-cycle', transitionChanged: true } }, transitionRecoveries: [{ accepted: false }] };
  const { report, summary } = buildSolverErrorContext({ input: rae(), event: {
    alphaContinuation, checkpoint: { restart: { input: { alpha: 2.1155 } } } },
    progress: { actualAlpha: 2.1156, residual: 5.1e-9 } });
  assert.equal(report.actualAlpha, 2.1155);
  assert.equal(report.targetAlpha, 2.68);
  assert.deepEqual(report.alphaContinuation, alphaContinuation);
  assert.match(summary, /Actual α 2.1155°/);
});
