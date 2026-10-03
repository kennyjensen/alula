import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { quadCoupledGridFailure, quadCoupledDisplayedStateLabel, quadCoupledStartingFlowLabel } from '../src/ui/quad-coupled-failure.js';

function fixture(tubes = [2, 3, 2]) {
  const input = { mach: .74, alpha: 2.68, reynolds: 2.7e6, ncrit: 4,
    gridIntervals: 32, gridTubes: 11, gridEllipticSmoothing: true,
    transitionMode: 'automatic', materialTrips: [[1, 1], [1, 1]],
    elements: [{ name: 'Main' }, { name: 'Flap' }], geometrySource: { id: 'test' } };
  const raw = { converged: false, stateConverged: false, physicalAcceptance: false,
    reason: 'Invalid final displacement grid', conditions: { mach: .2 },
    initialization: { attempts: [{ converged: false }] },
    machContinuation: { sourceMach: .2, actualMach: .2, targetMach: .74,
      coldBaselineUsed: true, reachedTarget: false },
    solverInput: { bodies: [{ element: 1, leadingIndex: 1, trailingIndex: 4 },
      { element: 0, leadingIndex: 2, trailingIndex: 6 }] },
    mesh: { topology: 'intrinsic-quadrilateral-streamtubes',
      cells: Array.from({ length: 8 * tubes.reduce((a, b) => a + b, 0) }, () => [0, 1, 2, 3]),
      initialization: { tubes, streamwiseSegments: 8 },
      quality: { valid: false, invalidCells: [28], minCornerSine: -.02, minArea: 1e-5 } } };
  return { raw, input };
}

test('cold baseline rejection identifies the failed Mach and first wake bank without calling it converged', () => {
  const { raw, input } = fixture(), before = structuredClone({ raw, input });
  const d = quadCoupledGridFailure(raw, input);
  assert.equal(d.stage, 'cold-coupled'); assert.equal(d.actualMach, .2); assert.equal(d.targetMach, .74);
  assert.match(d.message, /Cold coupled startup at Mach 0\.200 failed/);
  assert.match(d.message, /Flap upper first wake cell/);
  assert.doesNotMatch(d.message, /retained converged/);
  assert.deepEqual(d.locations[0], { cell: 28, resolved: true, group: 1, interval: 4, tube: 0,
    boundaries: [{ body: 0, element: 1, name: 'Flap', side: 'upper', region: 'wake', wakeInterval: 0 }] });
  assert.equal(d.solverReason, raw.reason); assert.equal(d.returnedStateConverged, false);
  assert.deepEqual({ raw, input }, before);
  d.invalidCells[0] = -1; d.case.materialTrips[0][0] = .1;
  assert.deepEqual({ raw, input }, before);
});

test('logical locations preserve unequal passage counts, both body mappings and interior ambiguity', () => {
  const { raw, input } = fixture(); raw.mesh.quality.invalidCells = [7, 24, 29];
  const d = quadCoupledGridFailure(raw, input);
  assert.deepEqual(d.locations.map(p => [p.group, p.interval, p.tube]), [[0, 3, 1], [1, 2, 2], [1, 4, 1]]);
  assert.deepEqual(d.locations[0].boundaries.map(b => [b.element, b.side, b.region]), [[1, 'lower', 'surface']]);
  assert.deepEqual(d.locations[1].boundaries.map(b => [b.element, b.side, b.region]), [[0, 'lower', 'surface']]);
  assert.deepEqual(d.locations[2].boundaries, []);
  const single = fixture([2, 1, 2]); single.raw.mesh.quality.invalidCells = [20];
  assert.deepEqual(quadCoupledGridFailure(single.raw, single.input).locations[0].boundaries
    .map(b => [b.body, b.side, b.region]), [[0, 'upper', 'wake'], [1, 'lower', 'surface']]);
});

test('missing or inconsistent topology never guesses a physical cell location', () => {
  for (const change of [r => delete r.mesh.initialization, r => r.mesh.cells.pop(),
    r => r.mesh.topology = 'other', r => r.mesh.quality.invalidCells = [1000]]) {
    const { raw, input } = fixture(); change(raw);
    const d = quadCoupledGridFailure(raw, input);
    assert.equal(d.locations[0].resolved, false); assert.match(d.message, /location unavailable/);
  }
  const { raw, input } = fixture(); delete raw.solverInput;
  const d = quadCoupledGridFailure(raw, input);
  assert.equal(d.locations[0].resolved, true); assert.deepEqual(d.locations[0].boundaries, []);
});

test('stage requires recorded route evidence; Mach mismatch alone is insufficient', () => {
  const { raw, input } = fixture(); delete raw.machContinuation; delete raw.initialization;
  assert.equal(quadCoupledGridFailure(raw, input).stage, 'coupled');
  raw.machContinuation = { sourceMach: .2, targetMach: .74, attempts: [{ mach: .4 }] };
  assert.equal(quadCoupledGridFailure(raw, input).stage, 'mach-continuation');
  raw.refinement = { parentUnknowns: 10 };
  assert.equal(quadCoupledGridFailure(raw, input).stage, 'refinement');
  delete raw.refinement; raw.failure = { stage: 'boundary-layer-initialization' };
  assert.equal(quadCoupledGridFailure(raw, input).stage, 'boundary-layer-initialization');
});

test('valid grids gain no failure fields; absent rejected-cell lists remain unresolved', () => {
  const { raw, input } = fixture(); raw.mesh.quality.valid = true; raw.converged = true;
  assert.equal(quadCoupledGridFailure(raw, input), null);
  raw.mesh.quality.valid = false; raw.converged = false; delete raw.mesh.quality.invalidCells;
  const d = quadCoupledGridFailure(raw, input);
  assert.deepEqual(d.locations, []); assert.doesNotMatch(d.message, /First rejected cell/);
  assert.equal('sourceRevision' in d, false);
});

test('target-failure state wording distinguishes failed startup from accepted retained flow', () => {
  assert.equal(quadCoupledDisplayedStateLabel({ stateConverged: true, mesh: { quality: { valid: true } } }), 'retained converged state');
  for (const stateConverged of [false, undefined, true]) {
    const raw = { stateConverged, mesh: { quality: { valid: false } }, failureDiagnostic: { stage: 'cold-coupled' } };
    assert.equal(quadCoupledDisplayedStateLabel(raw), 'failed startup state');
  }
  assert.equal(quadCoupledDisplayedStateLabel({ stateConverged: false, mesh: { quality: { valid: true } } }), 'displayed state');
  assert.equal(quadCoupledDisplayedStateLabel({}), 'displayed state');
  const app = fs.readFileSync(new URL('../src/ui/app.js', import.meta.url), 'utf8');
  assert.match(app, /Showing the \$\{quadCoupledDisplayedStateLabel\(next\)\}/);
});

test('starting-flow row cannot label a failed cold baseline as converged', () => {
  const failed = { stateConverged: false, mesh: { quality: { valid: false } }, failureDiagnostic: { stage: 'cold-coupled' } };
  assert.equal(quadCoupledStartingFlowLabel(failed, .2), 'Failed startup at Mach 0.200');
  assert.equal(quadCoupledStartingFlowLabel({ stateConverged: true, mesh: { quality: { valid: true } } }, .2), 'Converged Mach 0.200');
  assert.equal(quadCoupledStartingFlowLabel({}, .2), 'Source Mach 0.200');
  assert.equal(quadCoupledStartingFlowLabel({}, undefined), 'Source condition unavailable');
  const app = fs.readFileSync(new URL('../src/ui/app.js', import.meta.url), 'utf8');
  assert.match(app, /Starting flow<\/dt><dd>\$\{quadCoupledStartingFlowLabel\(next, continuation.sourceMach\)\}/);
});

// Exercise the real display composition using manufactured profiles and
// fixed algebraic presentation stubs. No physical kernel or solve is called.
async function displayModule(path) {
  const failureUrl = new URL('../src/ui/quad-coupled-failure.js', import.meta.url).href;
  const stubs = `
    import { quadCoupledGridFailure } from '${failureUrl}';
    import { quadCoupledGrid, quadCoupledGridResult } from '${new URL('../src/ui/quad-coupled-grid.js', import.meta.url).href}';
    import { quadCoupledNcrit, quadCoupledNcritResult } from '${new URL('../src/ui/quad-coupled-ncrit.js', import.meta.url).href}';
    const createContourCurve = () => ({ evaluate: x => ({ point: { x, y: 0 } }) });
    const createSurfaceContourCurve = createContourCurve;
    const prepareContour = x => x;
    const createIntegralKernel = () => ({ station: () => ({ hk: 2, cf: 0, rho: 1 }) });
    const isentropicState = () => ({ cp: 0 });
    const quadCoupledCoefficients = () => ({ cl: 1, cm: 0, cd: .01, warnings: [] });
    const quadCoupledTransonicCoefficientsFromResult = quadCoupledCoefficients;
    const createContourArc = () => { throw Error('Unused automatic marker'); };
    const streamtubeFlowSnapshot = () => ({ captured: true });
  `;
  const source = fs.readFileSync(path, 'utf8').replace(/^import .*;\n/gm, '');
  return import(`data:text/javascript;base64,${Buffer.from(stubs + source).toString('base64')}`);
}

function displayFixture() {
  const { raw, input } = fixture([2, 2]);
  input.elements = [{ name: 'Airfoil' }]; input.materialTrips = [[1, 1]];
  Object.assign(raw, { model: 'research-streamtube-euler-bl', solverLength: 1, referenceChord: 1,
    kernelReynolds: 1e6, mach: .2, alpha: 0, families: { euler: 1e-12, boundaryLayer: 1e-12, edgeMatching: 1e-12 },
    x: [1, 2], history: [{ iteration: 0 }], initialization: { euler: { gridSmoothing: {} }, attempts: [{}] },
    solverInput: { bodies: [{ element: 0, points: [], leadingIndex: 1, trailingIndex: 4 }] },
    boundaryLayer: { stations: [{ id: 0, i: 2, s: 1, theta: .001, deltaStar: .002, ue: 1, aux: .03, regime: 'turbulent' },
      { id: 1, i: 2, s: 1, theta: .001, deltaStar: .002, ue: 1, aux: .03, regime: 'turbulent' }],
      surfaces: [{ body: 0, side: 'upper', ids: [0], tripParameter: .5 },
        { body: 0, side: 'lower', ids: [1], tripParameter: .5 }], wakes: [] } });
  const nodes = Array.from({ length: 2 }, (_, g) => Array.from({ length: 9 }, (_, i) =>
    Array.from({ length: 3 }, (_, j) => ({ x: i, y: g * 3 + j }))));
  raw.flow = { nodes, undisplacedNodes: nodes, diagnostics: { maxMach: .3 } };
  raw.mesh.quality.invalidCells = [24];
  return { raw, input };
}

test('display preserves original reason, rejection, quality and exportable diagnostic', async () => {
  const { quadCoupledResultForDisplay } = await displayModule(new URL('../src/ui/quad-coupled-result.js', import.meta.url));
  const { raw, input } = displayFixture(), before = structuredClone({ raw, input });
  const shown = quadCoupledResultForDisplay(raw, input);
  assert.equal(shown.converged, false); assert.equal(shown.physicalAcceptance, false);
  assert.equal(shown.solverReason, 'Invalid final displacement grid');
  assert.match(shown.reason, /Cold coupled startup at Mach 0\.200/);
  assert.match(shown.reason, /Airfoil upper first wake cell/);
  assert.deepEqual(shown.mesh.quality, raw.mesh.quality);
  assert.deepEqual(shown.diagnostics.failure, shown.failureDiagnostic);
  assert.deepEqual(JSON.parse(JSON.stringify(shown)).failureDiagnostic, JSON.parse(JSON.stringify(shown.failureDiagnostic)));
  assert.deepEqual({ raw, input }, before);
});

test('converged display output remains exactly equal to archived adapter', async () => {
  const current = await displayModule(new URL('../src/ui/quad-coupled-result.js', import.meta.url));
  const previous = await displayModule(new URL('../docs/rae2822/failure-diagnostics/quad-coupled-result.before.js.txt', import.meta.url));
  const { raw, input } = displayFixture(); raw.converged = true; raw.stateConverged = true;
  raw.reason = 'residual'; raw.mesh.quality = { valid: true, invalidCells: [], minCornerSine: .2, minArea: 1e-5 };
  raw.machContinuation.reachedTarget = true; raw.machContinuation.targetMach = .2;
  input.mach = .2;
  assert.deepEqual(current.quadCoupledResultForDisplay(raw, input), previous.quadCoupledResultForDisplay(raw, input));
});
