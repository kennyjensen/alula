import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';

const max = a => Math.max(...Array.from(a, Math.abs));
const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const solved = new Map();
function bodyCase(elements) {
  if (!solved.has(elements)) {
    const input = intrinsicBodyFixture({ elements }), system = createStreamtubeBodySystem(input);
    const result = solveStreamtubeBody(system, { tolerance: 1e-11, maxIterations: 15 });
    assert.equal(result.converged, true, result.reason);
    solved.set(elements, { input, system, result });
  }
  return solved.get(elements);
}

test('intrinsic single/two-body solves close all equations and independently conserve physical-boundary flux and force', () => {
  for (const elements of [1, 2]) {
    const { input, system, result } = bodyCase(elements);
    assert.equal(result.residual.length, system.layout.n);
    assert.ok(max(result.residual) < 1e-11);
    const check = directBodyConservation(result, input.bodies, system.conditions);
    for (const block of check.blocks) for (const key of ['maxLocal', 'total', 'internalCancellation'])
      assert.ok(max(block[key]) < 2e-9, JSON.stringify(block));
    assert.ok(max(check.cutTraction) < 2e-9, JSON.stringify(check));
    assert.ok(max(check.balance) < 2e-9, JSON.stringify(check));
    assert.equal(result.surfaces.length, 2 * elements);
    assert.match(result.forceStatus, /Unvalidated/);
    assert.match(result.formulation, /no BL coupling/);
    // An intentionally wrong force must be detected by the independent CV.
    const corrupted = { ...result, diagnosticForces: result.diagnosticForces.map(f => ({ ...f, cy: f.cy + .01 })) };
    assert.ok(Math.abs(directBodyConservation(corrupted, input.bodies, system.conditions).balance[2]) > .004);
  }
});

test('single-body symmetry is recovered and two-body capture/stagnation/Kutta constraints participate in the solve', () => {
  const one = bodyCase(1), two = bodyCase(2);
  close(one.result.diagnosticForces[0].cl, 0, 2e-10);
  close(one.result.diagnosticForces[0].cm, 0, 2e-10);
  close(one.result.strengths.circulation, 0, 2e-10);
  close(one.result.stagnation[0], one.system.initialStagnation[0], 2e-10);
  const lower = one.result.surfaces.find(s => s.side === 'lower').points;
  const upper = one.result.surfaces.find(s => s.side === 'upper').points;
  lower.forEach((p, i) => { close(p.cp, upper[i].cp); close(p.x, upper[i].x); close(p.y, -upper[i].y); });
  const initial = createStreamtubeBodySystem(two.input).decode(new Float64Array(two.system.layout.n));
  assert.ok(Math.abs(two.result.captured[2] - initial.captured[2]) > .005, 'Free capture must change.');
  close(two.result.captured[1], initial.captured[1]);
  close(two.result.captured.at(-1), initial.captured.at(-1));
  assert.ok(two.result.stagnation.some((s, b) => Math.abs(s - two.system.initialStagnation[b]) > 1e-4));
  for (const kind of ['leadingKutta', 'trailingKutta']) {
    const rows = two.system.layout.rows.filter(r => r.kind === kind);
    assert.equal(rows.length, 2); rows.forEach(r => assert.ok(Math.abs(two.result.residual[r.index]) < 1e-11));
  }
});

test('normal-chart rebasing and physical-grid restart preserve flow and residual; rejected restarts are atomic', () => {
  const input = intrinsicBodyFixture({ elements: 2 }), system = createStreamtubeBodySystem(input);
  const state = system.initial.map((_, i) => 1e-5 * Math.sin(i + 1)), before = system.evaluate(state);
  const rebased = system.rebase(state), after = system.evaluate(rebased);
  before.residual.forEach((r, i) => close(r, after.residual[i], 2e-13));
  system.layout.positions.forEach(p => assert.equal(rebased[p.column], 0));
  assert.deepEqual(before.nodes, after.nodes);
  const fresh = createStreamtubeBodySystem(input), restored = fresh.adoptGeometry(rebased, after.nodes);
  const again = fresh.evaluate(restored);
  after.residual.forEach((r, i) => close(r, again.residual[i], 2e-13));
  const chart = fresh.geometryChart(), invalid = structuredClone(after.nodes);
  invalid[0][2][0] = { ...invalid[0][0][0] };
  assert.throws(() => fresh.adoptGeometry(restored, invalid), /Degenerate rebased body normal/);
  assert.deepEqual(fresh.geometryChart(), chart);
  assert.deepEqual(fresh.evaluate(restored).residual, again.residual);
  const disconnected = structuredClone(after.nodes); disconnected[1][0][0] = { ...disconnected[1][0][0], y: disconnected[1][0][0].y + .001 };
  assert.throws(() => fresh.adoptGeometry(restored, disconnected), /disconnected/);
  const wall = structuredClone(after.nodes); wall[0][input.bodies[0].leadingIndex].at(-1).y += .001;
  assert.throws(() => fresh.adoptGeometry(restored, wall), /wall does not match/);
  assert.deepEqual(fresh.geometryChart(), chart);
});

test('the complete body equations and diagnostic forces are invariant under a rigid coordinate change', () => {
  const { input, result } = bodyCase(2), angle = 27, radians = angle * Math.PI / 180, c = Math.cos(radians), s = Math.sin(radians);
  const move = p => ({ x: 3 + c * p.x - s * p.y, y: -2 + s * p.x + c * p.y });
  const movedInput = { ...input, alpha: input.alpha + angle, bodies: input.bodies.map(b => ({ ...b, points: b.points.map(move) })),
    outerLower: input.outerLower.map(move), outerUpper: input.outerUpper.map(move), cutPaths: input.cutPaths.map(path => path.map(move)) };
  const system = createStreamtubeBodySystem(movedInput), nodes = result.nodes.map(group => group.map(row => row.map(move)));
  const initial = system.adoptGeometry(result.x, nodes), moved = solveStreamtubeBody(system, { initial, maxIterations: 0, tolerance: 1e-10 });
  assert.equal(moved.converged, true);
  moved.diagnosticForces.forEach((f, b) => { close(f.cl, result.diagnosticForces[b].cl, 2e-9); close(f.cd, result.diagnosticForces[b].cd, 2e-9); close(f.cm, result.diagnosticForces[b].cm, 2e-9); });
  const flux = directBodyConservation(moved, movedInput.bodies, system.conditions);
  assert.ok(max(flux.balance) < 2e-9, JSON.stringify(flux.balance));
});

test('invalid body states are rejected and an unfinished solve remains explicitly unconverged', () => {
  const input = intrinsicBodyFixture({ elements: 2 }), system = createStreamtubeBodySystem(input);
  const invalidMass = system.initial.slice(); invalidMass[system.layout.globals.capture[1]] = 1;
  assert.equal(system.admissible(invalidMass), false);
  assert.throws(() => system.evaluate(invalidMass), /mass levels/);
  const folded = system.initial.slice(); folded[system.layout.nodes[0][4][0].column] = 10;
  assert.equal(system.admissible(folded), false);
  assert.throws(() => system.evaluate(folded), /Folded|area/);
  const stagnation = system.initial.slice(); stagnation[system.layout.globals.stagnation[0]] = 1;
  assert.equal(system.admissible(stagnation), false);
  for (const mach of [0, 1, NaN]) assert.throws(() => createStreamtubeBodySystem({ ...input, mach }), /subcritical/);
  const limited = solveStreamtubeBody(system, { maxIterations: 0 });
  assert.equal(limited.converged, false); assert.equal(limited.reason, 'iteration limit');
  assert.ok(max(limited.residual) > .001);
});

test('resolving the slat-normal direction permits a cold lifting two-body solve with independent conservation', () => {
  const input = intrinsicBodyFixture({ elements: 2, alpha: 2, tubes: 5, tubeGrowth: 3 });
  const system = createStreamtubeBodySystem(input), result = solveStreamtubeBody(system, { tolerance: 1e-11, maxIterations: 15 });
  assert.equal(result.converged, true, result.reason);
  assert.ok(max(result.residual) < 1e-11);
  const flux = directBodyConservation(result, input.bodies, system.conditions);
  assert.ok(max(flux.balance) < 2e-9, JSON.stringify(flux.balance));
  assert.ok(max(flux.cutTraction) < 2e-9);
  assert.match(result.forceStatus, /Unvalidated/);
  // The baseline's fold remains recorded. Convergence after normal
  // refinement does not certify the still-coarse surface discretization.
});
