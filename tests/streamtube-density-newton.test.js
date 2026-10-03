import test from 'node:test';
import assert from 'node:assert/strict';
import { proposeDensityNewton } from '../src/euler/streamtube-density-newton.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveLinear } from '../src/numerics/linear.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';

const near = (a, b, tol = 2e-14) => assert.ok(Math.abs(a - b) < tol * Math.max(1, Math.abs(b)), `${a} != ${b}`);
const mock = {
  layout: { n: 6, densityCount: 3, globals: { stagnation: [5] }, bodies: [{ leadingIndex: 1 }],
    positions: [{ column: 3 }, { column: 4 }] },
  curves: [{ length: 2 }],
  decode: () => ({ nodes: [Array(3).fill([{ x: 0, y: -.1 }]), Array(3).fill([{ x: 0, y: .1 }])] }),
};

test('additive physical density uses one scalar for all variables, preserving the factor-two bounds', () => {
  for (const [d, step, cause] of [
    [[4, -.5, .1, 3, -1, 0], .25, 'density-increase'],
    [[1, -4, .1, 3, -1, 0], .125, 'density-decrease'],
    [[.2, -.1, .1, 3, -1, .01], 1, 'full-step'],
  ]) {
    const rho = [1.2, .8, 2], state = Float64Array.of(...rho.map(Math.log), .1, .2, .3);
    const before = state.slice(), direction = Float64Array.from(d);
    const p = proposeDensityNewton(mock, state, direction);
    assert.equal(p.step, step); assert.equal(p.limiter.kind, cause);
    for (let c = 0; c < 3; c++) {
      near(Math.exp(p.x[c]), rho[c] + step * rho[c] * d[c]);
      const ratio = Math.exp(p.x[c]) / rho[c]; assert.ok(ratio >= .5 - 1e-14 && ratio <= 2 + 1e-14);
    }
    for (let c = 3; c < 6; c++) near(p.x[c], state[c] + step * d[c]);
    assert.deepEqual(state, before); assert.deepEqual([...direction], d);
    near(p.undampedUpdate.relativeDensity.maximum, Math.max(...d.slice(0, 3).map(Math.abs)));
  }
});

test('stagnation limiter uses the appendix step cap and applies it to density and normal motion too', () => {
  const d = Float64Array.of(2, -1, 0, 3, -1, .3);
  const p = proposeDensityNewton(mock, new Float64Array(6), d);
  near(p.step, .25); assert.equal(p.limiter.kind, 'stagnation-motion');
  near(p.stagnation[0].spacing, .1); near(p.stagnation[0].acceptedChange, .15);
  near(Math.exp(p.x[0]), 1.5); near(Math.exp(p.x[1]), .75); near(p.x[3], .75);
  assert.ok(p.undampedUpdate.relativeDensity.maximum > 1);
  const prose = proposeDensityNewton(mock, new Float64Array(6), d, { stagnationLimiter: 'prose' });
  near(prose.step, 1 / 12); near(prose.stagnation[0].acceptedChange, .05);
  near(Math.exp(prose.x[0]), 1 + 1 / 6); near(prose.x[3], .25);
  assert.equal(prose.stagnationLimiter, 'prose');
  const reduced = proposeDensityNewton(mock, new Float64Array(6), d, { maximumStep: .125 });
  near(reduced.step, .125); near(Math.exp(reduced.x[0]), 1.25); near(Math.exp(reduced.x[1]), .875);
  near(reduced.x[3], .375); near(reduced.stagnation[0].acceptedChange, .075);
  assert.deepEqual(reduced.undampedUpdate, p.undampedUpdate);
  assert.throws(() => proposeDensityNewton(mock, new Float64Array(6), d, { maximumStep: 0 }), /maximum/);
  assert.throws(() => proposeDensityNewton(mock, new Float64Array(6), d, { stagnationLimiter: 'unknown' }), /limiter/);
});

test('small physical-density corrections retain log1p accuracy and invalid inputs fail explicitly', () => {
  const x = new Float64Array(6), d = Float64Array.of(1e-18, -1e-18, 0, 0, 0, 0);
  const p = proposeDensityNewton(mock, x, d);
  assert.equal(p.x[0], 1e-18); assert.equal(p.x[1], -1e-18);
  assert.throws(() => proposeDensityNewton(mock, x, [NaN, 0, 0, 0, 0, 0]), /Invalid/);
  assert.throws(() => proposeDensityNewton(mock, x, [0]), /Invalid/);
});

test('log-coordinate Newton direction matches an independently differenced physical-density system', () => {
  for (const streamwiseMode of ['momentum', 'isentropic']) {
    const input = { ...intrinsicBodyFixture({ elements: 1, bodySegments: 4, tubes: 2, alpha: .25 }), streamwiseMode };
    const system = createStreamtubeBodySystem(input), n = system.layout.n, nd = system.layout.densityCount;
    const x = system.initial.map((v, i) => i < nd ? .001 * Math.sin(i + .3) : v);
    const physical = x.map((v, i) => i < nd ? Math.exp(v) : v);
    const evaluate = u => system.residual(u.map((v, i) => i < nd ? Math.log(v) : v));
    const jac = new Float64Array(n * n), h = 1e-7;
    for (let c = 0; c < n; c++) {
      const a = physical.slice(), b = physical.slice(); a[c] += h; b[c] -= h;
      const fa = evaluate(a), fb = evaluate(b);
      for (let r = 0; r < n; r++) jac[r * n + c] = (fa[r] - fb[r]) / (2 * h);
    }
    const rhs = evaluate(physical).map(v => -v), physicalDelta = solveLinear(jac, rhs);
    const logDelta = solveLinear(system.jacobian(x), rhs), proposal = proposeDensityNewton(system, x, logDelta);
    for (let c = 0; c < n; c++) {
      near(logDelta[c] * (c < nd ? physical[c] : 1), physicalDelta[c], 3e-7);
      near(c < nd ? Math.exp(proposal.x[c]) : proposal.x[c], physical[c] + proposal.step * physicalDelta[c], 3e-7);
    }
  }
});

test('controlled two-element density-Newton root retains independent momentum and energy balances', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const system = createStreamtubeBodySystem(input);
  const result = solveStreamtubeBody(system, { stepMethod: 'density-newton', maxIterations: 16, tolerance: 1e-11 });
  assert.equal(result.converged, true, `${result.reason}, R=${result.diagnostics.residual}`);
  assert.equal(result.projectedSteps, false); assert.equal(result.secondOrderSteps, false);
  assert.equal(result.surfaces.length, 4); assert.ok(result.history.slice(1).every(h => h.stepKind === 'density-newton'));
  assert.ok(result.linearDiagnostics.maxRelativeResidual < 1e-10);
  assert.equal(streamtubeMeshSnapshot({ system, nodes: result.nodes }).quality.valid, true);
  const balance = directBodyConservation(result, input.bodies, system.conditions);
  for (const v of balance.balance) assert.ok(Math.abs(v) < 2e-9);
  for (const v of balance.cutTraction) assert.ok(Math.abs(v) < 2e-9);
});

test('density-Newton iteration limits and unsupported backends cannot certify convergence', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ bodySegments: 4, tubes: 2 }));
  const result = solveStreamtubeBody(system, { stepMethod: 'density-newton', maxIterations: 0 });
  assert.equal(result.converged, false); assert.equal(result.reason, 'iteration limit');
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'density-newton', projectedSteps: true }), /Projected Euler/);
  assert.throws(() => solveStreamtubeBody(system, { stepMethod: 'density-newton', jacobianBackend: 'finite-difference' }), /compressible analytic/);
});

test('capture updates preserve positive passage masses for adjacent moving levels and either fixed primary', () => {
  for (const primaryBody of [0, 1]) {
    const system = createStreamtubeBodySystem({ ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), primaryBody });
    const state = system.initial, before = system.decode(state), scale = before.captured.at(-1) - before.captured[0];
    for (const sign of [-1, 1]) {
      const direction = new Float64Array(state.length);
      system.layout.globals.capture.forEach((col, b) => { if (col !== null) direction[col] = sign * (b % 2 ? -4 : 4); });
      const proposal = proposeDensityNewton(system, state, direction), after = system.decode(proposal.x);
      assert.equal(proposal.limiter.kind, 'passage-mass-decrease');
      assert.ok(proposal.step > 0 && proposal.step < 1);
      for (let g = 0; g < before.captured.length - 1; g++) {
        const mass = before.captured[g + 1] - before.captured[g];
        const next = after.captured[g + 1] - after.captured[g];
        assert.ok(next >= .5 * mass - 1e-14 * scale);
        assert.ok(after.allocation.groups[g].every(t => t.massFlow > 0));
      }
      assert.equal(after.captured[primaryBody + 1], before.captured[primaryBody + 1]);
      near(after.allocation.totalMass, before.allocation.totalMass);
      const small = proposeDensityNewton(system, state, direction, { maximumStep: proposal.step / 4 });
      assert.equal(small.step, proposal.step / 4);
      assert.equal(small.limiter.kind, 'maximum-step');
      assert.deepEqual(system.decode(state).captured, before.captured);
    }
  }
});


test('a passage between two moving capture levels limits their relative motion', () => {
  const captured = [0, .3, .4, .7, 1];
  const system = { ...mock, layout: { ...mock.layout, globals: { stagnation: [], capture: [3, 4, null] } },
    decode: () => ({ nodes: [], captured }) };
  const x = new Float64Array(6), direction = Float64Array.of(.2, 0, 0, 1, -1, 0);
  const p = proposeDensityNewton(system, x, direction);
  near(p.step, .025);
  assert.deepEqual(p.limiter, { kind: 'passage-mass-decrease', group: 1 });
  near(captured[2] + p.x[4] - captured[1] - p.x[3], .05);
  near(Math.exp(p.x[0]), 1 + .2 * p.step);
});
