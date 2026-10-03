import test from 'node:test';
import assert from 'node:assert/strict';
import { freestream, physicalFlux, hllc } from '../src/euler/gas.js';
import { buildMesh, channelMesh } from '../src/euler/mesh.js';
import { multielementMesh } from '../src/euler/multielement-mesh.js';
import { eulerResidual } from '../src/euler/residual.js';
import { solveEuler, encodeStates, decodeStates } from '../src/euler/solve.js';
import { eulerJacobian } from '../src/euler/jacobian.js';
import { solveStreamlineChannel, createStreamlineSystem } from '../src/euler/streamline.js';
import { naca4, transform, signedArea } from '../src/geometry/airfoil.js';
import { annularVortex, isentropicNozzle, weightedError } from './oracles/euler.js';

const close = (a, b, tol = 2e-12) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b} (tol ${tol})`);
const converged = r => { assert.equal(r.converged, true, r.reason); assert.ok(r.diagnostics.residual < 1e-9); assert.ok(r.diagnostics.wallLeakage < 1e-12); assert.ok(r.diagnostics.relativeMassImbalance < 1e-8); };

test('Euler flux consistency, contact preservation, rotation and owner/neighbor antisymmetry', () => {
  const left = { rho: 1.1, u: 0.8, v: -0.2, p: 1.3 }; const right = { rho: 0.6, u: -0.1, v: 0.3, p: 0.8 };
  for (const angle of [0, 0.3, 1.2, 2.6]) {
    const nx = Math.cos(angle); const ny = Math.sin(angle);
    const f = hllc(left, right, nx, ny); const reverse = hllc(right, left, -nx, -ny);
    f.forEach((v, i) => close(v, -reverse[i]));
    hllc(left, left, nx, ny).forEach((v, i) => close(v, physicalFlux(left, nx, ny)[i]));
    const rot = s => ({ ...s, u: s.u * nx - s.v * ny, v: s.u * ny + s.v * nx });
    const expected = hllc(left, right, 1, 0); const rotated = hllc(rot(left), rot(right), nx, ny);
    [expected[0], expected[1] * nx - expected[2] * ny, expected[1] * ny + expected[2] * nx, expected[3]].forEach((v, i) => close(v, rotated[i]));
  }
  const contact = hllc({ rho: 1, u: 0, v: 0.7, p: 1 }, { rho: 2, u: 0, v: -0.2, p: 1 }, 1, 0);
  contact.forEach((v, i) => close(v, [0, 1, 0, 0][i]));
  assert.throws(() => hllc({ ...left, p: -1 }, right, 1, 0), /Nonphysical/);
  assert.throws(() => hllc(left, right, 2, 0), /unit/);
});

test('Euler exact normal-shock states satisfy all four Rankine–Hugoniot flux balances', () => {
  const gamma = 1.4; const mach = 2;
  const left = { rho: 1, p: 1 / (gamma * mach * mach), u: 1, v: 0.2 };
  const densityRatio = (gamma + 1) * mach * mach / ((gamma - 1) * mach * mach + 2);
  const pressureRatio = 1 + 2 * gamma / (gamma + 1) * (mach * mach - 1);
  const right = { rho: densityRatio, p: left.p * pressureRatio, u: 1 / densityRatio, v: left.v };
  physicalFlux(left, 1, 0).forEach((v, i) => close(v, physicalFlux(right, 1, 0)[i]));
  // This verifies the conservation equations, not shock location/capture.
});

test('signed Euler mesh rejects inverted cells and missing/false boundary tags', () => {
  const mesh = channelMesh({ nx: 3, ny: 2 });
  const tags = mesh.faces.filter(f => f.boundary).map(f => ({ ...f.boundary, a: f.a, b: f.b }));
  const cells = mesh.cells.map(c => c.vertices);
  assert.throws(() => buildMesh(mesh.vertices, [cells[0].toReversed(), ...cells.slice(1)], tags), /Inverted/);
  assert.throws(() => buildMesh(mesh.vertices, cells, tags.slice(1)), /Untagged/);
  const interior = mesh.faces.find(f => f.neighbor !== null);
  assert.throws(() => buildMesh(mesh.vertices, cells, [...tags, { a: interior.a, b: interior.b, type: 'wall' }]), /Internal/);
});

test('uniform Euler flow is preserved on distorted grids and shared fluxes cancel for nonuniform states', () => {
  const mesh = channelMesh({ nx: 7, ny: 3, map: ({ x, y, i, j, nx, ny }) => ({
    x: x + 0.018 * Math.sin(Math.PI * i / nx) * Math.sin(Math.PI * j / ny),
    y: y + 0.035 * Math.sin(2 * Math.PI * i / nx) * Math.sin(Math.PI * j / ny) }) });
  const reference = freestream();
  const uniform = eulerResidual(mesh, mesh.cells.map(() => reference), reference);
  assert.ok(uniform.diagnostics.residual < 5e-15);
  const varying = mesh.cells.map(c => ({ rho: 1 + .03 * c.x, u: 1 + .02 * c.y, v: .01 * c.x, p: reference.p + .03 * c.y }));
  const r = eulerResidual(mesh, varying, reference);
  assert.ok(r.diagnostics.sharedFluxCancellation < 2e-14);
  assert.ok(r.diagnostics.residual > .001);
  assert.ok(mesh.diagnostics.metricClosure < 2e-16);
});

test('simultaneous Euler state/grid Newton closes every moving streamline mass equation', () => {
  const mesh = channelMesh({ nx: 10, ny: 3, upper: x => 1 - .08 * Math.sin(Math.PI * x / 2) ** 2 });
  const original = structuredClone(mesh.vertices); const progress = [], frames = [];
  const r = solveStreamlineChannel(mesh, { onIteration: it => progress.push(it), onMesh: frame => frames.push(frame) });
  converged(r); assert.ok(r.diagnostics.streamlineMassResidual < 1e-9);
  assert.ok(Math.max(...r.displacements.map(Math.abs)) > .005);
  assert.deepEqual(mesh.vertices, original);
  assert.equal(progress.length, r.history.length);
  assert.equal(frames.length, r.history.length);
  assert.deepEqual(frames[0].mesh.vertices, original);
  assert.deepEqual(frames.at(-1).mesh.vertices, r.mesh.vertices);
  for (let i = 1; i < frames.length; i++) {
    assert.equal(frames[i].iteration.iteration, r.history[i].iteration);
    assert.equal(frames[i].iteration.residual, r.history[i].residual);
    assert.notDeepEqual(frames[i].mesh.vertices, frames[i - 1].mesh.vertices);
    assert.ok(frames[i].mesh.diagnostics.minArea > 0);
    assert.equal(frames[i].iteration.maximumNodeMovement, Math.max(...frames[i].mesh.vertices.map((p, j) =>
      Math.hypot(p.x - frames[i - 1].mesh.vertices[j].x, p.y - frames[i - 1].mesh.vertices[j].y))));
  }
  assert.ok(r.mesh.diagnostics.minArea > 0);
  const retry = solveStreamlineChannel(mesh, { initial: mesh.cells.map(() => ({ rho: 1.002, u: .98, v: .001, p: freestream().p + .01 })) });
  converged(retry); r.states.forEach((s, i) => { for (const k of ['rho', 'u', 'v', 'p']) close(s[k], retry.states[i][k], 2e-8); });
});

test('face-local Euler Jacobian matches independent full-residual directional differences, including boundary conditions', () => {
  const reference = freestream();
  for (const mesh of [channelMesh({ nx: 4, ny: 2, upper: x => 1 - .05 * Math.sin(Math.PI * x / 2) ** 2 }), multielementMesh([naca4('0012',20)], { rows: 1 })]) {
    // Distinct states avoid the nondifferentiable equal-wave-speed ties in
    // min/max HLLC bounds; at a tie a unique classical Jacobian need not exist.
    const states = mesh.cells.map((_, i) => ({ rho: 1 + .02 * Math.sin(.77 * i + .1), u: 1 + .04 * Math.sin(1.13 * i + .4),
      v: .02 + .04 * Math.cos(.93 * i + .2), p: reference.p + .05 * Math.cos(.53 * i + .2) }));
    const x = encodeStates(states, reference); const n = x.length;
    const j = eulerJacobian(mesh, states, reference);
    for (const phase of [.4, 1.7, 3.5]) {
      const d = Float64Array.from(x, (_, i) => Math.sin((i + 1) * phase)); const h = 2e-6;
      const plus = eulerResidual(mesh, decodeStates(x.map((v, i) => v + h * d[i]), states.length, reference), reference).residual;
      const minus = eulerResidual(mesh, decodeStates(x.map((v, i) => v - h * d[i]), states.length, reference), reference).residual;
      for (let row = 0; row < n; row++) {
        let product = 0; for (let col = 0; col < n; col++) product += j[row * n + col] * d[col];
        close(product, (plus[row] - minus[row]) / (2 * h), 2e-7);
      }
    }
  }
});

test('flow/grid residual cross-blocks are nonzero and inverted candidate grids are inadmissible', () => {
  const mesh = channelMesh({ nx: 6, ny: 3, upper: x => 1 - .08 * Math.sin(Math.PI * x / 2) ** 2 });
  const sys = createStreamlineSystem(mesh, { initial: mesh.cells.map(c => ({ rho: 1, u: 1 + .02 * c.x, v: .02 * c.y, p: freestream().p + .03 * c.x })) });
  const { flow, grid } = sys.unknowns;
  assert.equal(sys.residual(sys.initial).length, flow + grid);
  const shifted = sys.initial.slice(); shifted[flow + 4] += .01;
  const base = sys.residual(sys.initial); const r = sys.residual(shifted);
  // A nonuniform state makes internal-face metrics affect the Euler balance;
  // a uniform state would cancel exactly by geometric conservation.
  assert.ok(Math.max(...r.slice(0, flow).map((v, i) => Math.abs(v - base[i]))) > 1e-5);
  const velocity = sys.initial.slice(); velocity[4 * mesh.faces[sys.faceIds[4]].owner + 2] += .01;
  const coupled = sys.residual(velocity);
  assert.ok(Math.max(...coupled.slice(flow).map((v, i) => Math.abs(v - base[flow + i]))) > 1e-3);
  const inverted = sys.initial.slice(); inverted[flow] = 10;
  assert.equal(sys.admissible(inverted), false);
});

test('independent exact 2-D compressible vortex has decreasing pressure/velocity errors under refinement', () => {
  const errors = [];
  for (const [nx, ny] of [[6, 2], [12, 4], [24, 8]]) {
    const { mesh, exact } = annularVortex(nx, ny);
    const r = solveEuler(mesh, { initial: mesh.cells.map(exact), alpha: 120 }); converged(r);
    errors.push({ p: weightedError(r, exact, 'p'), u: weightedError(r, exact, 'u') });
  }
  for (let i = 1; i < errors.length; i++) { assert.ok(errors[i].p < .6 * errors[i - 1].p); assert.ok(errors[i].u < .7 * errors[i - 1].u); }
  assert.ok(errors.at(-1).p < .014); assert.ok(errors.at(-1).u < .015);
});

test('slender nozzle tends toward the isentropic quasi-1D relation with grid refinement', () => {
  let previous = Infinity;
  for (const nx of [12, 24, 48]) {
    const r = solveEuler(channelMesh({ nx, ny: 2, upper: x => isentropicNozzle(x).area })); converged(r);
    const error = weightedError(r, c => isentropicNozzle(c.x), 'u');
    assert.ok(error < previous * .8); previous = error;
  }
  assert.ok(previous < .008); // Includes finite-height departure from quasi-1D.
});

test('multielement mesh preserves polygon area, hole count, walls, and shared inlet/wake cuts', () => {
  const contours = [naca4('0012',20), transform(naca4('0012',20), { chord: .4, x: 1.05, y: -.2, angle: -10 })];
  const mesh = multielementMesh(contours);
  const b = mesh.bounds;
  close(mesh.cells.reduce((s, c) => s + c.area, 0), (b.xmax - b.xmin) * (b.ymax - b.ymin) - contours.reduce((s, c) => s + signedArea(c), 0));
  assert.equal(mesh.vertices.length - mesh.faces.length + mesh.cells.length, 1 - contours.length);
  assert.ok(mesh.cuts.some(c => c.type === 'wake-cut'));
  for (const cut of mesh.cuts) { assert.notEqual(mesh.faces[cut.face].neighbor, null); assert.equal(mesh.faces[cut.face].boundary, null); }
  for (let e = 0; e < contours.length; e++) {
    const walls = mesh.faces.filter(f => f.boundary?.element === e);
    close(walls.reduce((s, f) => s + f.nx * f.length, 0), 0);
    close(walls.reduce((s, f) => s + f.ny * f.length, 0), 0);
  }
  assert.throws(() => multielementMesh([naca4('0012',20), transform(naca4('0012',20), { x: 1.1 })]), /cuts cross/);
  const flatNose = naca4('0012',20);
  flatNose.splice(10, 1, { x: -.001, y: .002 }, { x: -.001, y: -.002 });
  assert.throws(() => multielementMesh([flatNose]), /Vertical body edges/);
});

test('two bodies interact in one Euler solve; reflection and element ordering preserve forces', () => {
  const lower = transform(naca4('0012',20), { y: -.2 }); const upper = transform(naca4('0012',20), { y: .2 });
  const a = solveEuler(multielementMesh([lower, upper], { rows: 2 })); converged(a);
  close(a.wallForces[0].x, a.wallForces[1].x, 2e-10);
  close(a.wallForces[0].y, -a.wallForces[1].y, 2e-10);
  assert.ok(a.wallForces[0].y > .01); // Accelerated interelement gap attracts the bodies.
  const b = solveEuler(multielementMesh([upper, lower], { rows: 2 })); converged(b);
  close(a.wallForces[0].y, b.wallForces[1].y, 2e-10);
  const isolated = solveEuler(multielementMesh([lower], { rows: 2 })); converged(isolated);
  close(isolated.wallForces[0].y, 0, 2e-10);
});

test('Euler iteration limits remain explicit failures and zero Mach is rejected', () => {
  const mesh = channelMesh({ nx: 6, ny: 2, upper: x => 1 - .05 * Math.sin(Math.PI * x / 2) ** 2 });
  const failed = solveEuler(mesh, { maxIterations: 0 });
  assert.equal(failed.converged, false); assert.equal(failed.reason, 'iteration limit');
  assert.throws(() => solveEuler(mesh, { mach: 0 }), /positive Mach/);
  assert.throws(() => solveStreamlineChannel(channelMesh({ ny: 1 })), /internal rows/);
});
