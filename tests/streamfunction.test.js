import test from 'node:test';
import assert from 'node:assert/strict';
import { makePanel } from '../src/inviscid/panel.js';
import { solveInviscid, velocityAt, vortexBasis } from '../src/inviscid/linear-vortex.js';
import { vortexStreamfunctionBasis, streamfunctionAt, vortexPotentialBasis, potentialDifference } from '../src/inviscid/streamfunction.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('linear vortex streamfunction agrees with independent log-distance integration near and far from a panel', () => {
  for (const length of [1, 1e-6]) {
    const panel = makePanel({ x: 0, y: 0 }, { x: length, y: 0 });
    for (const point of [{ x: .3 * length, y: .2 * length }, { x: -.8 * length, y: -.1 * length }, { x: 10000, y: 3000 }]) {
      const n = 20000, expected = [0, 0];
      for (let i = 0; i < n; i++) {
        const t = (i + .5) / n, psi = -length * Math.log(Math.hypot(point.x - t * length, point.y)) / (2 * Math.PI * n);
        expected[0] += (1 - t) * psi; expected[1] += t * psi;
      }
      vortexStreamfunctionBasis(point, panel).forEach((v, k) => close(v, expected[k], 5e-10 * length));
    }
  }
});

test('panel streamfunction has exact finite endpoint limits and a continuous primitive/series boundary', () => {
  const panel = makePanel({ x: 0, y: 0 }, { x: 1, y: 0 });
  for (const [point, expected] of [[panel.a, [3 / (8 * Math.PI), 1 / (8 * Math.PI)]],
    [panel.b, [1 / (8 * Math.PI), 3 / (8 * Math.PI)]], [{ x: .5, y: 0 }, Array(2).fill((1 + Math.log(2)) / (4 * Math.PI))]]) {
    vortexStreamfunctionBasis(point, panel).forEach((v, k) => close(v, expected[k], 1e-15));
    const plus = vortexStreamfunctionBasis({ x: point.x, y: 1e-12 }, panel), minus = vortexStreamfunctionBasis({ x: point.x, y: -1e-12 }, panel);
    plus.forEach((v, k) => { close(v, expected[k], 1e-12); close(v, minus[k], 1e-15); });
  }
  for (const angle of [.1, 1, 2.5, 4.1]) {
    const basis = r => vortexStreamfunctionBasis({ x: .5 + r * Math.cos(angle), y: r * Math.sin(angle) }, panel);
    const inside = basis(1 - 1e-12), outside = basis(1 + 1e-12);
    inside.forEach((v, k) => close(v, outside[k], 1e-12));
  }
  assert.throws(() => vortexStreamfunctionBasis({ x: NaN, y: 0 }, panel), /Invalid/);
});

test('both streamfunction shape-function gradients reproduce independently implemented vortex velocities', () => {
  const panel = makePanel({ x: .2, y: -.1 }, { x: 1.1, y: .4 });
  for (const point of [{ x: .4, y: .1 }, { x: -.2, y: -.3 }, { x: 2, y: 4 }]) {
    const h = 1e-6, px = vortexStreamfunctionBasis({ x: point.x + h, y: point.y }, panel), mx = vortexStreamfunctionBasis({ x: point.x - h, y: point.y }, panel);
    const py = vortexStreamfunctionBasis({ x: point.x, y: point.y + h }, panel), my = vortexStreamfunctionBasis({ x: point.x, y: point.y - h }, panel);
    vortexBasis(point, panel).forEach((q, k) => { close(q.u, (py[k] - my[k]) / (2 * h), 2e-8); close(q.v, -(px[k] - mx[k]) / (2 * h), 2e-8); });
  }
});

test('multielement streamfunction differences equal integrated cross-line mass and preserve rigid-coordinate changes', () => {
  const elements = [{ points: naca4('0012', 80) }, { points: transform(naca4('0012', 60), { chord: .3, x: -.4, y: .2 }) }];
  const field = solveInviscid({ elements, alpha: 2 }).field;
  const a = { x: -1.2, y: -1.5 }, b = { x: -1.2, y: 1.5 }, n = 200;
  let flux = 0;
  for (let i = 0; i <= n; i++) {
    const q = velocityAt({ x: a.x + (b.x - a.x) * i / n, y: a.y + (b.y - a.y) * i / n }, field);
    flux += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * (q.u * (b.y - a.y) - q.v * (b.x - a.x)) / (3 * n);
  }
  const difference = streamfunctionAt(b, field) - streamfunctionAt(a, field);
  close(difference, flux, 2e-9);
  const movement = { x: 3, y: -2, angle: 27 }, moved = solveInviscid({ elements: elements.map(e => ({ points: transform(e.points, movement) })), alpha: 29 }).field;
  const [ma, mb] = transform([a, b], movement);
  close(streamfunctionAt(mb, moved) - streamfunctionAt(ma, moved), difference, 2e-12);
});

test('potential differences match independent velocity integrals across angle cuts and retain physical circulation', () => {
  const panel = { ...makePanel({ x: 0, y: 0 }, { x: 1, y: 0 }), node: 0 };
  const field = { panels: [panel], gamma: [2, -.5], u: .8, v: .3 };
  for (const [a, b] of [[{ x: -.2, y: .1 }, { x: -.2, y: -.1 }], [{ x: .2, y: .15 }, { x: 1.2, y: .4 }],
    [{ x: 2, y: -1 }, { x: 2, y: 1 }]]) {
    let integral = 0; const n = 400;
    for (let i = 0; i <= n; i++) {
      const q = velocityAt({ x: a.x + i * (b.x - a.x) / n, y: a.y + i * (b.y - a.y) / n }, field);
      integral += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * (q.u * (b.x - a.x) + q.v * (b.y - a.y)) / (3 * n);
    }
    close(potentialDifference(a, b, field), integral, 2e-10);
    close(potentialDifference(b, a, field), -integral, 2e-10);
  }
  const loop = [{ x: -1, y: -1 }, { x: 2, y: -1 }, { x: 2, y: 1 }, { x: -1, y: 1 }];
  close(loop.reduce((sum, a, i) => sum + potentialDifference(a, loop[(i + 1) % 4], field), 0), .75, 2e-14);
  assert.throws(() => potentialDifference({ x: .3, y: -.1 }, { x: .3, y: .1 }, field), /crosses a vortex sheet/);
});

test('both potential gradients agree with vortex velocity through near/far primitive branches', () => {
  const panel = makePanel({ x: 0, y: 0 }, { x: 1, y: 0 });
  for (const point of [{ x: .3, y: .2 }, { x: .5, y: 1 }, { x: -.2, y: -.1 }, { x: 10, y: -5 }]) {
    const h = 1e-6, px = vortexPotentialBasis({ x: point.x + h, y: point.y }, panel), mx = vortexPotentialBasis({ x: point.x - h, y: point.y }, panel);
    const py = vortexPotentialBasis({ x: point.x, y: point.y + h }, panel), my = vortexPotentialBasis({ x: point.x, y: point.y - h }, panel);
    vortexBasis(point, panel).forEach((q, k) => { close(q.u, (px[k] - mx[k]) / (2 * h), 2e-8); close(q.v, (py[k] - my[k]) / (2 * h), 2e-8); });
  }
});

test('oblique exact panel endpoints are admissible potential-path endpoints without allowing interior sheet crossings', () => {
  const contours = transform(naca4('0012', 120), { chord: .3, x: 1.05, y: -.15, angle: -5 });
  const panel = { ...makePanel(contours.at(-2), contours[0]), node: 0 };
  const field = { panels: [panel], gamma: [.7, -.2], u: 1, v: .1 };
  for (const a of [panel.a, panel.b]) for (const sign of [-1, 1]) {
    const b = { x: a.x + sign * .004 * panel.nx, y: a.y + sign * .004 * panel.ny };
    // Integrate velocity along t=u^2, so the logarithmic endpoint velocity
    // has the integrable limit 2u*log(u) -> 0. No singular velocity is sampled.
    let integral = 0; const n = 4000;
    for (let i = 1; i <= n; i++) {
      const u = i / n, t = u * u, q = velocityAt({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) }, field);
      integral += (i === n ? 1 : i % 2 ? 4 : 2) * 2 * u * (q.u * (b.x - a.x) + q.v * (b.y - a.y)) / (3 * n);
    }
    close(potentialDifference(a, b, field), integral, 3e-10);
    close(potentialDifference(b, a, field), -integral, 3e-10);
  }
  const p = { x: panel.x, y: panel.y };
  assert.throws(() => potentialDifference({ x: p.x - .01 * panel.nx, y: p.y - .01 * panel.ny },
    { x: p.x + .01 * panel.nx, y: p.y + .01 * panel.ny }, field), /crosses a vortex sheet/);
});
