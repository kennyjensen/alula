// SPDX-License-Identifier: GPL-2.0-or-later
// Independent stagnation-flow limits of the documented midpoint geometry.
// These are imposed-field consistency checks, not converged Euler/BL roots.
import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-cell.js';
import { streamtubeEdgeVelocity } from '../src/euler/streamtube-edge-velocity.js';

// Exact linear stagnation potential: phi=a(s²-n²)/2, psi=a*s*n.
// On the wall, phi=a*sWall²/2. The mass coordinate is rho*psi.
function stagnationCell(k, h, beta, a = 2, rho = 1) {
  const massFlow = beta * rho * h * h, b = massFlow / (rho * a);
  const lower = [k - 1, k, k + 1].map(i => ({ x: i * h, y: 0 }));
  const upper = lower.map(p => {
    const s = Math.sqrt(.5 * (p.x * p.x + Math.hypot(p.x * p.x, 2 * b)));
    return { x: s, y: b / s };
  });
  return evaluateIncompressibleStreamtubeCell({ lower, upper, massFlow, density: rho });
}

test('thin stagnation tubes approach the independently derived harmonic-position limit', () => {
  // For s=k*h, qmean/(a*s)=1-1/(4*k²-1). At k=1 the left
  // section tends to zero separately through the singular endpoint.
  for (const k of [1, 2, 3, 5]) for (const a of [2, 5]) {
    const h = .02, edge = streamtubeEdgeVelocity(stagnationCell(k, h, 1e-8, a));
    const ratio = edge.ue / (a * k * h), limit = 1 - 1 / (4 * k * k - 1);
    assert.ok(Math.abs(ratio - limit) < 1e-4, JSON.stringify({ k, a, ratio, limit }));
    assert.equal(edge.correction, 0);
  }
});

test('shrinking a self-similar first stagnation cell does not imply relative strain convergence', () => {
  for (const beta of [.5, .25, .125]) {
    const ratios = [.04, .02, .01, .005].map(h => streamtubeEdgeVelocity(stagnationCell(1, h, beta)).ue / (2 * h));
    assert.ok(Math.max(...ratios) - Math.min(...ratios) < 2e-14);
    assert.ok(ratios[0] < .99 && ratios[0] > 2 / 3);
  }
});

function inverseCylinder(phi, psi) {
  // Recover the smaller square-root component from 2*real*imag=psi.
  // Subtracting nearly equal hypot(x,psi) and |x| loses the tiny mass.
  const root = x => {
    const large = Math.sqrt(.5 * (Math.hypot(x, psi) + Math.abs(x))), small = psi / (2 * large);
    return x < 0 ? { x: small, y: large } : { x: large, y: small };
  };
  const a = root(phi - 2), b = root(phi + 2);
  return { x: .5 * (phi + a.x * b.x - a.y * b.y), y: .5 * (psi + a.x * b.y + a.y * b.x) };
}

test('cylinder wall speed converges quadratically at a fixed physical station', () => {
  const s = .04, exact = 2 * Math.sin(s);
  const errors = [.01, .005, .0025, .00125].map(h => {
    const massFlow = .5 * h * h, angles = [s - h, s, s + h];
    const lower = angles.map(t => ({ x: -Math.cos(t), y: Math.sin(t) }));
    const upper = angles.map(t => {
      const phi = -2 * Math.cos(t), p = inverseCylinder(phi, massFlow), r2 = p.x * p.x + p.y * p.y;
      assert.ok(r2 >= 1);
      assert.ok(Math.abs(p.x * (1 + 1 / r2) - phi) < 2e-13);
      assert.ok(Math.abs(p.y * (1 - 1 / r2) - massFlow) < 2e-13);
      return p;
    });
    const cell = evaluateIncompressibleStreamtubeCell({ lower, upper, massFlow });
    return Math.abs(streamtubeEdgeVelocity(cell).ue / exact - 1);
  });
  for (let i = 1; i < errors.length; i++) assert.ok(errors[i - 1] / errors[i] > 3.8 && errors[i - 1] / errors[i] < 4.2);
  assert.ok(errors.at(-1) < .0003);
});
