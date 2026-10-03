import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-cell.js';

// Independent exact cylinder flow: w=z+1/z, |z|>=1, U_inf=radius=1.
// Its inverse solves z^2-w*z+1=0; choose the exterior root. The incoming
// dividing streamline and the upper wall meet at z=-1, where p-p0=0.
function exterior(phi, psi) {
  const re = phi * phi - psi * psi - 4, im = 2 * phi * psi, magnitude = Math.hypot(re, im);
  const a = Math.sqrt((magnitude + re) / 2), b = Math.sign(im) * Math.sqrt((magnitude - re) / 2);
  const roots = [{ x: (phi + a) / 2, y: (psi + b) / 2 }, { x: (phi - a) / 2, y: (psi - b) / 2 }];
  return roots.reduce((a, b) => Math.hypot(a.x, a.y) > Math.hypot(b.x, b.y) ? a : b);
}

test('exact cylinder stagnation pressure has second-order local convergence when arc distances and normal distances are refined together', t => {
  const rows = [];
  for (const h of [.08, .04, .02, .01]) {
    // Near stagnation psi is quadratic in physical distance: quarter the
    // mass when halving distances. A doubled tube count alone does not do this.
    const mass = h * h, phi = [-1 - h - 1 / (1 + h), -2, -2 * Math.cos(h)];
    const wall = [{ x: -1 - h, y: 0 }, { x: -1, y: 0 }, { x: -Math.cos(h), y: Math.sin(h) }];
    const outer = phi.map(p => exterior(p, mass));
    outer.forEach((p, i) => {
      const r2 = p.x * p.x + p.y * p.y;
      assert.ok(r2 > 1); assert.ok(Math.abs(p.x * (1 + 1 / r2) - phi[i]) < 2e-13);
      assert.ok(Math.abs(p.y * (1 - 1 / r2) - mass) < 2e-13);
    });
    const upper = evaluateIncompressibleStreamtubeCell({ lower: wall, upper: outer, massFlow: mass });
    const mirror = row => row.map(p => ({ x: p.x, y: -p.y }));
    const lower = evaluateIncompressibleStreamtubeCell({ lower: mirror(outer), upper: mirror(wall), massFlow: mass });
    assert.ok(Math.abs(upper.interfacePressure.lower - lower.interfacePressure.upper) < 1e-14);
    const error = Math.abs(upper.interfacePressure.lower); // exact wall p-p0=0
    if (rows.length) {
      const ratio = error / rows.at(-1).error;
      assert.ok(ratio > .2 && ratio < .3, `pressure error ratio ${ratio}`);
    }
    rows.push({ h, mass, error });
  }
  assert.ok(rows.at(-1).error < 6e-5); t.diagnostic(JSON.stringify(rows));
});

test('stagnation-pressure refinement also converges on oblique exact-flow crosslines', t => {
  const cases = [];
  for (const skew of [-.5, .5, 1.2]) {
    const rows = [];
    for (const h of [.08, .04, .02, .01]) {
      const mass = h * h, phi = [-1 - h - 1 / (1 + h), -2, -2 * Math.cos(h)];
      const lower = [{ x: -1 - h, y: 0 }, { x: -1, y: 0 }, { x: -Math.cos(h), y: Math.sin(h) }];
      // phi = wall_phi + skew*psi tilts a crossline in the exact potential
      // plane. It changes geometry while preserving the same physical flow.
      const upper = phi.map(p => exterior(p + skew * mass, mass));
      const value = evaluateIncompressibleStreamtubeCell({ lower, upper, massFlow: mass });
      const error = Math.abs(value.interfacePressure.lower);
      assert.ok(error > 0);
      if (rows.length) assert.ok(error / rows.at(-1).error > .2 && error / rows.at(-1).error < .3);
      rows.push({ h, error });
    }
    assert.ok(rows.at(-1).error < 2e-4); cases.push({ skew, rows });
  }
  t.diagnostic(JSON.stringify(cases));
});
