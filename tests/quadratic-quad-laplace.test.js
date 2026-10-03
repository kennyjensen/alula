import test from 'node:test';
import assert from 'node:assert/strict';
import { quadraticQuadLaplaceMatrix } from '../src/numerics/tests/quadratic-quad-laplace.js';

test('Q2 rectangle stiffness matches exact polynomial mass/stiffness integrals', () => {
  const mass = [[4, 2, -1], [2, 16, 2], [-1, 2, 4]].map(r => r.map(v => v / 30));
  const stiffness = [[7, -8, 1], [-8, 16, -8], [1, -8, 7]].map(r => r.map(v => v / 3));
  for (const [width, height] of [[1, 1], [2, .3], [.1, 4]]) for (const quadratureOrder of [3, 5]) {
    const k = quadraticQuadLaplaceMatrix([{ x: 0, y: 0 }, { x: width, y: 0 }, { x: width, y: height }, { x: 0, y: height }], { quadratureOrder });
    for (let a = 0; a < 9; a++) for (let b = 0; b < 9; b++) {
      const i = Math.floor(a / 3), j = a % 3, u = Math.floor(b / 3), v = b % 3;
      const exact = height / width * stiffness[i][u] * mass[j][v] + width / height * mass[i][u] * stiffness[j][v];
      assert.ok(Math.abs(k[9 * a + b] - exact) < 5e-13 * Math.max(1, Math.abs(exact)));
    }
  }
});

test('Q2 distorted-cell energy, constant nullspace and harmonic bubble satisfy physical identities', () => {
  const p = [{ x: 0, y: 0 }, { x: 1.3, y: .15 }, { x: 1.1, y: 1.4 }, { x: -.2, y: .8 }];
  const map = (s, t) => Object.fromEntries(['x', 'y'].map(key => [key, (1-s)*(1-t)*p[0][key]+s*(1-t)*p[1][key]+s*t*p[2][key]+(1-s)*t*p[3][key]]));
  const nodes = Array.from({ length: 9 }, (_, k) => map(Math.floor(k / 3) / 2, k % 3 / 2));
  const k = quadraticQuadLaplaceMatrix(p), area = .5 * p.reduce((s, a, i) => { const b = p[(i+1)%4]; return s+a.x*b.y-a.y*b.x; }, 0);
  const energy = f => k.reduce((s, v, ij) => s+v*f(nodes[Math.floor(ij/9)])*f(nodes[ij%9]), 0);
  // The mapped physical coordinates and both harmonic quadratics are in Q2.
  for (const f of [p => 1, p => p.x, p => p.y, p => p.x*p.x-p.y*p.y, p => 2*p.x*p.y])
    assert.ok(Math.abs(nodes.reduce((s, p, j) => s+k[4*9+j]*f(p), 0)) < 1e-12);
  assert.ok(Math.abs(energy(p => 2*p.x-3*p.y) - 13*area) < 2e-12);
  for (let i = 0; i < 9; i++) {
    assert.ok(Math.abs(k.slice(9*i,9*i+9).reduce((s,v)=>s+v,0)) < 2e-14);
    for (let j = 0; j < 9; j++) assert.equal(k[9*i+j], k[9*j+i]);
  }
  const moved = p.map(p => ({ x: 500+3*(.8*p.x-.6*p.y), y: -300+3*(.6*p.x+.8*p.y) }));
  const transformed = quadraticQuadLaplaceMatrix(moved);
  k.forEach((v,i)=>assert.ok(Math.abs(v-transformed[i]) < 2e-12));
});

test('Q2 reference rejects degenerate, concave, reversed and nonfinite cells', () => {
  const p = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 1 }, { x: 0, y: 1 }];
  for (const bad of [p.slice(1), p.toReversed(), [p[0],p[1],p[1],p[3]], [p[0],p[1],{x:.1,y:.1},p[3]], p.map(p=>({...p,x:NaN}))])
    assert.throws(()=>quadraticQuadLaplaceMatrix(bad));
  assert.throws(()=>quadraticQuadLaplaceMatrix(p,{quadratureOrder:2}),/Gauss/);
});
