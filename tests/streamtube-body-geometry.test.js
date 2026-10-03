import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { multipoleBasis } from '../src/potential/farfield.js';
import { multipoleGeometryDerivatives } from '../src/potential/farfield-geometry.js';

test('body geometry derivatives include stagnation/cut interpolation and refresh after normal rebasing', () => {
  for (const elements of [1, 2]) {
    const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements })), { layout } = system;
    let state = system.initial.map((_, i) => i < layout.densityCount ? 0 : 1e-4 * Math.sin(i));
    for (let chart = 0; chart < 2; chart++) {
      const derivatives = system.geometryDerivatives(state);
      for (let col = layout.densityCount; col < layout.n; col++) {
        const h = 1e-7, plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
        const a = system.decode(plus).nodes, b = system.decode(minus).nodes;
        for (let g = 0; g < a.length; g++) for (let i = 0; i <= layout.nx; i++) for (let j = 0; j <= layout.tubes[g]; j++) {
          const d = derivatives[g][i][j].get(col) ?? { x: 0, y: 0 };
          for (const key of ['x', 'y']) assert.ok(Math.abs(d[key] - (a[g][i][j][key] - b[g][i][j][key]) / (2 * h)) < 1e-7);
        }
      }
      state = system.rebase(state);
    }
  }
});

test('moving farfield velocity derivatives recover independent differences and irrotational symmetry', () => {
  for (const mach of [0, .2, .8]) for (const alpha of [-37, 0, 63]) for (const point of [{ x: 2, y: .3 }, { x: -1, y: -2 }]) {
    const settings = { center: { x: .25, y: -.1 }, alpha, mach, gamma: 1.4 }, r = multipoleGeometryDerivatives(point, settings);
    assert.deepEqual(r.velocity, multipoleBasis(point, settings).velocity);
    for (let axis = 0; axis < 2; axis++) {
      const key = ['x', 'y'][axis], h = 2e-6;
      const plus = multipoleBasis({ ...point, [key]: point[key] + h }, settings), minus = multipoleBasis({ ...point, [key]: point[key] - h }, settings);
      for (let k = 0; k < 5; k++) for (let component = 0; component < 2; component++)
        assert.ok(Math.abs(r.velocityDerivatives[k][component][axis] - (plus.velocity[k][component] - minus.velocity[k][component]) / (2 * h)) < 2e-9);
    }
    for (const d of r.velocityDerivatives) assert.ok(Math.abs(d[0][1] - d[1][0]) < 2e-14);
  }
});
