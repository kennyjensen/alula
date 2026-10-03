import test from 'node:test';
import assert from 'node:assert/strict';
import { isentropicSectionDensity, initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('isentropic section initialization recovers prescribed subsonic gas states and rejects choking', () => {
  // Construct expected states from local Mach and stagnation relations,
  // independently of the implementation's speed/mass inversion.
  for (const gamma of [1.2, 1.4, 1.67]) for (const mach of [.05, .2, .4, .8]) {
    for (const localMach of [.001, .1, .6, .99]) {
      const temperature = (1 + .5 * (gamma - 1) * mach ** 2) / (1 + .5 * (gamma - 1) * localMach ** 2);
      const rho = temperature ** (1 / (gamma - 1)), q = localMach * Math.sqrt(temperature) / mach;
      const recovered = isentropicSectionDensity({ massFlux: rho * q, mach, gamma });
      assert.ok(Math.abs(recovered / rho - 1) < 2e-13);
    }
    const sonicTemperature = (1 + .5 * (gamma - 1) * mach ** 2) / (1 + .5 * (gamma - 1));
    const sonicFlux = sonicTemperature ** (1 / (gamma - 1) + .5) / mach;
    assert.throws(() => isentropicSectionDensity({ massFlux: sonicFlux * (1 + 1e-10), mach, gamma }), error => {
      assert.equal(error.code, 'streamtube-sonic-capacity');
      assert.ok(Math.abs(error.diagnostics.criticalMassFlux / sonicFlux - 1) < 2e-15);
      assert.ok(Math.abs(error.diagnostics.capacityRatio - (1 + 1e-10)) < 2e-15);
      return true;
    });
  }
  for (const massFlux of [-1, 0, NaN, Infinity]) assert.throws(() => isentropicSectionDensity({ massFlux, mach: .2 }));
  for (const mach of [0, 1, NaN]) assert.throws(() => isentropicSectionDensity({ massFlux: 1, mach }));
});

test('sonic startup reports the most restrictive section without modifying the supplied state', () => {
  // Straight centerline with width decreasing 2 -> 1 -> 0.5. The second
  // section's normal area is 0.75 and has the larger capacity excess.
  const nodes = [[[{x:0,y:-1},{x:0,y:1}], [{x:1,y:-.5},{x:1,y:.5}], [{x:2,y:-.25},{x:2,y:.25}]]];
  const state = new Float64Array([.01,.02]), system = {
    conditions: { gamma: 1.4, mach: .4 }, layout: { nx: 2, tubes: [1], densityIndex: i => i },
    decode: () => ({ nodes, allocation: { groups: [[{ massFlow: 3 }]] } }),
  };
  assert.throws(() => initializeStreamtubeDensities(system,state), error => {
    assert.equal(error.code,'streamtube-sonic-capacity');
    assert.deepEqual(error.diagnostics.section,{i:1,group:0,tube:0});
    assert.equal(error.diagnostics.stage,'gas-initialization');
    assert.equal(error.diagnostics.massFlux,4);
    assert.ok(Math.abs(error.diagnostics.criticalMassFlux-1.59014)<1e-14);
    return true;
  });
  assert.deepEqual(state,new Float64Array([.01,.02]));
});

test('section initialization gives uniform stagnation entropy throughout a nonuniform body grid', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, alpha: 2, tubes: 5, tubeGrowth: 3 }));
  const state = initializeStreamtubeDensities(system, system.initial), r = system.evaluate(state);
  assert.ok(system.initial.every(x => x === 0), 'input state was mutated');
  for (const row of r.sections) for (const group of row) for (const s of group) {
    assert.ok(Math.abs(s.p / s.rho ** system.conditions.gamma / system.conditions.pInf - 1) < 2e-14);
    assert.ok(s.machSquared < 1);
  }
  assert.ok(r.diagnostics.residualByFamily.inletDensity < 3e-15);
});

test('geometry preflight rejects a downstream concave cell before an upstream choking condition', () => {
  // Minimal two-section reproducer: mean areas are positive, but the final
  // strip has a reversed corner. No panel tracing or Newton solve is needed.
  const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }],
    [{ x: 2, y: 0 }, { x: .8, y: .2 }]]];
  const state = new Float64Array(2), system = {
    conditions: { gamma: 1.4, mach: .2 }, layout: { nx: 2, tubes: [1], densityIndex: i => i },
    decode: () => ({ nodes, allocation: { groups: [[{ massFlow: 10 }]] } }),
  };
  assert.throws(() => initializeStreamtubeDensities(system, state), /Initial grid i=1, group=0, tube=0: Folded/);
  assert.deepEqual(state, new Float64Array(2));
  nodes[0][2][1] = { x: 2, y: 1 };
  assert.throws(() => initializeStreamtubeDensities(system, state), /Initial section i=0.*sonic/);
});
