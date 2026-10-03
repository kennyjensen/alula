import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createPhysicalBLDomain } from '../src/viscous/physical-domain.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

test('compressible BL margins and derivatives recover original Fortran BLKIN/HKIN', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/physical-bl-domain.json', import.meta.url)));
  for (const { parameters, station, expected } of fixture.cases) {
    const domain = createPhysicalBLDomain(parameters)(station);
    const denominator = station.theta * domain.enthalpy * (1 + .113 * expected.machSquared);
    assert.ok(Math.abs(domain.shape / denominator - (expected.hk - 1)) < 2e-13);
    assert.equal(domain.shape > 0, expected.hk > 1);
    if (expected.hk > 1) {
      // Preserve BLKIN's kinematic shape separately from BLVAR's closure
      // bounds, for both wall and wake properties and nonzero wake gaps.
      for (const regime of ['laminar', 'wake']) {
        const p = createIntegralKernel(parameters).station({ s: .1, aux: .03, ...station }, regime);
        assert.ok(Math.abs(p.rawHk - expected.hk) < 2e-13);
        assert.equal(p.hk, Math.max(p.rawHk, regime === 'wake' ? 1.00005 : 1.05));
      }
    }
    for (const units of [1e-3, 1, 1e3]) {
      const scaled = { ...station, theta: station.theta * units, deltaStar: station.deltaStar * units, wakeGap: station.wakeGap * units };
      const check = createPhysicalBLDomain(parameters)(scaled);
      assert.ok(Math.abs(check.shape / units - domain.shape) < 1e-15);
      // Differentiate the independent native Hk expression, including its
      // positive denominator, to check all three polynomial derivatives.
      const h = station.theta, e = domain.enthalpy, m = expected.machSquared, k = expected.hk - 1;
      const denominatorGradient = [e * (1 + .113 * m), 0,
        h * (domain.enthalpyGradient[2] * (1 + .113 * m) + e * .113 * expected.machSquaredUe)];
      for (let j = 0; j < 3; j++) {
        const native = expected.hkGradient[j] * denominator + k * denominatorGradient[j];
        assert.ok(Math.abs(native - domain.shapeGradient[j]) < 2e-12 * Math.max(1, Math.abs(native)));
      }
    }
  }
});

test('physical domain derivatives match independent differences, with thermal failures remaining evaluable', () => {
  for (const mach of [0, .2, .75]) {
    const domain = createPhysicalBLDomain({ mach }), station = { theta: .02, deltaStar: .027, ue: 1.3, wakeGap: .001 };
    const base = domain(station);
    for (const [j, key] of ['theta', 'deltaStar', 'ue'].entries()) {
      const h = 1e-4 * station[key], f = [-2, -1, 1, 2].map(m => domain({ ...station, [key]: station[key] + m * h }));
      for (const name of ['shape', 'enthalpy']) {
        const difference = (f[0][name] - 8 * f[1][name] + 8 * f[2][name] - f[3][name]) / (12 * h);
        assert.ok(Math.abs(difference - base[`${name}Gradient`][j]) < 2e-9);
      }
    }
    if (mach) {
      const outside = domain({ ...station, ue: 20 });
      assert.ok(outside.enthalpy < 0); assert.ok(Number.isFinite(outside.shape));
    } else assert.equal(base.shape, station.deltaStar - station.wakeGap - station.theta);
  }
});

test('a state with deltaStar greater than theta can still violate the compressible kernel domain', () => {
  const station = { s: .5, theta: .001, deltaStar: .00101, ue: 1, aux: .03 };
  assert.ok(station.deltaStar > station.theta);
  assert.ok(createPhysicalBLDomain({ mach: .2 })(station).shape < 0);
  assert.throws(() => createIntegralKernel({ mach: .2 }).station(station, 'wake'), /compressible BL edge state/);
  const valid = { ...station, deltaStar: .0011 };
  assert.ok(createPhysicalBLDomain({ mach: .2 })(valid).shape > 0);
  assert.ok(createIntegralKernel({ mach: .2 }).station(valid, 'wake').hk > 1);
});
