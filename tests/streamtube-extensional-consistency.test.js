import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-cell.js';
import { streamtubeEdgeVelocity } from '../src/euler/streamtube-edge-velocity.js';

// Exact potential phi=x+a*(x*x-y*y)/2 gives (u,v)=(1+a*x,-a*y),
// psi=y*(1+a*x). The wall is psi=0, the other streamline is psi=mass.
// These are truncation checks on imposed analytic streamlines, not Euler roots.
const a = 4, center = .1, wallSpeed = 1 + a * center;
const cellAtSpacing = (left, right, mass) => {
  const lower = [center - left, center, center + right].map(x => ({ x, y: 0 }));
  const upper = lower.map(p => ({ x: p.x, y: mass / (1 + a * p.x) }));
  return evaluateIncompressibleStreamtubeCell({ lower, upper, massFlow: mass });
};
const cellAt = (h, mass) => cellAtSpacing(h, h, mass);

test('accelerating potential-flow wall velocity and pressure approach the analytic values under paired refinement', t => {
  const rows = [.1, .05, .025, .0125].map(h => {
    const c = cellAt(h, .2 * h * h);
    return { h, velocityError: Math.abs(streamtubeEdgeVelocity(c).ue - wallSpeed),
      pressureError: Math.abs(c.interfacePressure.lower + .5 * wallSpeed * wallSpeed),
      pressureCorrection: c.pressureCorrection, curvature: c.geometry.pressureCurvature };
  });
  for (let i = 1; i < rows.length; i++) {
    assert.ok(rows[i].velocityError / rows[i - 1].velocityError < .26);
    assert.ok(rows[i].pressureError / rows[i - 1].pressureError < .31);
    assert.ok(Math.abs(rows[i].curvature / rows[i - 1].curvature) < .26);
  }
  assert.ok(rows.at(-1).velocityError < 5e-4 && rows.at(-1).pressureError < 7e-5);
  t.diagnostic(JSON.stringify(rows));
});

test('normal refinement alone retains the independently derived finite-streamwise wall-speed error', t => {
  const h = .1, d = a * h, q = wallSpeed;
  // As mass -> 0, each section speed is the harmonic mean of the two
  // analytic endpoint speeds. Averaging the adjacent sections gives this
  // closed form. The pressure geometry factor is independent of mass here.
  const limit = q - q * d * d / (4 * q * q - d * d);
  const curvature = -2 * d * d / (2 * q * q - d * d);
  const rows = [.002, .0005, .000125, .00003125].map(mass => {
    const c = cellAt(h, mass), velocity = streamtubeEdgeVelocity(c).ue;
    assert.ok(Math.abs(c.geometry.pressureCurvature - curvature) < 1e-13);
    return { mass, velocity, distanceFromLimit: Math.abs(velocity - limit), error: velocity - q };
  });
  for (let i = 1; i < rows.length; i++) assert.ok(rows[i].distanceFromLimit < .064 * rows[i - 1].distanceFromLimit);
  assert.ok(rows.at(-1).distanceFromLimit < 1e-9);
  assert.ok(Math.abs(rows.at(-1).error) > .029);
  t.diagnostic(JSON.stringify({ limit, curvature, rows }));
});

test('arithmetic wall-speed reconstruction loses an order at a fixed spacing jump', t => {
  const rows = [.04, .02, .01, .005].map(h => {
    // Smooth grading: right-left = O(h^2). Abrupt grading: right/left=1/3.
    const error = right => Math.abs(streamtubeEdgeVelocity(cellAtSpacing(h, right, .2 * h * h)).ue - wallSpeed);
    const right = h / 3, q = wallSpeed;
    // Independently derived thin-tube section speeds are harmonic means
    // of the exact endpoint velocities, for unequal intervals as well.
    const limit = .5 * (2 * q * (q - a * h) / (2 * q - a * h)
      + 2 * q * (q + a * right) / (2 * q + a * right));
    const thin = streamtubeEdgeVelocity(cellAtSpacing(h, right, 1e-7)).ue;
    assert.ok(Math.abs(thin - limit) < 1e-12);
    return { h, smoothError: error(h * (1 + 2 * h)), jumpError: error(right) };
  });
  for (let i = 1; i < rows.length; i++) {
    const smoothRatio = rows[i].smoothError / rows[i - 1].smoothError;
    const jumpRatio = rows[i].jumpError / rows[i - 1].jumpError;
    assert.ok(smoothRatio > .2 && smoothRatio < .3);
    assert.ok(jumpRatio > .43 && jumpRatio < .53);
  }
  t.diagnostic(JSON.stringify(rows));
});
