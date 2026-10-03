import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-cell.js';
import { linearizeIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-linearization.js';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { isentropicSectionDensity } from '../src/euler/streamtube-initial-state.js';
import { createIncompressibleStreamtubePatch, solveIncompressibleStreamtubePatch } from '../src/euler/tests/incompressible-streamtube-patch.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const cell = { lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
  upper: [{ x: .1, y: .25 }, { x: .58, y: .3 }, { x: 1.18, y: .4 }], massFlow: .27 };

test('incompressible analytic chain rule checks all coordinates, density, mass and curvature coefficient', () => {
  const parameters = { ...cell, density: 1.2, pressureCorrectionFactor: .1 }, linearization = linearizeIncompressibleStreamtubeCell(parameters);
  const tangents = [{ density: 1 }, { massFlow: 1 }, { pressureCorrectionFactor: 1 }];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const row = Array.from({ length: 3 }, () => ({ x: 0, y: 0 })); row[i][key] = 1; tangents.push({ [side]: row });
  }
  const flatten = v => [v.streamwiseResidual, v.pressureCorrection, ...Object.values(v.interfacePressure), ...v.states.flatMap(Object.values)];
  for (const tangent of tangents) {
    const h = 1e-6, shifted = sign => {
      const p = structuredClone(parameters);
      for (const side of ['lower', 'upper']) if (tangent[side]) p[side].forEach((v, i) => { for (const key of ['x', 'y']) v[key] += sign * h * tangent[side][i][key]; });
      for (const key of ['density', 'massFlow', 'pressureCorrectionFactor']) p[key] += sign * h * (tangent[key] ?? 0);
      return flatten(evaluateIncompressibleStreamtubeCell(p));
    };
    const a = shifted(1), b = shifted(-1);
    flatten(linearization.apply(tangent)).forEach((v, i) => assert.ok(Math.abs(v - (a[i] - b[i]) / (2 * h)) < 2e-8 * Math.max(1, Math.abs(v))));
  }
  assert.throws(() => linearization.apply({ massFlow: NaN }), /Invalid/);
});

test('incompressible cell satisfies Bernoulli, section mass and independent normal momentum balance', () => {
  const c = evaluateIncompressibleStreamtubeCell(cell), g = c.geometry, [a, b] = c.states;
  for (const [i, s] of c.states.entries()) {
    close(s.rho * s.q * g.normalAreas[i], cell.massFlow);
    close(s.p + .5 * s.rho * s.q ** 2, 0);
  }
  const flux = { x: cell.massFlow * (b.q * g.directions[1].x - a.q * g.directions[0].x),
    y: cell.massFlow * (b.q * g.directions[1].y - a.q * g.directions[0].y) };
  for (const [p, edge, sign] of [[a.p, g.sections[0], -1], [b.p, g.sections[1], 1],
    [c.interfacePressure.lower, g.sides.lower, 1], [c.interfacePressure.upper, g.sides.upper, -1]]) {
    flux.x += sign * p * edge.y; flux.y -= sign * p * edge.x;
  }
  close((flux.x * g.transverse.x + flux.y * g.transverse.y) / g.area, 0);
  close((flux.x * g.streamwise.x + flux.y * g.streamwise.y) / g.area, c.streamwiseResidual);
  assert.match(c.pressureReference, /relative/);
  assert.throws(() => evaluateIncompressibleStreamtubeCell({ ...cell, massFlow: 0 }), /Invalid/);
  assert.throws(() => evaluateIncompressibleStreamtubeCell({ ...cell, lower: cell.upper, upper: cell.lower }), /Folded/);
});

test('finite-Mach constant-total-pressure Euler cells approach the incompressible pressure and curvature limits quadratically', () => {
  const c = evaluateIncompressibleStreamtubeCell(cell), gamma = 1.4; let previous = Infinity;
  for (const mach of [.04, .02, .01]) {
    const densities = c.geometry.normalAreas.map(area => isentropicSectionDensity({ massFlux: cell.massFlow / area, mach, gamma }));
    const e = evaluateStreamtubeCell({ ...cell, densities, gamma, stagnationEnthalpy: 1 / ((gamma - 1) * mach ** 2) + .5 });
    const p0 = (1 + .5 * (gamma - 1) * mach ** 2) ** (gamma / (gamma - 1)) / (gamma * mach ** 2);
    const error = Math.max(...['lower', 'upper'].map(side => Math.abs(e.interfacePressure[side] - p0 - c.interfacePressure[side])),
      Math.abs(e.pressureCorrection - c.pressureCorrection));
    assert.ok(error < .26 * previous); previous = error;
  }
  assert.ok(previous < 1.6e-5);
});

test('incompressible cell retains rotation, translation, density and length scaling', () => {
  const original = evaluateIncompressibleStreamtubeCell(cell), angle = .6, scale = 3, density = 2;
  const move = row => row.map(p => ({ x: 4 + scale * (p.x * Math.cos(angle) - p.y * Math.sin(angle)),
    y: -2 + scale * (p.x * Math.sin(angle) + p.y * Math.cos(angle)) }));
  const moved = evaluateIncompressibleStreamtubeCell({ ...cell, lower: move(cell.lower), upper: move(cell.upper), density,
    massFlow: density * scale * cell.massFlow });
  for (const side of ['lower', 'upper']) close(moved.interfacePressure[side] / density, original.interfacePressure[side]);
  close(moved.pressureCorrection / density, original.pressureCorrection);
  moved.states.forEach((s, i) => close(s.q, original.states[i].q));
});

test('small incompressible relaxation patches refine streamline radius and physical pressure toward an exact circular vortex', () => {
  let previous = { radius: Infinity, pressure: Infinity };
  for (const [nx, nt] of [[4, 2], [8, 4], [16, 8]]) {
    // Exact q=1/r, psi=-log(r/2), p-p0=-1/(2*r^2).
    const exact = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
      const r = 2 * Math.exp(-.5 * j / nt), t = -.4 + .8 * i / nx;
      return { x: r * Math.cos(t), y: r * Math.sin(t) };
    }));
    // Straight interpolation of inner/outer curves displaces every interior
    // streamline. End nodes retain the prescribed exact mass allocation.
    const nodes = exact.map((row, i) => row.map((p, j) => i === 0 || i === nx ? p
      : { x: (1 - j / nt) * row[0].x + j / nt * row[nt].x, y: (1 - j / nt) * row[0].y + j / nt * row[nt].y }));
    const before = structuredClone(nodes), system = createIncompressibleStreamtubePatch({ nodes, massFlows: Array(nt).fill(.5 / nt) });
    const r = solveIncompressibleStreamtubePatch(system);
    assert.equal(r.converged, true, r.reason); assert.ok(r.history.at(-1).residual <= 1e-10);
    assert.deepEqual(nodes, before);
    const error = { radius: 0, pressure: 0 };
    r.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Normal movement can also change the angular parameter. Compare the
      // physical streamline, not the original node's angular coordinate.
      error.radius = Math.max(error.radius, Math.abs(Math.hypot(p.x, p.y) - 2 * Math.exp(-.5 * j / nt)));
      if (i === 0 || i === nx || j === 0 || j === nt) assert.deepEqual(p, before[i][j]);
    }));
    r.cells.forEach((row, i) => row.forEach((c, j) => {
      for (const [side, k] of [['lower', j], ['upper', j + 1]]) {
        const p = r.nodes[i + 1][k]; error.pressure = Math.max(error.pressure, Math.abs(c.interfacePressure[side] + .5 / (p.x * p.x + p.y * p.y)));
      }
      c.states.forEach((s, k) => close(s.rho * s.q * c.geometry.normalAreas[k], .5 / nt));
      if (j) close(row[j - 1].interfacePressure.upper, c.interfacePressure.lower, 1e-10);
    }));
    for (const key of Object.keys(error)) assert.ok(error[key] < .3 * previous[key], `${key}: ${JSON.stringify(error)} after ${JSON.stringify(previous)}`);
    previous = error;
  }
  assert.ok(previous.radius < 5.1e-5 && previous.pressure < 6.6e-4);
});
