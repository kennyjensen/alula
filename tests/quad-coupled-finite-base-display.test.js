// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';
import { quadCoupledCoefficients } from '../src/ui/quad-coupled-coefficients.js';
import { quadCoupledIterationCoefficients } from '../src/ui/quad-coupled-iteration-coefficients.js';
import { quadCoupledWakeGap } from '../src/ui/quad-coupled-wake-gap.js';
import { createSurfaceContourCurve, createContourTopology } from '../src/geometry/contour-topology.js';
import { isentropicState } from '../src/potential/isentropic.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { nlrFrozenQuadDisplay } from './fixtures/nlr-frozen-quad-display.js';

const near = (a, b, tol = 3e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const fixture = nlrFrozenQuadDisplay();
const staticBL = raw => ({ ...raw.boundaryLayer, scale: 1 / Math.sqrt(raw.kernelReynolds),
  stations: raw.boundaryLayer.stations.map(({ id, kind, body, side, i, k, regime }) => ({ id, kind, body, side, i, k, regime })) });
const liveInput = raw => ({ checkpoint: raw.checkpoint, flow: raw.flow, bl: staticBL(raw), bodies: raw.solverInput.bodies,
  solverLength: raw.solverLength, referenceChord: raw.referenceChord, momentReference: { x: .25, y: 0 }, alpha: raw.alpha, mach: raw.mach });

test('saved finite NLR automatic seed displays all original base segments, transition branches and dimensional gaps without promotion', () => {
  const { raw, input } = fixture, before = JSON.stringify(raw), view = quadCoupledResultForDisplay(raw, input);
  assert.equal(view.converged, false); assert.equal(view.physicalAcceptance, false);
  assert.equal(view.coefficientStatus, 'unconverged'); assert.equal(view.diagnostics.iterations, 0);
  assert.equal(view.boundaryLayer.surfaces.length, 4); assert.equal(view.boundaryLayer.wakes.length, 2);
  assert.deepEqual(view.elements.map(e => e.points), input.elements.map(e => e.points));
  assert.deepEqual(view.coefficients.elements.map(e => e.basePressureModel.panelCount), [32, 32]);
  for (const wake of view.boundaryLayer.wakes) for (const p of wake.stations) {
    const source = raw.boundaryLayer.stations[p.index];
    assert.equal(p.wakeGap, source.wakeGap * raw.solverLength);
    near(p.fluidDeltaStar + p.wakeGap, p.deltaStar, 2e-17);
    assert.equal(p.fluidH, (source.deltaStar - source.wakeGap) / source.theta);
    assert.equal(p.h, source.deltaStar / source.theta);
  }
  const terminal = view.boundaryLayer.surfaces.find(s => s.element === 1 && s.side === 'lower');
  assert.equal(terminal.transitionKind, 'laminar-to-te');
  const body = raw.solverInput.bodies.find(b => b.element === 1);
  assert.deepEqual(terminal.transitionPoint, body.points[body.trailingEdge.lowerIndex]);
  assert.ok([view.cl, view.cd, view.cm].every(Number.isFinite));
  assert.equal(JSON.stringify(raw), before);
});

test('fixed-trip display uses only the wetted finite surface and does not include the base in its material coordinate', () => {
  // A display-only branch control, not a fixed-trip flow solution.
  const raw = structuredClone(fixture.raw);
  delete raw.boundaryLayer.transitions;
  const expected = [];
  for (const surface of raw.boundaryLayer.surfaces) {
    const body = raw.solverInput.bodies[surface.body], curve = createSurfaceContourCurve(body.points, body);
    const trip = curve.branch(surface.side, .3, raw.flow.stagnation[surface.body]);
    surface.tripParameter = trip.parameter; expected.push({ element: body.element, side: surface.side, point: trip.point });
  }
  const view = quadCoupledResultForDisplay(raw, fixture.input);
  for (const surface of view.boundaryLayer.surfaces) {
    assert.equal(surface.transitionKind, 'forced');
    assert.deepEqual(surface.transitionPoint, expected.find(e => e.element === surface.element && e.side === surface.side).point);
  }
});

test('finite-base display and live loads preserve physical length scaling and agree at the same saved iterate', () => {
  const { raw, input } = fixture, a = quadCoupledResultForDisplay(raw, input);
  assert.deepEqual(quadCoupledIterationCoefficients(liveInput(raw)), a.coefficients);
  const k = 3, b = structuredClone(raw), bi = structuredClone(input), scale = p => ({ x: k * p.x, y: k * p.y });
  b.solverLength *= k; b.referenceChord *= k; bi.referenceChord *= k;
  bi.elements.forEach(e => { e.points = e.points.map(scale); });
  b.solverInput.bodies.forEach(body => { body.points = body.points.map(scale); });
  b.flow.stagnation = b.flow.stagnation.map(s => k * s);
  b.boundaryLayer.surfaces.forEach(s => { s.tripParameter *= k; });
  for (const key of ['nodes', 'undisplacedNodes']) b.flow[key] = b.flow[key].map(g => g.map(row => row.map(scale)));
  b.mesh.vertices = b.mesh.vertices.map(scale);
  const c = quadCoupledResultForDisplay(b, bi);
  for (const key of ['cl', 'cm', 'cd']) near(c[key], a[key]);
  for (const wake of a.boundaryLayer.wakes) for (const p of wake.stations) {
    const q = c.boundaryLayer.stations[p.index];
    for (const key of ['s', 'theta', 'deltaStar', 'wakeGap', 'fluidDeltaStar']) near(q[key], k * p[key]);
    for (const key of ['h', 'fluidH', 'ue']) assert.equal(q[key], p[key]);
  }
});

test('live gap reconstruction matches every saved decoded NLR gap and the exact current-bank arc convention', () => {
  const source = fixture.raw;
  for (const wake of source.boundaryLayer.wakes) for (const id of wake.ids) {
    const station = source.boundaryLayer.stations[id], { wakeGap, ...metadata } = station;
    const args = { body: source.solverInput.bodies[wake.body], bodyIndex: wake.body,
      flow: source.flow, station: metadata, solverLength: source.solverLength };
    near(quadCoupledWakeGap(args), wakeGap * source.solverLength, 3e-17);
    // Decoded current states take priority even when live geometry is absent.
    assert.equal(quadCoupledWakeGap({ ...args, flow: null, station }), wakeGap * source.solverLength);
  }
  const raw = structuredClone(fixture.raw), wake = raw.boundaryLayer.wakes[0], body = raw.solverInput.bodies[0];
  const topology = createContourTopology(body.points, body), center = topology.trailingEdge.center;
  for (let i = body.trailingIndex + 1; i < raw.flow.nodes[0].length; i++) {
    const x = center.x + .0001 * (i - body.trailingIndex) / (raw.flow.nodes[0].length - 1 - body.trailingIndex);
    raw.flow.nodes[0][i][raw.flow.nodes[0][i].length - 1] = { x, y: center.y - .0002 };
    raw.flow.nodes[1][i][0] = { x, y: center.y + .0002 };
  }
  const result = quadCoupledIterationCoefficients(liveInput(raw));
  assert.equal(result.cd, null); assert.match(result.warnings.join(' '), /dead-air gap.*extend the wake/);
  assert.ok([result.cl, result.cm].every(Number.isFinite));
});

test('default and historical final drag withhold a premature nonzero-gap exit instead of subtracting its gap', () => {
  const raw = structuredClone(fixture.raw), id = raw.boundaryLayer.wakes[0].ids.at(-1);
  raw.boundaryLayer.stations[id].wakeGap = .0001;
  const before = raw.boundaryLayer.stations[id].deltaStar;
  const ordinary = quadCoupledResultForDisplay(raw, fixture.input);
  assert.equal(ordinary.cd, null); assert.match(ordinary.warnings.join(' '), /dead-air gap/);
  raw.conditions.blThermodynamics = raw.checkpoint.restart.options.blThermodynamics = 'historical-common-isentrope';
  raw.solverInput.streamwiseMode = 'hybrid';
  raw.solverInput.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
  const historical = quadCoupledResultForDisplay(raw, fixture.input);
  assert.equal(historical.cd, null); assert.match(historical.coefficients.errors.viscousExit, /dead-air gap/);
  assert.ok([historical.cl, historical.cm, historical.coefficients.eulerWaveDragCoefficient].every(Number.isFinite));
  assert.equal(raw.boundaryLayer.stations[id].deltaStar, before);
});

test('manufactured retained-base pressure closes constant loads and matches independent edge Gauss integration with unequal TE pressures', () => {
  const wetted = naca4('0012', 20).map((p, i) => ({ x: p.x, y: p.y + (i < 10 ? 1 : -1) * .01 * p.x }));
  const body = { element: 0, points: [...wetted, { x: 1.02, y: 0 }, { ...wetted[0] }],
    trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 20 } };
  const base = createContourTopology(body.points, body).base;
  const mach = .2, p0 = isentropicState(0, 0, { mach }).cp;
  for (const varying of [false, true]) {
    const cp = p => p0 + (varying ? .3 * p.x + .4 * p.y : 0);
    const node = p => ({ ...p, cp: cp(p) });
    const surfaces = [{ element: 0, side: 'upper', stations: wetted.slice(0, 10).reverse().map(node) },
      { element: 0, side: 'lower', stations: wetted.slice(11).map(node) }];
    const result = quadCoupledCoefficients({ surfaces, wakes: [{ element: 0, stations: [{ theta: .001, deltaStar: .002, ue: 1, wakeGap: 0 }] }],
      stagnation: [wetted[10]], bodies: [body], alpha: 0, mach, referenceChord: 1, momentReference: { x: .25, y: 0 } });
    const points = wetted.map(node), lo = points.at(-1).cp, up = points[0].cp;
    let arc = 0;
    for (let i = 1; i < base.points.length; i++) {
      arc += base.panels[i - 1].length;
      const f = i === base.points.length - 1 ? 1 : arc / base.length;
      points.push({ ...base.points[i], cp: (1 - f) * lo + f * up });
    }
    let cx = 0, cy = 0, cm = 0;
    for (let i = 1; i < points.length; i++) for (const f of [.5 - .5 / Math.sqrt(3), .5 + .5 / Math.sqrt(3)]) {
      const a = points[i - 1], b = points[i], p = a.cp * (1 - f) + b.cp * f;
      const fx = -.5 * p * (b.y - a.y), fy = .5 * p * (b.x - a.x);
      cx += fx; cy += fy;
      cm -= ((1 - f) * a.x + f * b.x - .25) * fy - ((1 - f) * a.y + f * b.y) * fx;
    }
    near(result.cx, cx, 3e-16); near(result.cy, cy, 3e-16); near(result.cm, cm, 3e-16);
    if (!varying) for (const value of [cx, cy, cm]) near(value, 0, 3e-16);
    assert.equal(result.elements[0].basePressureModel.panelCount, 2);
  }
});

test('complete sharp display remains byte-identical to the immediately archived adapter', async () => {
  const file = new URL('../docs/nlr-finite-base/quad-display/before/quad-coupled-result.js.txt', import.meta.url);
  const base = new URL('../src/ui/quad-coupled-result.js', import.meta.url);
  const code = fs.readFileSync(file, 'utf8').replace(/from '(\.[^']+)'/g, (_, p) => `from '${new URL(p, base).href}'`);
  const old = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-coupled-startup-default-browser.json', import.meta.url)));
  assert.equal(JSON.stringify(quadCoupledResultForDisplay(saved.result, saved.input)),
    JSON.stringify(old.quadCoupledResultForDisplay(saved.result, saved.input)));
});
