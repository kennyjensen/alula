import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { prepareContour } from '../src/geometry/airfoil.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

function fixture() {
  // Presentation uses the current captured public result, including both
  // natural transitions and laminar-to-TE surfaces. Exact numerical restart
  // has its own strict tests; this check performs no flow or Newton solve.
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
  assert.equal(saved.result.converged, true);
  const raw = { ...saved.result, boundaryLayer: saved.result.numericalBoundaryLayer, physicalAcceptance: false };
  assert.equal(raw.boundaryLayer.transitions.filter(t => t.kind === 'natural').length, 2);
  assert.equal(raw.boundaryLayer.transitions.filter(t => t.kind === 'trailing-edge').length, 2);
  return { input: saved.input, raw };
}

test('automatic display maps solved arc roots onto the correct solid element and keeps terminal skin friction laminar', () => {
  const { raw, input } = fixture(), before = JSON.stringify(raw), view = quadCoupledResultForDisplay(raw, input);
  assert.equal(JSON.stringify(raw), before); assert.equal(view.physicalAcceptance, false);
  assert.deepEqual(view.boundaryLayer.surfaces.map(s => [s.element, s.side]), [[0, 'upper'], [0, 'lower'], [1, 'upper'], [1, 'lower']]);
  for (const s of view.boundaryLayer.surfaces) {
    const body = raw.solverInput.bodies.findIndex(b => b.element === s.element);
    const transition = raw.boundaryLayer.transitions.find(t => t.body === body && t.side === s.side);
    assert.equal(s.transitionDistance, transition.s * raw.solverLength); assert.equal(s.forced, false);
    const curve = createContourCurve(prepareContour(raw.solverInput.bodies[body].points));
    const stag = raw.flow.stagnation[body], end = s.side === 'upper' ? 0 : curve.length;
    if (s.transitionKind === 'laminar-to-te') {
      const te = curve.evaluate(end).point;
      assert.ok(Math.hypot(s.transitionPoint.x - te.x, s.transitionPoint.y - te.y) < 1e-12);
      const node = s.stations.at(-1), native = raw.boundaryLayer.stations[node.id];
      assert.equal(node.regime, 'laminar'); assert.ok(node.amplification < raw.conditions.ncrit);
      const kernel = createIntegralKernel({ reynolds: raw.kernelReynolds, mach: raw.mach, ncrit: raw.conditions.ncrit });
      const properties = kernel.station(native, 'laminar');
      assert.equal(node.cf, properties.cf * properties.rho * native.ue ** 2);
    } else {
      // Independent fine polyline arclength, rather than the production arc
      // quadrature/inversion, checks orientation, dimensional scale and body.
      let previous = curve.evaluate(stag).point, distance = 0, point;
      for (let i = 1; i <= 8192; i++) {
        const next = curve.evaluate(stag + i / 8192 * (end - stag)).point, ds = Math.hypot(next.x - previous.x, next.y - previous.y);
        if (distance + ds >= s.transitionDistance) {
          const f = (s.transitionDistance - distance) / ds;
          point = { x: previous.x + f * (next.x - previous.x), y: previous.y + f * (next.y - previous.y) }; break;
        }
        distance += ds; previous = next;
      }
      assert.ok(point); assert.ok(Math.hypot(s.transitionPoint.x - point.x, s.transitionPoint.y - point.y) < 2e-6);
      assert.equal(s.transitionKind, 'natural');
    }
  }
});

test('automatic transition markers scale with physical geometry while surface fractions stay dimensionless', () => {
  const { raw, input } = fixture(), a = quadCoupledResultForDisplay(raw, input), scale = p => ({ x: 2 * p.x, y: 2 * p.y });
  const scaled = structuredClone(raw), scaledInput = structuredClone(input);
  scaled.solverLength *= 2; scaled.referenceChord *= 2; scaledInput.referenceChord *= 2;
  scaledInput.elements.forEach(e => { e.points = e.points.map(scale); });
  scaled.solverInput.bodies.forEach(b => { b.points = b.points.map(scale); });
  scaled.flow.stagnation = scaled.flow.stagnation.map(v => 2 * v);
  for (const key of ['nodes', 'undisplacedNodes']) scaled.flow[key] = scaled.flow[key].map(g => g.map(row => row.map(scale)));
  scaled.mesh.vertices = scaled.mesh.vertices.map(scale);
  const b = quadCoupledResultForDisplay(scaled, scaledInput);
  a.boundaryLayer.surfaces.forEach((s, i) => {
    const t = b.boundaryLayer.surfaces[i];
    assert.ok(Math.hypot(t.transitionPoint.x - 2 * s.transitionPoint.x, t.transitionPoint.y - 2 * s.transitionPoint.y) < 1e-12);
    assert.ok(Math.abs(t.transitionFraction - s.transitionFraction) < 1e-12);
  });
});
