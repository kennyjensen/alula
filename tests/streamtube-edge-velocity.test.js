import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { evaluateIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-cell.js';
import { linearizeIncompressibleStreamtubeCell } from '../src/euler/incompressible-streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from '../src/euler/streamtube-edge-velocity.js';

const close = (a, b, tol = 3e-8) => assert.ok(Math.abs(a - b) < tol * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const parameters = { lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
  upper: [{ x: .1, y: .25 }, { x: .58, y: .3 }, { x: 1.18, y: .4 }],
  densities: [1.05, 1.02], massFlow: .27, stagnationEnthalpy: 8, gamma: 1.4 };

test('section-speed correction responds with opposite signs to constricting and expanding sawtooth geometry', () => {
  for (const epsilon of [-.01, .01]) {
    // Both section areas remain 1, both centerline speeds remain 1.
    // Successive lower-edge cross product is -4e, upper is +4e;
    // S=(1,0), N=(0,1), so the signed geometry factor is -4e.
    const lower = [{ x: 0, y: -epsilon }, { x: 1, y: epsilon }, { x: 2, y: -epsilon }];
    const upper = lower.map(p => ({ x: p.x, y: 1 - p.y }));
    const cell = evaluateStreamtubeCell({ lower, upper, densities: [1, 1], massFlow: 1,
      stagnationEnthalpy: .5 + 1 / (.4 * .25) });
    const edge = streamtubeEdgeVelocity(cell);
    close(edge.meanSpeed, 1, 1e-13); close(edge.meanMachSquared, .25, 1e-13);
    close(edge.correction, .15 * epsilon, 1e-13);
    assert.equal(Math.sign(edge.ue - 1), Math.sign(epsilon));
    const reflected = { ...cell, geometry: { ...cell.geometry, pressureCurvature: -cell.geometry.pressureCurvature } };
    close(streamtubeEdgeVelocity(reflected).correction, -edge.correction, 1e-13);
    close(streamtubeEdgeVelocity(cell, 0).ue, 1, 1e-13);
  }
});

test('section-speed chain differentiates every node coordinate, both densities, mass and enthalpy on each smooth branch', () => {
  const changes = [];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const row = Array.from({ length: 3 }, () => ({ x: 0, y: 0 })); row[i][key] = 1;
    changes.push({ tangent: { [side]: row }, shift: (p, h) => { p[side][i][key] += h; } });
  }
  for (const i of [0, 1]) changes.push({ tangent: { densities: [Number(i === 0), Number(i === 1)] }, shift: (p, h) => { p.densities[i] += h; } });
  for (const key of ['massFlow', 'stagnationEnthalpy']) changes.push({ tangent: { [key]: 1 }, shift: (p, h) => { p[key] += h; } });
  for (const interpolation of ['arithmetic','distance-weighted']) for (const massFlow of [.27, .65]) {
    const p = { ...parameters, massFlow }, cell = linearizeStreamtubeCell(p);
    for (const { tangent, shift } of changes) {
      const plus = structuredClone(p), minus = structuredClone(p), h = 1e-6;
      shift(plus, h); shift(minus, -h);
      const fd = (streamtubeEdgeVelocity(evaluateStreamtubeCell(plus), .2, interpolation).ue - streamtubeEdgeVelocity(evaluateStreamtubeCell(minus), .2, interpolation).ue) / (2 * h);
      close(streamtubeEdgeVelocityTangent(cell.value, cell.apply(tangent), .2, interpolation), fd);
    }
  }
  const cell = linearizeIncompressibleStreamtubeCell(parameters);
  for (const interpolation of ['arithmetic','distance-weighted']) for (const {tangent,shift} of changes.slice(0,12)) {
    const plus=structuredClone(parameters),minus=structuredClone(parameters),h=1e-6;shift(plus,h);shift(minus,-h);
    const fd=(streamtubeEdgeVelocity(evaluateIncompressibleStreamtubeCell(plus),.2,interpolation).ue
      -streamtubeEdgeVelocity(evaluateIncompressibleStreamtubeCell(minus),.2,interpolation).ue)/(2*h);
    close(streamtubeEdgeVelocityTangent(cell.value,cell.apply(tangent),.2,interpolation),fd);
  }
});

test('section-speed matching preserves rotations, translations, length scaling and its stated incompressible limit', () => {
  for (const interpolation of ['arithmetic','distance-weighted']) {
  const cell = evaluateStreamtubeCell(parameters), expected = streamtubeEdgeVelocity(cell,.2,interpolation).ue;
  for (const length of [.1, 3]) for (const angle of [.7, -2]) {
    const map = p => ({ x: 2 + length * (p.x * Math.cos(angle) - p.y * Math.sin(angle)), y: -3 + length * (p.x * Math.sin(angle) + p.y * Math.cos(angle)) });
    const next = evaluateStreamtubeCell({ ...parameters, lower: parameters.lower.map(map), upper: parameters.upper.map(map), massFlow: parameters.massFlow * length });
    close(streamtubeEdgeVelocity(next,.2,interpolation).ue, expected, 1e-12);
  }
  }
  const cell=evaluateStreamtubeCell(parameters);
  const incompressible = evaluateIncompressibleStreamtubeCell(parameters);
  const edge = streamtubeEdgeVelocity(incompressible);
  close(edge.ue, .5 * (incompressible.states[0].q + incompressible.states[1].q), 1e-13);
  assert.equal(edge.correction, 0);
  assert.throws(() => streamtubeEdgeVelocity(cell, -1), /Invalid/);
  assert.throws(() => streamtubeEdgeVelocity({ ...cell, geometry: { pressureCurvature: NaN } }), /Invalid/);
});

test('distance weighting is linear-exact at the station on unequal intervals and retains second-order scalar interpolation', () => {
  // Scalar interpolation oracle: section samples are halfway along each
  // adjacent centerline segment. This does not assert an exact Euler state.
  for (const ratio of [.2,1,2.5,7]) for(const h of [.1,.01,.001]) {
    const left=h,right=ratio*h,slope=.4,curvature=.7;
    const at=s=>1+slope*s+curvature*s*s;
    const cell={states:[{q:at(-left/2)},{q:at(right/2)}],geometry:{pressureCurvature:0,streamwiseLengths:[left,right]}};
    const weighted=streamtubeEdgeVelocity(cell,.2,'distance-weighted');
    close(weighted.ue,1+curvature*left*right/4,1e-14);
    const linear={...cell,states:[{q:1-slope*left/2},{q:1+slope*right/2}]};
    close(streamtubeEdgeVelocity(linear,.2,'distance-weighted').ue,1,1e-14);
    close(streamtubeEdgeVelocity(linear).ue-1,slope*(right-left)/4,1e-14);
    assert.ok(weighted.weights.every(w=>w>0&&w<1));
    close(weighted.weights[0]+weighted.weights[1],1,1e-14);
    if(ratio===1)assert.equal(weighted.ue,streamtubeEdgeVelocity(cell).ue);
  }
});

test('distance weighting preserves the printed correction and exposes a geometry-dependent weight derivative', () => {
  const cell=linearizeStreamtubeCell(parameters),weighted=streamtubeEdgeVelocity(cell.value,.2,'distance-weighted');
  assert.equal(weighted.correction,streamtubeEdgeVelocity(cell.value).correction);
  const tangent=cell.apply({lower:[{x:0,y:0},{x:1,y:0},{x:0,y:0}]});
  const frozen={...tangent,geometry:{...tangent.geometry,streamwiseLengths:[0,0]}};
  assert.ok(Math.abs(streamtubeEdgeVelocityTangent(cell.value,tangent,.2,'distance-weighted')
    -streamtubeEdgeVelocityTangent(cell.value,frozen,.2,'distance-weighted'))>1e-4);
  assert.throws(()=>streamtubeEdgeVelocity(cell.value,.2,'unknown'),/interpolation/);
  assert.throws(()=>streamtubeEdgeVelocity({...cell.value,geometry:{...cell.value.geometry,streamwiseLengths:[0,1]}},.2,'distance-weighted'),/positive/);
});

test('surface-streamtube velocity tends to the independent vortex wall speed as the wall-adjacent tube is refined', t => {
  const rows = [];
  for (const h of [.08, .04, .02, .01]) {
    // Exact irrotational incompressible vortex: q(r)=1/r and mass
    // integral from wall r=1 to outer streamline r=1+h is log(1+h).
    const at = r => [-h, 0, h].map(a => ({ x: r * Math.cos(a), y: r * Math.sin(a) }));
    const cell = evaluateIncompressibleStreamtubeCell({ lower: at(1 + h), upper: at(1), massFlow: Math.log1p(h) });
    const error = Math.abs(streamtubeEdgeVelocity(cell).ue - 1);
    if (rows.length) assert.ok(error / rows.at(-1).error < .6);
    rows.push({ h, error });
  }
  assert.ok(rows.at(-1).error < .005); t.diagnostic(JSON.stringify(rows));
});
