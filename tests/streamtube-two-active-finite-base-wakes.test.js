import test from 'node:test';
import assert from 'node:assert/strict';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';

const relativeError = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
const close = (a, b, tolerance, label) => assert.ok(relativeError(a, b) <= tolerance,
  `${label}: analytic ${a}, independent ${b}, error ${relativeError(a, b)}`);
const fd4 = (f, x, d, h) => {
  const at = k => f(x.map((v, j) => v + k * h * d[j]));
  const pp = at(2), p = at(1), m = at(-1), mm = at(-2);
  return p.map((_, i) => (-pp[i] + 8 * p[i] - 8 * m[i] + mm[i]) / (12 * h));
};
const dotMap = (map, d) => [...map].reduce((v, [col, a]) => v + a * d[col], 0);

// Independent expanded cubic with g(0)=ANTE, g'(0)=kappa and double zero
// at L=2.5*ANTE. Do not call the production dead-air-gap implementation.
function cubic(base, s) {
  const u = base.upperDerivative, l = base.lowerDerivative;
  const cross = (u.x * l.y - u.y * l.x) / (Math.hypot(u.x, u.y) * Math.hypot(l.x, l.y));
  const kappa = Math.max(-1.2, Math.min(1.2, cross / Math.sqrt(1 - cross * cross)));
  const g = base.width, length = 2.5 * g, q = s / length;
  return s >= length ? 0 : g + kappa * s - (3 * g + 2 * kappa * length) * q ** 2
    + (2 * g + kappa * length) * q ** 3;
}

function independentWakeData(system, x) {
  system.euler.setDisplacement(system.bl.thicknesses(x.subarray(system.ne)));
  const nodes = system.euler.decode(x.subarray(0, system.ne)).nodes, length = system.euler.conditions.lengthScale;
  return system.bl.wakes.map(wake => {
    const base = system.euler.baseGeometry[wake.body], te = system.euler.layout.bodies[wake.body].trailingIndex;
    const lower = nodes[wake.body][te].at(-1), upper = nodes[wake.body + 1][te][0];
    let previous = { x: .5 * (lower.x + upper.x), y: .5 * (lower.y + upper.y) }, arc = 0;
    return wake.ids.map((id, k) => {
      const i = system.bl.stations[id].i;
      if (k) {
        const a = nodes[wake.body][i].at(-1), b = nodes[wake.body + 1][i][0];
        const center = { x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) };
        arc += Math.hypot(center.x - previous.x, center.y - previous.y); previous = center;
      }
      return { id, i, arc, gap: cubic(base, arc) / length };
    });
  });
}

test('two distinct finite-base cubic tails share a convex unequal-mass passage and retain physical fluid BL shape', t => {
  const { system, x } = twoActiveFiniteBaseWakes(), value = system.evaluate(x);
  const quality = streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes }).quality;
  assert.equal(system.admissible(x), true); assert.equal(quality.valid, true);
  assert.equal(system.bl.surfaces.length, 4); assert.equal(system.bl.wakes.length, 2);
  const masses = value.outer.allocation.groups.map(g => g.map(t => t.massFlow));
  assert.equal(masses.length, 3); masses.forEach(g => assert.notEqual(g[0], g[1]));
  assert.notEqual(system.euler.baseGeometry[0].width, system.euler.baseGeometry[1].width);
  const geo = system.bl.geometry(x.subarray(0, system.ne), true), oracle = independentWakeData(system, x);
  for (const [b, rows] of oracle.entries()) {
    const active = rows.filter(p => p.arc > 0 && p.gap > 0);
    assert.equal(active.length, 3); assert.equal(rows.at(-1).gap, 0);
    for (const p of rows) {
      close(geo.coordinates[p.id].wakeGap, p.gap, 3e-16, `body${b} station${p.id} cubic`);
      close(geo.coordinates[p.id].wakeDistance, p.arc / system.euler.conditions.lengthScale, 3e-16, 'physical TE-center arc');
      const state = value.layers.states[p.id];
      assert.ok(state.deltaStar - state.wakeGap > state.theta);
    }
    const te = rows[0].id, surfaces = system.bl.surfaces.filter(s => s.body === b);
    const total = surfaces.reduce((v, s) => v + value.layers.states[s.ids.at(-1)].deltaStar, 0);
    close(value.layers.states[te].deltaStar - total, system.euler.baseGeometry[b].width / system.euler.conditions.lengthScale, 2e-16, 'TE gap added once');
  }
  assert.deepEqual(oracle[0].slice(1,4).map(p => p.i), oracle[1].slice(1,4).map(p => p.i));
  t.diagnostic(JSON.stringify({ unknowns: system.n, eulerUnknowns: system.ne, cells: system.euler.layout.nx * system.euler.layout.nt,
    quality, masses, activeTailStations: oracle.map(rows => rows.slice(1,4)), families: value.families,
    scope: 'Admissible manufactured off-root state; no BL initialization, Newton or linear solve.' }));
});

test('each moving mean-bank arc differentiates its own active tail, with no opposite-body, interior or stagnation contamination', t => {
  const { system, x } = twoActiveFiniteBaseWakes(); system.evaluate(x);
  const geo = system.bl.geometry(x.subarray(0, system.ne), true), layout = system.euler.layout;
  const first = system.bl.wakes.map(w => system.bl.stations[w.ids[1]].i), cases = [];
  // Four independent physical banks, including both sides of passage 1.
  for (let body = 0; body < 2; body++) for (const side of ['upper','lower']) {
    const p = layout.positions.find(p => p.kind === 'cut' && p.body === body && p.side === side && p.i === first[body]);
    assert.ok(p); const d = new Float64Array(system.n); d[p.column] = .01;
    cases.push({ label: `body${body}-${side}`, body, d });
  }
  for (const [label, col] of [
    ['shared-interior', layout.nodes[1][first[0]][1].column],
    ['stagnation0', layout.globals.stagnation[0]], ['stagnation1', layout.globals.stagnation[1]],
  ]) { const d = new Float64Array(system.n); d[col] = .01; cases.push({ label, body: null, d }); }
  let maximum = 0; const changes = [];
  for (const { label, body, d } of cases) {
    const numerical = fd4(z => independentWakeData(system, z).flat().map(p => p.gap), x, d, 2e-6);
    let row = 0; const byBody = [];
    for (const wake of system.bl.wakes) {
      let norm = 0;
      for (const id of wake.ids) {
        const analytic = dotMap(geo.coordinates[id].wakeGapDerivatives, d);
        close(analytic, numerical[row], 2e-10, `${label} gap station${id}`);
        maximum = Math.max(maximum, relativeError(analytic, numerical[row++])); norm = Math.max(norm, Math.abs(analytic));
      }
      byBody.push(norm);
    }
    if (body === null) assert.deepEqual(byBody, [0,0]);
    else { assert.ok(byBody[body] > 1e-5, 'Active cubic sensitivity must be measurable.'); assert.equal(byBody[1-body], 0); }
    changes.push({ label, byBody });
  }
  t.diagnostic(JSON.stringify({ maximum, independentBankCases: changes }));
});

test('complete two-tail coupled and fluid-domain Jacobians cover both bodies, wakes, all edge columns and unequal mass motion', t => {
  const { system, x } = twoActiveFiniteBaseWakes(), before = x.slice();
  const j = system.jacobian(x, { sparse: false }), constraints = system.stepConstraints(x);
  const directions = [.41,.83,1.31].map(phase => ({ name: `full-mixed-${phase}`, d: x.map((v,k) =>
    Math.sin((k+1)*phase) * (k < system.ne ? .01 : .2 * Math.max(.03, Math.abs(v)))) }));
  for (let body=0; body<2; body++) {
    const d = new Float64Array(system.n);
    for (const s of system.bl.stations.filter(s => s.body === body)) for (let k=0;k<4;k++) {
      const col = system.ne + 4*s.id+k; d[col] = Math.cos(.57*(col+1)) * .2 * Math.max(.03, Math.abs(x[col]));
    }
    directions.push({name:`body${body}-all-BL-and-edge`,d});
  }
  const capture = new Float64Array(system.n); capture[system.euler.layout.globals.capture[1]] = .01;
  directions.push({name:'unequal-passage-mass-capture',d:capture});
  const maxima = {}, comparisons = { full: 0, shape: 0 }; let shapeMaximum = 0;
  const rowName = row => row < system.ne ? `Euler:${system.euler.layout.rows[row].kind}`
    : `body${system.bl.stations[Math.floor((row-system.ne)/4)].body}:${system.bl.stations[Math.floor((row-system.ne)/4)].kind}:${(row-system.ne)%4===3?'edge':'BL'}`;
  for (const {name,d} of directions) for (const h of [2e-6,1e-6]) {
    const numerical = fd4(z => system.residual(z), x, d, h);
    for (let row=0;row<system.n;row++) {
      let exact=0;for(let col=0;col<system.n;col++) exact+=j[row*system.n+col]*d[col];
      const error=relativeError(exact,numerical[row]), family=rowName(row);
      maxima[family]=Math.max(maxima[family]??0,error); comparisons.full++;
      close(exact,numerical[row],5e-6,`${name} h${h} row${row} ${family}`);
    }
    const numericalConstraints=fd4(z=>system.constraintValues(z),x,d,h);
    constraints.forEach((c,row)=>{
      if(c.kind!=='kinematic-shape')return;
      const exact=dotMap(c.gradient,d); shapeMaximum=Math.max(shapeMaximum,relativeError(exact,numericalConstraints[row]));
      comparisons.shape++;close(exact,numericalConstraints[row],5e-7,`${name} physical shape${row}`);
    });
  }
  for(let body=0;body<2;body++)for(const part of ['surface','wake'])for(const family of ['BL','edge'])
    assert.ok(Object.hasOwn(maxima,`body${body}:${part}:${family}`));
  assert.deepEqual(x,before); assert.equal(system.admissible(x),true);
  t.diagnostic(JSON.stringify({unknowns:system.n,directions:directions.map(d=>d.name),steps:[2e-6,1e-6],comparisons,
    familyMaximumRelativeErrors:maxima,shapeMaximum,globalNewtonUpdates:0,linearSolves:0}));
});
