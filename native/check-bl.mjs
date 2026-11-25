// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { checkAutomaticTransition, evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
import { evaluateLeadingTransitionInterval } from '../src/viscous/leading-transition-interval.js';

const binary = fileURLToPath(new URL('./check-bl', import.meta.url));
const fixture = name => JSON.parse(readFileSync(new URL(`../tests/fixtures/fortran/${name}.json`, import.meta.url)));
const regimes = ['similarity', 'laminar', 'turbulent', 'wake', 'transition'];
const props = ['rho', 'viscosity', 'machSquared', 'h', 'rawHk', 'hk', 'hStar', 'reTheta', 'cf',
  'dissipation', 'equilibriumShear', 'transitionShear'];
const cases = [];
for (const { input, expected } of fixture('integral').cases) {
  if (input.parameters.velocityConvention !== 'physical') continue;
  const [a,b] = input.states;
  if (input.regime === 'te') {
    const m = input.matched;
    const half = { aux: m.aux, theta: .5*m.theta, deltaStar: .5*m.deltaStar };
    const [upper,lower] = input.surfaceStates ?? [half,half];
    cases.push({ ...input, upstream: { ...a,...upper }, downstream: { ...a,...lower }, wake: b, expected });
  } else cases.push({ ...input, upstream: a, downstream: b, expected });
}
for (const name of ['transonic-integral','bl-active-limits'])
  for (const { name: label, parameters, input, expected } of fixture(name).cases)
    cases.push({ name: label, ...input, parameters, expected });
for (const { input, expected } of fixture('automatic-transition').cases)
  cases.push({ ...input, regime: 'transition', expected });
const criticalCycles = JSON.parse(readFileSync(new URL('../tests/fixtures/transition-root-failures.json',import.meta.url)));
for (const [i,{ payload:{ input,parameters } }] of criticalCycles.cases.entries())
  cases.push({ ...input,name:`critical amplification cycle ${i}`,parameters,regime:'transition' });

const referenceLaminar = { s: .24, aux: .3, theta: .00045, deltaStar: .00115, ue: 1.06 };
for (const tripS of [.27, Number.MAX_VALUE]) cases.push({ name: 'transition check', onlyTransition: true,
  parameters: { reynolds: 1e6, mach: .3 }, tripS, upstream: referenceLaminar,
  downstream: { s: .31, theta: .00056, deltaStar: .0013, ue: 1.08, aux: .025 } });
cases.push({ name: 'inadmissible supersonic shape factor', regime: 'similarity',
  parameters: { reynolds: 1e6, mach: .74 },
  downstream: { s: .25, ue: 1.6, theta: .001, deltaStar: .00101 }, domain: true });

function encode(s) {
  return [s.s,s.theta,s.deltaStar,s.ue,s.aux ?? 0,s.wakeGap ?? 0].join(' ');
}
function serialize(c) {
  const { reynolds = 1e6, mach = 0, gamma = 1.4, ncrit = 9 } = c.parameters;
  const op = c.onlyTransition ? 5 : c.regime === 'te' ? 6 : regimes.indexOf(c.regime);
  assert.ok(op >= 0, c.name);
  const rows = [[op,reynolds,mach,gamma,ncrit,op === 6 ? c.gap ?? 0 : c.tripS ?? Number.MAX_VALUE,
    c.similarityExponent ?? 1].join(' '),
    encode(c.upstream ?? c.downstream), encode(c.downstream)];
  if (op === 6) rows.push(encode(c.wake));
  return rows.join('\n');
}

const output = execFileSync(binary, { input: cases.map(serialize).join('\n') + '\n', encoding: 'utf8' })
  .trim().split('\n').map(row => row.split(' ').map(Number));
assert.equal(output.length, cases.length);
const near = (actual,expected,label,tolerance=5e-11) => {
  assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) <= tolerance*Math.max(1,Math.abs(expected)),
    `${label}: ${actual} vs ${expected}`);
};
for (const [i,c] of cases.entries()) {
  const [status,r1,r2,r3,active,forced,s,amp,...properties] = output[i];
  if (c.domain) { assert.equal(status,2); continue; }
  assert.equal(status,0,`${c.name}, status ${status}`);
  const kernel = createIntegralKernel({ ...c.parameters, exactJacobian: true, transitionTolerance: 1e-12 });
  const tripS = c.tripS ?? Number.MAX_VALUE;
  if (c.onlyTransition) {
    const r = checkAutomaticTransition(kernel, { ...c, tripS });
    assert.equal(active,+r.transition);
    assert.equal(forced,+r.forced);
    near(amp,r.amplification,c.name);
    near(s,r.s,c.name);
    continue;
  }
  const reference = c.regime === 'te' ? kernel.trailingEdge(c.upstream,c.downstream,c.wake,c.gap) :
    c.regime === 'transition' ? evaluateTransitionInterval(kernel,{ ...c, tripS }) : kernel.interval(c);
  [r1,r2,r3].forEach((r,j) => near(r,reference.residual[j],`${c.name}, residual ${j}`));
  // The archived natural-transition values use TRCHEK's original looser
  // stopping tolerance. Compare the native C roots to the strict JS kernel.
  if (c.regime !== 'transition' && c.expected?.residual)
    [r1,r2,r3].forEach((r,j) => near(r,c.expected.residual[j],`${c.name}, Fortran residual ${j}`,3e-9));
  if (reference.transition) {
    assert.equal(active,1);
    assert.equal(forced,+reference.transition.forced);
    near(s,reference.transition.s,`${c.name}, transition coordinate`);
    near(amp,reference.transition.amplification,`${c.name}, transition amplification`);
  }
  for (const [j,key] of props.entries()) near(properties[j],reference.properties[key],`${c.name}, ${key}`);
}
console.log(`Native boundary-layer parity passed on ${cases.length} compressible, transition, and wake blocks.`);

function march({ edges, parameters, tripS = Number.MAX_VALUE, similarityExponent = 0, junction }) {
  const { reynolds = 1e6, mach = 0, gamma = 1.4, ncrit = 9 } = parameters;
  const input = [[junction ? 8 : 7,reynolds,mach,gamma,ncrit,tripS,similarityExponent].join(' '),
    ...(junction ?? []).map(encode), edges.length, ...edges.map(e => [e.s,e.ue,e.wakeGap ?? 0].join(' '))].join('\n') + '\n';
  const rows = execFileSync(binary,{ input,encoding:'utf8' }).trim().split('\n').map(row => row.split(' ').map(Number));
  const [status,complete,iterations,,maximumResidual,active,forced,s] = rows[0];
  assert.equal(status,0,`BL profile stopped after ${complete}/${edges.length} stations`);
  assert.equal(complete,edges.length);
  assert.ok(iterations > 0 && maximumResidual < 1e-9);
  const kernel = createIntegralKernel({ ...parameters,exactJacobian:true,transitionTolerance:1e-12 });
  let upstream;
  for (const [i,row] of rows.slice(1).entries()) {
    const [kind,s,theta,deltaStar,ue,aux,wakeGap] = row;
    assert.equal(ue,edges[i].ue);
    assert.equal(s,edges[i].s);
    const downstream = { s,theta,deltaStar,ue,aux,wakeGap }, regime = regimes[kind];
    const reference = junction && i === 0 ? kernel.trailingEdge(...junction,downstream,wakeGap) :
      regime === 'transition' ? (i ? evaluateTransitionInterval(kernel,{ upstream,downstream,tripS }) :
        evaluateLeadingTransitionInterval(kernel,{ downstream,tripS })) :
      kernel.interval({ upstream,downstream,regime,similarityExponent });
    for (const v of reference.residual) near(v,0,`profile station ${i}, ${regime}`,1e-9);
    upstream = downstream;
  }
  return { last: upstream, active,forced,s };
}

const edges = Array.from({ length:61 },(_,i) => ({ s:10**(-4+4*i/60),ue:1 }));
const laminar = march({ edges,parameters:{ ncrit:25 } });
assert.equal(laminar.active,0);
assert.ok(Math.abs(laminar.last.theta/.000664-1)<.01);
const natural = march({ edges,parameters:{ mach:.2,ncrit:2 } });
assert.equal(natural.active,1);
assert.equal(natural.forced,0);
assert.ok(natural.s>0 && natural.s<1);
const forced = march({ edges,parameters:{ mach:.2 },tripS:.43 });
assert.equal(forced.active,1);
assert.equal(forced.forced,1);
assert.equal(forced.s,.43);
const leading = march({ edges,parameters:{ mach:.2 },tripS:.00004,similarityExponent:1 });
assert.equal(leading.s,.00004);
const wakeEdges = Array.from({length:25},(_,i) => ({s:1+.02*i,ue:1-.03*(1-Math.exp(-i/5)),
  wakeGap:.0004*Math.max(0,1-i/5)**3}));
const wake = march({ edges:wakeEdges,parameters:{mach:.2},junction:[forced.last,forced.last] });
assert.ok(wake.last.deltaStar>wake.last.theta && wake.last.aux>0);
console.log('Native BL profiles passed for laminar, natural and forced transition, a leading trip, and a finite-base wake.');
