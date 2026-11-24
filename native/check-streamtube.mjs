// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';

import { linearizeStreamtubeSpeedUpwind } from '../src/euler/streamtube-speed-upwind.js';
import { linearizeStreamtubeMomentumBlend } from '../src/euler/streamtube-momentum-blend.js';

const binary = fileURLToPath(new URL('./check-streamtube', import.meta.url));
const zero = () => [{x: 0, y: 0}, {x: 0, y: 0}, {x: 0, y: 0}];
const blank = () => ({ massFlow: 0, stagnationEnthalpy: 0, pressureCorrectionFactor: 0,
  densities: [0, 0], transportSpeeds: [0, 0], lower: zero(), upper: zero() });
const nodes = (lower, upper) => lower.flatMap((p,i) => [p.x,p.y,upper[i].x,upper[i].y]);
const encode = ({p,d,epsilon,depsilon}) => [
  ...(epsilon === undefined ? ['C'] : ['B',epsilon,depsilon]),
  +(p.geometryDomain === 'positive-simple'), +(p.transportSpeeds !== undefined),
  p.massFlow,p.stagnationEnthalpy,p.gamma,p.pressureCorrectionFactor,...p.densities,
  ...(p.transportSpeeds ?? [0,0]),...nodes(p.lower,p.upper),
  d.massFlow,d.stagnationEnthalpy,d.pressureCorrectionFactor,...d.densities,...d.transportSpeeds,
  ...nodes(d.lower,d.upper),
].join(' ');
const flatten = v => {
  const g = v.geometry;
  return [g.area,g.pressureCurvature,...g.normalAreas,...g.streamwiseLengths,
    v.interfacePressure.lower,v.interfacePressure.upper,v.streamwiseResidual,v.entropyJump,
    v.isentropicResidual,v.pressureCorrection,
    ...v.states.flatMap(s => [s.rho,s.q,s.p,s.enthalpy,s.machSquared])];
};

const base = {
  lower: [{x:0,y:0},{x:.7,y:.02},{x:1.1,y:.08}],
  upper: [{x:.04,y:.36},{x:.73,y:.43},{x:1.18,y:.43}],
  densities: [1,.98], gamma: 1.4, massFlow: .13, stagnationEnthalpy: 1.2,
  pressureCorrectionFactor: .1,
};
const configurations = [
  base,
  {...base, transportSpeeds: [.39,.36], massFlow: .25},
  {...base, massFlow: .32, stagnationEnthalpy: .95},
  {...base, geometryDomain: 'positive-simple', massFlow: .015,
    lower: [{x:0,y:0},{x:1,y:0},{x:2,y:0}],
    upper: [{x:0,y:1},{x:.49,y:.49},{x:2,y:1}]},
  {...base, massFlow: .0001, lower: base.lower.map(p => ({x:p.x+32,y:p.y+16})),
    upper: base.lower.map((p,i) => ({x:p.x+32+.00001*i,y:p.y+16+.0003+.00002*i}))},
];
const cases = [];
for (const p of configurations) {
  const l = linearizeStreamtubeCell(p);
  const scalar = ['massFlow','stagnationEnthalpy','pressureCorrectionFactor'];
  const pairs = ['densities',...(p.transportSpeeds ? ['transportSpeeds'] : [])];
  for (const key of scalar) {
    const d = blank(); d[key] = 1;
    cases.push({p,d,value:flatten(l.value),derivative:flatten(l.apply(d))});
  }
  for (const key of pairs) for (let i=0;i<2;i++) {
    const d = blank(); d[key][i] = 1;
    cases.push({p,d,value:flatten(l.value),derivative:flatten(l.apply(d))});
  }
  for (const key of ['lower','upper']) for (let i=0;i<3;i++) for (const axis of ['x','y']) {
    const d = blank(); d[key][i][axis] = 1;
    cases.push({p,d,value:flatten(l.value),derivative:flatten(l.apply(d))});
  }
}
for (const [status,p] of [
  [1,{...base, densities: [1,-1]}],
  [2,{...base, upper: base.lower}],
  [2,{...configurations[3], geometryDomain: 'convex'}],
  [2,{...base, geometryDomain: 'positive-simple', upper: base.lower.toReversed()}],
  [3,{...base, stagnationEnthalpy: .01}],
  [5,{...base, transportSpeeds: [1.6,1.6]}],
]) {
  assert.throws(() => evaluateStreamtubeCell(p));
  cases.push({p,d:blank(),status});
}
const blendVector = v => [v.fraction,v.residual,...v.lossIndicators,v.compression];
const straight = { ...base,
  lower: [{x:0,y:0},{x:1,y:0},{x:2,y:0}], upper: [{x:0,y:.4},{x:1,y:.4},{x:2,y:.4}],
};
const fractions = [];
for (const density of [1.0008,1.005,1.12,.95]) {
  const p = { ...straight, densities:[1,density], massFlow:.3, transportSpeeds:[.76,.75] };
  const d = blank(); d.densities = [.002,.001]; d.transportSpeeds = [.0001,.0004]; d.upper[2].y=.002;
  const epsilon=1e-5, depsilon=1e-7;
  const l=linearizeStreamtubeCell(p), dl=l.apply(d);
  const blend=linearizeStreamtubeMomentumBlend({...l.value,epsilonP:epsilon});
  const derivative=blend.apply({...dl,epsilonP:depsilon});
  fractions.push(blend.value.fraction);
  cases.push({p,d,value:flatten(l.value),derivative:flatten(dl),epsilon,depsilon,
    blend:[0,...blendVector(blend.value),...blendVector(derivative)]});
}
assert.ok(fractions.includes(0) && fractions.includes(1) && fractions.some(x => x>0 && x<1), String(fractions));

const input = cases.map(encode).join('\n')+'\n';
const output = execFileSync(binary, {input, encoding:'utf8'}).trim().split('\n');
assert.equal(output.length, cases.length);
function close(actual, expected, label) {
  assert.equal(actual.length,expected.length);
  for (const [i,x] of actual.entries())
    assert.ok(Math.abs(x-expected[i]) <= 2e-11*Math.max(1,Math.abs(expected[i])),
      `${label} ${i}: ${x} != ${expected[i]}`);
}
for (const [i,line] of output.entries()) {
  const [status,...values] = line.split(' ').map(Number), c=cases[i];
  assert.equal(status,c.status ?? 0,`case ${i}`);
  if (c.status) continue;
  const n=c.value.length;
  close(values.slice(0,n),c.value,`case ${i} value`);
  assert.equal(values[n],0);
  close(values.slice(n+1,2*n+1),c.derivative,`case ${i} derivative`);
  close(values.slice(2*n+1),c.blend ?? [],`case ${i} blend`);
}
console.log(`Native Euler cell parity passed for ${cases.length} stencils and directions.`);

const speedBase = { speeds:[.94,.82,.65], machSquared:[1.05,1.14],
  spacing:[.14,.22], mucon:1, mcrit:.99, gamma:1.4 };
const upwindCases = [
  speedBase,
  {...speedBase, mucon:-1.4},
  {...speedBase, machSquared:[.01,.02]},
  {...speedBase, machSquared:[0,0]},
  {...speedBase, machSquared:[.99,1.01]},
  {...speedBase, mcrit:1},
  {...speedBase, mucon:0},
];
const upwindDir = p => ({ speeds:[.24,-.17,.05], machSquared:[.02,.04], spacing:[-.03,.08],
  mucon:p.mucon===0 ? 0 : -.07, mcrit:p.mcrit===1 ? 0 : .001, gamma:.03 });
const upwindInput = x => [...x.speeds,...x.machSquared,...x.spacing,x.mucon,x.mcrit,x.gamma];
const upwindOutput = v => [v.speed,v.coefficient,v.correction];
const upwindResults = execFileSync(binary, {encoding:'utf8',input:upwindCases.map(p =>
  ['U',...upwindInput(p),...upwindInput(upwindDir(p))].join(' ')).join('\n')+'\n'}).trim().split('\n');
assert.equal(upwindResults.length,upwindCases.length);
for (const [i,p] of upwindCases.entries()) {
  const l=linearizeStreamtubeSpeedUpwind(p);
  close(upwindResults[i].split(' ').map(Number),
    [0,...upwindOutput(l.value),...upwindOutput(l.apply(upwindDir(p)))],`upwind case ${i}`);
}
console.log(`Native upwind parity passed for ${upwindCases.length} stencils.`);
