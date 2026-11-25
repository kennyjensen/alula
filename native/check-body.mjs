// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { createContourTopology, createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { initialStreamtubeDisplacement, streamtubeBaseGeometry } from '../src/euler/streamtube-geometry.js';
import { streamtubeEquationOrder } from '../src/euler/streamtube-body-layout.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { solveNewton } from '../src/numerics/newton.js';
import { solveSparseDirectAligned } from '../src/numerics/klu.js';
import { intrinsicBodyFixture } from '../tests/fixtures/intrinsic-body.js';
import { twoActiveFiniteBaseWakes } from '../tests/fixtures/two-active-finite-base-wakes.js';

const binary=fileURLToPath(new URL('./check-body',import.meta.url));
execFileSync(binary,['--check-sparse']);
const families=['streamwise','inletDensity','internalPressure','farfieldPressure','cutPressure',
  'endTangency','wakeGap','trailingKutta','leadingKutta','farfieldMatch'];
const upwind={ mucon:1, mcrit:.5, boundary:{kind:'unfiltered-first-two'} };
const R=process.env.ALULA_BODY_BENCH === '1' ? 100 : 1;

function linearization(system, state, value) {
  const {layout:l,conditions:c}=system;
  const dir=state.map((_, i)=>.1*Math.sin(.27*i+.73));
  const includeDisplacement=l.displacedBoundaries;
  const widths=system.displacementParameters.map((_, i)=>.0003*Math.cos(.7*i));
  const packed=[...dir,...widths];
  const dn=system.geometryDerivatives(state,{includeDisplacement}).flat(2).map(map=>{
    const d=[0,0];
    for (const [col,v] of map) { d[0]+=packed[col]*v.x; d[1]+=packed[col]*v.y; }
    return d;
  });
  const drho=value.sections.flat(2).map((v,i)=>v.rho*dir[i]);
  const dm=value.allocation.groups.flatMap(group=>group.map(({derivatives})=>{
    let d=0;
    for (const [body, v] of derivatives) d+=c.massScale*dir[l.globals.capture[body]]*v;
    return d;
  }));
  const dw=l.independentWakeBanks?widths.filter((_,i)=>system.displacementParameters[i].kind==='wake'):[];
  const ds=['circulation','source','doubletX','doubletY'].map((k,i)=>
    c.lengthScale**(i<2?1:2)*dir[l.globals[k]]);
  const block=system.jacobian(state,{sparse:true,includeDisplacement});
  const expected=sparseProduct(includeDisplacement?block.state:block,dir);
  if (includeDisplacement)
    for (const [i,row] of block.displacement.entries())
      for (const [col,v] of row) expected[i]+=v*widths[col];
  else {
    const h=1e-6;
    const rp=system.residual(state.map((v,i)=>v+h*dir[i]));
    const rm=system.residual(state.map((v,i)=>v-h*dir[i]));
    for (const [i,v] of expected.entries())
      assert.ok(Math.abs((rp[i]-rm[i])/(2*h)-v)<2e-6*Math.max(1,Math.abs(v)), `difference at row ${i}`);
  }
  return {input:[...dn,drho,dm,dw,ds],expected};
}

function newtonChart(system,state,value,iterations) {
  const {layout:l,conditions:c}=system;
  const rows=system.geometryDerivatives(state).flat(2).flatMap(map=>[
    [...map].map(([k,p])=>[k,p.x]), [...map].map(([k,p])=>[k,p.y]),
  ]);
  for(let i=0;i<l.densityCount;i++) rows.push([[i,1]]);
  for(const group of value.allocation.groups)
    for(const {derivatives} of group) rows.push([...derivatives].map(([k,v])=>[l.globals.capture[k],c.massScale*v]));
  ['circulation','source','doubletX','doubletY'].forEach((key,i)=>
    rows.push([[l.globals[key],c.lengthScale**(i<2?1:2)]]));
  if (l.independentWakeBanks)
    for (const _ of system.displacement.wakes.flat()) rows.push([]);
  const compressed=rows.map(row=>row.filter(([k,v])=>Number.isInteger(k)&&v!==0).sort((a,b)=>a[0]-b[0]));
  const count=compressed.reduce((a,row)=>a+row.length,0);
  const columnToRow=streamtubeEquationOrder(l), order=new Int32Array(l.n);
  columnToRow.forEach((row,col)=>order[row]=col);
  const d=state.map((_,i)=>Math.cos(.17*i+.4));
  return {input:[[count,iterations,1e-10],...compressed.map(row=>[row.length,...row.flat()]),Array.from(order),
    Array.from(state),Array.from(d)], expected:sparseProduct(system.jacobian(state,{sparse:true}),d)};
}

function physical(system,value) {
  const {layout:l}=system;
  return [
    ...value.nodes.flat(2).flatMap(({x,y})=>[x,y]),...value.sections.flat(2).map(s=>s.rho),
    ...value.allocation.groups.flatMap(group=>group.map(m=>m.massFlow)),
    ...['circulation','source','doubletX','doubletY'].map(k=>value.strengths[k]),
    ...(l.independentWakeBanks?system.displacement.wakes.flat():[]),
  ];
}

function coordinateFixture(system,state,iterations) {
  const {layout:l,conditions:c}=system, base=system.decode(system.initial);
  const input=[
    [l.primaryBody,+(c.stagnationMotion==='walls-only'),+(c.normalStencil==='body-stations'),
      +(c.wakeDisplacementMotion==='te-center'),iterations,1e-10],
    ...system.originalNodes.flat(2).map(p=>[p.x,p.y]),
    base.allocation.groups.flatMap(group=>group.map(m=>m.massFlow)),
    base.captured,
  ];
  for(const [g,b] of l.bodies.entries()) {
    const points=b.trailingEdge?.kind==='finite-base' ?
      createContourTopology(b.points,b).surface.points : b.points;
    input.push([points.length,system.initialStagnation[g],+(b.wakeTangentialReference==='material-te')],
      ...points.map(p=>[p.x,p.y]),system.fractions[g].lower,system.fractions[g].upper);
    if(l.displacedBoundaries) input.push(system.displacement.surfaces[g].lower,
      system.displacement.surfaces[g].upper,system.displacement.wakes[g]);
  }
  const d=state.map((_,i)=>.1*Math.cos(.312*i+.54));
  return {input:[...input,Array.from(state),Array.from(d)],
    expected:sparseProduct(system.jacobian(state,{sparse:true}),d)};
}

function call(system, state, value, settings={}) {
  const {layout:l,conditions:c}=system;
  const mode=['momentum','isentropic','hybrid'].indexOf(c.streamwiseMode);
  assert.ok(mode>=0);
  const points=value.nodes.flat(2), densities=value.sections.flat(2).map(s=>s.rho);
  const nodes=points.map(p=>[p.x,p.y]);
  const mass=value.allocation.groups.flatMap(g=>g.map(s=>s.massFlow));
  const widths=l.independentWakeBanks ? system.displacement.wakes.flat() : [];
  const strengths=value.strengths;
  settings.edit?.({nodes,densities,mass,widths,strengths});
  const repeats=settings.repeat ?? R;
  const derivative=settings.derivative ? linearization(system,state,value) : undefined;
  const newton=settings.newton===undefined?undefined:newtonChart(system,state,value,settings.newton);
  const chart=settings.chart===undefined?undefined:coordinateFixture(system,state,settings.chart);
  const input=[
    [l.nx,l.elements,mode,repeats,chart?3:newton?2:+Boolean(derivative)],
    [c.gamma,c.mach,c.alpha,c.lengthScale,c.pressureCorrectionFactor,c.center.x,c.center.y],
    [+Boolean(c.upwind),c.hybrid?.ismom ?? 4,+(c.geometryDomain==='positive-simple'),
      +l.displacedBoundaries,+l.independentWakeBanks,+(c.wakeOutlet==='banks'),
      c.upwind?.mucon ?? 1,c.upwind?.mcrit ?? .99,c.hybrid?.epsilonP ?? 1e-5],
    l.tubes,...l.bodies.map((b,i)=>{
      const r=c.hybrid?.entropyRegions?.regions[i];
      return [b.leadingIndex,b.trailingIndex,+b.roundedLeadingEdge,r?.throughRow ?? 0,r?.lowerTube ?? 0,r?.upperTube ?? 0];
    }),
    ...(chart?chart.input:[
      ...nodes,densities,mass,widths,[strengths.circulation,strengths.source,strengths.doubletX,strengths.doubletY],
      ...(derivative?.input ?? newton?.input ?? []),
    ]),
  ].map(row=>row.join(' ')).join('\n')+'\n';
  const rows=execFileSync(binary,{input,encoding:'utf8'}).trim().split('\n').map(line=>line.split(' ').map(Number));
  if(chart) {
    const expected=physical(system,value);
    assert.equal(rows[0].length,expected.length);
    for(const [i,x] of rows[0].entries())
      assert.ok(Math.abs(x-expected[i])<8e-13*Math.max(1,Math.abs(expected[i])),`physical chart scalar ${i}: ${x} vs ${expected[i]}`);
    assert.equal(rows[1].length,l.n);
    for(const [i,x] of rows[1].entries())
      assert.ok(Math.abs(x-chart.expected[i])<3e-10*Math.max(1,Math.abs(chart.expected[i])),`chart ${l.rows[i].kind} at ${i}: ${x} vs ${chart.expected[i]}`);
    const [status,updates,backtracks,colors,norm,error,ms]=rows[2];
    assert.ok(error<3e-14);
    return {status,updates,backtracks,colors,norm,error,ms,result:rows[3],x:rows[4]};
  }
  if (newton) {
    assert.equal(rows[0].length,l.n);
    for (const [i,x] of rows[0].entries())
      assert.ok(Math.abs(x-newton.expected[i])<=3e-10*Math.max(1,Math.abs(newton.expected[i])),
        `sparse ${l.rows[i].kind} at row ${i}: ${x} vs ${newton.expected[i]}`);
    const [status,iterations,backtracks,colors,norm,mach,ms,factors,refactors,linearResidual,reuses]=rows[1];
    return {status,iterations,backtracks,colors,norm,mach,ms,factors,refactors,linearResidual,reuses,x:Float64Array.from(rows[2])};
  }
  const [status,n,ms]=rows[0];
  assert.equal(n,l.n);
  if (settings.status !== undefined) {
    assert.equal(status,settings.status);
    return;
  }
  assert.equal(status,0);
  const [mach,...diagnostics]=rows[1];
  assert.ok(Math.abs(mach-value.diagnostics.maxMach)<2e-12);
  assert.equal(rows[2].length,l.n);
  let err=0;
  for (const [i,x] of rows[2].entries()) {
    const d=Math.abs(x-value.residual[i]);
    err=Math.max(err,d);
    assert.ok(d<=2e-11*Math.max(1,Math.abs(value.residual[i])),
      `${l.rows[i].kind} at row ${i}: ${x} vs ${value.residual[i]}`);
  }
  for (const [i,f] of families.entries()) assert.ok(Math.abs(diagnostics[i]-(value.diagnostics.residualByFamily[f] ?? 0))<2e-11);
  if (derivative) {
    assert.equal(rows[3].length,l.n);
    for (const [i,x] of rows[3].entries())
      assert.ok(Math.abs(x-derivative.expected[i])<=3e-10*Math.max(1,Math.abs(derivative.expected[i])),
        `${l.rows[i].kind} derivative at row ${i}: ${x} vs ${derivative.expected[i]}`);
  }
  if (repeats>1 && !derivative) {
    system.evaluate(state);
    const t=performance.now();
    for (let i=0;i<repeats;i++) system.evaluate(state);
    console.log(`${l.n} rows: C ${ms.toFixed(3)} ms, JS ${((performance.now()-t)/repeats).toFixed(3)} ms, max difference ${err}`);
  }
}

function moved(system) {
  const {layout:l}=system;
  const x=system.initial.map((v,i)=>v+(i<l.densityCount?.002:.000005)*Math.sin(.7*(i+1)));
  for (const [k,v] of [['circulation',.03],['source',.005],['doubletX',.001],['doubletY',-.002]]) x[l.globals[k]]=v;
  return x;
}

let count=0;
for (const elements of [1,2]) for (const mode of ['momentum','isentropic',1,2,3,4]) {
  const input=intrinsicBodyFixture({elements,bodySegments:8,tubes:3,mach:.5,alpha:.35});
  const hybrid=Number.isInteger(mode);
  const system=createStreamtubeBodySystem({...input,streamwiseMode:hybrid?'hybrid':mode,
    ...(mode==='isentropic'?{}:{upwind}),...(hybrid?{hybrid:{ismom:mode,epsilonP:1e-5}}:{})});
  const state=moved(system), value=system.evaluate(state);
  call(system,state,value);
  call(system,state,value,{repeat:1,derivative:true});
  count++;
  if (elements===2 && mode===4) {
    call(system,state,value,{status:2,edit:({mass})=>{mass[0]=0;}});
    const section=value.sections[0][0][0], gm=system.conditions.gamma-1;
    const sonicSpeed=Math.sqrt(gm*system.conditions.h0*1.1/(1+.5*gm*1.1));
    call(system,state,value,{status:3,edit:({densities})=>{densities[0]*=section.q/sonicSpeed;}});
  }
}

const {system:coupled,x,input:wakeInput}=twoActiveFiniteBaseWakes();
const euler=coupled.euler, xe=x.slice(0,coupled.ne);
coupled.evaluate(x); // installs the two BL-driven wake widths in the outer system
const wakeValue=euler.evaluate(xe);
call(euler,xe,wakeValue);
call(euler,xe,wakeValue,{repeat:1,derivative:true});
call(euler,xe,wakeValue,{status:1,edit:({widths})=>{widths[0]=-1;}});
const separated={...wakeValue,nodes:structuredClone(wakeValue.nodes)};
const layout=euler.layout, cut=layout.elements-1, i=layout.nx, j=layout.tubes[cut];
const previous=separated.nodes[cut+1][i][0].y;
separated.nodes[cut+1][i][0].y=separated.nodes[cut][i][j].y-.25*(previous-separated.nodes[cut][i][j].y);
call(euler,xe,separated,{status:4});
count++;

const input=intrinsicBodyFixture({elements:2,alpha:2,tubes:5,tubeGrowth:3});
const system=createStreamtubeBodySystem(input);
const root=solveStreamtubeBody(system,{tolerance:1e-11,maxIterations:15});
assert.equal(root.converged,true,root.reason);
call(system,root.x,system.evaluate(root.x));
count++;

const large=createStreamtubeBodySystem({
  ...intrinsicBodyFixture({elements:2,bodySegments:32,tubes:20,tubeGrowth:1.12,mach:.3}),
  streamwiseMode:'hybrid',upwind,hybrid:{ismom:4,epsilonP:1e-5},
});
const xl=large.initial, vl=large.evaluate(xl);
call(large,xl,vl,{repeat:R>1?30:1});
call(large,xl,vl,{repeat:1,derivative:true});
console.log(`Native residual parity passed for ${count+1} airfoil states and Jacobian-product parity for ${count}.`);


for(const options of [{elements:1,tubes:3,mach:.2},{elements:2,tubes:5,tubeGrowth:3,alpha:2}]) {
  const input=intrinsicBodyFixture(options);
  input.bodies.forEach(body=>{body.roundedLeadingEdge=false;});
  const body=createStreamtubeBodySystem(input);
  const initial=body.initial;
  const solved=call(body,initial,body.evaluate(initial),{newton:20});
  assert.equal(solved.status,0,JSON.stringify({...solved,x:undefined}));
  assert.ok(solved.iterations>1 && solved.refactors+solved.reuses>0 && solved.norm<=1e-10 && solved.linearResidual<=1e-10);
  assert.ok(solved.colors<50);
  assert.ok(body.evaluate(solved.x).diagnostics.residual<1e-10);
  const order=streamtubeEquationOrder(body.layout);
  const t=performance.now();
  const reference=solveNewton({initial,maxIterations:20,tolerance:1e-10,
    residual:body.residual,admissible:body.admissible,jacobian:x=>body.jacobian(x,{sparse:true}),
    linearSolve:(a,b)=>solveSparseDirectAligned(a,b,order).x});
  const js=performance.now()-t;
  assert.equal(reference.converged,true,reference.reason);
  for(const [i,x] of solved.x.entries())
    assert.ok(Math.abs(x-reference.x[i])<2e-8,`fixed-edge root at column ${i}: ${x} vs ${reference.x[i]}`);
  console.log(`Native sparse Euler ${options.elements}-element solve: ${initial.length} rows, ${solved.colors} colors, ${solved.iterations} updates, C ${solved.ms.toFixed(1)} ms / JS ${js.toFixed(1)} ms`);
}
const largeMatrix=call(large,xl,vl,{newton:0});
assert.equal(largeMatrix.status,1);
assert.ok(largeMatrix.colors<80);
const wakeMatrix=call(euler,xe,wakeValue,{newton:0});
assert.equal(wakeMatrix.status,1);
console.log('Sparse matrix parity passed for the large rounded-body and independent-wake grids.');


for (const opts of [
  {elements:1,alpha:1,tubes:3},
  {elements:2,alpha:2,tubes:5,tubeGrowth:3},
]) {
  const body=createStreamtubeBodySystem(intrinsicBodyFixture(opts));
  const initial=moved(body), v=body.evaluate(initial);
  const actual=call(body,initial,v,{chart:20});
  assert.equal(actual.status,0,JSON.stringify({...actual,result:undefined,x:undefined}));
  assert.ok(actual.updates>0 && actual.norm<=1e-10);
  const t=performance.now();
  const reference=solveStreamtubeBody(body,{initial:body.rebase(initial),maxIterations:20,tolerance:1e-10});
  const js=performance.now()-t;
  assert.equal(reference.converged,true,reference.reason);
  const target=physical(body,reference);
  for(const [i,actualValue] of actual.result.entries())
    assert.ok(Math.abs(actualValue-target[i])<3e-8*Math.max(1,Math.abs(target[i])),
      `rebased root scalar ${i}: ${actualValue} vs ${target[i]}`);
  console.log(`Native rounded ${opts.elements}-element Euler: ${actual.updates} updates, residual ${actual.norm.toExponential(2)}, C ${actual.ms.toFixed(1)} ms / JS ${js.toFixed(1)} ms`);
}

const wakeCurves=wakeInput.bodies.map(body=>createSurfaceContourCurve(body.points,body));
const displacement=initialStreamtubeDisplacement({bodies:wakeInput.bodies,nx:wakeInput.outerLower.length-1},
  streamtubeBaseGeometry(wakeInput.bodies,wakeCurves));
const wakeChart=createStreamtubeBodySystem({...wakeInput,displacement});
const separatedChart=call(wakeChart,wakeChart.initial,wakeChart.evaluate(wakeChart.initial),{chart:0});
assert.equal(separatedChart.status,1);
