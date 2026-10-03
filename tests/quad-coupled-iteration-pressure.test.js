import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledIterationPressure } from '../src/ui/quad-coupled-iteration-pressure.js';

const fixture = () => {
  const mach=.5, gamma=1.4, pInf=1/(gamma*mach*mach);
  const bodies=[{element:1,leadingIndex:1,trailingIndex:3},{element:0,leadingIndex:1,trailingIndex:3,
    trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:3}}];
  const nodes=Array.from({length:3},(_,g)=>Array.from({length:5},(_,i)=>[{x:i/10,y:g},{x:i/10,y:g+.5}]));
  const bl={stations:[],surfaces:[]};
  for(let body=0;body<2;body++)for(const side of ['upper','lower']){
    const ids=[];for(let i=2;i<=3;i++){const id=bl.stations.length;ids.push(id);bl.stations.push({id,kind:'surface',body,side,i});}
    bl.surfaces.push({body,side,ids});
  }
  const checkpoint={version:1,restart:{input:{mach,gamma,bodies},options:{blThermodynamics:'historical-common-isentrope'},
    initialBL:bl.stations.flatMap(()=>[.03,1,2,1]),initialEuler:{nodes:structuredClone(nodes),undisplacedNodes:structuredClone(nodes)}}};
  const flow={nodes:structuredClone(nodes),undisplacedNodes:structuredClone(nodes),
    cells:Array.from({length:4},(_,i)=>Array.from({length:3},(_,g)=>[{interfacePressure:{lower:pInf+(10*g+i+1)/2,upper:pInf+(10*g+i+1.5)/2}}]))};
  return {checkpoint,flow,bl,bodies,elements:[{name:'Main'},{name:'Flap'}],referenceChord:2,mach};
};

test('live Cp matches the saved browser terminal plot at floating-point precision with exact coordinates',()=>{
  const saved=JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json',import.meta.url)));
  const r=saved.result, input={checkpoint:saved.checkpoint,flow:r.flow,bl:r.numericalBoundaryLayer,
    bodies:saved.checkpoint.restart.input.bodies,elements:r.elements,referenceChord:r.referenceChord,mach:r.mach};
  const before=JSON.stringify(input), p=quadCoupledIterationPressure(input);
  assert.deepEqual(p.elements.map(e=>({name:e.name,points:e.cp.map(({x,y})=>({x,y}))})),
    r.elements.map(e=>({name:e.name,points:e.cp.map(({x,y})=>({x,y}))})));
  // Browser-captured and Node-recomputed isentropic pressure differ by one
  // last bit at one station despite identical saved Ue/Mach/gamma.
  p.elements.forEach((e,k)=>e.cp.forEach((point,j)=>assert.ok(
    Math.abs(point.cp-r.elements[k].cp[j].cp)<=4*Number.EPSILON*Math.max(1,Math.abs(point.cp)))));
  assert.equal(JSON.stringify(input),before);
  assert.equal(p.pressureKind,'isentropic BL-edge pressure from solved edge speed');
});

test('historical Cp uses actual interface pressure, original element order and both finite TE endpoints',()=>{
  const f=fixture(),before=JSON.stringify(f),p=quadCoupledIterationPressure(f);
  assert.deepEqual(p.elements.map(e=>e.name),['Main','Flap']);
  // Body1 is the user's first element. Upper is reversed, then lower.
  assert.deepEqual(p.elements[0].cp.map(({x,y})=>({x,y})),[
    {x:.3,y:2},{x:.2,y:2},{x:.2,y:1.5},{x:.3,y:1.5}]);
  [23,22,12.5,13.5].forEach((cp,i)=>assert.ok(Math.abs(p.elements[0].cp[i].cp-cp)<1e-13));
  assert.equal(p.elements[0].cp.length,4,'No force-only stagnation or modeled base vertices.');
  assert.equal(p.actualMach,.5);assert.equal(p.referenceChord,2);assert.equal(p.physicalAcceptance,false);
  assert.equal(JSON.stringify(f),before);
});

test('pressure frames reject requested-Mach substitution, stale geometry/topology and invalid pressure',()=>{
  for(const mutate of [f=>{f.mach=.74;},f=>{f.bl.stations[0].i++;},
    f=>{f.flow.undisplacedNodes[1][2][0].x+=.001;},f=>{f.flow.nodes[1][2][0].x+=.001;},
    f=>{f.flow.cells[1][1][0].interfacePressure.lower=NaN;},f=>{f.bodies[1].element=1;}]){
    const f=fixture();mutate(f);assert.throws(()=>quadCoupledIterationPressure(f));
  }
});

test('changing BL Ue does not replace historical Euler pressure with an isentropic estimate',()=>{
  const f=fixture(),expected=quadCoupledIterationPressure(f);
  f.checkpoint.restart.initialBL[3]=1.4;
  assert.deepEqual(quadCoupledIterationPressure(f),expected);
});
