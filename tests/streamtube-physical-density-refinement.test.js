import test from 'node:test';
import assert from 'node:assert/strict';
import { prolongPhysicalLogDensity as prolong } from '../src/euler/streamtube-density-prolongation.js';
const close = (a,b,tol=3e-14) => assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
const centers = rows => rows.map(row => {let sum=0;return row.map(v=>{const answer=sum+v/2;sum+=v;return answer;});});
function fixture({ q=1, factor=2, logRho=()=>0, childShape=1 }={}) {
  const nx=3,masses=[[1,3],[2,1,5]],subdivisions=[[3,2],[1,3,2]],mc=centers(masses),h0=3.5,gamma=1.4;
  const source={masses,logDensity:Array.from({length:nx},(_,i)=>masses.map((row,g)=>row.map((_,j)=>logRho(i+.5,mc[g][j],g))))};
  source.normalAreas=source.logDensity.map(row=>row.map((r,g)=>r.map((v,j)=>masses[g][j]/(Math.exp(v)*q))));
  const childMasses=masses.map((row,g)=>row.flatMap((v,j)=>Array(subdivisions[g][j]).fill(v/subdivisions[g][j])));
  const cm=centers(childMasses),nodes=Array.from({length:nx*factor+1},(_,i)=>i/factor);
  const target={nodeCoordinates:nodes,subdivisions,masses:childMasses,normalAreas:nodes.slice(1).map((_,i)=>masses.map((row,g)=>childMasses[g].map((m,j)=>{
    const s=Math.max(.5,Math.min(nx-.5,(nodes[i]+nodes[i+1])/2)),p=Math.max(mc[g][0],Math.min(mc[g].at(-1),cm[g][j]));
    return childShape*m/(Math.exp(logRho(s,p,g))*q);
  })))};
  return {source,target,h0,gamma};
}
for(const q of [.5,2]) test(`constant physical state q=${q} preserves the ${q<1?'subsonic':'supersonic'} branch and mass/enthalpy`,()=>{
  const f=fixture({q,logRho:()=>.125}),r=prolong(f);
  assert.ok(r.logDensity.every(row=>row.every(g=>g.every(v=>v===.125))));
  assert.equal(r.diagnostics.interpolatedMachBranchChanges,0);
  close(r.diagnostics.maximumEntropyInterpolationDeparture,0);
  close(r.diagnostics.maximumMach,q/Math.sqrt((f.gamma-1)*(f.h0-q*q/2)));
  for(let g=0;g<f.source.masses.length;g++) close(f.source.masses[g].reduce((a,b)=>a+b)*f.h0,f.target.masses[g].reduce((a,b)=>a+b)*f.h0);
});

test('affine log density is exact on interval-section and unequal physical mass centers',()=>{
  const value=(s,m,g)=>.03125*s+.015625*m+.125*g,f=fixture({logRho:value}),r=prolong(f);
  const old=centers(f.source.masses),next=centers(f.target.masses);
  r.logDensity.forEach((row,i)=>row.forEach((group,g)=>group.forEach((v,j)=>{
    const s=Math.max(.5,Math.min(2.5,(f.target.nodeCoordinates[i]+f.target.nodeCoordinates[i+1])/2));
    const m=Math.max(old[g][0],Math.min(old[g].at(-1),next[g][j]));close(v,value(s,m,g));
  })));
  close(r.diagnostics.maximumEntropyInterpolationDeparture,0);
});

test('odd subdivision retains exactly collocated density unknowns; even subdivision has new section centers',()=>{
  const f=fixture({factor:3,logRho:(s,m)=>.125*s+.0625*m}),r=prolong(f);
  assert.ok(r.diagnostics.collocatedValues>0);assert.equal(r.diagnostics.collocatedValuesExact,true);
  assert.equal(r.targetSectionCoordinates[1],.5);assert.equal(r.targetSectionCoordinates[4],1.5);
  const even=prolong(fixture({factor:2}));assert.equal(even.diagnostics.collocatedValues,0);
  assert.deepEqual(even.targetSectionCoordinates,[.25,.75,1.25,1.75,2.25,2.75]);
});

test('entropy change from a changed target area is measured instead of reset by an isentropic inverse',()=>{
  const f=fixture({q:1,childShape:1.2}),r=prolong(f);
  assert.ok(r.logDensity.every(row=>row.every(group=>group.every(v=>v===0))));
  const expected=Math.abs(Math.log((3.5-.5/(1.2*1.2))/(3.5-.5))/(1.4-1));
  close(r.diagnostics.maximumEntropyInterpolationDeparture,expected);
  assert.equal(r.diagnostics.commonIsentropeInversions,0);
});

test('interpolating a physical density jump is bounded and does not claim exact entropy preservation',()=>{
  const f=fixture(); f.source.logDensity[0]=f.source.logDensity[0].map(row=>row.map(()=>-.2));
  f.source.logDensity[2]=f.source.logDensity[2].map(row=>row.map(()=>.3));
  const r=prolong(f),values=r.logDensity.flat(2);
  assert.ok(values.every(v=>v>=-.2&&v<=.3));assert.ok(values.some(v=>v>-.2&&v<0));
  assert.ok(r.diagnostics.maximumEntropyInterpolationDeparture>0);
});

test('thermal impossibility, lost parent mass and malformed section maps reject',()=>{
  for(const mutate of [
    f=>f.target.normalAreas[0][0][0]*=.01,
    f=>f.target.normalAreas[0][0][0]=-1,
    f=>f.target.masses[0][0]*=1.01,
    f=>f.target.nodeCoordinates[1]=0,
    f=>{ f.target.nodeCoordinates[f.target.nodeCoordinates.length-1]+=.01; },
    f=>f.source.logDensity[0][0][0]=NaN,
    f=>f.target.subdivisions[0][0]=0,
  ]) {const f=fixture();mutate(f);assert.throws(()=>prolong(f));}
});

test('source and target data remain untouched and results are detached',()=>{
  const f=fixture(),before=JSON.stringify(f),r=prolong(f);assert.equal(JSON.stringify(f),before);
  r.logDensity[0][0][0]=99;r.nodeCoordinates[0]=99;assert.equal(JSON.stringify(f),before);
});
