import test from 'node:test';
import assert from 'node:assert/strict';
import { createOrthogonalBoundaryControl, createOrthogonalBoundaryFeedback } from '../src/geometry/orthogonal-boundary-control.js';

const fixture = () => {
  const xi=[0,.2,.4,.6,.8,1],eta=[0,.08,.3,.7,1];
  return {xi,eta,nodes:xi.map(u=>eta.map(e=>({x:u+.08*Math.sin(Math.PI*u)*e*(1-e),y:e}))),
    background:xi.map(u=>eta.map(e=>.1+u+.2*e)),decay:{lower:3,upper:2}};
};

for(const sourceForm of ['poisson','metric-stretch'])test(`inactive ${sourceForm} boundary stations retain background and have zero feedback sensitivity`,()=>{
  const args=fixture(),activeStations={lower:args.xi.map((_,i)=>i===2),upper:args.xi.map(()=>false)};
  const feedback=createOrthogonalBoundaryFeedback({...args,sourceForm,activeStations});
  const field=sourceForm==='poisson'?'poisson':'stretch';
  const reference=createOrthogonalBoundaryFeedback({...args,sourceForm,sides:['lower']});
  for(let i=1;i<args.xi.length-1;i++)for(let j=1;j<args.eta.length-1;j++){
    const actual=feedback.evaluate(args.nodes,i,j);
    if(i===2)assert.deepEqual(actual,reference.evaluate(args.nodes,i,j));
    else {assert.equal(actual[field],args.background[i][j]);assert.deepEqual(actual.derivative,{x:0,y:0});}
  }
  const before=feedback.evaluate(args.nodes,2,1);activeStations.lower.fill(false);
  assert.deepEqual(feedback.evaluate(args.nodes,2,1),before);
});

test('inactive corner and cut nodes impose no normal-speed branch condition',()=>{
  const args=fixture(),activeStations={lower:args.xi.map((_,i)=>i===3),upper:args.xi.map(()=>false)};
  // The control component must not reconstruct a normal on this inactive
  // node; global cell validity is separately enforced by the grid solver.
  args.nodes[1][1].y=-1;
  const result=createOrthogonalBoundaryControl({...args,activeStations,corners:{lower:[2]}});
  assert.equal(result.lower[1],null);assert.equal(result.lower[2],null);assert.ok(result.lower[3]);
  assert.ok(result.upper.every(c=>c===null));
});

test('station support rejects malformed masks and incomplete active-corner neighborhoods',()=>{
  const args=fixture();
  for(const activeStations of [null,[],{wrong:[]},{lower:[true]},{lower:Array(6)},{lower:Array(6).fill(1)}])
    assert.throws(()=>createOrthogonalBoundaryFeedback({...args,activeStations}),/Active boundary/);
  assert.throws(()=>createOrthogonalBoundaryFeedback({...args,corners:{lower:[2]},
    activeStations:{lower:args.xi.map((_,i)=>i===2)}}),/both neighboring/);
});
