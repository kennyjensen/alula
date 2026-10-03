import test from 'node:test';import assert from 'node:assert/strict';
import {coupledAssemblyConditions,solveCoupledStreamtubeAssembly} from '../src/euler/streamtube-coupled-assembly.js';
test('reference Reynolds scaling preserves dimensional viscosity and trips follow original elements',()=>{
 const input={elements:[{},{}],referenceChord:2,reynolds:3e6,ncrit:7,materialTrips:[[.1,.2],[.3,.4]]},before=structuredClone(input);
 const a=coupledAssemblyConditions(input,[{element:1},{element:0}],1.5);
 assert.equal(a.options.reynolds,2.25e6);assert.equal(1.5/a.options.reynolds,2/input.reynolds);
 assert.deepEqual(a.options.tripFractions,[[.3,.4],[.1,.2]]);assert.deepEqual(a.elementOrder,[1,0]);
 a.options.tripFractions[0][0]=.8;assert.deepEqual(input,before);
 const b=coupledAssemblyConditions({...input,referenceChord:4},[{element:1},{element:0}],3);
 assert.equal(b.options.reynolds,a.options.reynolds);
});
test('ambiguous trip/topology data and invalid physical controls fail before grid work',()=>{
 const input={elements:[{},{}]};
 for(const data of [{materialTrips:[[.1,.2]]},{materialTrips:[[.1,.2],[.3,1]]},{reynolds:0},{referenceChord:0},{ncrit:NaN}])
  assert.throws(()=>solveCoupledStreamtubeAssembly({...input,...data}));
 assert.throws(()=>coupledAssemblyConditions(input,[{element:0},{element:0}],1));
 assert.throws(()=>coupledAssemblyConditions(input,[{element:0},{element:2}],1));
 assert.throws(()=>solveCoupledStreamtubeAssembly(input,{maxIterations:-1}));
 for(const maxStartupAttempts of [0,3,1.5,NaN])assert.throws(()=>solveCoupledStreamtubeAssembly(input,{maxStartupAttempts}));
});
test('a caller can cancel at a stage boundary without constructing a grid',()=>{
 const stop=new Error('cancelled by caller'),stages=[];
 assert.throws(()=>solveCoupledStreamtubeAssembly({elements:[{}]}, {onStage:s=>{stages.push(s.stage);throw stop;}}),e=>e===stop&&e.stage==='euler');
 assert.deepEqual(stages,['euler']);
});
test('automatic assembly mode maps per-element trips and permits natural transition to the TE',()=>{
 const input={elements:[{},{}],transitionMode:'automatic'},bodies=[{element:1},{element:0}];
 const a=coupledAssemblyConditions(input,bodies,1);
 assert.equal(a.options.transitionMode,'automatic');assert.deepEqual(a.options.tripFractions,[[1,1],[1,1]]);
 const b=coupledAssemblyConditions({...input,materialTrips:[[1,.2],[.3,1]]},bodies,1);
 assert.deepEqual(b.options.tripFractions,[[.3,1],[1,.2]]);
 assert.throws(()=>coupledAssemblyConditions({...input,transitionMode:'unknown'},bodies,1),/Unknown/);
});
