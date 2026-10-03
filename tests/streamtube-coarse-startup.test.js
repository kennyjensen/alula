import test from 'node:test';
import assert from 'node:assert/strict';
import { coupledCoarseStartupCase, solveCoupledStreamtubeAssembly } from '../src/euler/streamtube-coupled-assembly.js';
import { planCoupledGridLevels } from '../src/euler/streamtube-coupled-grid-levels.js';
const prepared = (factor=.25) => ({ initialization:{thicknessFactor:factor}, system:{bl:{transitionMode:'automatic',hasFiniteBase:false,trips:[[1,1]]}} });
const input = () => ({elements:[{points:[{x:1,y:0},{x:0,y:0},{x:1,y:0}]}],gridIntervals:32,gridTubes:11,
  mach:.2,alpha:2.68,reynolds:2.7e6,ncrit:4,referenceChord:1,transitionMode:'automatic',materialTrips:[[1,1]],gridEllipticSmoothing:true});
test('coarse startup retains all requested physical controls and returns a separate bounded initial grid',()=>{
 const c=input(),before=structuredClone(c),p=prepared(),source=coupledCoarseStartupCase(c,p,{maxIterations:40});
 assert.deepEqual(source,{...before,gridIntervals:16,gridTubes:7}); assert.deepEqual(c,before);
 source.elements[0].points[0].x=99;source.materialTrips[0][0]=.1;assert.deepEqual(c,before);
 assert.equal(coupledCoarseStartupCase({...c,gridIntervals:16},p,{maxIterations:40}),null);
 assert.deepEqual(coupledCoarseStartupCase(c,prepared(1),{maxIterations:40}),{...before,gridIntervals:16,gridTubes:7});
 assert.equal(coupledCoarseStartupCase(c,p,{maxIterations:0}),null);
 assert.equal(coupledCoarseStartupCase(c,p,{maxIterations:40,enabled:false}),null);
 assert.equal(coupledCoarseStartupCase({...c,gridTubes:5},p,{maxIterations:40}).gridTubes,5);
});
test('fine automatic sharp terminal-trip startup is independent of a valid full-thickness guess',()=>{
 const c=input(),before=structuredClone(c),full=prepared(1);
 full.mesh={quality:{valid:true}};full.families={euler:.01,boundaryLayer:.39,edgeMatching:.1};
 for(const factor of [1,.5,.25]){
  const p={...full,initialization:{thicknessFactor:factor}},saved=structuredClone(p);
  assert.deepEqual(coupledCoarseStartupCase(c,p,{maxIterations:40}),{...before,gridIntervals:16,gridTubes:7});
  assert.deepEqual(p,saved);
 }
 assert.deepEqual(c,before);
});
test('unsupported physical branch and caller limits retain the ordinary startup path',()=>{
 const c=input();
 const finite=prepared();finite.system.bl.hasFiniteBase=true;
 assert.deepEqual(coupledCoarseStartupCase(c,finite,{maxIterations:40}),{...c,gridIntervals:16,gridTubes:7});
 for(const change of [bl=>{bl.transitionMode='fixed-trip'},bl=>{bl.trips=[[.03,1]]}]){
  const p=prepared();change(p.system.bl);assert.equal(coupledCoarseStartupCase(c,p,{maxIterations:40}),null);
 }
 assert.throws(()=>solveCoupledStreamtubeAssembly(c,{coarseStartup:'yes'}),/iteration controls/);
});

test('explicit inlet and wake counts belong to the target grid, not every nested level', () => {
 for (const intervals of [32,64,128]) {
  const c={...input(),gridIntervals:intervals,gridInletIntervals:128,gridOutletIntervals:128},before=structuredClone(c);
  const coarse=coupledCoarseStartupCase(c,prepared(),{maxIterations:40});
  assert.equal(coarse.gridInletIntervals,128/(intervals/16));
  assert.equal(coarse.gridOutletIntervals,128/(intervals/16));
  assert.deepEqual(c,before);
 }
 const c={...input(),gridIntervals:128,gridInletIntervals:4,gridOutletIntervals:256};
 const coarse=coupledCoarseStartupCase(c,prepared(),{maxIterations:40});
 assert.equal(coarse.gridInletIntervals,4);assert.equal(coarse.gridOutletIntervals,32);
 const automatic=coupledCoarseStartupCase(input(),prepared(),{maxIterations:40});
 assert.equal(Object.hasOwn(automatic,'gridInletIntervals'),false);
 assert.equal(Object.hasOwn(automatic,'gridOutletIntervals'),false);
});

test('reported 292-interval source exceeds the budget only when explicit farfield counts are duplicated', () => {
 const c={...input(),gridIntervals:128,gridInletIntervals:128,gridOutletIntervals:128};
 const coarse=coupledCoarseStartupCase(c,prepared(),{maxIterations:40});
 // Dimension-only planning; 292 intervals and [10,10] source tubes are
 // from the diagnostic. No synthetic geometry is passed into a flow solve.
 const shape=(nx,tubes)=>({bodies:[{}],outerLower:Array(nx+1).fill(null),weights:[0,1].map(()=>Array(tubes).fill(1))});
 const args={sourceInput:shape(292,10),input:shape(512,14),sourceGridIntervals:16,requestedGridIntervals:128};
 assert.throws(()=>planCoupledGridLevels(args),/70110 nodes.*50000/);
 const coarseNx=292-128-128+coarse.gridInletIntervals+coarse.gridOutletIntervals;
 const plan=planCoupledGridLevels({...args,sourceInput:shape(coarseNx,10)});
 assert.equal(plan.finalNominalGridIntervals,128);
 assert.equal(plan.maximumNodes,16350);
 assert.ok(plan.maximumNodes<plan.maxNodes);
});
