import test from 'node:test';
import assert from 'node:assert/strict';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { captureStreamtubeInletFractions } from '../src/geometry/streamtube-grid-maintenance.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { nestedRefinementCheckpoint } from '../src/euler/streamtube-nested-checkpoint.js';
import { nestedRefinementCheckpoint as nodeCertificate } from '../scripts/validation/nested-refinement-checkpoint.js';
const serialize=x=>JSON.parse(JSON.stringify(x,(_,v)=>ArrayBuffer.isView(v)?[...v]:v));
const create=f=>createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:f.initialBL});

test('public continuation certificate preserves fluid finite-wake displacement and real zero-update history',t=>{
 const fixture=twoActiveFiniteBaseWakes(),before=fixture.system.evaluate(fixture.x);
 // Explicit manufactured UNCONVERGED continuation data, never an accepted
 // root. The current physical chart and inlet coordinates define its history.
 // An actual zero-update resume verifies every family before certification.
 const sourceInput={...fixture.input,stagnationMotion:'walls-only',normalStencil:'body-stations',geometryDomain:'convex'};
 const options={reynolds:1e6,ncrit:9,edgeMatching:'section-velocity'};
 const source=createCoupledStreamtubeBody(sourceInput,{...options,initialEuler:{x:fixture.x.slice(0,fixture.system.ne),nodes:before.outer.nodes,undisplacedNodes:before.outer.undisplacedNodes},initialBL:fixture.x.slice(fixture.system.ne)});
 const value=source.evaluate(source.initial);assert.ok(source.admissible(source.initial));
 const provisional={version:1,families:value.families,
  restart:{input:sourceInput,options,initialEuler:{x:source.initial.slice(0,source.ne),nodes:value.outer.nodes,undisplacedNodes:value.outer.undisplacedNodes},initialBL:source.initial.slice(source.ne)},
  continuation:{fractions:captureStreamtubeInletFractions(value.outer.nodes,source.euler.layout.bodies),lastRedistributedStagnation:value.outer.stagnation,
   preferredOrdering:'amd',pivotTolerance:.001,iterationGeometry:'convex',stepAcceptance:'admissible',stagnationLimiter:'listing'}};
 const zero=solveCoupledStreamtubeIses(undefined,{resume:serialize(provisional),maxIterations:0,tolerance:1e-10,stepAcceptance:'admissible'});
 assert.equal(zero.converged,false);assert.equal(zero.history.length,1);assert.equal(zero.linearDiagnostics.solves,0);
 assert.ok(zero.checkpoint,zero.reason);
 const checkpoint=serialize(zero.checkpoint),parent=create(checkpoint.restart),original=parent.evaluate(parent.initial);
 const controls={streamwiseFactor:2,normalFactor:2},r=refineCoupledStreamtubeBody(checkpoint.restart.input,parent,controls);
 const config={streamwiseSubdivisions:Array(parent.euler.layout.nx).fill(2),normalSubdivisions:r.diagnostics.normalSubdivisions};
 const restart=serialize({input:r.input,options:r.options,initialEuler:r.initialEuler,initialBL:r.initialBL});
 const frozen=serialize({checkpoint,restart,config}),cert=nestedRefinementCheckpoint(checkpoint,restart,config);
 const mirrored=nodeCertificate(checkpoint,restart,config);
 assert.deepEqual(mirrored,cert,'Independent Node and browser/runtime certificates must agree completely.');
 assert.ok(cert.diagnostics.finiteWakeDisplacement.length>0);
 assert.ok(cert.diagnostics.finiteWakeDisplacement.every(d=>d.fluidError<1e-14));
 const changed=cert.diagnostics.finiteWakeDisplacement.find(d=>Math.abs(d.totalChange)>1e-12);
 assert.ok(changed,'A retained total thickness must change measurably when the physical cubic-gap arc changes.');
 assert.ok(Math.abs(changed.totalChange-changed.gapChange)<1e-14);
 assert.equal(cert.diagnostics.exactSerializedReplay,true);assert.equal(cert.diagnostics.initialSMOVERepeated,false);
 const resumed=solveCoupledStreamtubeIses(undefined,{resume:cert.checkpoint,maxIterations:0,tolerance:1e-10,stepAcceptance:'admissible'});
 assert.equal(resumed.converged,false);assert.equal(resumed.history.length,1);assert.equal(resumed.linearDiagnostics.solves,0);
 assert.equal(resumed.initialRedistribution.resumed,true);assert.deepEqual(resumed.families,cert.checkpoint.families);
 const corrupted=structuredClone(restart),station=r.system.bl.stations.find(s=>s.kind==='wake'&&s.body===changed.body&&s.i===changed.refinedIndex);
 const oldStation=parent.bl.stations.find(s=>s.kind==='wake'&&s.body===changed.body&&s.i===changed.parentIndex);
 // Restore the old total without its updated geometric gap. A total-only
 // comparison would pass this station; the fluid comparison must reject it.
 corrupted.initialBL[4*station.id+2]=original.layers.states[oldStation.id].deltaStar/r.system.bl.scale;
 for(const certify of [nestedRefinementCheckpoint,nodeCertificate])
  assert.throws(()=>certify(checkpoint,corrupted,config),/retained physical BL thickness/);
 const changedModel=structuredClone(restart);changedModel.options.ncrit=10;
 for(const certify of [nestedRefinementCheckpoint,nodeCertificate])
  assert.throws(()=>certify(checkpoint,changedModel,config),/effective coupled physics/);
 assert.deepEqual(serialize({checkpoint,restart,config}),frozen);
 t.diagnostic(JSON.stringify({sourceKind:'explicitly manufactured unconverged checkpoint, verified by zero-update resume',sourceConverged:zero.converged,newtonUpdates:0,globalLinearSolves:0,independentNodeCertificateExact:true,finiteWakeStations:cert.diagnostics.finiteWakeDisplacement.length,changed,
  sourceFamilies:zero.families,childFamilies:cert.checkpoint.families}));
});
