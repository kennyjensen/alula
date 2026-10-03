import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createCoupledStreamtubeBody} from '../src/euler/streamtube-coupled.js';
import {sparseProduct} from '../src/numerics/sparse.js';
import {solveSparseDirect} from '../src/numerics/klu.js';
import {gridReplayDeparture} from '../src/euler/streamtube-coupled-grid-levels.js';

test('retained RAE wake contact retains its exact residual and no longer generates the 2880x thickness direction', t => {
  const cp=JSON.parse(fs.readFileSync(new URL('../docs/solver-reliability/rae-recovery-regression/before-checkpoint.json',import.meta.url)));
  const f=cp.restart, s=createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:f.initialBL});
  const original=s.initial.slice(),value=s.evaluate(original);
  assert.deepEqual(value.families,cp.families);
  const matrix=s.jacobian(original),linear=solveSparseDirect(matrix,value.residual.map(v=>-v));
  let maxWake=0, worstStation;
  for(const station of s.bl.stations.filter(t=>t.kind==='wake')) {
    const k=s.ne+4*station.id,relative=Math.abs(linear.x[k+2]/original[k+2]);
    if(relative>maxWake){maxWake=relative;worstStation=station.id;}
  }
  assert.ok(linear.relativeResidual<1e-10);
  assert.ok(maxWake<40,`Wake relative thickness direction ${maxWake}`);
  // Independently check the assembled contact column from the admissible
  // side; centered differences would cross the nondifferentiable floor.
  const column=s.ne+4*185+2, basis=new Float64Array(s.n); basis[column]=1;
  const analytic=sparseProduct(matrix,basis), h=original[column]*1e-7;
  const samples=[0,1,2,3,4].map(m=>{const x=original.slice();x[column]+=m*h;return s.residual(x);});
  let maximumDerivativeError=0;
  for(let row=0;row<s.n;row++) {
    const fd=(-25*samples[0][row]+48*samples[1][row]-36*samples[2][row]+16*samples[3][row]-3*samples[4][row])/(12*h);
    maximumDerivativeError=Math.max(maximumDerivativeError,Math.abs(fd-analytic[row])/Math.max(1,Math.abs(fd),Math.abs(analytic[row])));
  }
  assert.ok(maximumDerivativeError<5e-5,`Assembled outward derivative error ${maximumDerivativeError}`);
  assert.deepEqual(s.initial,original);
  assert.deepEqual(s.evaluate(original).residual,value.residual);
  t.diagnostic(JSON.stringify({maxWakeRelativeThickness:maxWake,worstStation,maximumDerivativeError,linearRelativeResidual:linear.relativeResidual}));
});

test('converged RAE wake root transfers between Mach conditions despite coordinate roundoff, and rejects edited nodes', async () => {
  const {initializeCoupledStreamtubeFromFlow}=await import('../src/euler/tests/streamtube-coupled-flow-restart.js');
  const cp=JSON.parse(fs.readFileSync(new URL('../docs/solver-reliability/rae-wake-contact/cold-checkpoint.json',import.meta.url)));
  const original=structuredClone(cp),same=initializeCoupledStreamtubeFromFlow(.2,cp);
  const reconstructed = createCoupledStreamtubeBody(cp.restart.input, { ...cp.restart.options,
    initialEuler: cp.restart.initialEuler, initialBL: cp.restart.initialBL });
  const replay = reconstructed.evaluate(reconstructed.initial);
  const coordinateReplay = gridReplayDeparture(replay.outer.nodes, cp.restart.initialEuler.nodes);
  assert.ok(coordinateReplay.maximum > 0 && coordinateReplay.equivalent);
  assert.deepEqual(replay.families, cp.families);
  assert.equal(same.diagnostics.targetConverged,true);
  assert.ok(same.diagnostics.sourceNodeReconstructionDeparture>0 && same.diagnostics.sourceNodeReconstructionDeparture<1e-15);
  assert.deepEqual(same.value.families,cp.families);
  const target=initializeCoupledStreamtubeFromFlow(.201,cp);
  assert.equal(target.diagnostics.physicalDensityPreserved,true);
  assert.equal(target.diagnostics.packedBLPreserved,true);
  assert.deepEqual(target.value.outer.nodes,same.value.outer.nodes);
  assert.deepEqual(cp,original);
  const bad=structuredClone(cp);bad.restart.initialEuler.nodes[0][0][0].x+=1e-8;
  assert.throws(()=>initializeCoupledStreamtubeFromFlow(.201,bad),/physical nodes within roundoff/);
});
