import test from 'node:test';
import assert from 'node:assert/strict';
import { auditStreamtubeShocks, auditDissipationDamping } from '../src/euler/streamtube-shock-audit.js';
import { coupledMachGrowth, compareShockAudits } from '../src/euler/tests/streamtube-shock-continuation.js';
function flow(shear=0) {
  const a={rho:1,p:1,q:1.2,machSquared:1.44},b={rho:1.3,p:1.5,q:.8,machSquared:.64};
  return {nodes:[Array.from({length:3},(_,i)=>[{x:i,y:0},{x:i+shear,y:1}])],
    sections:[[[a]],[[b]]],cells:[[[{states:[a,b],streamwiseResidual:.12,isentropicResidual:.06}]]],
    hybridCells:[[[{fraction:.2,lossIndicators:[.01,.02]}]]]};
}
test('shock audit independently reports compression, blend, pressure loss and grid shear without mutation',()=>{
  const f=flow(),before=structuredClone(f),a=auditStreamtubeShocks(f);
  assert.deepEqual(f,before); assert.equal(a.candidateCount,1);assert.equal(a.weakMomentumCandidates,1);
  assert.equal(a.maximumSkewDegrees,0);assert.equal(a.supersonicBoundarySections,1);
  assert.equal(a.candidates[0].momentumResidualOverPressure,.12/1.25);
  assert.ok(Math.abs(auditStreamtubeShocks(flow(1)).maximumSkewDegrees-45)<1e-12);
  assert.equal(auditStreamtubeShocks({...f,hybridCells:[[[{fraction:1}]]]}).weakMomentumCandidates,0);
  const expansion=flow();expansion.cells[0][0][0].states.reverse();
  assert.equal(auditStreamtubeShocks(expansion).candidateCount,0);
});
test('small accepted updates can suppress the unchanged sixth-order broadening law',()=>{
  const a=auditDissipationDamping(.99,.5,.001);
  assert.ok(a.dampingSuppressesBroadening);assert.ok(a.acceptedMcrit>.989);assert.ok(a.undampedMcrit<.76);
  assert.equal(auditDissipationDamping(.99,.5,1).dampingSuppressesBroadening,false);
  assert.equal(auditDissipationDamping(.99,1e-7,1).acceptedMcrit,.99);
  assert.throws(()=>auditDissipationDamping(.99,1,2));
});
test('Mach increments grow only after inexpensive stages, and never grow immediately after subdivision',()=>{
  const easy={iterations:3,budget:20};assert.equal(coupledMachGrowth(easy).factor,1.5);
  for(const update of [{failedLargerStep:true},{iterations:16},{backtracks:4},
    {shockAudit:{maximumSkewDegrees:60}},{shockAudit:{weakMomentumCandidates:1}}])
    assert.equal(coupledMachGrowth({...easy,...update}).factor,1);
});

test('compression movement is measured in local cell widths and unmatched candidates are explicit',()=>{
  const a=auditStreamtubeShocks(flow()),b=structuredClone(a);
  b.candidates[0].x+=3;
  const m=compareShockAudits(a,b);
  assert.equal(m.matched,1);assert.equal(m.maximumCellWidths,3);
  assert.equal(coupledMachGrowth({iterations:2,budget:20,shockMovement:m}).factor,1);
  assert.equal(compareShockAudits(undefined,b).maximumCellWidths,null);
});

test('saved transonic RAE reporting audits the real flow without altering its equations or state',async()=>{
  const fs=await import('node:fs');
  const {solveCoupledStreamtubeIses}=await import('../src/euler/streamtube-coupled-ises.js');
  const cp=JSON.parse(fs.readFileSync('docs/temporary-shock-broadening/fine07146875-step7/result.json')).checkpoint;
  const before=JSON.stringify(cp),c=cp.continuation;
  let report;
  const r=solveCoupledStreamtubeIses(undefined,{resume:cp,maxIterations:0,
    iterationGeometry:c.iterationGeometry,stepAcceptance:c.stepAcceptance,stagnationLimiter:c.stagnationLimiter,
    onIteration:h=>{report=h;}});
  assert.equal(JSON.stringify(cp),before);
  assert.ok(report.shockAudit.candidateCount>0);
  assert.ok(report.shockAudit.weakMomentumCandidates>0);
  assert.deepEqual(r.families,cp.families);
  assert.equal(r.checkpoint.restart.input.upwind.mucon,cp.restart.input.upwind.mucon);
  assert.equal(r.checkpoint.restart.input.upwind.mcrit,cp.restart.input.upwind.mcrit);
});
