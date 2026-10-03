// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { shockRefinementRequest } from '../src/euler/tests/streamtube-shock-refinement.js';
import { solveCoupledStreamtubeAutomatic } from '../src/euler/tests/streamtube-coupled-automatic.js';
import { solveCoupledStreamtubeAlpha } from '../src/euler/tests/streamtube-coupled-alpha.js';
import { createQuadSolveProgress } from '../src/ui/quad-solve-progress.js';
test('real transonic RAE root yields for refinement before another Mach or alpha attempt', () => {
 const {checkpoint}=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-target-grid-continuation/warm-rae64/accepted-0.7205468526482599-198.json'));
 const before=structuredClone(checkpoint);
 const m=solveCoupledStreamtubeAutomatic(.74,{initialCheckpoint:checkpoint,stopAtShock:true});
 assert.equal(m.converged,false);assert.equal(m.stateConverged,true);
 assert.equal(m.continuation.attempts.length,0);assert.equal(m.refinementRequested.trigger,'sonic-compression');
 assert.deepEqual(m.checkpoint,checkpoint);
 const a=solveCoupledStreamtubeAlpha(2.8,{initialCheckpoint:checkpoint,targetMach:.74,stopAtShock:true});
 assert.equal(a.operatingPointContinuation.attempts.length,0);
 assert.equal(a.refinementRequested.trigger,'sonic-compression');
 assert.deepEqual(a.checkpoint,checkpoint);assert.deepEqual(checkpoint,before);
 assert.equal(shockRefinementRequest({...m,converged:false}),null,'unconverged trials cannot request a root transfer');
 assert.equal(shockRefinementRequest({converged:true,mesh:{quality:{valid:true}},flow:{diagnostics:{maxMach:.9}}}),null);
 const p=createQuadSolveProgress({mach:.74,alpha:2.8});
 assert.match(p.update({stage:'coupled-grid-refinement',gridStrategy:'shock-triggered-refinement',gridLevel:64}, {stageChange:true}).current,/shock-triggered refinement.*grid 64/);
});
