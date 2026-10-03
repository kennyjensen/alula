// SPDX-License-Identifier: GPL-2.0-or-later
// Reinitialize only a wake's streamwise station correspondence. The new
// same-condition state must solve every Euler/BL row before it can replace
// the accepted source. No tangential-gap constraint is added to those rows.
import { initializeStreamtubeWakeCorrespondence } from '../streamtube-wake-correspondence.js';
import { initializeStreamtubeWakeArcCorrespondence } from './streamtube-wake-arc-correspondence.js';
import { initializeCoupledStreamtubeFromFlow } from './streamtube-coupled-flow-restart.js';
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

// Kept as a re-export for archived reproduction scripts. The app uses the
// same implementation from the production startup library.
export { correctCoupledWakeCoordinates } from '../streamtube-coupled-startup.js';

export function coupledWakeRecoveryPlan(trial, checkpoint, tolerance) {
  const input = checkpoint?.restart?.input, invalid = trial?.mesh?.quality?.invalidCells;
  const residual = trial?.residual;
  // A first-wake fold can stop a target iteration before its equations close.
  // Eligibility only permits reinitializing the accepted SOURCE's coordinates;
  // the failed target is never accepted or used as that source.
  if (!Number.isFinite(tolerance) || tolerance <= 0
    || input?.wakeGeometry !== 'independent-banks' || !invalid?.length || trial.mesh.quality.valid
    || !(Array.isArray(residual) || ArrayBuffer.isView(residual) && !(residual instanceof DataView))
    || !residual.length || !residual.every(Number.isFinite)
    || !trial.families || !['euler','boundaryLayer','edgeMatching'].every(k =>
      Number.isFinite(trial.families[k]) && trial.families[k] >= 0)) return null;
  const nodes = trial.flow?.nodes;
  if (!Array.isArray(nodes) || nodes.length !== input.bodies.length + 1) return null;
  const nx = input.outerLower.length - 1, tubes = nodes.map(group => group[0].length - 1);
  const starts = [0]; for (const count of tubes) starts.push(starts.at(-1) + nx * count);
  const firstWakeCells = new Map();
  input.bodies.forEach((body,b) => {
    firstWakeCells.set(starts[b] + body.trailingIndex * tubes[b] + tubes[b] - 1,b);
    firstWakeCells.set(starts[b+1] + body.trailingIndex * tubes[b+1],b);
  });
  if (!invalid.every(id => firstWakeCells.has(id))) return null;
  return { bodies: [...new Set(invalid.map(id => firstWakeCells.get(id)))], invalidCells: invalid.slice(),
    targetFamilies: { ...trial.families }, residualConverged: residual.every(v => Math.abs(v) <= tolerance),
    reason: 'The finite failed target has nonconvex cells only at first wake banks; reinitialize the accepted source at its unchanged Mach.' };
}

export function recoverCoupledWakeCorrespondence(checkpoint, { tolerance = 1e-10, maxIterations = 12,
  maxBacktracks = 12, onIteration, onMesh, onCheckpoint } = {}) {
  const mach = checkpoint.restart.input.mach;
  const prepared = initializeCoupledStreamtubeFromFlow(mach,checkpoint,{tolerance});
  const { system, initial, value } = prepared, { layout } = system.euler;
  const massFractions = value.outer.allocation.groups.map(group => {
    const total = group.reduce((s,t) => s+t.massFlow,0), fractions = [0];
    for (const tube of group) fractions.push(fractions.at(-1)+tube.massFlow/total);
    fractions[fractions.length-1] = 1; return fractions;
  });
  // Original reference station fractions may restore a crowded mean wake.
  // Old checkpoints without a reference path keep their original behavior;
  // malformed explicitly supplied paths must fail rather than be ignored.
  const referenceCutPaths = checkpoint.restart.input.cutPaths;
  const initialized = referenceCutPaths === undefined
    ? initializeStreamtubeWakeCorrespondence({ nodes:value.outer.nodes,layout,massFractions })
    : initializeStreamtubeWakeArcCorrespondence({ nodes:value.outer.nodes,layout,massFractions,referenceCutPaths });
  const quality = streamtubeMeshSnapshot({system:system.euler,nodes:initialized.nodes}).quality;
  if (!quality.valid) throw new Error('Wake correspondence initializer has invalid cells.');
  const x = initial.slice();
  x.set(system.euler.adoptGeometry(x.subarray(0,system.ne),initialized.nodes));
  const adopted = system.evaluate(x), f = checkpoint.restart;
  const initialEuler = { x:x.slice(0,system.ne), nodes:adopted.outer.nodes, undisplacedNodes:adopted.outer.undisplacedNodes };
  const candidate = createCoupledStreamtubeBody(f.input,{...f.options,initialEuler,initialBL:f.initialBL});
  if (!candidate.admissible(candidate.initial)) throw new Error('Wake correspondence initializer fails physical admissibility.');
  const replay = candidate.evaluate(candidate.initial);
  const resume = structuredClone({...checkpoint,families:replay.families,restart:{...f,initialEuler}});
  const result = solveCoupledStreamtubeIses(undefined,{resume,tolerance,maxIterations,maxBacktracks,
    iterationGeometry:checkpoint.continuation.iterationGeometry,stepAcceptance:checkpoint.continuation.stepAcceptance,
    stagnationLimiter:checkpoint.continuation.stagnationLimiter,onIteration,onMesh,onCheckpoint});
  return { result, diagnostics: { method:'one-time-wake-station-reinitialization', mach,
    geometry:initialized.diagnostics, initialQuality:quality, initialFamilies:replay.families,
    converged:result.converged && result.mesh.quality.valid, finalFamilies:result.families,
    finalQuality:result.mesh.quality, iterations:result.history.length-1,
    equationsChanged:false, physicalConditionsChanged:false, tangentialGapConstrainedDuringNewton:false,
    sourceRelationship:'An explicit coordinate initialization; not a recovered ISES/MSES wake redistribution formula.' } };
}
