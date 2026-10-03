// SPDX-License-Identifier: GPL-2.0-or-later
// Opt-in starting-profile preparation before a discrete transition event.
// This is a separately tested event boundary: it prepares a warm physical
// profile but does not evaluate or alter the coupled residual.
import { prepareCoupledMrchduProfiles } from './streamtube-coupled-mrchdu-predictor.js';
import { incrementIndependentWakeWidths } from './streamtube-wake-geometry.js';
import { extendWarmBoundaryIncrements } from './streamtube-displacement.js';

export function prepareCoupledTransitionProfileTrial(system, proposedPacked, { decoded, baseNodes = decoded?.nodes } = {}) {
  const { bl, euler, ne } = system;
  const originalPhase = bl.snapshotActive();
  if (bl.transitionMode !== 'automatic' || !bl.trips.every(pair => pair.every(x => x === 1))
    || originalPhase.some(j => !Number.isInteger(j) || j <= 0)
    || !euler.layout.independentWakeBanks || euler.layout.wakeDisplacementMotion !== 'te-center')
    throw new Error('Native event-profile trial requires automatic terminal trips, resolved old transitions and independent TE-center wake banks.');
  const oldBL = proposedPacked.subarray(ne), eulerState = proposedPacked.subarray(0, ne);
  const targets = bl.activeTargets(eulerState, oldBL);
  if (targets.every(t => t.from === t.to)) return { x: proposedPacked, nodes: baseNodes, decoded,
    diagnostics: { method: 'xfoil-mrchdu-event-profile', active: false, initialGuessOnly: true, equationsChanged: false } };
  const originalThicknesses = bl.thicknesses(oldBL), candidate = proposedPacked.slice();
  let diagnostics;
  try {
    const geometry = bl.geometry(eulerState);
    const states = bl.stations.map(({ id }) => ({ s: geometry.coordinates[id].s,
      aux: oldBL[4 * id], theta: bl.scale * oldBL[4 * id + 1], deltaStar: bl.scale * oldBL[4 * id + 2], ue: oldBL[4 * id + 3],
      ...(geometry.coordinates[id].wakeGap === undefined ? {} : { wakeGap: geometry.coordinates[id].wakeGap }) }));
    const prediction = prepareCoupledMrchduProfiles({ bl, states, initialBL: oldBL, targetMach: euler.conditions.mach });
    diagnostics = { method: 'xfoil-mrchdu-event-profile', active: true, initialGuessOnly: true, equationsChanged: false,
      oldPhase: originalPhase, rawTargets: targets.map(({ body, side, from, to, kind }) => ({ body, side, from, to, kind })),
      nativePhase: prediction.transitionState.slice(), prediction: prediction.diagnostics };
    candidate.set(prediction.initialBL, ne);
    const predictedBL = candidate.subarray(ne);
    euler.setDisplacement(bl.thicknesses(predictedBL));
    bl.restoreActive(prediction.transitionState);
    diagnostics.auxiliaryReconciliation = bl.updateActive(predictedBL, eulerState, { reinitializeAmplification: true });
    diagnostics.tightPhase = bl.snapshotActive();
    bl.restoreActive(originalPhase);
    const afterThicknesses = bl.thicknesses(predictedBL);
    euler.setDisplacement(afterThicknesses);
    const after = euler.decode(eulerState);
    const wake = incrementIndependentWakeWidths({ layout: euler.layout, nodes: after.nodes,
      beforeWidths: originalThicknesses.wakes, afterWidths: afterThicknesses.wakes });
    const nodes = extendWarmBoundaryIncrements({ sourceNodes: baseNodes, targetNodes: wake.nodes,
      masses: decoded.allocation.groups.map(group => group.map(tube => tube.massFlow)) });
    diagnostics.wakeWidthIncrement = wake.diagnostics;
    return { x: candidate, nodes, decoded: after, diagnostics };
  } catch (error) {
    if (error && typeof error === 'object') error.eventProfile = diagnostics ?? {
      method: 'xfoil-mrchdu-event-profile', active: true, oldPhase: originalPhase, initialGuessOnly: true, equationsChanged: false };
    throw error;
  } finally {
    bl.restoreActive(originalPhase);
    euler.setDisplacement(originalThicknesses);
  }
}
