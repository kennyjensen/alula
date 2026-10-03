// SPDX-License-Identifier: GPL-2.0-or-later
// Target-gas initial profile only. The original mixed march changes Ue and
// thickness together with phase; it does not solve the coupled Euler field.
import { relaxXfoilProfiles } from '../viscous/xfoil-profile-relaxation.js';

const require = (ok,message) => { if (!ok) throw new Error(message); };
export function prepareCoupledMrchduProfiles({ bl, states, initialBL, targetMach }) {
  require(bl?.transitionMode === 'automatic' && Array.isArray(bl.surfaces) && Array.isArray(bl.wakes)
    && bl.surfaces.length === 2 * bl.wakes.length && bl.wakes.length > 0,
  'MRCHDU Mach prediction requires complete automatic surface and wake profiles.');
  require(Array.isArray(bl.trips) && bl.trips.length === bl.wakes.length
    && bl.trips.every(p => Array.isArray(p) && p.length === 2 && p.every(x => x === 1)),
  'MRCHDU Mach prediction currently supports terminal material trips only.');
  require(Number.isFinite(targetMach) && targetMach > 0 && targetMach < 1 && Number.isFinite(bl.scale) && bl.scale > 0
    && Array.isArray(states) && initialBL?.length === 4 * states.length,
  'Invalid MRCHDU target gas, normalization or complete BL state.');
  require(states.every((p,id) => p.aux === initialBL[4*id] && p.theta === bl.scale*initialBL[4*id+1]
    && p.deltaStar === bl.scale*initialBL[4*id+2] && p.ue === initialBL[4*id+3]),
  'MRCHDU source profiles disagree with their packed physical BL values.');
  const parameters = { ...bl.kernel.parameters, mach: targetMach };
  const predicted = Float64Array.from(initialBL), phases = bl.surfaces.map(s => s.transition), bodies = [], covered = new Set();
  const put = (ids,profile) => profile.forEach((p,k) => {
    const id = ids[k]; require(Number.isInteger(id) && id >= 0 && id < states.length && !covered.has(id),
      'MRCHDU body profiles overlap or have invalid station ids.');
    covered.add(id); predicted.set([p.aux,p.theta/bl.scale,p.deltaStar/bl.scale,p.ue],4*id);
  });
  for (const wake of bl.wakes) {
    const surfaceIndices = ['upper','lower'].map(side => bl.surfaces.findIndex(s => s.body === wake.body && s.side === side));
    require(surfaceIndices.every(i => i >= 0), 'MRCHDU body requires an upper and lower surface.');
    const surfaces = surfaceIndices.map(i => bl.surfaces[i]);
    const input = { surfaces: surfaces.map(s => s.ids.map(id => states[id])), wake: wake.ids.map(id => states[id]),
      phases: surfaces.map(s => s.transition), tripS: surfaces.map(s => states[s.ids.at(-1)].s),
      normalGap: states[wake.ids[0]].wakeGap ?? 0 };
    const result = relaxXfoilProfiles(input,parameters);
    surfaces.forEach((s,k) => { put(s.ids,result.surfaces[k].states); phases[surfaceIndices[k]] = result.surfaces[k].transition; });
    put(wake.ids,result.wake);
    const output = [...result.surfaces.map(s => s.states),result.wake], before = [...input.surfaces,input.wake];
    bodies.push({ body:wake.body, before:input.phases, after:result.surfaces.map(s => s.transition),
      ...(result.inputShapeRecovery ? { inputShapeRecovery: structuredClone(result.inputShapeRecovery) } : {}),
      transitionS:result.surfaces.map(s => s.s), forced:result.surfaces.map(s => s.forced),
      localConvergenceWarnings:result.localConvergenceWarnings.slice(), messages:result.messages.slice(),
      profileChanges:output.map((profile,k) => ({ part:['upper','lower','wake'][k],stations:profile.length,
        maxRelative:Object.fromEntries(['theta','deltaStar','ue'].map(name => [name,
          Math.max(...profile.map((p,i) => Math.abs(p[name]/before[k][i][name]-1)))])) })) });
  }
  require(covered.size === states.length, 'MRCHDU prediction omitted one or more BL stations.');
  return { initialBL:predicted, transitionState:phases, diagnostics:{ method:'xfoil-mrchdu', parameters,
    sourceTransitionMap:bl.surfaces.map(s => s.transition),targetTransitionMap:phases.slice(),bodies,
    physicalBLPreserved:predicted.every((x,i) => i%4 === 0 || x === initialBL[i]),
    interpretation:'One target-gas mixed current-profile march, using source phases and source wake distances. All coupled equations and the current displaced geometry must still be evaluated and solved. Native local warnings remain explicit.',
    operations:{ translatedMRCHDUCalls:bl.wakes.length,globalNewtonUpdates:0,globalLinearSolves:0 } } };
}
