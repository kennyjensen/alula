// SPDX-License-Identifier: GPL-2.0-or-later
// Research nonlinear elimination of N/Ctau only. Euler coordinates, theta,
// delta* and Ue are untouched; their equations remain globally coupled.
import { selectSurfaceTransition } from '../../src/viscous/transition-selection.js';
import { evaluateTransitionInterval } from '../../src/viscous/transition-interval.js';
import { evaluateLeadingTransitionInterval } from '../../src/viscous/leading-transition-interval.js';

// At fixed station geometry, thicknesses and Ue, the native turbulent lag
// equation has the form a-b*C-c*log(C), b>=0,c>0. Its positive root is
// unique. Bracket it, then use safeguarded scalar Newton without clipping
// the solution. A failed local solve returns no partial state to the caller.
function shearRoot(evaluate, guess, tolerance) {
  let x = Number.isFinite(guess) && guess > 0 ? guess : .03, calls = 0;
  const at = x => {
    const v = evaluate(x); calls++;
    const r = 20 * v.residual[0], d = 20 * v.downstream[0][0];
    if (!Number.isFinite(r) || !Number.isFinite(d) || !(d < 0)) throw new Error('Unresolved monotone shear equation.');
    return { r, d };
  };
  let v = at(x);
  if (Math.abs(v.r) <= tolerance) return { x, calls, residual: v.r };
  let lo = x, hi = x, vl = v, vh = v;
  for (let i = 0; (vl.r < 0 || vh.r > 0) && i < 64; i++) {
    if (vl.r < 0) { lo *= .5; vl = at(lo); }
    if (vh.r > 0) { hi *= 2; vh = at(hi); }
  }
  if (!(vl.r >= 0 && vh.r <= 0)) throw new Error('Unable to bracket positive shear root.');
  for (let i = 0; i < 50; i++) {
    const trial = x - v.r / v.d;
    x = trial > lo && trial < hi ? trial : .5 * (lo + hi);
    v = at(x);
    if (Math.abs(v.r) <= tolerance) return { x, calls, residual: v.r };
    if (v.r > 0) { lo = x; vl = v; } else { hi = x; vh = v; }
  }
  throw new Error('Shear root did not reach the requested accuracy.');
}

export function prepareCoupledAuxiliaries(system, state, { tolerance = 2e-13, captureIntervals = false } = {}) {
  if (system.bl.transitionMode !== 'automatic' || state.length !== system.n || !state.every(Number.isFinite)
    || !Number.isFinite(tolerance) || tolerance <= 0) throw new Error('Invalid auxiliary-elimination input.');
  const { bl, ne } = system, kernel = bl.kernel;
  const geometry = bl.geometry(state.subarray(0, ne));
  const x = state.slice(), states = bl.stations.map(({ id }) => ({ s: geometry.coordinates[id].s,
    aux: x[ne + 4 * id], theta: bl.scale * x[ne + 4 * id + 1], deltaStar: bl.scale * x[ne + 4 * id + 2], ue: x[ne + 4 * id + 3] }));
  const intervals = [], summary = { scalarCalls: 0, maximumShearResidual: 0, maximumRelativeShearChange: 0,
    maximumAmplificationChange: 0, surfaces: bl.surfaces.length, wakes: bl.wakes.length };
  const put = (id, aux) => { states[id] = { ...states[id], aux }; x[ne + 4 * id] = aux; };
  const solve = (id, previous, regime, tripS) => {
    const downstream = states[id], upstream = states[previous];
    const input = aux => ({ upstream, downstream: { ...downstream, aux }, regime, tripS });
    const evaluate = aux => regime === 'leading-transition' ? evaluateLeadingTransitionInterval(kernel, input(aux))
      : regime === 'transition' ? evaluateTransitionInterval(kernel, input(aux)) : kernel.interval(input(aux));
    const result = shearRoot(evaluate, downstream.aux, tolerance); put(id, result.x);
    summary.scalarCalls += result.calls;
    summary.maximumShearResidual = Math.max(summary.maximumShearResidual, Math.abs(result.residual));
    summary.maximumRelativeShearChange = Math.max(summary.maximumRelativeShearChange,
      Math.abs(result.x - downstream.aux) / Math.max(.01, Math.abs(downstream.aux)));
    if (captureIntervals) intervals.push({ id, previous, input: input(result.x), beforeAux: downstream.aux, scalar: result });
  };
  for (const [i, surface] of bl.surfaces.entries()) {
    const tripS = geometry.surfaceData[i].tripS;
    const selected = selectSurfaceTransition(kernel, surface.ids.map(id => states[id]), { tripS });
    if (selected.index !== surface.transition) throw new Error('Prepare the transition interval before eliminating its auxiliaries.');
    for (let j = 0; j < surface.transition; j++) {
      const id = surface.ids[j], aux = selected.amplification[j];
      summary.maximumAmplificationChange = Math.max(summary.maximumAmplificationChange, Math.abs(aux - states[id].aux)); put(id, aux);
      if (captureIntervals) intervals.push({ id, previous: j ? surface.ids[j - 1] : id,
        input: { upstream: states[j ? surface.ids[j - 1] : id], downstream: states[id], regime: j ? 'laminar' : 'similarity', tripS } });
    }
    for (let j = surface.transition; j < surface.ids.length; j++) solve(surface.ids[j], j ? surface.ids[j - 1] : surface.ids[j],
      j === 0 ? 'leading-transition' : j === surface.transition ? 'transition' : 'turbulent', tripS);
  }
  for (const wake of bl.wakes) {
    const [upper, lower] = bl.surfaces.filter(s => s.body === wake.body).map(s => states[s.ids.at(-1)]), id = wake.ids[0];
    const theta = upper.theta + lower.theta, aux = (upper.theta * upper.aux + lower.theta * lower.aux) / theta;
    if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Invalid matched wake shear.');
    put(id, aux);
    if (captureIntervals) intervals.push({ id, input: { regime: 'te', upper, lower, downstream: states[id],
      matched: { aux, theta, deltaStar: upper.deltaStar + lower.deltaStar } } });
    for (const id of wake.ids.slice(1)) solve(id, id - 1, 'wake', Number.MAX_VALUE);
  }
  return { x, summary, ...(captureIntervals ? { intervals } : {}) };
}
