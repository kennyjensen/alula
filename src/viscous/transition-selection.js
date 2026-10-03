// SPDX-License-Identifier: GPL-2.0-or-later
// Determine a surface's active transition interval from the native e^N
// criterion. This prepares discrete regime changes; it does not solve the
// surface momentum/thickness equations or alter any caller-owned state.
import { checkAutomaticTransition, evaluateTransitionInterval } from './transition-interval.js';
import { initializeMixedTransitionShear } from './transition-shear-initialization.js';

export function selectSurfaceTransition(kernel, states, { tripS = Number.MAX_VALUE, initialAmplification = 0 } = {}) {
  if (!Array.isArray(states) || !states.length || !(tripS > 0) || !Number.isFinite(tripS)
    || !Number.isFinite(initialAmplification) || initialAmplification < 0 || initialAmplification >= kernel.parameters.ncrit)
    throw new Error('Invalid automatic-transition surface controls.');
  // Validate the whole surface, including the currently turbulent tail. A
  // later interval change must not discover an invalid hidden station.
  states.forEach((state, i) => {
    if (!(state.s > 0) || (i && !(state.s > states[i - 1].s))) throw new Error('Transition surface stations must increase downstream.');
    try { kernel.station(state, 'laminar'); }
    catch (error) {
      error.stationIndex = i;
      error.diagnostics = { ...error.diagnostics, surfaceStationIndex: i };
      throw error;
    }
  });
  const amplification = [initialAmplification];
  if (tripS <= states[0].s) return { index: 0, kind: 'forced', s: tripS, amplification };
  for (let j = 1; j < states.length; j++) {
    const upstream = { ...states[j - 1], aux: amplification[j - 1] };
    const result = checkAutomaticTransition(kernel, { upstream, downstream: states[j], tripS });
    amplification[j] = result.amplification;
    if (result.transition) {
      if (!(result.s > upstream.s && result.s <= states[j].s)) throw new Error('Automatic transition root left its selected interval.');
      return { index: j, kind: result.forced ? 'forced' : 'natural', s: result.s, amplification };
    }
    if (!Number.isFinite(result.amplification) || result.amplification >= kernel.parameters.ncrit)
      throw new Error('Unresolved laminar amplification during transition selection.');
  }
  // Distinguish a laminar surface from an artificial turbulent surface at TE.
  // Its conversion to a turbulent wake is a separate trailing-edge condition.
  return { index: null, kind: 'laminar', s: null, amplification };
}

// Prepare all auxiliary conversions before returning anything to commit.
// A station's aux is N on the laminar side and Ctau on the turbulent side.
// The caller commits all surfaces and their active metadata atomically.
export function prepareSurfaceTransition(kernel, states, { previousIndex, reinitializeAmplification = false, ...controls } = {}) {
  if (previousIndex !== null && (!Number.isInteger(previousIndex) || previousIndex < 0 || previousIndex >= states.length))
    throw new Error('Invalid previous transition interval.');
  const target = selectSurfaceTransition(kernel, states, controls);
  const oldIndex = previousIndex ?? states.length, nextIndex = target.index ?? states.length;
  const auxiliary = states.map(s => s.aux), converted = [];
  if (oldIndex === nextIndex && !reinitializeAmplification) return { ...target, changed: false, auxiliary, converted };
  const changedMixedInterval = oldIndex !== nextIndex && nextIndex > 0 && nextIndex < states.length;
  for (let j = 0; j < states.length; j++) {
    const wasTurbulent = j >= oldIndex, turbulent = j >= nextIndex;
    // A crossing or a new surface station distribution needs upstream N
    // consistent with the criterion used to select the interval. This is
    // initial-guess transfer, not a projection during unchanged Newton steps.
    // Thickness and Ue remain global Newton unknowns.
    if (!turbulent) auxiliary[j] = target.amplification[j];
    else if (!wasTurbulent && !(changedMixedInterval && j === nextIndex)) {
      auxiliary[j] = kernel.station({ ...states[j], aux: .03 }, 'turbulent').transitionShear;
      if (!(auxiliary[j] > 0) || !Number.isFinite(auxiliary[j])) throw new Error('Invalid newly turbulent shear.');
    }
    if (wasTurbulent !== turbulent) converted.push({ index: j, from: wasTurbulent ? 'turbulent' : 'laminar',
      to: turbulent ? 'turbulent' : 'laminar', oldAux: states[j].aux, aux: auxiliary[j] });
  }
  if (changedMixedInterval) {
    // Prepare a positive shear candidate for the newly mixed interval.
    // Closing this row alone need not reduce the complete BL residual;
    // a coupled caller can select shear using both affected intervals.
    const initialized = initializeMixedTransitionShear(kernel, {
      upstream: { ...states[nextIndex - 1], aux: auxiliary[nextIndex - 1] },
      downstream: { ...states[nextIndex], aux: nextIndex >= oldIndex ? states[nextIndex].aux : .03 },
      tripS: controls.tripS, expected: target });
    auxiliary[nextIndex] = initialized.aux;
    let conversion = converted.find(c => c.index === nextIndex);
    if (!conversion) {
      conversion = { index: nextIndex, from: 'turbulent', to: 'transition', oldAux: states[nextIndex].aux };
      converted.push(conversion);
    }
    Object.assign(conversion, { aux: initialized.aux, transitionShearInitialization: true,
      shearInitialization: initialized.diagnostics,
      ...(target.kind === 'forced' && target.s === states.at(-1).s && nextIndex === states.length - 1
        ? { terminalShearInitialization: true } : {}) });
  }
  return { ...target, changed: oldIndex !== nextIndex, auxiliary, converted };
}

// Ctau at an interior mixed endpoint affects this interval and the following
// turbulent interval. Compare all six affected equations with the caller's
// shear-row weight; selecting on the mixed shear equation alone can worsen
// downstream transport. This selects an event guess, not a nonlinear root.
function mixedTransitionShearScore(kernel, { upstream, downstream, following, tripS, shearWeight = 1 }) {
  if (!(shearWeight > 0) || !Number.isFinite(shearWeight)) throw new Error('Invalid transition shear residual weight.');
  return aux => {
    try {
      if (!(aux > 0) || !Number.isFinite(aux)) throw new Error('Invalid candidate shear.');
      const endpoint = { ...downstream, aux };
      const mixed = evaluateTransitionInterval(kernel, { upstream, downstream: endpoint, tripS }).residual;
      const tail = kernel.interval({ upstream: endpoint, downstream: following, regime: 'turbulent' }).residual;
      const squaredResidual = [mixed, tail].reduce((sum, row) => sum + row.reduce((s, v, i) => s + (v * (i === 0 ? shearWeight : 1)) ** 2, 0), 0);
      if (!Number.isFinite(squaredResidual)) throw new Error('Nonfinite transition shear residual.');
      return { aux, squaredResidual };
    } catch (error) { return { aux, squaredResidual: Infinity, reason: error.message }; }
  };
}

export function chooseMixedTransitionShear(kernel, { initializedAux, ...intervals }) {
  const score = mixedTransitionShearScore(kernel, intervals), shearWeight = intervals.shearWeight ?? 1;
  const existing = score(intervals.downstream.aux), initialized = score(initializedAux);
  const selected = existing.squaredResidual <= initialized.squaredResidual ? 'existing' : 'initialized';
  const choice = selected === 'existing' ? existing : initialized;
  if (!Number.isFinite(choice.squaredResidual)) throw new Error('No admissible mixed-transition shear candidate.');
  return { aux: choice.aux, diagnostics: { method: 'mixed-and-downstream-residual', selected, shearWeight, existing, initialized } };
}

// An upstream crossing converts N into Ctau, so there is no existing shear
// to preserve at the new mixed endpoint. Closing its shear row alone can
// introduce a large downstream transport error. Minimize the six affected
// rows instead, keeping thickness, edge speed, transition and the tail fixed.
// This bounded scalar search only prepares an event guess; Newton still
// solves every original equation. Retaining the initial candidate guarantees
// that an unsuccessful search cannot worsen the local squared residual.
export function minimizeMixedTransitionShear(kernel, { initializedAux, ...intervals }) {
  if (!(initializedAux > 0) || !Number.isFinite(initializedAux)) throw new Error('Invalid initial transition shear.');
  const score = mixedTransitionShearScore(kernel, intervals);
  const seed = Math.min(initializedAux, .25);
  let evaluations = 0;
  const sample = z => {
    evaluations++;
    const aux = seed * Math.exp(z);
    // The same Ctau ceiling is enforced by the coupled XFOIL update.
    return { z, ...(aux > 0 && aux <= .25 ? score(aux) : { aux, squaredResidual: Infinity }) };
  };
  const initial = sample(0);
  let best = initial, lo = -.5, hi = .5, left = sample(lo), right = sample(hi);
  for (let k = 0; k < 16 && (left.squaredResidual < best.squaredResidual || right.squaredResidual < best.squaredResidual); k++) {
    const lower = left.squaredResidual < right.squaredResidual;
    best = lower ? left : right;
    if (lower) { hi = right.z; right = best; lo -= 2 ** k; left = sample(lo); }
    else { lo = left.z; left = best; hi += 2 ** k; right = sample(hi); }
  }
  for (let k = 0; k < 30; k++) {
    const a = sample(lo + (hi - lo) / 3), b = sample(hi - (hi - lo) / 3);
    if (a.squaredResidual < best.squaredResidual) best = a;
    if (b.squaredResidual < best.squaredResidual) best = b;
    if (a.squaredResidual < b.squaredResidual) hi = b.z; else lo = a.z;
  }
  if (!Number.isFinite(best.squaredResidual)) throw new Error('No admissible mixed-transition shear candidate.');
  return { aux: best.aux, diagnostics: { method: 'mixed-and-following-shear-minimum',
    selected: best.aux === initializedAux ? 'initialized' : 'minimum', initializedAux, initial, minimum: best,
    shearWeight: intervals.shearWeight ?? 1, evaluations, equationsChanged: false, initialGuessOnly: true } };
}
