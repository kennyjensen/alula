// SPDX-License-Identifier: GPL-2.0-or-later
import { viscousObservables } from './observables.js';

// Reuse the normal load calculation on an unfinished iterate. Invalid states
// can still have no usable loads; preserve that distinction from a zero load.
export function provisionalViscousCoefficients(state, conditions) {
  try {
    const o = viscousObservables(typeof state === 'function' ? state() : state, conditions);
    const values = Object.fromEntries(['cl', 'cm', 'cd', 'cdf'].map(k => [k, o[k]]));
    if (!Object.values(values).every(Number.isFinite)) throw new Error('Nonfinite pressure or wake load.');
    return { ...values, coefficientStatus: 'unconverged',
      coefficientConditions: { mach: conditions.mach, referenceChord: conditions.referenceChord, momentReference: conditions.momentReference } };
  } catch (error) {
    return { cl: null, cm: null, cd: null, cdf: null, coefficientStatus: 'unavailable', coefficientReason: error.message,
      ...(error.code?{coefficientCode:error.code}:{}) };
  }
}
