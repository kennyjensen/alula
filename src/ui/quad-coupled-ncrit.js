// SPDX-License-Identifier: GPL-2.0-or-later
// Condition labels for the displayed state. No flow or transition evaluation.
const positive = value => Number.isFinite(value) && value > 0;
const first = (...values) => values.find(positive);

export function quadCoupledNcrit(value, requestedNcrit) {
  const actualNcrit = first(value?.checkpoint?.restart?.options?.ncrit,
    value?.conditions?.ncrit, value?.actualNcrit, value?.ncritContinuation?.actualNcrit);
  const targetNcrit = first(requestedNcrit, value?.targetNcrit, value?.ncritContinuation?.targetNcrit);
  const active = value?.ncritContinuation !== undefined || value?.actualNcrit !== undefined || value?.targetNcrit !== undefined
    || actualNcrit !== undefined && targetNcrit !== undefined && actualNcrit !== targetNcrit;
  if (!active) return {};
  return { ...(actualNcrit === undefined ? {} : { actualNcrit }), ...(targetNcrit === undefined ? {} : { targetNcrit }) };
}

export function quadCoupledNcritLabel(value) {
  const { actualNcrit, targetNcrit } = quadCoupledNcrit(value);
  return actualNcrit === undefined ? '' : ` · Ncrit ${actualNcrit}`
    + (targetNcrit !== undefined && targetNcrit !== actualNcrit ? ` → target ${targetNcrit}` : '');
}

// A converged intermediate root does not satisfy the submitted condition.
// Preserve its physical checkpoint, coefficients and actual-case provenance.
export function quadCoupledNcritResult(raw, requested) {
  const point = quadCoupledNcrit(raw, requested?.ncrit ?? 9);
  if (!Object.keys(point).length) return raw;
  const missed = point.actualNcrit !== point.targetNcrit || raw.ncritContinuation?.reachedTarget === false;
  const sourceCase = { ...(raw.sourceCase ?? requested),
    ...(point.actualNcrit === undefined ? {} : { ncrit: point.actualNcrit }),
    ...((raw.conditions?.mach ?? raw.mach) === undefined ? {} : { mach: raw.conditions?.mach ?? raw.mach }) };
  return { ...raw, ...point, sourceCase,
    ...(missed ? { requestedCase: structuredClone(requested), converged: false,
      stateConverged: raw.stateConverged ?? (raw.converged === true && raw.mesh?.quality?.valid === true),
      status: 'research-coupled-target-not-reached',
      ...(raw.coefficientStatus === undefined || raw.coefficientStatus === 'unavailable' ? {} : { coefficientStatus: 'unconverged' }) } : {}) };
}
