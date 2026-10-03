// SPDX-License-Identifier: GPL-2.0-or-later
// A retained grid level is useful evidence, but is not the requested solve.
// These guards only label existing states; they never rebuild or evaluate one.
const count = value => Number.isInteger(value) && value > 0;

export function quadCoupledGrid(value, requested) {
  const sequence = value?.gridSequence
    ?? (value?.automaticRefinement?.kind === 'coarse-to-fine' ? value.automaticRefinement : undefined)
    ?? (count(value?.gridLevel) && count(value?.requestedGridIntervals)
      ? { actualGridIntervals: value.gridLevel, requestedGridIntervals: value.requestedGridIntervals } : undefined);
  if (!sequence) return {};
  const actualGridIntervals = sequence.actualGridIntervals;
  const requestedGridIntervals = requested?.gridIntervals ?? sequence.requestedGridIntervals;
  const differentGrid = count(actualGridIntervals) && count(requestedGridIntervals)
    && actualGridIntervals !== requestedGridIntervals;
  return { ...(count(actualGridIntervals) ? { actualGridIntervals } : {}),
    ...(count(requestedGridIntervals) ? { requestedGridIntervals } : {}),
    differentGrid, targetNotReached: sequence.reachedTarget === false || differentGrid };
}

export function quadCoupledGridLabel(value) {
  const grid = quadCoupledGrid(value);
  return grid.differentGrid
    ? ` · Grid ${grid.actualGridIntervals} → target ${grid.requestedGridIntervals} intervals/side` : '';
}

export function quadCoupledGridResult(raw, requested) {
  const grid = quadCoupledGrid(raw, requested);
  if (!grid.targetNotReached) return raw;
  return { ...raw, converged: false,
    stateConverged: raw.stateConverged ?? (raw.converged === true && raw.mesh?.quality?.valid === true),
    requestedCase: structuredClone(requested),
    ...(grid.differentGrid ? { status: 'research-coupled-target-grid-not-reached',
      sourceCase: { ...(raw.sourceCase ?? requested), gridIntervals: grid.actualGridIntervals } } : {}),
    ...(raw.coefficientStatus === undefined || raw.coefficientStatus === 'unavailable'
      ? {} : { coefficientStatus: 'unconverged' }) };
}
