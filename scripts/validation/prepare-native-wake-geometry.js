// SPDX-License-Identifier: GPL-2.0-or-later
// Validation-only coordinate initialization. No gas, BL, residual or solve.
import { extendWarmBoundaryIncrements } from '../../src/euler/streamtube-displacement.js';

const require = (ok, message) => { if (!ok) throw new Error(message); };
const point = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const add = (a, b) => ({ x: a.x + b.x, y: a.y + b.y });
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const mean = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
const dot = (a, b) => a.x * b.x + a.y * b.y;
const norm = a => Math.hypot(a.x, a.y);
const unit = a => {
  const length = norm(a);
  require(length > 0 && Number.isFinite(length), 'Degenerate retained wake centerline or bank segment.');
  return { x: a.x / length, y: a.y / length };
};
const center = (nodes, b, i) => mean(nodes[b][i].at(-1), nodes[b + 1][i][0]);
const frame = (nodes, b, i, nx) => {
  const tangent = unit(sub(center(nodes, b, Math.min(nx, i + 1)), center(nodes, b, i - 1)));
  return { tangent, normal: { x: -tangent.y, y: tangent.x } };
};

/**
 * rawNodes/targetNodes: physical [passage][station][crossline point].
 * targetNodes contains the prepared solid walls. Exterior chart roundoff is
 * checked explicitly; output exterior banks always copy raw values exactly.
 * trailingIndices[b], wakeGaps[b][i-TE-1]: TE index and desired physical TOTAL gap.
 * masses[g][j]: positive physical tube mass. Only cumulative fractions are used.
 * solidTrailingCenters[b]: unchanged bare solid TE / finite-base center.
 * lengthScale: physical geometry length per BL-kernel unit, for arc diagnostics.
 * Caller must qualify geometry, gas, phase and complete unchanged equations.
 */
export function prepareNativeWakeGeometry({ rawNodes, targetNodes, trailingIndices,
  wakeGaps, masses, solidTrailingCenters, lengthScale = 1 }) {
  require(Number.isFinite(lengthScale) && lengthScale > 0, 'Positive finite length scale required.');
  require(Array.isArray(rawNodes) && rawNodes.length >= 2 && Array.isArray(targetNodes)
    && targetNodes.length === rawNodes.length, 'Matching body passage groups required.');
  const bodies = rawNodes.length - 1, nx = rawNodes[0]?.length - 1;
  require(Number.isInteger(nx) && nx >= 2 && Array.isArray(masses) && masses.length === bodies + 1,
    'Matching passage stations and masses required.');
  rawNodes.forEach((group, g) => {
    const target = targetNodes[g], mass = masses[g];
    require(Array.isArray(mass) && mass.length > 0 && mass.every(m => Number.isFinite(m) && m > 0)
      && Number.isFinite(mass.reduce((a, b) => a + b, 0)), 'Positive finite physical tube masses required.');
    require(Array.isArray(group) && group.length === nx + 1 && Array.isArray(target) && target.length === group.length
      && group.every((row, i) => Array.isArray(row) && row.length === mass.length + 1 && row.every(point)
        && Array.isArray(target[i]) && target[i].length === row.length && target[i].every(point)),
    'Matching finite physical node arrays required.');
  });
  require(Array.isArray(trailingIndices) && trailingIndices.length === bodies
    && trailingIndices.every(te => Number.isInteger(te) && te >= 1 && te < nx), 'Interior body TE indices required.');
  require(Array.isArray(wakeGaps) && wakeGaps.length === bodies
    && wakeGaps.every((row, b) => Array.isArray(row) && row.length === nx - trailingIndices[b]
      && row.every(g => Number.isFinite(g) && g >= 0)), 'Finite nonnegative physical wake gaps required.');
  require(Array.isArray(solidTrailingCenters) && solidTrailingCenters.length === bodies
    && solidTrailingCenters.every(point), 'Bare solid trailing-edge centers required for BL arc diagnostics.');
  let coordinateMagnitude = lengthScale, exteriorTargetDeparture = 0;
  for (const groups of [rawNodes, targetNodes]) for (const g of groups) for (const row of g) for (const p of row)
    coordinateMagnitude = Math.max(coordinateMagnitude, Math.abs(p.x), Math.abs(p.y));
  const exteriorRoundoffTolerance = 64 * Number.EPSILON * coordinateMagnitude;
  for (let i = 0; i <= nx; i++) for (const [g, j] of [[0, 0], [bodies, masses[bodies].length]])
    exteriorTargetDeparture = Math.max(exteriorTargetDeparture, norm(sub(rawNodes[g][i][j], targetNodes[g][i][j])));
  require(exteriorTargetDeparture <= exteriorRoundoffTolerance, 'Prepared target changed an exterior bank beyond coordinate roundoff.');

  const target = targetNodes.map(g => g.map(row => row.map(p => ({ ...p }))));
  for (let i = 0; i <= nx; i++) {
    target[0][i][0] = { ...rawNodes[0][i][0] };
    target[bodies][i][masses[bodies].length] = { ...rawNodes[bodies][i].at(-1) };
  }
  // Every body reads immutable raw/prepared arrays. Shared passages get both
  // boundary increments before the single cumulative-mass extension below.
  const translations = trailingIndices.map((te, b) => sub(center(targetNodes, b, te), center(rawNodes, b, te)));
  trailingIndices.forEach((te, b) => {
    const translation = translations[b];
    for (let i = te + 1; i <= nx; i++) {
      const lower = rawNodes[b][i].at(-1), upper = rawNodes[b + 1][i][0];
      const { normal } = frame(rawNodes, b, i, nx), separation = sub(upper, lower);
      const change = wakeGaps[b][i - te - 1] - dot(separation, normal);
      const shift = { x: .5 * change * normal.x, y: .5 * change * normal.y };
      // Add opposite normal increments to the retained pair, never discard
      // its existing tangential offset by replacing it with +/- gap*n/2.
      target[b][i][masses[b].length] = sub(add(lower, translation), shift);
      target[b + 1][i][0] = add(add(upper, translation), shift);
    }
  });
  const nodes = extendWarmBoundaryIncrements({ sourceNodes: rawNodes, targetNodes: target, masses });
  const diagnostics = { validationOnly: true, initialGuessOnly: true, geometryQualified: false,
    gasQualified: false, converged: false, physicalAcceptance: false,
    method: 'uniform TE-mean translation and retained tangential-offset normal-gap replacement',
    exteriorBanks: { targetDeparture: exteriorTargetDeparture, coordinateRoundoffTolerance: exteriorRoundoffTolerance,
      outputMatchesRawExactly: true, interpretation: 'Only prepared-input chart roundoff is allowed; no exterior displacement is applied.' },
    operations: { boundaryExtensions: 1, gasEvaluations: 0, boundaryLayerEvaluations: 0,
      jacobians: 0, linearSolves: 0, newtonUpdates: 0 }, bodies: [] };
  trailingIndices.forEach((te, b) => {
    const translation = translations[b], firstOld = center(rawNodes, b, te + 1), firstNew = center(nodes, b, te + 1);
    const oldBareLength = norm(sub(firstOld, solidTrailingCenters[b])), newBareLength = norm(sub(firstNew, solidTrailingCenters[b]));
    require(oldBareLength > 0 && newBareLength > 0 && Number.isFinite(oldBareLength + newBareLength),
      'Collapsed or nonfinite bare-TE BL wake interval.');
    let maxGapError = 0, maxTangentialOffsetError = 0, maxCenterTranslationError = 0, maxCenterIntervalLengthError = 0;
    const stations = [];
    for (let i = te + 1; i <= nx; i++) {
      const oldFrame = frame(rawNodes, b, i, nx), newFrame = frame(nodes, b, i, nx);
      const oldD = sub(rawNodes[b + 1][i][0], rawNodes[b][i].at(-1));
      const newD = sub(nodes[b + 1][i][0], nodes[b][i].at(-1));
      const desired = wakeGaps[b][i - te - 1], actual = dot(newD, newFrame.normal);
      const oldTau = dot(oldD, oldFrame.tangent), newTau = dot(newD, newFrame.tangent);
      const cOld = center(rawNodes, b, i), cNew = center(nodes, b, i);
      const before = norm(sub(cOld, center(rawNodes, b, i - 1))), after = norm(sub(cNew, center(nodes, b, i - 1)));
      require(before > 0 && after > 0 && [desired, actual, oldTau, newTau, before, after].every(Number.isFinite),
        'Degenerate or nonfinite initialized wake geometry.');
      maxGapError = Math.max(maxGapError, Math.abs(actual - desired));
      maxTangentialOffsetError = Math.max(maxTangentialOffsetError, Math.abs(newTau - oldTau));
      maxCenterTranslationError = Math.max(maxCenterTranslationError, norm(sub(sub(cNew, cOld), translation)));
      maxCenterIntervalLengthError = Math.max(maxCenterIntervalLengthError, Math.abs(after - before));
      stations.push({ i, desiredGap: desired, actualGap: actual, gapError: actual - desired,
        oldTangentialOffset: oldTau, newTangentialOffset: newTau });
    }
    const bankOutlet = ['lower', 'upper'].map(side => {
      const g = side === 'lower' ? b : b + 1, j = side === 'lower' ? masses[b].length : 0;
      const oldT = unit(sub(rawNodes[g][nx][j], rawNodes[g][nx - 1][j]));
      const newT = unit(sub(nodes[g][nx][j], nodes[g][nx - 1][j]));
      return { side, directionChangeRadians: Math.atan2(oldT.x * newT.y - oldT.y * newT.x, dot(oldT, newT)),
        endpointTranslation: sub(nodes[g][nx][j], rawNodes[g][nx][j]) };
    });
    diagnostics.bodies.push({ body: b, trailingIndex: te, translation, maxGapError, maxTangentialOffsetError,
      maxCenterTranslationError, maxCenterIntervalLengthError, stations,
      bareTeArc: { oldFirstLength: oldBareLength, newFirstLength: newBareLength,
        physicalShift: newBareLength - oldBareLength, kernelShift: (newBareLength - oldBareLength) / lengthScale,
        interpretation: 'All downstream BL arclengths shift by this amount; the fixed bare TE is the runtime anchor. Finite-base dead-air gaps must be recomputed by the caller.' },
      outlet: { qualified: false, banks: bankOutlet,
        interpretation: 'Bank directions and positions may change. Tangency to the position-dependent farfield flow is unevaluated.' } });
  });
  return { nodes, diagnostics };
}
