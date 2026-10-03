// SPDX-License-Identifier: GPL-2.0-or-later
// Internal alternative objective: retain physical locations of fine
// density clusters while enforcing the original linear spacing constraints.
import { stationGrowthFeasibility, fitConstrainedStations as feasibleSpacing } from './constrained-stations.js';
import { projectPositionKkt, accurateRowResidual } from './position-kkt-projection.js';
import { projectPositionChains } from './position-chain-projection.js';

export function fitConstrainedStationPositions({ positions, firstSpacing, lastSpacing, maximumGrowth = 1.5,
  tolerance = 1e-10, maxSweeps = 20000 }) {
  if (!Array.isArray(positions) || positions.length < 4 || positions[0] !== 0 || positions.at(-1) !== 1
    || positions.some((s, i) => !Number.isFinite(s) || i && !(s > positions[i - 1]))
    || ![firstSpacing, lastSpacing].every(v => Number.isFinite(v) && v > 0) || !(firstSpacing + lastSpacing < 1)
    || !Number.isFinite(tolerance) || tolerance <= 0 || !Number.isInteger(maxSweeps) || maxSweeps < 1)
    throw new Error('Position-constrained stations require finite ordered normalized nodes and positive endpoint lengths.');
  const count = positions.length - 1;
  const bounds = stationGrowthFeasibility({ intervals: count, firstSpacing, lastSpacing, maximumGrowth });
  if (!bounds.feasible) throw new Error('Infeasible endpoint/growth constraints.');
  if (count === 3) {
    const middle = 1 - firstSpacing - lastSpacing;
    if (Math.max(firstSpacing / middle, middle / firstSpacing, lastSpacing / middle, middle / lastSpacing) > maximumGrowth + tolerance)
      throw new Error('Fixed station positions fail physical growth.');
    return { positions: [0, firstSpacing, 1 - lastSpacing, 1], sweeps: 0, primalResidual: 0, stationarityResidual: 0,
      complementarity: 0, totalResidual: 0, objective: 0, objectiveScale: 1, ...bounds, maximumGrowth,
      rawStationarityResidual: 0, rawComplementarity: 0, rawStationarityBackwardError: 0, rawComplementarityBackwardError: 0,
      certificate: 'unique fixed station positions; physical growth',
      forwardCertificate: { physicalPositionErrorBound: 4 * Number.EPSILON, inverseDefectBound: 0,
        maximumInactiveViolationBound: 0, activeConstraintCount: 0, minimumActiveMultiplierBound: null },
      method: 'fixed unique station positions' };
  }
  const natural = positions.slice(1).map((v, i) => v - positions[i]);
  const reference = [...positions]; reference[1] = firstSpacing; reference[count - 1] = 1 - lastSpacing;
  const scales = Array.from({ length: count - 3 }, (_, k) => .5 * (natural[k + 1] + natural[k + 2]));
  const objectiveScale = Math.min(...scales) ** 2;
  if (!(objectiveScale > 0 && Number.isFinite(objectiveScale)) || scales.some(v => !(v > 0 && Number.isFinite(v))))
    throw new Error('Unrepresentable position objective weights.');
  const feasible = feasibleSpacing({ positions, firstSpacing, lastSpacing, maximumGrowth, tolerance, method: 'primal-dual' });
  const x = scales.map((scale, k) => (feasible.positions[k + 2] - reference[k + 2]) / scale), rows = [];
  const add = (terms, intervalTerms) => {
    let rhs = 0; const free = [];
    for (const [i, a] of terms) {
      rhs -= a * reference[i];
      if (i >= 2 && i <= count - 2) free.push([i - 2, a * scales[i - 2]]);
    }
    const norm = Math.hypot(...free.map(t => t[1]));
    if (!norm) { if (rhs < -tolerance) throw new Error('Infeasible fixed position row.'); return; }
    let intervalRhs = 0; const intervals = [];
    for (const [i, a] of intervalTerms) {
      if (i === 0) intervalRhs -= a * firstSpacing;
      else if (i === count - 1) intervalRhs -= a * lastSpacing;
      else intervals.push([i - 1, a]);
    }
    intervals.sort((a, b) => a[0] - b[0]);
    rows.push({ terms: free.map(([i, a]) => [i, a / norm]), rhs: rhs / norm, correction: 0, norm, intervals, intervalRhs });
  };
  for (let i = 1; i < count; i++) {
    add([[i - 1, maximumGrowth], [i, -1 - maximumGrowth], [i + 1, 1]], [[i, 1], [i - 1, -maximumGrowth]]);
    add([[i - 1, -1], [i, 1 + maximumGrowth], [i + 1, -maximumGrowth]], [[i - 1, 1], [i, -maximumGrowth]]);
  }
  const dot = (row, value) => row.terms.reduce((s, [i, a]) => s + a * value[i], 0);
  const project = active => {
    try { return projectPositionChains({ active, rows, reference, scales, natural }); }
    catch { return projectPositionKkt({ active, rows, scales }); }
  };
  let sweeps = 0, change = Infinity, primalResidual = Infinity, complete = false, forwardCertificate;
  const active = [], limits = Math.min(maxSweeps, 10 * (rows.length + 1));
  for (; sweeps < limits; sweeps++) {
    let { target, multipliers } = project(active);
    if (![...target, ...multipliers].every(Number.isFinite)) throw new Error('Nonfinite position projection.');
    primalResidual = Math.max(0, ...rows.map(row => dot(row, target) - row.rhs));
    if (primalResidual <= tolerance && multipliers.every(v => v >= 0)) {
      ({ target, multipliers, forwardCertificate } = projectPositionKkt({ active, rows, scales, certify: true }));
      if (![...target, ...multipliers].every(Number.isFinite)) throw new Error('Nonfinite position KKT projection.');
      primalResidual = Math.max(0, ...rows.map(row => dot(row, target) - row.rhs));
    }
    const direction = target.map((v, i) => v - x[i]);
    change = Math.max(0, ...direction.map(Math.abs));
    if (primalResidual <= tolerance && multipliers.every(v => v >= 0)) {
      target.forEach((v, i) => { x[i] = v; }); active.forEach((k, j) => { rows[k].correction = multipliers[j]; });
      complete = true; sweeps++; break;
    }
    if (change <= tolerance) {
      let drop = -1;
      for (let j = 0; j < active.length; j++) if (multipliers[j] < 0 && (drop < 0 || multipliers[j] < multipliers[drop])) drop = j;
      if (drop < 0) throw new Error('Position active set stopped without a certificate.');
      active.splice(drop, 1); continue;
    }
    let fraction = 1, blocking = -1;
    rows.forEach((row, k) => {
      if (active.includes(k)) return;
      const speed = dot(row, direction);
      if (!(speed > 0)) return;
      const candidate = Math.max(0, row.rhs - dot(row, x)) / speed;
      if (candidate < fraction) { fraction = candidate; blocking = k; }
    });
    x.forEach((_, i) => { x[i] += fraction * direction[i]; });
    if (blocking >= 0) active.push(blocking);
  }
  if (!complete) throw new Error(`Position projection did not converge in ${limits} active-set iterations (${change}, ${primalResidual}).`);
  const output = reference.map((s, i) => i >= 2 && i <= count - 2 ? s + x[i - 2] * scales[i - 2] : s);
  // The fixed target weights are normalized to a maximum of one in physical
  // position coordinates. This positive common factor preserves the minimizer.
  // Primal rows and physical growth/density requirements are not rescaled.
  const stationarity = x.map((v, i) => [[i, objectiveScale]]);
  const packed = [...x, ...rows.map(row => objectiveScale * row.correction)];
  let complementarity = 0;
  for (let k = 0; k < rows.length; k++) {
    const row = rows[k];
    for (const [i, a] of row.terms) stationarity[i].push([x.length + k, a]);
    complementarity = Math.max(complementarity, Math.abs(objectiveScale * row.correction * accurateRowResidual(row.terms, x, row.rhs)));
  }
  const stationarityResidual = Math.max(0, ...stationarity.map(terms => Math.abs(accurateRowResidual(terms, packed))));
  primalResidual = Math.max(0, ...rows.map(row => -accurateRowResidual(row.terms, x, row.rhs)));
  const h = output.slice(1).map((v, i) => v - output[i]);
  const growth = Math.max(...h.slice(1).map((v, i) => Math.max(v / h[i], h[i] / v)));
  const rawPacked = [...x, ...rows.map(row => row.correction)];
  const rawStationarityResidual = Math.max(0, ...stationarity.map(terms => Math.abs(accurateRowResidual(
    terms.map(([i, a], k) => [i, k ? a : 1]), rawPacked))));
  const rawComplementarity = Math.max(0, ...rows.map(row => Math.abs(row.correction * accurateRowResidual(row.terms, x, row.rhs))));
  const rawStationarityBackwardError = Math.max(0, ...stationarity.map(terms => {
    const rawTerms = terms.map(([i, a], k) => [i, k ? a : 1]);
    const scale = Math.max(1, rawTerms.reduce((s, [i, a]) => s + Math.abs(a * rawPacked[i]), 0));
    return Math.abs(accurateRowResidual(rawTerms, rawPacked)) / scale;
  }));
  const rawComplementarityBackwardError = Math.max(0, ...rows.map(row => {
    const scale = Math.max(1, Math.abs(row.correction) * (Math.abs(row.rhs)
      + row.terms.reduce((s, [i, a]) => s + Math.abs(a * x[i]), 0)));
    return Math.abs(row.correction * accurateRowResidual(row.terms, x, row.rhs)) / scale;
  }));
  const reconstructionRoundingBound = Math.max(...scales.map((s, i) => 4 * Number.EPSILON
    * (Math.abs(reference[i + 2]) + Math.abs(s * x[i])) + 4 * Number.MIN_VALUE));
  forwardCertificate = { ...forwardCertificate, reconstructionRoundingBound,
    physicalPositionErrorBound: forwardCertificate.physicalPositionErrorBound + reconstructionRoundingBound };
  if (![primalResidual, stationarityResidual, complementarity, growth, rawStationarityResidual, rawComplementarity].every(Number.isFinite)
    || h.some(v => !(v > 0)) || growth > maximumGrowth + tolerance
    || Math.max(primalResidual, stationarityResidual, complementarity) > tolerance
    || !Number.isFinite(rawStationarityBackwardError) || !Number.isFinite(rawComplementarityBackwardError)
    || Math.max(rawStationarityBackwardError, rawComplementarityBackwardError) > 256 * Number.EPSILON
    || !(forwardCertificate.physicalPositionErrorBound <= tolerance)
    || forwardCertificate.maximumInactiveViolationBound > 0 || forwardCertificate.minimumActiveMultiplierBound < 0)
    throw new Error(`Position KKT or physical growth certificate failed: ${primalResidual}, ${stationarityResidual}, ${complementarity}, growth ${growth}.`);
  return { positions: output, sweeps, primalResidual, stationarityResidual, complementarity, totalResidual: 0,
    forwardCertificate,
    rawStationarityResidual, rawComplementarity, certificate: 'normalized-objective KKT; original primal rows and physical growth',
    rawStationarityBackwardError, rawComplementarityBackwardError,
    objectiveScale, objective: objectiveScale * .5 * x.reduce((s, v) => s + v * v, 0), ...bounds, maximumGrowth,
    method: 'position-space active set with chain endpoint projection', policy: 'minimum squared station displacement weighted by minimum/local target spacing; unchanged endpoint positions and growth constraints' };
}
