// SPDX-License-Identifier: GPL-2.0-or-later
// Validation-only alternative objective: retain physical locations of fine
// density clusters while enforcing the original linear spacing constraints.
import { stationGrowthFeasibility, fitConstrainedStations as feasibleSpacing } from '../../src/numerics/constrained-stations.js';
import { projectPositionKkt, accurateRowResidual } from './position-kkt-project-prototype.js';
import { projectPositionChains } from './position-chain-project-prototype.js';
export { stationGrowthFeasibility };

export function fitConstrainedStations({ positions, firstSpacing, lastSpacing, maximumGrowth = 1.5,
  tolerance = 1e-10, maxSweeps = 20000 }) {
  const count = positions.length - 1;
  const bounds = stationGrowthFeasibility({ intervals: count, firstSpacing, lastSpacing, maximumGrowth });
  if (!bounds.feasible) throw new Error('Infeasible endpoint/growth constraints.');
  const natural = positions.slice(1).map((v, i) => v - positions[i]);
  const reference = [...positions]; reference[1] = firstSpacing; reference[count - 1] = 1 - lastSpacing;
  const scales = Array.from({ length: count - 3 }, (_, k) => .5 * (natural[k + 1] + natural[k + 2]));
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
  const project = active => projectPositionChains({ active, rows, reference, scales, natural });
  let sweeps = 0, change = Infinity, primalResidual = Infinity, complete = false;
  const active = [], limits = Math.min(maxSweeps, 10 * (rows.length + 1));
  for (; sweeps < limits; sweeps++) {
    let { target, multipliers } = project(active);
    primalResidual = Math.max(0, ...rows.map(row => dot(row, target) - row.rhs));
    if (primalResidual <= tolerance && multipliers.every(v => v >= 0)) {
      ({ target, multipliers } = projectPositionKkt({ active, rows, scales }));
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
  const objectiveScale = Math.min(...scales) ** 2;
  const stationarity = x.map((v, i) => [[i, objectiveScale]]);
  const packed = [...x, ...rows.map(row => objectiveScale * row.correction)];
  let complementarity = 0;
  for (let k = 0; k < rows.length; k++) {
    const row = rows[k];
    for (const [i, a] of row.terms) stationarity[i].push([x.length + k, a]);
    complementarity = Math.max(complementarity, Math.abs(objectiveScale * row.correction * accurateRowResidual(row.terms, x, row.rhs)));
  }
  const stationarityResidual = Math.max(0, ...stationarity.map(terms => Math.abs(accurateRowResidual(terms, packed))));
  if (Math.max(primalResidual, stationarityResidual, complementarity) > tolerance) throw new Error(`Position KKT certificate failed: ${primalResidual}, ${stationarityResidual}, ${complementarity}.`);
  return { positions: output, sweeps, primalResidual, stationarityResidual, complementarity, totalResidual: 0,
    objectiveScale, objective: objectiveScale * .5 * x.reduce((s, v) => s + v * v, 0), ...bounds, maximumGrowth,
    method: 'position-space active set with chain endpoint projection', policy: 'minimum squared displacement of interior stations in units of their local natural spacing; unchanged endpoint positions and growth constraints' };
}
