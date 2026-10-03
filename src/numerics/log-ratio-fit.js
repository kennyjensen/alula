// SPDX-License-Identifier: GPL-2.0-or-later
// Least squares of log(value[right]/value[left]) - log(ratio).
// First minimize edge mismatch. The free scale in each connected component
// then minimizes squared log changes from the supplied natural values.
import { solveLinear } from './linear.js';

export function fitLogRatios({ natural, connections, maximumValues }) {
  if (!Array.isArray(natural) || !natural.length || natural.some(v => !Number.isFinite(v) || !(v > 0))
    || !Array.isArray(connections) || connections.some(e => !e || ![e.left, e.right].every(i => Number.isInteger(i) && i >= 0 && i < natural.length)
      || e.left === e.right || !Number.isFinite(e.ratio) || !(e.ratio > 0)))
    throw new Error('Log-ratio fitting requires positive values and valid connections.');
  if (maximumValues !== undefined && (!Array.isArray(maximumValues) || maximumValues.length !== natural.length
    || Array.from(maximumValues).some(v => v !== Infinity && (!Number.isFinite(v) || !(v > 0)))))
    throw new Error('Endpoint upper bounds must be positive and match the fitted values.');
  const neighbors = natural.map(() => []);
  connections.forEach((e, k) => { neighbors[e.left].push(k); neighbors[e.right].push(k); });
  const logs = natural.map(Math.log), visited = new Set(), groups = [], limitedGroups = [];
  for (let first = 0; first < natural.length; first++) {
    if (visited.has(first)) continue;
    const stack = [first], members = [], edges = new Set();
    while (stack.length) {
      const i = stack.pop(); if (visited.has(i)) continue;
      visited.add(i); members.push(i);
      for (const k of neighbors[i]) { edges.add(k); const e = connections[k]; stack.push(e.left === i ? e.right : e.left); }
    }
    members.sort((a, b) => a - b);
    if (members.length > 1) {
      // Fix one temporary gauge and solve the positive-definite reduced
      // graph Laplacian. Parallel edges and cycles retain all constraints.
      const indices = new Map(members.map((v, i) => [v, i - 1])), n = members.length - 1;
      const matrix = new Float64Array(n * n), rhs = new Float64Array(n);
      for (const k of edges) {
        const e = connections[k], a = indices.get(e.left), b = indices.get(e.right), target = Math.log(e.ratio);
        if (a >= 0) { matrix[a * n + a]++; rhs[a] -= target; }
        if (b >= 0) { matrix[b * n + b]++; rhs[b] += target; }
        if (a >= 0 && b >= 0) { matrix[a * n + b]--; matrix[b * n + a]--; }
      }
      const relative = [0, ...solveLinear(matrix, rhs)];
      const shift = members.reduce((sum, i, k) => sum + Math.log(natural[i]) - relative[k], 0) / members.length;
      // A common log shift leaves every edge residual unchanged. Therefore
      // one-sided endpoint limits can retain the unconstrained best ratio
      // fit exactly; only the minimum-change gauge is constrained.
      const permitted = maximumValues === undefined ? shift : Math.min(shift,
        ...members.map((i, k) => Math.log(maximumValues[i]) - relative[k]));
      if (permitted < shift) limitedGroups.push(groups.length);
      members.forEach((i, k) => { logs[i] = relative[k] + permitted; });
    } else if (maximumValues && natural[first] > maximumValues[first]) {
      logs[first] = Math.log(maximumValues[first]); limitedGroups.push(groups.length);
    }
    groups.push(members);
  }
  const values = logs.map(Math.exp), gradient = natural.map(() => 0);
  if (values.some(v => !Number.isFinite(v) || !(v > 0))) throw new Error('Unresolved fitted log-ratio values.');
  const residuals = connections.map(e => {
    const residual = logs[e.right] - logs[e.left] - Math.log(e.ratio);
    gradient[e.left] -= residual; gradient[e.right] += residual;
    return residual;
  });
  return { values, residuals, groups,
    maximumMismatch: Math.exp(Math.max(0, ...residuals.map(Math.abs))),
    stationarityResidual: Math.max(...gradient.map(Math.abs)),
    objective: `unweighted least squares of log-ratio residuals; then minimum squared log change per endpoint variable${maximumValues ? ' subject to endpoint upper bounds' : ''}`,
    ...(maximumValues ? { maximumValues: maximumValues.map(v => v === Infinity ? null : v), limitedGroups } : {}),
    exactEquality: residuals.every(r => Math.abs(r) < 1e-12) };
}
