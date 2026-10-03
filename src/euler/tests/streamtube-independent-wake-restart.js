// SPDX-License-Identifier: GPL-2.0-or-later
// Lift a complete centerline Euler/BL restart to independent wake banks.
// No remeshing, BL marching, density reset, or extra flow iteration occurs.
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { transferStreamtubeGeometry } from '../streamtube-geometry.js';

export function independentWakeBankRestart(input, options) {
  if (input.wakeGeometry === 'independent-banks' || !options?.initialEuler || !options?.initialBL)
    throw new Error('Supply a complete centerline Euler/BL restart for independent wake conversion.');
  const source = createCoupledStreamtubeBody(input, options), before = source.evaluate(source.initial);
  const nextInput = { ...input, wakeGeometry: 'independent-banks' };
  const euler = createStreamtubeBodySystem({ ...nextInput, displacement: source.euler.displacement });
  const x = transferStreamtubeGeometry(source.euler, source.initial.subarray(0, source.ne), euler);
  const decoded = euler.decode(x), initialEuler = { x, ...decoded }, initialBL = source.initial.slice(source.ne);
  const { initialEuler: unusedEuler, initialBL: unusedBL, ...settings } = options;
  const system = createCoupledStreamtubeBody(nextInput, { ...settings, initialEuler, initialBL });
  const after = system.evaluate(system.initial), oldRows = source.euler.layout.rows, newRows = system.euler.layout.rows.filter(r => r.kind !== 'wakeGap');
  if (oldRows.length !== newRows.length || oldRows.some((r, i) => r.kind !== newRows[i].kind)) throw new Error('Incompatible wake conversion equation ordering.');
  const maximum = a => a.reduce((peak, value) => Math.max(peak, Math.abs(value)), 0);
  const retainedEulerResidualChange = maximum(newRows.map((r, i) => after.outer.residual[r.index] - before.outer.residual[oldRows[i].index]));
  const boundaryLayerResidualChange = maximum(Array.from(before.residual.subarray(source.ne), (v, i) => v - after.residual[system.ne + i]));
  const gapResidual = after.outer.diagnostics.residualByFamily.wakeGap;
  const maxGeometryChange = maximum(before.outer.nodes.flatMap((group, g) => group.flatMap((row, i) => row.map((p, j) =>
    Math.hypot(p.x - after.outer.nodes[g][i][j].x, p.y - after.outer.nodes[g][i][j].y)))));
  if (retainedEulerResidualChange > 1e-10 || boundaryLayerResidualChange > 1e-10 || gapResidual > 1e-10
    || maxGeometryChange > 1e-12 * source.euler.conditions.lengthScale) throw new Error('Wake conversion did not preserve the complete physical restart.');
  return { input: nextInput, options: settings, system, initialEuler: { ...after.outer, x: system.initial.slice(0, system.ne) },
    initialBL: system.initial.slice(system.ne), diagnostics: { oldUnknowns: source.n, newUnknowns: system.n,
      addedWakeEquations: system.euler.layout.rowCounts.wakeGap, retainedEulerResidualChange, boundaryLayerResidualChange, gapResidual, maxGeometryChange } };
}
