// SPDX-License-Identifier: GPL-2.0-or-later
// Inverse finite-element harmonic streamfunction with prescribed transverse
// node guides. Solve K(z)*eta=0 directly; mesh-energy stationarity is different.
// Browser promotion requires separate geometry, field and refinement gates.
import { createTransverseStreamfunctionGrid } from './transverse-streamfunction-grid.js';
import { createTransverseHarmonicResidual } from './transverse-harmonic-residual.js';
import { createTransverseNodeChart } from './transverse-node-chart.js';
import { isPolynomialGridGeometry } from './polynomial-grid-geometry.js';
import { solveStreamwiseScalarLines } from '../numerics/streamwise-scalar-lines.js';
import { sparseProduct } from '../numerics/sparse.js';

const copy = nodes => nodes.map(row => row.map(p => ({ x: p.x, y: p.y })));
const norm = values => {
  let scale = 0;
  for (const v of values) scale = Math.max(scale, Math.abs(v));
  return scale ? scale * Math.sqrt(values.reduce((sum, v) => sum + (v / scale) ** 2, 0)) : 0;
};

export function createTransverseHarmonicGrid(input) {
  if (input.geometryCorrection !== undefined) throw new Error('Curved smoothing requires a certified polynomial geometry, not a raw correction callback.');
  if (input.curvedGeometry !== undefined) {
    const geometry = input.curvedGeometry;
    if (!isPolynomialGridGeometry(geometry)) throw new Error('Curved smoothing requires a polynomial grid geometry with whole-cell certificates.');
    const chart = createTransverseNodeChart(input);
    if (geometry.nx !== chart.nx || geometry.nt !== chart.nt) throw new Error('Curved geometry and nodal chart dimensions differ.');
    for (let i = 0; i <= chart.nx; i++) for (let j = 0; j <= chart.nt; j++)
      if (!i || !j || i === chart.nx || j === chart.nt) {
        const p = chart.initial[i][j], q = geometry.initial[i][j];
        if (p.x !== q.x || p.y !== q.y) throw new Error('Curved smoothing boundaries must match the fixed reference curves.');
      }
    const flow = createTransverseHarmonicResidual({ ...input, geometryCorrection: geometry.correction });
    const quality = nodes => { chart.validate(nodes); return geometry.quality(nodes, chart.directions); };
    if (!quality(chart.initial).valid) throw new Error('Curved initialization requires certified positive cells and transverse guides.');
    return { ...chart, quality, geometryModel: flow.geometryModel, quadratureOrder: flow.quadratureOrder,
      evaluate: (nodes, controls) => {
        if (!quality(nodes).valid) throw new Error('Curved harmonic state lacks a whole-cell geometry certificate.');
        return flow.evaluate(nodes, controls);
      } };
  }
  const chart = createTransverseStreamfunctionGrid(input), flow = createTransverseHarmonicResidual(input);
  return { ...chart, evaluate: (nodes, controls) => { chart.validate(nodes); return flow.evaluate(nodes, controls); }, energy: chart.evaluate,
    geometryModel: flow.geometryModel, quadratureOrder: flow.quadratureOrder };
}

export function smoothTransverseHarmonicGrid(system, { maxIterations = 40, tolerance = 1e-9,
  linearSolve, linearSolverName, onIteration } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || !(tolerance > 0)
    || (linearSolve !== undefined && typeof linearSolve !== 'function')) throw new Error('Invalid transverse harmonic controls.');
  const solve = linearSolve ?? ((matrix, rhs) => solveStreamwiseScalarLines(matrix, rhs, system));
  const linearBackend = linearSolve ? (linearSolverName ?? 'custom certified linear solver') : 'scalar streamwise SLOR';
  let nodes = copy(system.initial), reason = 'iteration limit'; const history = [];
  for (let iteration = 0; iteration <= maxIterations; iteration++) {
    const e = system.evaluate(nodes, { linearize: true });
    const scales = e.rowScale.map(v => v * system.lengthScale);
    if (!scales.every(v => Number.isFinite(v) && v > 0)) { reason = 'singular harmonic Jacobian'; break; }
    const scaled = e.residual.map((v, k) => v / scales[k]);
    const residual = scaled.reduce((m, v) => Math.max(m, Math.abs(v)), 0), merit = .5 * norm(scaled) ** 2;
    // The negative Jacobian has positive pivots at a regular affine flow.
    // This sign change leaves the exact Newton system unchanged. No row
    // scaling or diagonal shifts enter the linear system or its certificate.
    const matrix = { ...e.matrix, values: e.matrix.values.map(v => -v) };
    let linear;
    try { linear = solve(matrix, e.residual); }
    catch (error) { reason = 'linear solve failed'; history.push({ iteration, residual, linearError: error.message }); break; }
    if (linear.converged === false || linear.x?.length !== system.n || !linear.x.every(Number.isFinite)) {
      reason = 'linear solve failed'; history.push({ iteration, residual, linearResidual: linear.relativeResidual }); break;
    }
    const rawNorm = norm(e.residual), product = sparseProduct(e.matrix, linear.x);
    const defect = norm(product.map((v, k) => v + e.residual[k]));
    const linearResidual = rawNorm ? defect / rawNorm : defect === 0 ? 0 : Infinity;
    if (!(linearResidual <= 1e-10)) { reason = 'linear solve failed'; history.push({ iteration, residual, linearResidual }); break; }
    const maxUpdate = linear.x.reduce((m, v) => Math.max(m, Math.abs(v) / system.lengthScale), 0);
    const h = { iteration, residual, rawResidualNorm: rawNorm, merit, maxUpdate, linearResidual, linearSweeps: linear.sweeps };
    history.push(h); onIteration?.(h, nodes);
    if (maxUpdate <= tolerance && residual <= tolerance) { reason = 'harmonic residual converged'; break; }
    if (iteration === maxIterations) break;
    // Freeze the merit scaling for this search. Its directional derivative
    // uses the analytic residual Jacobian, not a mesh-energy surrogate.
    const descent = scaled.reduce((s, v, k) => s + v * product[k] / scales[k], 0);
    if (!(descent < 0)) { reason = 'Newton direction is not descent'; break; }
    let step = 1, next;
    for (; step >= 2 ** -24; step *= .5) {
      next = system.move(nodes, linear.x, step);
      if (!system.quality(next).valid) continue;
      const trial = system.evaluate(next, { linearize: false }).residual;
      const nextMerit = .5 * norm(trial.map((v, k) => v / scales[k])) ** 2;
      if (nextMerit <= merit + 1e-4 * step * descent) break;
    }
    h.step = step;
    if (step < 2 ** -24) { reason = 'line search stalled'; break; }
    nodes = next;
  }
  return { nodes, converged: reason === 'harmonic residual converged', reason, history, quality: system.quality(nodes),
    physicsValidated: false, linearBackend, geometryModel: system.geometryModel, quadratureOrder: system.quadratureOrder,
    formulation: 'Q1 physical-space harmonic streamfunction residual with fixed mass labels, physical boundaries and transverse guide lines; analytic inverse-map Newton Jacobian. Research only.' };
}
