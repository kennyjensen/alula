// SPDX-License-Identifier: GPL-2.0-or-later
// Iterate Sorenson-style damped boundary controls and an implicit Giles
// line sweep. The damping is an iteration device; convergence is measured
// against the undamped Q=0 equations and the undamped boundary targets.
import { createEllipticStreamtubeGrid } from '../elliptic-streamtube-grid.js';
import { createOrthogonalBoundaryControl } from '../orthogonal-boundary-control.js';

export function smoothOrthogonalBoundaryGrid(options, { maxSweeps = 600, tolerance = 1e-10, controlTolerance = 1e-9,
  omega = 1.3, relaxation = .3, changeLimit = 1, onSweep } = {}) {
  const config = options?.orthogonalBoundaryControl;
  if (!config || config.previous !== undefined || config.relaxation !== undefined || config.changeLimit !== undefined
    || !Number.isInteger(maxSweeps) || maxSweeps < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(controlTolerance) || controlTolerance <= 0 || !Number.isFinite(omega) || !(omega > 0 && omega < 2)
    || !Number.isFinite(relaxation) || !(relaxation > 0 && relaxation <= 1) || !Number.isFinite(changeLimit) || changeLimit <= 0)
    throw new Error('Invalid damped boundary-angle SLOR inputs.');
  const rawSystem = createEllipticStreamtubeGrid(options), sides = config.sides ?? ['lower', 'upper'];
  const field = config.sourceForm === 'poisson' ? 'poisson' : 'stretch';
  const boundaryOptions = { xi: rawSystem.xi, eta: rawSystem.eta, sides, corners: config.corners ?? {}, discretization: 'giles-1985', sourceForm: config.sourceForm };
  let nodes = structuredClone(rawSystem.initial), previous = Object.fromEntries(sides.map(side => [side,
    config.background.map(row => row[side === 'lower' ? 0 : rawSystem.nt])]));
  let converged = false, reason = 'sweep limit', sweeps = 0; const history = [];
  for (let iteration = 0; iteration <= maxSweeps; iteration++) {
    try {
      const target = createOrthogonalBoundaryControl({ ...boundaryOptions, nodes });
      let controlResidual = 0;
      for (const side of sides) for (let i = 1; i < rawSystem.nx; i++) controlResidual = Math.max(controlResidual,
        Math.abs(target[side][i][field] - previous[side][i]) / Math.max(1, Math.abs(target[side][i][field])));
      const { residual, interiorResidual, boundaryResidual } = rawSystem.residuals(nodes), quality = rawSystem.quality(nodes);
      history.push({ iteration, residual, interiorResidual, boundaryResidual, controlResidual,
        minCornerSine: quality.minCornerSine, invalidCells: quality.invalidCells.length });
      onSweep?.(history.at(-1), nodes);
      if (!quality.valid) { reason = 'folded grid'; break; }
      if (residual <= tolerance && controlResidual <= controlTolerance) { converged = true; reason = 'converged'; break; }
      if (iteration === maxSweeps) break;
      const blended = createEllipticStreamtubeGrid({ ...options, orthogonalBoundaryControl: { ...config, previous, relaxation, changeLimit } });
      const next = blended.sweep(nodes, omega); nodes = next.nodes; sweeps++;
      history.at(-1).maxUpdate = next.maxUpdate;
      const committed = createOrthogonalBoundaryControl({ ...boundaryOptions, nodes, previous, relaxation, changeLimit });
      previous = Object.fromEntries(sides.map(side => [side, committed[side].map((row, i) => row ? row[field] : previous[side][i])]));
    } catch (error) { reason = error.message; break; }
  }
  return { nodes, converged, reason, sweeps, history, quality: rawSystem.quality(nodes), boundaryControls: previous,
    coordinateEquations: rawSystem.coordinateEquations, controls: { omega, relaxation, changeLimit, tolerance, controlTolerance, maxSweeps },
    formulation: `Giles secants and implicit block-line SLOR with damped boundary ${field === 'poisson' ? 'P' : 'F'} updates; acceptance uses undamped equations and targets; harmonic mass eta; fixed boundaries` };
}
