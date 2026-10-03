// SPDX-License-Identifier: GPL-2.0-or-later
// Wall, cut and end nodes are prescribed. The default also fixes farfields;
// explicit modes implement Giles's x-copy rule with fixed indexed y, with
// either a horizontal-only restriction or NORLIN's nonconstant y values.
// This prepares geometry; compressible gas and Euler/BL rows need separate
// initialization and verification. It does not mutate/adopt the input chart.
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid, gilesStreamwiseCoordinates } from '../geometry/elliptic-streamtube-grid.js';
import { createBoundaryStretchControl } from '../geometry/boundary-stretch-control.js';
import { smoothPairedBoundaryGrid } from '../geometry/paired-boundary-slor.js';
import { createPartialSlorSeedRecorder, isRecoverableSlorTermination } from '../geometry/partial-slor-seed.js';
import { slorNumericalError, slorObserver } from '../geometry/slor-termination.js';

export function relaxStreamtubeInitialGrid({ system, initial = system.initial, nodes = system.decode(initial).nodes, guideField }, {
  seed = 'linear', discretization = 'giles-1985', boundaryControl = 'none', farfieldBoundary = 'fixed', farfieldCurves, streamwiseCoordinates, onSweep, allowPartialInitialGuess = false, ...controls } = {}) {
  if (typeof allowPartialInitialGuess !== 'boolean') throw new Error('Invalid partial elliptic initial-guess control.');
  if (!['linear', 'supplied'].includes(seed)) throw new Error('Unknown elliptic grid seed.');
  if (!['fixed', 'giles-vertical', 'giles-indexed-y', 'normal-curve'].includes(farfieldBoundary)) throw new Error('Unknown elliptic farfield boundary condition.');
  if (!['none', 'wall-angle'].includes(boundaryControl) || boundaryControl === 'wall-angle'
    && (discretization !== 'giles-1985' || farfieldBoundary !== 'fixed' || seed !== 'supplied'))
    throw new Error('Wall-angle SLOR requires Giles differences, fixed boundaries and the supplied initial grid.');
  const { allocation } = system.decode(initial);
  if (!Array.isArray(nodes) || nodes.length !== allocation.groups.length) throw new Error('Invalid elliptic fluid regions.');
  if (farfieldBoundary === 'normal-curve' && farfieldCurves === undefined) {
    if (typeof guideField?.velocityAt !== 'function') throw new Error('Normal farfields need the prepared panel velocity field or explicit boundary curves.');
    // A declared immutable C1 approximation of the traced farfield: preserve
    // every actual sampled node, with tangent dy/dx from the same panel field.
    // This does not assert exact between-node panel streamfunction values.
    farfieldCurves = [nodes[0].map(row => row[0]), nodes.at(-1).map(row => row.at(-1))].map(boundary => ({
      points: boundary.map(p => ({ x: p.x, y: p.y })),
      slopes: boundary.map(p => {
        const q = guideField.velocityAt(p);
        if (!Number.isFinite(q?.u) || !Number.isFinite(q?.v) || !(q.u > 0))
          throw new Error('Unresolved streamwise parameterization of a farfield curve.');
        return q.v / q.u;
      }),
    }));
  }
  if (farfieldBoundary === 'normal-curve' && (!Array.isArray(farfieldCurves) || farfieldCurves.length !== 2))
    throw new Error('Supply lower and upper normal farfield curves.');
  // Giles OUTLIN uses one suction body/stagnation outline for both exterior
  // passages. Here the explicitly selected primary body supplies that common
  // coordinate for all passages. This multielement choice is a reconstruction;
  // the later MSET X law is not reproduced. Giles's indexed-y outer-row
  // rule retains each NORLIN station's y, not a fixed physical curve.
  const primaryBody = system.layout.primaryBody;
  const xi = streamwiseCoordinates ?? (discretization === 'giles-1985'
    ? gilesStreamwiseCoordinates(nodes[primaryBody + 1].map(row => row[0])) : undefined);
  const stationCoordinate = { source: streamwiseCoordinates ? 'prescribed station coordinates'
    : xi ? 'primary upper body/stagnation outline' : 'uniform station index',
    primaryBody, sharedAcrossPassages: true, exactModernMsetGauge: false,
    xi: xi ? [...xi] : Array.from({ length: nodes[0].length }, (_, i) => i / (nodes[0].length - 1)) };
  const result = structuredClone(nodes), regions = [], seeds = [], partial = [];
  for (let g = 0; g < nodes.length; g++) {
    const massFlows = allocation.groups[g].map(tube => tube.massFlow);
    const gridControls = { massFlows, streamwiseCoordinates: xi, discretization,
      boundaryConditions: { lower: g === 0 ? farfieldBoundary : 'fixed',
        upper: g === nodes.length - 1 ? farfieldBoundary : 'fixed' },
      ...(farfieldBoundary === 'normal-curve' ? { boundaryCurves: {
        ...(g === 0 ? { lower: farfieldCurves[0] } : {}),
        ...(g === nodes.length - 1 ? { upper: farfieldCurves[1] } : {}),
      } } : {}) };
    const grid = createEllipticStreamtubeGrid({ nodes: nodes[g], ...gridControls });
    const starting = nodes[g].map((row, i) => row.map((p, j) => seed === 'supplied' || !i || i === grid.nx || !j || j === grid.nt ? { ...p }
      : { x: (1 - grid.eta[j]) * row[0].x + grid.eta[j] * row[grid.nt].x,
        y: (1 - grid.eta[j]) * row[0].y + grid.eta[j] * row[grid.nt].y }));
    seeds.push(starting);
    let smoothingSystem = createEllipticStreamtubeGrid({ nodes: starting, ...gridControls });
    let unavailableAngleTarget = null;
    if (boundaryControl === 'wall-angle') {
      const { eta, nx, nt } = smoothingSystem;
      if (nt < 3) throw new Error('Wall-angle SLOR needs at least three streamtubes per region.');
      const activeStations = Object.fromEntries(['lower', 'upper'].map(side => {
        const body = system.layout.bodies[side === 'lower' ? g - 1 : g];
        return [side, Array.from({ length: nx + 1 }, (_, i) => Boolean(i && i < nx
          && (!body || i > body.leadingIndex && i < body.trailingIndex)))];
      }));
      try {
        smoothingSystem = createEllipticStreamtubeGrid({ nodes: starting, ...gridControls,
          lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
          orthogonalBoundaryControl: { sourceForm: 'metric-stretch', activeStations,
            background: createBoundaryStretchControl({ nodes: starting, xi, eta, metric: 'polygon-arc' }).values,
            decay: { lower: .45 / eta[1], upper: .45 / (1 - eta[nt - 1]) } } });
      } catch (error) {
        if (error.code !== 'orthogonal-normal-branch') throw error;
        unavailableAngleTarget = { converged: false, reason: error.message, code: error.code,
          diagnostics: error.diagnostics, stage: 'boundary-control-initialization', history: [] };
      }
    }
    const stage = 'requested-boundary-angle';
    const recorder = allowPartialInitialGuess && boundaryControl === 'wall-angle' && !unavailableAngleTarget
      ? createPartialSlorSeedRecorder({ nodes: starting, stage }) : null;
    const relax = boundaryControl === 'wall-angle' ? smoothPairedBoundaryGrid : smoothEllipticStreamtubeGrid;
    // An additional pointwise angle target can have no resolved positive
    // normal-speed root on a valid coarse grid. The base inverse-Laplace
    // problem still exists; certify that problem with the same boundaries,
    // mass coordinates and original seed, without the extra angle control.
    let relaxed = unavailableAngleTarget
      ? smoothEllipticStreamtubeGrid(smoothingSystem, { ...controls, omega: 1, requireConvex: true,
        onSweep: onSweep === undefined ? undefined : (history, row) => slorObserver(onSweep, { region: g, ...history, harmonicFallback: true }, row) })
      : relax(smoothingSystem, { ...controls,
        onSweep: recorder || onSweep !== undefined ? (history, row) => {
          recorder?.observe({ ...history, stage }, row);
          slorObserver(onSweep, { region: g, ...history }, row);
        } : undefined });
    if (unavailableAngleTarget) relaxed.harmonicFallback = { attempted: true, converged: relaxed.converged,
      reason: unavailableAngleTarget.reason, originalAttempt: unavailableAngleTarget,
      initialGuess: 'original supplied grid', streamwiseControl: 'none; harmonic xi',
      harmonicMass: true, fixedBoundaries: true, pointwiseWallAngleEnforced: false,
      retainedAdditionalPoissonSpacingControl: false, source: 'Giles ELLIP inverse-Laplace equations, P=Q=0' };
    if (!unavailableAngleTarget && boundaryControl === 'wall-angle' && !relaxed.converged
      && relaxed.termination?.origin === 'solver' && relaxed.termination.termination === 'no-admissible-decreasing-line-step') {
      // Near a prescribed wake/solid junction the additional angle
      // feedback can drive a cell toward a fold. Retry the original seed
      // with its existing fixed spacing control and harmonic mass equation.
      // This is ordinary elliptic SLOR, not acceptance of a failed angle solve.
      const angleAttempt = relaxed;
      const spacingSystem = createEllipticStreamtubeGrid({ nodes: starting, ...gridControls,
        streamwiseStretch: createBoundaryStretchControl({ nodes: starting, xi, eta: smoothingSystem.eta,
          metric: 'polygon-arc' }).values });
      relaxed = smoothEllipticStreamtubeGrid(spacingSystem, { ...controls, omega: 1,
        onSweep: (history, row) => {
          if (history.invalidCells) throw slorNumericalError('Fixed-spacing SLOR trial has a nonconvex cell.', 'spacing-nonconvex');
          slorObserver(onSweep, { region: g, ...history, spacingFallback: true }, row);
        } });
      const { nodes: unusedAngleNodes, ...angleReport } = angleAttempt;
      relaxed.spacingFallback = { attempted: true, converged: relaxed.converged,
        reason: angleAttempt.reason, originalAttempt: angleReport,
        initialGuess: 'original supplied grid', streamwiseControl: 'fixed boundary polygon-arc stretch',
        harmonicMass: true, fixedBoundaries: true, pointwiseWallAngleEnforced: false };
      if (!relaxed.converged && (isRecoverableSlorTermination(relaxed.termination)
        || relaxed.termination?.origin === 'solver' && relaxed.termination.termination === 'spacing-nonconvex')) {
        // The additional Poisson spacing target need not admit a convex
        // discrete grid around a strongly turning dividing streamline.
        // Return to Giles ELLIP's inverse-Laplace equations (P=Q=0), with
        // the SAME prescribed boundaries, xi labels, mass eta and seed.
        // This still performs and certifies SLOR; it does not accept an
        // unsmoothed seed or a failed angle/spacing target as converged.
        const spacingAttempt = relaxed, spacingFallback = relaxed.spacingFallback;
        const harmonicSystem = createEllipticStreamtubeGrid({ nodes: starting, ...gridControls });
        relaxed = smoothEllipticStreamtubeGrid(harmonicSystem, { ...controls, omega: 1, requireConvex: true,
          onSweep: onSweep === undefined ? undefined : (history, row) => slorObserver(onSweep, { region: g, ...history, harmonicFallback: true }, row) });
        const { nodes: unusedSpacingNodes, spacingFallback: unusedFallback, ...spacingReport } = spacingAttempt;
        relaxed.spacingFallback = spacingFallback;
        relaxed.harmonicFallback = { attempted: true, converged: relaxed.converged,
          reason: spacingAttempt.reason, originalAttempt: spacingReport,
          initialGuess: 'original supplied grid', streamwiseControl: 'none; harmonic xi',
          harmonicMass: true, fixedBoundaries: true, pointwiseWallAngleEnforced: false,
          retainedAdditionalPoissonSpacingControl: false, source: 'Giles ELLIP inverse-Laplace equations, P=Q=0' };
      }
    }
    partial[g] = !relaxed.converged ? recorder?.afterRequestedFailure({ ...relaxed.termination,
      completedRequestedStages: true, converged: false, reason: relaxed.reason,
      requestedTolerance: controls.tolerance ?? 1e-9 }) : null;
    result[g] = relaxed.nodes;
    const { nodes: unused, ...report } = relaxed; regions.push(report);
  }
  const converged = regions.every(r => r.converged);
  const partialInitialGuess = !converged && regions.every((r, g) => r.converged || partial[g])
    ? { nodes: result.map((row, g) => partial[g]?.nodes ?? row),
      passages: partial.flatMap((r, g) => r ? [g] : []),
      reports: partial.flatMap((r, g) => r ? [{ region: g, ...r.report }] : []) } : null;
  return { nodes: result, seeds, regions, converged,
    ...(partialInitialGuess ? { partialInitialGuess } : {}), seed, boundaryControl, stationCoordinate, farfieldBoundary,
    status: `${farfieldBoundary === 'fixed' ? 'Fixed-boundary' : farfieldBoundary === 'normal-curve'
      ? 'Prescribed-curve normal-farfield' : farfieldBoundary === 'giles-indexed-y'
      ? 'Giles indexed-y farfield' : 'Giles horizontal-farfield'} elliptic initial grid; Euler/BL equations not solved.` };
}
