// SPDX-License-Identifier: GPL-2.0-or-later
// Public workbench adapter for the research moving-quadrilateral equations.
// The initialized quad grid (optionally SLOR smoothed) is displayed and solved.
import { createInitialStreamtubeTopology } from '../geometry/streamtube-topology.js';
import { preparePanelStreamtubeTopology } from './streamtube-panel-topology.js';
import { createPanelStreamtubeGrid } from './streamtube-body-initializer.js';
import { sampleStreamtubePanelVelocity, initializeStreamtubePanelVelocity } from './streamtube-panel-topology.js';
import { initializeStreamtubeStartup } from './streamtube-startup.js';
import { recoverStreamtubePassageGrid } from './streamtube-passage-startup.js';
import { recoverStreamtubeStagnationMotion } from './streamtube-stagnation-startup.js';
import { relaxStreamtubeInitialGrid } from './streamtube-elliptic-initializer.js';
import { solveStreamtubeIses } from './streamtube-ises-update.js';
import { refineStreamtubeBody } from './streamtube-refinement.js';
import { streamtubeMeshSnapshot, prepareStreamtubeMesh } from './streamtube-mesh-preview.js';
import { prepareCurvedStreamtubePreview } from './streamtube-curved-geometry.js';
import { streamtubeForceCoefficients } from './streamtube-forces.js';
import { streamtubeSolidPressureForces } from './streamtube-forces.js';
import { streamtubeSurfaceDiagnostics } from './streamtube-body.js';
import { streamtubeEquationControls } from './streamtube-equation-selection.js';
import { streamtubeFlowSnapshot, compareStreamtubeFlow } from './streamtube-flow-preview.js';
import { capturePreparedStreamtubeAssembly, restorePreparedStreamtubeAssembly } from './streamtube-prepared-assembly.js';

export function streamtubeConvergenceDiagnostics(flow, quality) {
  const gridReason = quality.valid ? '' : `Final grid has ${quality.invalidCells.length} non-convex or degenerate cell${quality.invalidCells.length === 1 ? '' : 's'}; minimum corner sine ${quality.minCornerSine.toExponential(3)}.`;
  return { residualConverged: flow.residualConverged, gridValid: quality.valid, solverStopReason: flow.reason,
    reason: quality.valid ? flow.reason : flow.residualConverged ? gridReason : `${flow.reason} ${gridReason}` };
}

// Standard startup selects the validated final-grid smoothing/shock startup
// for transonic ISMOM4 precursors and fine inviscid grids. Resolve once in preparation
// and solve so the geometry, dissipation policy and iteration budget agree.
// Coupled startup retains its own policy; no operating condition is changed.
export function resolveStreamtubeStartup(input) {
  const requested = input.eulerStartup ?? 'standard';
  if (!['standard', 'harmonic-shock'].includes(requested)) throw new Error('Unknown Euler startup method.');
  if (requested === 'harmonic-shock') {
    if (input.eulerIsmom === undefined)
      throw new Error('Harmonic shock startup requires an explicit ISMOM formulation.');
    return requested;
  }
  // The viscous precursor needs the same healthy final-grid shock seed.
  return typeof input.quadBoundaryLayers === 'boolean' && input.eulerIsmom === 4
    && input.mach >= .7 && (input.quadBoundaryLayers || input.gridIntervals >= 64) && input.gridEllipticSmoothing !== false
    ? 'harmonic-shock' : 'standard';
}

// Automatic fine-grid preparation avoids applying a second nose-clustering
// law to already clustered source coordinates. Keep explicit distributions
// and end counts, and retain the established coarser-grid preparation.
export function resolveStreamtubeGridControls(input) {
  const automatic = (input.gridSurfaceSpacing ?? 'automatic') === 'automatic';
  const balanced = input.quadBoundaryLayers === false && automatic && input.gridIntervals >= 128 && resolveStreamtubeStartup(input) === 'harmonic-shock';
  return { surfaceSpacing: automatic ? (balanced ? 'source' : 'supplied') : input.gridSurfaceSpacing,
    normalSpacing: input.gridStagnationAspectRatio !== undefined ? 'stagnation' : balanced ? 'automatic' : 'supplied',
    inletIntervals: input.gridInletIntervals ?? (balanced ? Math.ceil(input.gridIntervals / 2) : undefined),
    outletIntervals: input.gridOutletIntervals ?? (balanced ? Math.ceil(input.gridIntervals / 2) : undefined),
    balanced };
}

// Return a reusable geometry/state description before any gas inversion or
// Euler iteration. The GUI solve and mesh-only preview share this path.
export function prepareStreamtubeAssembly(input, { meshOnly = false, onMesh,
  gridRepair = false, refinementInterpolation = 'streamfunction-quadratic', wallSubdivisions = 4,
  interiorInitialization = 'traced', panelVelocityInitialization = false } = {}) {
  if (typeof panelVelocityInitialization !== 'boolean') throw new Error('Invalid panel velocity initialization control.');
  const startup = resolveStreamtubeStartup(input);
  const gridSelection = resolveStreamtubeGridControls(input);
  const equationControls = streamtubeEquationControls(input.eulerIsmom);
  if (!['linear', 'streamfunction-quadratic'].includes(refinementInterpolation))
    throw new Error('Unknown Euler refinement interpolation.');
  if (!Number.isInteger(wallSubdivisions) || wallSubdivisions < 2 || wallSubdivisions > 4)
    throw new Error('Euler wall refinement requires two to four subdivisions.');
  if (!['traced', 'linear'].includes(interiorInitialization)) throw new Error('Unknown initial interior-grid construction.');
  const mach = input.mach ?? .2;
  const referenceChord = input.referenceChord ?? 1;
  const momentReference = input.momentReference ?? { x: referenceChord / 4, y: 0 };
  if (!(referenceChord > 0) || ![referenceChord, momentReference?.x, momentReference?.y].every(Number.isFinite))
    throw new Error('Invalid reference chord or moment reference.');
  const smoothingMethod = input.gridSmoothingMethod ?? 'elliptic';
  if (!['elliptic', 'curved-harmonic'].includes(smoothingMethod)) throw new Error('Unknown grid smoothing method.');
  const curvedPreview = Boolean(input.gridEllipticSmoothing && smoothingMethod === 'curved-harmonic');
  if (curvedPreview && !meshOnly) throw new Error('Curved harmonic smoothing is an experimental mesh preview. Use Build mesh to inspect it; physical validation for Euler/BL is unfinished.');
  const crosslinePlacement = input.gridCrosslinePlacement ?? 'potential';
  if (!['supplied-x', 'potential'].includes(crosslinePlacement)) throw new Error('Unknown grid cross-line placement.');
  if (!meshOnly && (!Number.isFinite(mach) || mach <= 0 || mach >= 1))
    throw new Error(input.eulerIsmom === undefined
      ? 'Quad Euler requires a freestream Mach number between 0 and 1. All local flow must remain subsonic.'
      : 'Quad Euler requires a freestream Mach number between 0 and 1.');
  const topologyInput = { ...input, mach: meshOnly ? .2 : mach };
  const topologyControls = {
    surfaceIntervals: input.gridIntervals ?? 16, surfaceChordExponent: input.gridChordExponent ?? 0, tubes: input.gridTubes ?? 7,
    upperTubes: input.gridUpperTubes, lowerTubes: input.gridLowerTubes, gapTubes: input.gridGapTubes,
    inletIntervals: gridSelection.inletIntervals, outletIntervals: gridSelection.outletIntervals };
  let topology, panelSolution;
  try { topology = createInitialStreamtubeTopology(topologyInput, topologyControls); }
  catch (error) {
    if (crosslinePlacement !== 'potential' || error.code !== 'STREAMTUBE_NONMONOTONE_BODY' || error.diagnostics?.side !== 'lower') throw error;
    ({ topology, panelSolution } = preparePanelStreamtubeTopology(topologyInput, topologyControls, { crosslinePlacement }));
  }
  // The smoothing checkbox must not select a different surface schedule,
  // cut geometry or farfield layout. Construct the same unsmoothed mesh first.
  const panelInput = { ...topology, flowModel: 'compressible', streamwiseMode: 'isentropic',
    normalStencil: 'body-stations', stagnationMotion: 'walls-only', geometryDomain: 'positive-simple',
    ...equationControls };
  const panelControls = { ...(panelSolution ? { panelSolution, panelBoundaryCondition: 'streamfunction' } : {}),
    surfaceSpacing: gridSelection.surfaceSpacing, curvatureSpacing: input.gridCurvatureSpacing ?? {},
    normalSpacing: gridSelection.normalSpacing,
    stagnationAspectRatio: input.gridStagnationAspectRatio ?? 2.5,
    interiorInitialization, crosslinePlacement,
    // Potential locates the multielement block intersections. Station spacing
    // must use physical distance: potential is quadratic at stagnation, so
    // interpolating it directly creates large cut/wall spacing jumps even
    // when SLOR is disabled. Reconcile both banks before tracing interiors.
    ...(crosslinePlacement === 'potential' ? { surfaceCrosslineMetric: 'arc', cutStationSpacing: 'physical-x',
      outerCrosslineSpread: .85, resolveSurfaceTurning: true,
      surfaceStationPlacement: 'passage-density', passageCountPlanning: true } : {}),
  };
  let panelGrid;
  try { panelGrid = createPanelStreamtubeGrid(panelInput, panelControls); }
  catch (error) {
    if (error.code !== 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE') throw error;
    // The smooth wall can lie inside its straight-panel approximation.
    // Rebuild the complete common reference field with the existing XFOIL-
    // style streamfunction BIE before using nodal gamma as a surface speed.
    // Never mix gamma from a normal-velocity field with a zero-interior-
    // velocity surface convention, or change the wall to bypass the guard.
    panelGrid = createPanelStreamtubeGrid(panelInput, { ...panelControls, panelBoundaryCondition: 'streamfunction' });
    panelGrid.diagnostics.panelInitializationRecovery = { reason: error.message, code: error.code,
      from: 'normal-velocity', to: 'streamfunction', attempts: 2 };
  }
  const base = prepareStreamtubeMesh(panelGrid, { onMesh, gridRepair });
  let prepared = input.gridEllipticSmoothing
    ? curvedPreview ? prepareCurvedStreamtubePreview(base, { onMesh })
      : prepareStreamtubeMesh(base, { onMesh, gridRepair, ellipticSmoothing: {
        boundaryControl: 'wall-angle', maxSweeps: 200, tolerance: 1e-9,
        ...(typeof input.gridEllipticSmoothing === 'object' ? input.gridEllipticSmoothing : {}),
      } }) : base;
  if (curvedPreview) return { status: 'mesh-ready', mesh: prepared.mesh, experimental: true, physicalAcceptance: false };
  if (prepared.diagnostics.gridSmoothing && !prepared.diagnostics.gridSmoothing.converged
    && prepared.diagnostics.gridSmoothing.initialGuessAccepted !== true)
    throw new Error(`Elliptic grid initialization failed: ${prepared.diagnostics.gridSmoothing.reason}. The original unsmoothed mesh is retained; Euler was not started.`);
  if (!prepared.mesh.quality.valid)
    throw new Error(`Initial quadrilateral grid contains ${prepared.mesh.quality.invalidCells.length} folded or degenerate cells. ${gridRepair === false
      ? 'The mesh is shown unchanged; flow was not started.'
      : `Bounded repair failed: ${prepared.diagnostics.gridRepair?.reason ?? 'unknown reason'}`}`);
  // Resolve the wall-adjacent tubes without changing streamwise stations,
  // material banks or captured flow. Do this for Build mesh too, before gas
  // inversion. SLOR has already finished on its original logical grid.
  const normalSubdivisions = prepared.system.layout.tubes.map((n, g) => Array.from({ length: n }, (_, j) =>
    g > 0 && j === 0 || g < prepared.system.layout.elements && j === n - 1 ? wallSubdivisions : 1));
  const refined = refineStreamtubeBody(prepared.input, prepared.system, { initial: prepared.initial,
    streamwiseFactor: 1, normalSubdivisions, initializeFlow: false, normalInterpolation: refinementInterpolation });
  const diagnostics = { ...prepared.diagnostics, gridSelection, gridRefinement: refined.diagnostics,
    parentGridDiagnostics: true };
  prepared = { ...refined, nodes: refined.initialEuler.nodes, diagnostics };
  if (startup === 'harmonic-shock') {
    // Wall-normal subdivision changes the interior stencil. Smooth that final
    // grid, keeping every wall, wake bank, inlet, outlet and farfield node fixed.
    const preview = onMesh ? structuredClone(prepared.nodes) : null;
    const relaxed = relaxStreamtubeInitialGrid(prepared, {
      seed: 'supplied', maxSweeps: 1000, tolerance: 1e-8, omega: 1, requireConvex: true, detectCoordinateRoundoff: true,
      onSweep: onMesh ? (h, region) => {
        preview[h.region] = region;
        if (onMesh && h.iteration % 20 === 0) onMesh(streamtubeMeshSnapshot({ system: prepared.system,
          nodes: preview, diagnostics: { ...prepared.diagnostics,
            gridSmoothing: { attempted: true, method: 'final-grid harmonic', converged: false, current: h } } }), 'smoothing');
      } : undefined,
    });
    const precisionLimited = !relaxed.converged && relaxed.regions.every(r => r.converged || r.coordinatePrecision?.limited);
    if (!relaxed.converged && !precisionLimited) throw Object.assign(new Error('Final-grid harmonic smoothing did not converge; Euler was not started.'), {
      code: 'STREAMTUBE_FINAL_GRID_SMOOTHING',
      diagnostics: { regions: relaxed.regions.map((r, region) => ({ region, converged: r.converged,
        reason: r.reason, termination: r.termination, lastSweep: r.history.at(-1) })) },
    });
    const initial = prepared.system.adoptGeometry(prepared.initial, relaxed.nodes);
    const nodes = prepared.system.decode(initial).nodes;
    for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length; i++)
      for (let j = 0; j < nodes[g][i].length; j++) {
        if (i && i < nodes[g].length - 1 && j && j < nodes[g][i].length - 1) continue;
        const before = prepared.nodes[g][i][j], after = nodes[g][i][j];
        if (before.x !== after.x || before.y !== after.y) throw new Error('Final-grid smoothing changed a prescribed boundary.');
      }
    prepared = { ...prepared, initial, nodes, initialEuler: { x: initial, nodes },
      diagnostics: { ...prepared.diagnostics, finalGridSmoothing: {
        method: 'harmonic', converged: relaxed.converged, initialGuessAccepted: true, precisionLimited, fixedBoundaries: true,
        regions: relaxed.regions.map(r => ({ converged: r.converged,
          ...(r.coordinatePrecision ? { coordinatePrecision: r.coordinatePrecision } : {}),
          iterations: r.history.at(-1)?.iteration, residual: r.history.at(-1)?.residual })),
      } } };
  }
  const initialMesh = streamtubeMeshSnapshot(prepared);
  onMesh?.(initialMesh, 'initial');
  const panelVelocitySeed = panelVelocityInitialization ? sampleStreamtubePanelVelocity(prepared.system, prepared.initial, panelGrid.guideField.velocityAt) : undefined;
  return { ...prepared, ...(panelVelocitySeed ? { panelVelocitySeed } : {}), mesh: initialMesh, referenceChord, momentReference, mach };
}

export function solveStreamtubeAssembly(input, { meshOnly = false, onMesh, onIteration,
  maxIterations = resolveStreamtubeStartup(input) === 'harmonic-shock' ? 80 : 40, tolerance = 1e-10, preparedEuler, onEulerPrepared, retainBestCheckpoint = false, stepAcceptance,
  adaptiveMcrit = input.eulerIsmom !== undefined,
  firstOrderStartup = input.eulerIsmom !== undefined && input.quadBoundaryLayers === false, ...gridControls } = {}) {
  // The standalone transonic solve needs the same maintained-residual check
  // as the direct coupled precursor. Physical admissibility alone can admit
  // a diverging shock update until the outlet reaches its sonic boundary.
  stepAcceptance ??= input.eulerIsmom === undefined ? 'admissible' : 'armijo';
  if (!['admissible', 'armijo'].includes(stepAcceptance)) throw new Error('Invalid Euler step acceptance.');
  if (typeof adaptiveMcrit !== 'boolean') throw new Error('Invalid Euler adaptive MCRIT control.');
  if (typeof firstOrderStartup !== 'boolean') throw new Error('Invalid Euler first-order startup control.');
  if (onEulerPrepared !== undefined && typeof onEulerPrepared !== 'function') throw new Error('Invalid prepared Euler observer.');
  if (typeof retainBestCheckpoint !== 'boolean') throw new Error('Invalid Euler checkpoint retention control.');
  if (preparedEuler !== undefined && Object.keys(gridControls).length)
    throw new Error('A prepared Euler mesh cannot be combined with grid-construction options.');
  const eulerIsmom = input.eulerIsmom, explicitEquations = eulerIsmom !== undefined;
  const prepared = preparedEuler === undefined ? prepareStreamtubeAssembly(input, { ...gridControls, meshOnly, onMesh })
    : restorePreparedStreamtubeAssembly(preparedEuler, input);
  if (prepared.experimental) return prepared;
  if (preparedEuler !== undefined) onMesh?.(prepared.mesh, 'initial');
  let initialMesh = prepared.mesh;
  if (meshOnly) return { status: 'mesh-ready', mesh: initialMesh };
  if (onEulerPrepared !== undefined) onEulerPrepared(capturePreparedStreamtubeAssembly(prepared, input));
  const { system, referenceChord, momentReference, mach } = prepared;
  const coefficientsAt = (flow, currentSystem) => {
    const coefficients = explicitEquations
      ? streamtubeSolidPressureForces({ flow, layout: currentSystem.layout, conditions: currentSystem.conditions, referenceChord, momentReference })
      : streamtubeForceCoefficients(flow.diagnosticForces ? flow : { ...flow, ...streamtubeSurfaceDiagnostics(currentSystem, flow) },
        currentSystem.conditions, { referenceChord, momentReference });
    if (explicitEquations) coefficients.cd = coefficients.pressureIntegralDrag;
    return coefficients;
  };
  const publishMesh = (state, mesh, phase) => {
    if (!onMesh) return;
    // Read the accepted flow already used for this mesh. No extra residual
    // evaluation, checkpoint, or force calculation during line-search trials.
    const coefficientProgress = { kind: 'euler-pressure', iteration: state.iteration.iteration,
      mach: state.system.conditions.mach, actualAlpha: state.system.conditions.alpha,
      gridLevel: input.gridIntervals };
    try { coefficientProgress.coefficients = coefficientsAt(state.flow, state.system); }
    catch (error) { coefficientProgress.message = error.message; }
    onMesh({ ...mesh, coefficientProgress }, phase);
  };
  let initial, initialFlow, gasInitialization;
  try {
    const startup = prepared.panelVelocitySeed
      ? initializeStreamtubePanelVelocity(system, prepared.initial, prepared.panelVelocitySeed)
      : initializeStreamtubeStartup(system, prepared.initial);
    initial = startup.initial; initialFlow = startup.flow; gasInitialization = startup.diagnostics;
  }
  catch (error) {
    if (error.code === 'streamtube-sonic-capacity') {
      const excess = 100 * (error.diagnostics.capacityRatio - 1);
      const fallback = error.diagnostics.stagnationDensityFallback;
      if (explicitEquations)
        throw Object.assign(new Error(`Cannot initialize ISMOM ${eulerIsmom} at Mach ${mach.toFixed(2)}: the exact-isentropic starting guess reaches or exceeds the mesh's sonic flow capacity (${excess.toFixed(1)}% excess at the most restrictive section).${fallback ? ` Uniform stagnation-density initialization was also rejected: ${fallback.reason}` : ''} Euler iterations have not started. The selected equations and initial mesh are retained; no alternate equation mode was attempted.`, { cause: error }), {
          code: error.code, diagnostics: error.diagnostics, stage: 'gas-initialization',
        });
      throw Object.assign(new Error(`Cannot initialize subsonic flow at Mach ${mach.toFixed(2)}: the exact-isentropic starting guess reaches or exceeds the mesh's sonic flow capacity (${excess.toFixed(1)}% excess at the most restrictive section).${fallback ? ` Uniform stagnation-density initialization was also rejected: ${fallback.reason}` : ''} Euler iterations have not started. This failure concerns the starting guess at the reported Mach; no higher-Mach flow has been solved. The initial mesh remains available for inspection.`, { cause: error }), {
        code: error.code, diagnostics: error.diagnostics, stage: 'gas-initialization',
      });
    }
    throw Object.assign(new Error(`Quad Euler initialization failed: ${error.message} The initial mesh remains available for inspection.`, { cause: error }), {
      stage: 'gas-initialization', ...(error.code === undefined ? {} : { code: error.code }),
      ...(error.diagnostics === undefined ? {} : { diagnostics: error.diagnostics }),
    });
  }
  prepared.diagnostics = { ...prepared.diagnostics, gasInitialization };
  initialMesh.initialization.gasInitialization = gasInitialization;
  // Publish the initialized flow before the first linear solve. Keep the
  // original mesh-only snapshot distinct from this provisional flow estimate.
  let initialFlowSnapshot = streamtubeFlowSnapshot(initialFlow);
  let previousMesh = { ...initialMesh, flow: compareStreamtubeFlow(initialFlowSnapshot, initialFlowSnapshot) };
  onMesh?.(previousMesh, 'solving');
  // The transonic fallback gas seed has no resolved shock. Temporarily omit
  // second-order anti-dissipation, then restore it on this same grid. Smooth
  // isentropic starts retain their existing path.
  const broadShockStartup = resolveStreamtubeStartup(input) === 'harmonic-shock';
  const useFirstOrderStartup = broadShockStartup || firstOrderStartup && explicitEquations
    && gasInitialization.method === 'stagnation-density';
  const startupAttempts = [];
  let flow;
  for (let attempt = 0; attempt < (broadShockStartup ? 2 : 1); attempt++) {
    const recoverStagnation = explicitEquations && mach <= .3 && input.elements.length > 1 && !broadShockStartup;
    const observers = {
      onIteration: h => onIteration?.({ ...h, stage: 'quad-euler', mach, startupAttempt: attempt + 1,
        ...(attempt ? { startupStrategy: 'broader-shock-retry' } : {}) }),
      onMesh: state => {
        const mesh = streamtubeMeshSnapshot({ ...state, diagnostics: prepared.diagnostics });
        // ISES startup redistribution is initialization, not a Newton update.
        // Compare subsequent motion and speed to the grid/flow it actually uses.
        if (state.iteration.iteration === 0) {
          initialFlowSnapshot = mesh.flow; initialMesh = mesh;
          mesh.flow = compareStreamtubeFlow(mesh.flow, initialFlowSnapshot);
          previousMesh = mesh; publishMesh(state, mesh, 'initial'); return;
        }
        mesh.flow = compareStreamtubeFlow(mesh.flow, initialFlowSnapshot, previousMesh.flow);
        mesh.iteration.maximumNodeMovement = mesh.vertices.reduce((maximum, p, i) => Math.max(maximum,
          Math.hypot(p.x - previousMesh.vertices[i].x, p.y - previousMesh.vertices[i].y)), 0);
        mesh.iteration.maximumNodeMovementFromInitial = mesh.vertices.reduce((maximum, p, i) => Math.max(maximum,
          Math.hypot(p.x - initialMesh.vertices[i].x, p.y - initialMesh.vertices[i].y)), 0);
        mesh.iteration.movementReferenceChord = referenceChord;
        previousMesh = mesh; publishMesh(state, mesh, 'solving');
      } };
    flow = solveStreamtubeIses(prepared.input, { initialEuler: { x: initial, nodes: initialFlow.nodes },
      maxIterations, tolerance, iterationGeometry: 'ises-sampled', stepAcceptance, retainBestCheckpoint,
      adaptiveMcrit: broadShockStartup || adaptiveMcrit,
      // A normal-only Newton direction can exhaust the convexity margin at
      // an inlet cut. Allow the compensated tangential direction when that
      // happens; it still has to pass the full grid and Armijo checks.
      ...(stepAcceptance === 'armijo' ? { gridCorrectionBacktracking: 'before-damping', tangentialGridRecovery: true } : {}),
      ...(useFirstOrderStartup ? { firstOrderStartup: true } : {}),
      ...(broadShockStartup ? { broadShockStartup: true, broadShockMucon: attempt ? 4 : 2,
        stopOnGridStagnation: true, gridCorrectionBacktracking: 'before-damping' } : {}),
      ...(recoverStagnation ? { stopOnGridStagnation: true } : {}),
      ...(explicitEquations ? { retainCheckpoint: true } : {}), ...observers });
    if (recoverStagnation && flow.lastRejectedStep?.code === 'EULER_GRID_STAGNATION')
      flow = recoverStreamtubeStagnationMotion(flow, { maxIterations, tolerance, ...observers });
    if (recoverStagnation && input.gridEllipticSmoothing && ['EULER_GRID_STAGNATION', 'EULER_PASSAGE_STAGNATION'].includes(flow.lastRejectedStep?.code))
      flow = recoverStreamtubePassageGrid(flow, prepared, { maxIterations, tolerance, ...observers });
    startupAttempts.push({ attempt: attempt + 1, mucon: flow.history[0].dissipation?.mucon,
      iterations: flow.history.length - 1, residual: flow.diagnostics.residual, reason: flow.reason });
    // Retry only a diagnosed temporary-law grid stall. A fresh gas state on
    // the SAME prepared grid avoids attempting to rescue an almost folded
    // accepted mesh by changing viscosity in place. Never alter the target.
    if (flow.lastRejectedStep?.code !== 'EULER_GRID_STAGNATION' || !(flow.solverInput.upwind?.mucon < 0)) break;
  }
  if (explicitEquations)
    flow.formulation = `Research simultaneous quadrilateral Euler with explicit ISMOM ${eulerIsmom}, shared-speed upwinding, additive-density Newton and ISES grid maintenance; final convexity and equation residual gates required. No boundary-layer coupling; shock and force accuracy remain unvalidated.`;
  const mesh = streamtubeMeshSnapshot({ system, nodes: flow.nodes, flow, diagnostics: prepared.diagnostics, iteration: flow.history.at(-1) });
  if (previousMesh.iteration?.iteration === mesh.iteration.iteration) {
    mesh.iteration = { ...previousMesh.iteration }; mesh.flow = previousMesh.flow;
  } else mesh.flow = compareStreamtubeFlow(mesh.flow, initialFlowSnapshot, previousMesh.flow);
  const converged = flow.converged && mesh.quality.valid;
  const coefficients = coefficientsAt(flow, system);
  mesh.initialization.flowSolved = converged;
  return { model: 'research-streamtube-euler', status: converged ? 'research-converged' : 'unconverged',
    mach, alpha: system.conditions.alpha, referenceChord, momentReference,
    cl: coefficients.cl, cm: coefficients.cm, cd: coefficients.cd, coefficients,
    coefficientStatus: converged ? 'research-unvalidated' : 'unconverged', boundaryLayer: null, mesh, flow,
    conditions: system.conditions, bodies: prepared.input.bodies, solverInput: flow.solverInput,
    diagnostics: { ...flow.diagnostics, equationResidual: flow.diagnostics.residual,
      ...streamtubeConvergenceDiagnostics(flow, mesh.quality),
      startupAttempts, totalIterations: startupAttempts.reduce((n, a) => n + a.iterations, 0), iterations: flow.history.length - 1, unknowns: system.layout.n, cells: mesh.cells.length },
    solverSettings: { stepMethod: 'density-newton', streamwiseMode: flow.streamwiseMode,
      ...(explicitEquations ? { eulerIsmom,
        hybrid: { ...system.conditions.hybrid },
        upwind: structuredClone(system.conditions.upwind) } : {}),
      jacobianBackend: flow.jacobianBackend, linearBackend: flow.linearBackend,
      ...(explicitEquations ? { adaptiveMcrit: broadShockStartup || adaptiveMcrit } : {}),
      ...(explicitEquations ? { requestedEulerStartup: input.eulerStartup ?? 'standard', eulerStartup: resolveStreamtubeStartup(input) } : {}),
      ...(broadShockStartup ? { broadShockStartup: true } : {}),
      ...(useFirstOrderStartup ? { firstOrderStartup: true } : {}),
      stagnationLimiter: flow.stagnationLimiter, iterationGeometry: flow.iterationGeometry,
      stepAcceptance: flow.stepAcceptance, ...(flow.gridAcceptance ? { gridAcceptance: flow.gridAcceptance } : {}), tolerance, maxIterations },
    forceStatus: `Unvalidated pressure coefficients from the ${converged ? 'converged research state' : 'last unconverged iterate'}. No viscous drag.`,
    limitations: explicitEquations
      ? `Research quadrilateral Euler with explicit ISMOM ${eulerIsmom} and speed upwinding. No boundary layers or viscous wakes. Local equation selection, shock structure, conservation defects and force/refinement accuracy require independent validation.`
      : 'Research subcritical isentropic quadrilateral Euler only. No boundary layers, viscous wakes or shocks. Streamwise isentropy replaces streamwise momentum; finite momentum defects and force/refinement accuracy remain unvalidated.',
    warnings: [...(!converged ? ['Unconverged coefficients are provisional values from the last iterate and may change substantially.'] : []),
      'Euler equation convergence does not establish physical or force accuracy. Drag is pressure-only; this workbench mode has no boundary-layer coupling.'] };
}
