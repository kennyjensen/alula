// SPDX-License-Identifier: GPL-2.0-or-later
// Continue a converged public coupled result on the tested paired nested grid.
// No cold mesh/BL initialization or initial SMOVE is part of this path.
import { coupledResultShearCoordinate } from './streamtube-coupled-shear-policy.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from './streamtube-coupled-refinement.js';
import { nestedRefinementCheckpoint, checkpointDataEqual } from './streamtube-nested-checkpoint.js';
import { coupledCheckpointHkPolicy, coupledAssemblyConditions } from './streamtube-coupled-assembly.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { coupledConvergenceSatisfied } from './streamtube-coupled-convergence.js';

export const COUPLED_REFINEMENT_MAX_NODES = 50000;
const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const require = (condition, message) => { if (!condition) throw new Error(message); };

// Cheap guards and the node budget precede construction or evaluation of a
// numerical system. sourceCase travels with the result across worker restarts.
export function coupledRefinementPlan(caseData, parentResult, { maxNodes = COUPLED_REFINEMENT_MAX_NODES,
  streamwiseInterpolation = 'linear' } = {}) {
  require(['linear', 'surface-pchip'].includes(streamwiseInterpolation), 'Unknown coupled streamwise refinement interpolation.');
  require(caseData?.flowModel === 'streamtube-grid' && caseData.quadBoundaryLayers === true,
    'Refinement requires the coupled quadrilateral flow model.');
  require(parentResult?.model === 'research-streamtube-euler-bl' && parentResult.converged === true
    && parentResult.mesh?.quality?.valid === true, 'Refinement requires a converged coupled solution with a valid grid.');
  require(parentResult.sourceCase && checkpointDataEqual(parentResult.sourceCase, caseData),
    'The coupled solution is stale: solve the current unchanged case before refining.');
  const checkpoint = parentResult.checkpoint, input = checkpoint?.restart?.input;
  require(checkpoint?.version === 1 && input && checkpoint.continuation && checkpoint.restart.options
    && checkpoint.restart.initialEuler?.x?.length && checkpoint.restart.initialBL?.length,
    'The coupled solution has no complete continuation checkpoint; solve it again before refining.');
  coupledCheckpointHkPolicy(caseData, checkpoint.restart.options);
  require(checkpointDataEqual(checkpoint.families, parentResult.families), 'The coupled result and checkpoint residuals differ.');
  const tolerance = parentResult.solverSettings?.tolerance;
  require(Number.isFinite(tolerance) && tolerance > 0 && coupledConvergenceSatisfied(parentResult, tolerance),
    'The parent checkpoint does not satisfy its coupled convergence tolerance.');
  require(Array.isArray(input.bodies) && input.bodies.length > 0 && Array.isArray(input.weights)
    && input.weights.length === input.bodies.length + 1 && input.weights.every(row => Array.isArray(row) && row.length >= 2)
    && Array.isArray(input.outerLower) && input.outerLower.length >= 3, 'The coupled checkpoint has invalid grid dimensions.');
  require(Number.isInteger(maxNodes) && maxNodes > 0 && maxNodes <= COUPLED_REFINEMENT_MAX_NODES,
    'Invalid coupled refinement node budget.');
  const nx = input.outerLower.length - 1, tubes = input.weights.map(row => row.length);
  const normalSubdivisions = tubes.map(n => Array(n).fill(1));
  for (let b = 0; b < input.bodies.length; b++) {
    normalSubdivisions[b][tubes[b] - 1] = 4;
    normalSubdivisions[b + 1][0] = 4;
  }
  const targetTubes = normalSubdivisions.map(row => row.reduce((sum, n) => sum + n, 0));
  const nodeCount = (2 * nx + 1) * targetTubes.reduce((sum, n) => sum + n + 1, 0);
  require(nodeCount <= maxNodes, `Coupled refinement needs ${nodeCount.toLocaleString()} nodes, exceeding the ${maxNodes.toLocaleString()} node budget.`);
  return { streamwiseFactor: 2, streamwiseSubdivisions: Array(nx).fill(2), normalSubdivisions,
    normalInterpolation: 'streamfunction-quadratic', maxNodes, nodeCount, tolerance,
    ...(streamwiseInterpolation === 'surface-pchip' ? { streamwiseInterpolation } : {}),
    parent: { nx, tubes }, child: { nx: 2 * nx, tubes: targetTubes } };
}

export function prepareCoupledStreamtubeRefinement(caseData, parentResult, { maxNodes, onMesh, onStage, streamwiseInterpolation } = {}) {
  const plan = coupledRefinementPlan(caseData, parentResult, { maxNodes, streamwiseInterpolation });
  onStage?.({ stage: 'coupled-refinement', startupAttempt: 1 });
  const checkpoint = serialize(parentResult.checkpoint), f = checkpoint.restart;
  const parent = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const settings = coupledAssemblyConditions(caseData, f.input.bodies, parent.euler.conditions.lengthScale, { checkpointOptions: f.options });
  require(checkpointDataEqual(settings.normalization, {
    referenceChord: parentResult.referenceChord, referenceReynolds: parentResult.referenceReynolds,
    solverLength: parentResult.solverLength, kernelReynolds: parentResult.kernelReynolds,
  }), 'Refinement changed the public physical normalization.');
  for (const key of ['reynolds', 'ncrit', 'edgeMatching', 'tripFractions'])
    require(checkpointDataEqual(settings.options[key], f.options[key]), `Refinement case and checkpoint disagree on ${key}.`);
  require((settings.options.transitionMode ?? 'fixed-trip') === parent.bl.transitionMode,
    'Refinement case and checkpoint disagree on transition mode.');
  require((caseData.alpha ?? 0) === parent.euler.conditions.alpha && (caseData.mach ?? .2) === parent.euler.conditions.mach,
    'Refinement case and checkpoint disagree on incidence or Mach number.');
  const mapped = refineCoupledStreamtubeBody(f.input, parent, {
    streamwiseFactor: 2, normalSubdivisions: plan.normalSubdivisions,
    normalInterpolation: plan.normalInterpolation, maxNodes: plan.maxNodes,
    ...(plan.streamwiseInterpolation ? { streamwiseInterpolation: plan.streamwiseInterpolation } : {}),
  });
  const child = mapped.system, value = child.evaluate(child.initial);
  // Export the final phase prepared by the refiner, exactly as in the validated
  // paired seed producer, rather than its pre-preparation option snapshot.
  const restart = serialize({ input: mapped.input, options: { ...mapped.options, transitionState: child.bl.snapshotActive() },
    initialEuler: { x: child.initial.slice(0, child.ne), nodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes },
    initialBL: child.initial.slice(child.ne) });
  const mesh = streamtubeMeshSnapshot({ system: child.euler, nodes: value.outer.nodes, flow: value.outer,
    iteration: { iteration: 0, ...value.families, residual: Math.max(...Object.values(value.families)) } });
  mesh.initialization.gridSmoothing = parentResult.initialization?.euler?.gridSmoothing;
  mesh.initialization.gridRefinement = mapped.diagnostics;
  // This first transfer preview precedes history certification and Newton.
  onMesh?.(mesh, 'initial');
  const transfer = nestedRefinementCheckpoint(checkpoint, restart, {
    streamwiseSubdivisions: plan.streamwiseSubdivisions, normalSubdivisions: plan.normalSubdivisions,
    ...(plan.streamwiseInterpolation ? { streamwiseInterpolation: plan.streamwiseInterpolation } : {}),
  });
  return { checkpoint: transfer.checkpoint, mesh, settings, plan, system: child,
    refinement: { kind: 'paired-nested', level: (parentResult.refinement?.level ?? 0) + 1,
      streamwiseFactor: 2, wallTubeSubdivisions: 4, normalInterpolation: plan.normalInterpolation,
      ...(plan.streamwiseInterpolation ? { streamwiseInterpolation: plan.streamwiseInterpolation } : {}),
      nodeCount: plan.nodeCount, parent: { ...plan.parent, unknowns: parent.n, cells: parentResult.mesh.cells.length },
      child: { ...plan.child, unknowns: child.n, cells: mesh.cells.length },
      parentCoefficients: { cl: parentResult.cl, cd: parentResult.cd, cm: parentResult.cm },
      initialization: mapped.diagnostics, transfer: transfer.diagnostics } };
}

export function solveCoupledStreamtubeRefinement(caseData, parentResult, {
  maxIterations = caseData?.maxIterations ?? (caseData?.transitionMode === 'automatic' ? 40 : 20),
  onMesh, onIteration, onStage, onCheckpoint, onPrepared,
  streamwiseInterpolation,
} = {}) {
  require(Number.isInteger(maxIterations) && maxIterations >= 0, 'Invalid coupled refinement iteration limit.');
  const prepared = prepareCoupledStreamtubeRefinement(caseData, parentResult, { onMesh, onStage, streamwiseInterpolation });
  const { checkpoint, settings, refinement, plan } = prepared;
  onPrepared?.(prepared);
  onStage?.({ stage: 'coupled', startupAttempt: 1, refinement: true });
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  const result = solveCoupledStreamtubeIses(undefined, { resume: checkpoint, maxIterations, tolerance: plan.tolerance,
    convergence: parentResult.solverSettings?.convergence ?? 'residual',
    iterationGeometry, stepAcceptance, stagnationLimiter, onCheckpoint,
    onIteration: h => onIteration?.({ ...h, stage: 'coupled', startupAttempt: 1, refinement: true }),
    onMesh: state => onMesh?.(streamtubeMeshSnapshot({ ...state, iteration: { ...state.iteration, startupAttempt: 1 } }),
      state.iteration.iteration === 0 ? 'initial' : 'solving', state) });
  const attempt = { kind: 'paired-nested-refinement', startupAttempt: 1, thicknessFactor: 1,
    converged: result.converged, reason: result.reason, families: result.families, quality: result.mesh.quality,
    iterations: result.history.length - 1, history: result.history, lastRejectedStep: result.lastRejectedStep };
  const solverSettings = { ...parentResult.solverSettings, maxIterations, tolerance: plan.tolerance, refinement: 'paired-nested',
    shearCoordinate: coupledResultShearCoordinate(result),
    hkFloorLinearization: coupledCheckpointHkPolicy(caseData, result.checkpoint?.restart.options
      ?? result.coupledOptions ?? result.conditions ?? checkpoint.restart.options) };
  const projectionGeometry = result.checkpoint
    ? result.checkpoint.continuation?.projectionGeometry ?? 'fixed' : result.projectionGeometry ?? 'fixed';
  if (projectionGeometry === 'boundary-increment' || solverSettings.projectionGeometry !== undefined)
    solverSettings.projectionGeometry = projectionGeometry;
  // This label describes this transfer, not the parent's previous transfer.
  if (plan.streamwiseInterpolation) solverSettings.streamwiseInterpolation = plan.streamwiseInterpolation;
  else delete solverSettings.streamwiseInterpolation;
  return { ...result, model: 'research-streamtube-euler-bl', physicalAcceptance: false,
    alpha: caseData.alpha ?? 0, mach: caseData.mach ?? .2, ...settings.normalization,
    materialTrips: settings.materialTrips, elementOrder: settings.elementOrder, bodies: result.solverInput.bodies,
    ...(result.checkpoint ? { restart: result.checkpoint.restart } : {}), sourceCase: structuredClone(caseData), refinement,
    initialization: { euler: structuredClone(parentResult.initialization.euler),
      boundaryLayer: { method: 'paired-nested-refinement', thicknessFactor: 1, source: 'Converged parent boundary layers and wakes.' },
      attempts: [attempt] },
    solverSettings,
    coefficientStatus: 'unavailable', forceStatus: 'Solid-wall forces and viscous drag are not validated for this model.' };
}
